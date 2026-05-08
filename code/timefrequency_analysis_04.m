% SF_TF_longshort.m
%
% Time-frequency (TF) analysis for the sensory feedback study.
% Part 1: Morlet wavelet convolution on epoched EEG → saves [nSub x nFrex x nTimes x nCh] .mat per group.
% Part 2: Loads TF .mat files and applies dB baseline normalisation (-1000 to -200 ms).
% Parts 3–6: Topoplots, band-power maps, and time-series figures for
%   alpha (8–13 Hz) and beta (13–30 Hz) bands, comparing PF vs VF groups.
%
% Expected folder layout (relative to code/):
%   ../derivatives/epoched/PF/  — preprocessed PF epochs (.set)
%   ../derivatives/epoched/VF/  — preprocessed VF epochs (.set)
%   ../derivatives/TF/          — TF .mat output (Part 1 write / Part 2 read)
%   ../figs/                    — figure output
%
% Dependencies: EEGLAB

clear all;
close all;
clc;

%% Shared parameters

tf_path_save     = '../derivatives/TF/';
tf_path_load     = './';
figs_path = '../figs/';

nSubject         = 21;
nChannels        = 60;
min_freq         = 2;
max_freq         = 80;
num_frex         = 50;
frex             = logspace(log10(min_freq), log10(max_freq), num_frex);
range_cycles     = [2 13];
SR               = 1000;
s                = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex) ./ (2*pi*frex);
wavtime          = -2:1/SR:2;
half_wave        = (length(wavtime)-1)/2;
times_full       = -3000:1:4999;   % full resolution (Part 1)
times            = -3000:2:4999;   % downsampled 2:1  (Parts 2-3)
nTimes           = length(times);
baseline_windows = [-1000 -200];

% Group configs: path, condition labels, output .mat name
% Each group has one merged epoch file per subject (all articulations).
% Produced by SPF_pre_longshort.m / SVF_pre_longshort.m.
groups(1).name       = 'PF';
groups(1).epoched    = '../derivatives/epoched/PF/';
groups(1).matfile    = 'tf_pf_merged_ica';
groups(1).events     = {'PF'};
groups(1).cond_names = {'pf_exp'};

groups(2).name       = 'VF';
groups(2).epoched    = '../derivatives/epoched/VF/';
groups(2).matfile    = 'tf_vf_merged_ica';
groups(2).events     = {'VF'};
groups(2).cond_names = {'vf_exp'};

nGroups     = numel(groups);
nTimes_full = length(times_full);


%% Part 1 — Morlet wavelet convolution

for g = 1:nGroups

    disp(['=== Part 1: Morlet convolution — ' groups(g).name ' ===']);

    nEvents = numel(groups(g).events);
    tf = zeros(nSubject, num_frex, nTimes_full, nChannels, nEvents, 'single');
    fileID = fopen(['loop_indices_' groups(g).name '.log'], 'w');

    for sub = 1:nSubject
        disp(sub)
        for evt = 1:nEvents

            eeg_file = ['SF' int2str(sub) '_' groups(g).events{evt} '.set'];
            EEG      = pop_loadset('filename', eeg_file, 'filepath', groups(g).epoched);

            for ch = 1:nChannels

                nWave   = length(wavtime);
                nData   = EEG.pnts * EEG.trials;
                nConv   = nWave + nData - 1;
                alldata = reshape(EEG.data(ch,:,:), 1, []);
                dataX   = fft(alldata, nConv);

                for fi = 1:num_frex
                    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2 ./ (2*s(fi)^2));
                    waveletX = fft(wavelet, nConv);
                    waveletX = waveletX ./ max(waveletX);

                    as = ifft(waveletX .* dataX);
                    as = as(half_wave+1:end-half_wave);
                    as = reshape(as, EEG.pnts, EEG.trials);

                    tf(sub, fi, :, ch, evt) = single(mean(abs(as).^2, 2));
                end

            end
        end
        fprintf(fileID, 'sub: %d, evt: %d\n', sub, evt);
    end

    % Downsample 2:1 in time
    tf = tf(:,:,1:2:end,:,:);
    fclose(fileID);

    filename = [tf_path_save groups(g).matfile];
    save(filename, 'tf', '-v7.3');
    disp(['Saved: ' filename]);

end


%% Part 2 — Load .mat files
% Run this section directly if .mat files already exist (skip Part 1)
tf_path_load     = '../derivatives/TF/';

for g = 1:nGroups
    disp(['Loading: ' groups(g).matfile '.mat']);
    data = load([ tf_path_load groups(g).matfile '.mat'], 'tf');
    groups(g).tf = data.tf;
    disp(['Done: ' groups(g).name]);
end


%% Part 2 — Baseline normalisation (dB)

tf_db = struct();

baseidx = reshape(dsearchn(times', baseline_windows(:)), [], 2);

for g = 1:nGroups

    disp(['=== Part 2: Baseline normalisation — ' groups(g).name ' ===']);

    tf = groups(g).tf;

    nEvents = numel(groups(g).events);
    for evt = 1:nEvents

        tf_cond = tf(:,:,:,:,evt);   % [nSub x nFrex x nTimes x nCh]
        tf_db_cond = zeros(nSubject, num_frex, nTimes, nChannels);

        for sub = 1:nSubject
            for ch = 1:nChannels
                activity = squeeze(tf_cond(sub,:,:,ch));
                baseline = mean(tf_cond(sub, :, baseidx(1):baseidx(2), ch), 3)';
                tf_db_cond(sub,:,:,ch) = 10*log10(bsxfun(@rdivide, activity, baseline));
            end
        end

        % Store in struct with meaningful name, e.g. tf_db.pf_exp_long
        tf_db.(groups(g).cond_names{evt}) = tf_db_cond;

    end

end

% Workspace now contains:
%   tf_db.pf_exp  [21 x 50 x 4000 x 60]  PF, all feedback trials (long + short merged)
%   tf_db.vf_exp  [21 x 50 x 4000 x 60]  VF, all feedback trials (long + short merged)


%% Part 3 — Topoplot diagnostics (avg over subjects after normalisation)

load EEG_chlocs.mat
close all;
fList = find(frex >= 13 & frex <= 30);  % alpha band
z1   = -3;
z2   =  3;
t0   = find(times == 0);
tbs  = find(times == -1000);
n    = 12;       % number of time windows per row
step = 125;      % samples per time window (= 250 ms at 500 Hz after 2:1 downsample)
c    = jet;
figs_path = '../figs/';

for g = 1:nGroups

    disp(['=== Part 3: Topoplots — ' groups(g).name ' ===']);

    all_conds = cellfun(@(name) tf_db.(name), groups(g).cond_names, 'UniformOutput', false);
    num_conds = numel(all_conds);

    figure(g);
    clf;

    for k = 1:num_conds

        current = all_conds{k};

        % Baseline column
        subplot(num_conds, n+1, (k-1)*(n+1) + 1);
        temp_base = mean(mean(squeeze(mean(current(:, fList, tbs:tbs+400, :), 1))));
        topoplot(squeeze(temp_base)', EEG_chlocs, 'maplimits', [z1 z2], 'electrodes', 'off');
        title([groups(g).name ' Base']);

        % Time-window columns
        for i = 1:n
            time_window = t0 + ((i-1)*step) : t0 + (i*step) - 1;
            subplot(num_conds, n+1, (k-1)*(n+1) + (i+1));
            temp_event = mean(mean(squeeze(mean(current(:, fList, time_window, :), 1))));
            topoplot(squeeze(temp_event)', EEG_chlocs, 'maplimits', [z1 z2], 'electrodes', 'off');
            if k == 1
                title([num2str(i * 250) ' ms'], 'FontSize', 7);
            end
        end

    end

    colormap(c);
    sgtitle(sprintf('%s — All Subjects Average: Color Range %g to %g dB', groups(g).name, z1, z2));
    set(gcf, 'Position', [100, 50, 1450, 550]);
    set(gcf, 'Color', 'w');
    drawnow;
    print(gcf, [figs_path 'supp3_' groups(g).name '_beta.svg'], '-dsvg');

end


%% Part 4 — Figure 4: alpha & beta topoplots, PF vs VF
% Two figures: active window (0–750ms) and late window (750–2000ms)
% Each: 2×2 grid — rows = alpha/beta, columns = PF/VF
load EEG_chlocs.mat
close all;

active_idx4 = find(times >= 0    & times <= 750);
late_idx4   = find(times >= 2000  & times <= 3000);
alpha_idx   = find(frex >= 8     & frex <= 13);
beta_idx    = find(frex >= 13    & frex <= 30);

% Average across subjects → [nFrex x nTimes x nCh]
pf_grand = squeeze(mean(tf_db.pf_exp, 1));
vf_grand = squeeze(mean(tf_db.vf_exp, 1));


% Topography vectors per band × window → [nCh x 1]
topo = @(grand, fidx, tidx) squeeze(mean(mean(grand(fidx, tidx, :), 1), 2));

pf_alpha_active = topo(pf_grand, alpha_idx, active_idx4);
pf_beta_active  = topo(pf_grand, beta_idx,  active_idx4);
vf_alpha_active = topo(vf_grand, alpha_idx, active_idx4);
vf_beta_active  = topo(vf_grand, beta_idx,  active_idx4);

pf_alpha_late   = topo(pf_grand, alpha_idx, late_idx4);
pf_beta_late    = topo(pf_grand, beta_idx,  late_idx4);
vf_alpha_late   = topo(vf_grand, alpha_idx, late_idx4);
vf_beta_late    = topo(vf_grand, beta_idx,  late_idx4);

% Electrode cluster markers
contra_marker_labels = {'C1', 'C3', 'CP1', 'CP3'};
ipsi_marker_labels   = {'C2', 'C4', 'CP2', 'CP4'};
contra_marker_chs = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), contra_marker_labels);
ipsi_marker_chs   = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), ipsi_marker_labels);
all_marker_chs    = [contra_marker_chs, ipsi_marker_chs];

clim_val   = 3;
col_titles = {'PF', 'VF'};
row_labels4 = {'Alpha  8–13 Hz', 'Beta  13–30 Hz'};

windows = { ...
    {pf_alpha_active, vf_alpha_active; pf_beta_active, vf_beta_active}, 'Stimulation', 'figure4a'; ...
    {pf_alpha_late,   vf_alpha_late;   pf_beta_late,   vf_beta_late},   'Rebound',   'figure4b'  ...
};

for w = 1:2
    plot_data4 = windows{w, 1};
    win_label  = windows{w, 2};
    fname      = windows{w, 3};

    fig4 = figure('Units', 'centimeters', 'Position', [5 5 15 9], 'Color', 'w');
    clf;
    ax4 = gobjects(2, 2);

    for row = 1:2
        for col = 1:2
            ax4(row, col) = subplot(2, 2, (row-1)*2 + col);
            topoplot(plot_data4{row, col}, EEG_chlocs, ...
                'maplimits', [-clim_val clim_val], ...
                'electrodes', 'off');
            if row == 1
                title(col_titles{col}, 'FontSize', 10, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(row_labels4{row}, 'FontSize', 9, 'FontWeight', 'bold');
            end
        end
    end

    sgtitle(sprintf('Window: %s', win_label), 'FontSize', 10);
    colormap(fig4, jet);

    for row = 1:2
        for col = 1:2
            pos = get(ax4(row, col), 'Position');
            set(ax4(row, col), 'Position', [pos(1), pos(2), pos(3)*0.82, pos(4)]);
        end
    end

    ax_cb4 = axes(fig4, 'Position', [0.86 0.1 0.01 0.8], 'Visible', 'off');
    caxis(ax_cb4, [-clim_val clim_val]);
    cb4 = colorbar(ax_cb4, 'eastoutside');
    cb4.Position   = [0.85 0.1 0.03 0.8];
    cb4.Ticks      = [-clim_val, 0, clim_val];
    cb4.FontSize   = 14;
    cb4.Label.FontSize = 14;

    drawnow;
    set(fig4, 'Color', 'none');
    print(fig4, [figs_path fname '.svg'], '-dsvg');
    set(fig4, 'Color', 'w');
end
% copygraphics(gcf, 'BackgroundColor', 'white');


%% Part 5 — Figure 4b: Time-frequency maps, CP3 & CP4, PF vs VF
% 2×2 grid: rows = contralateral, ipsilateral; columns = PF, VF
% Grand average across subjects; time -1000 to 3000 ms; frequency 2-50 Hz

tf_time_idx = find(times >= -1000 & times <= 3000);
tf_freq_idx = find(frex >= 2 & frex <= 50);
tf_times_ms = times(tf_time_idx);
tf_frex_hz  = frex(tf_freq_idx);
n_tf_freq   = numel(tf_frex_hz);

contra_labels = {'C1', 'C3', 'CP1', 'CP3'};
ipsi_labels   = {'C2', 'C4', 'CP2', 'CP4'};
contra_chs = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), contra_labels);
ipsi_chs   = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), ipsi_labels);

ftick_hz  = [2, 4, 8, 13, 20, 30, 50];
ftick_idx = dsearchn(tf_frex_hz', ftick_hz');
idx_8hz   = dsearchn(tf_frex_hz', 8);
idx_13hz  = dsearchn(tf_frex_hz', 13);

clim_tf     = 3;
col_titles5 = {'PF', 'VF'};

% Grand averages — [50 x 4000 x 60], reused by Part 6
pf_avg5 = squeeze(mean(tf_db.pf_exp, 1));
vf_avg5 = squeeze(mean(tf_db.vf_exp, 1));

pf_cp3_tf = squeeze(mean(pf_avg5(tf_freq_idx, tf_time_idx, contra_chs), 3));
pf_cp4_tf = squeeze(mean(pf_avg5(tf_freq_idx, tf_time_idx, ipsi_chs),   3));
vf_cp3_tf = squeeze(mean(vf_avg5(tf_freq_idx, tf_time_idx, contra_chs), 3));
vf_cp4_tf = squeeze(mean(vf_avg5(tf_freq_idx, tf_time_idx, ipsi_chs),   3));
plot_data5 = {pf_cp3_tf, vf_cp3_tf; pf_cp4_tf, vf_cp4_tf};

fig5 = figure('Units', 'centimeters', 'Position', [5 5 18 12], 'Color', 'w');
clf;
ax5 = gobjects(2, 2);

for row = 1:2
    for col = 1:2
        ax5(row, col) = subplot(2, 2, (row-1)*2 + col);
        imagesc(tf_times_ms, 1:n_tf_freq, plot_data5{row, col});
        axis xy;
        caxis([-clim_tf clim_tf]);
        set(gca, 'YTick', ftick_idx, ...
                 'YTickLabel', arrayfun(@num2str, ftick_hz, 'UniformOutput', false));
        set(gca, 'XTick', [-1000 -500 0 1000 2000 3000]);
        xlabel('Time (ms)', 'FontSize', 9);
        if col == 1, ylabel('Frequency (Hz)', 'FontSize', 9); end
        if row == 1, title(col_titles5{col}, 'FontSize', 10, 'FontWeight', 'bold'); end
        hold on;
        if col == 1
            text(tf_times_ms(1)-50, (idx_8hz+idx_13hz)/2, '\alpha', ...
                'Color','w','FontSize',11,'FontWeight','bold','HorizontalAlignment','right','Clipping','off');
            text(tf_times_ms(1)-50, (idx_13hz+n_tf_freq)/2, '\beta', ...
                'Color','w','FontSize',11,'FontWeight','bold','HorizontalAlignment','right','Clipping','off');
        end
        % xline(0,   '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
        % xline(750, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
        hold off;
    end
end

colormap(fig5, jet);
for row = 1:2
    for col = 1:2
        pos = get(ax5(row,col), 'Position');
        set(ax5(row,col), 'Position', [pos(1), pos(2), pos(3)*0.85, pos(4)]);
    end
end
ax_cb5 = axes(fig5, 'Position', [0.88 0.1 0.01 0.8], 'Visible', 'off');
caxis(ax_cb5, [-clim_tf clim_tf]);
cb5 = colorbar(ax_cb5, 'eastoutside');
cb5.Position       = [0.90 0.1 0.03 0.8];
cb5.Ticks          = [-clim_tf, 0, clim_tf];
cb5.FontSize       = 12;
cb5.Label.String   = 'Power [dB]';
cb5.Label.FontSize = 12;

set(fig5, 'Color', 'w');
drawnow;
print(fig5, [figs_path 'figure4b.pdf'], '-dpdf', '-bestfit');


%% Part 6 — Figure 4c: Band power timeseries, CP3 & CP4, PF vs VF
% Two panels (contralateral, ipsilateral), 4 lines each:
% alpha PF, alpha VF (solid); beta PF, beta VF (dashed)

col_pf6 = [0.92, 0.46, 0.32];
col_vf6 = [0.28, 0.60, 0.42];

alpha_idx6 = find(frex >= 8  & frex <= 13);
beta_idx6  = find(frex >= 13 & frex <= 30);

% pf_avg5 / vf_avg5 from Part 5 — [50 x 4000 x 60]
pf_alpha_cp3 = squeeze(mean(mean(pf_avg5(alpha_idx6, tf_time_idx, contra_chs), 3), 1));
pf_beta_cp3  = squeeze(mean(mean(pf_avg5(beta_idx6,  tf_time_idx, contra_chs), 3), 1));
vf_alpha_cp3 = squeeze(mean(mean(vf_avg5(alpha_idx6, tf_time_idx, contra_chs), 3), 1));
vf_beta_cp3  = squeeze(mean(mean(vf_avg5(beta_idx6,  tf_time_idx, contra_chs), 3), 1));

pf_alpha_cp4 = squeeze(mean(mean(pf_avg5(alpha_idx6, tf_time_idx, ipsi_chs), 3), 1));
pf_beta_cp4  = squeeze(mean(mean(pf_avg5(beta_idx6,  tf_time_idx, ipsi_chs), 3), 1));
vf_alpha_cp4 = squeeze(mean(mean(vf_avg5(alpha_idx6, tf_time_idx, ipsi_chs), 3), 1));
vf_beta_cp4  = squeeze(mean(mean(vf_avg5(beta_idx6,  tf_time_idx, ipsi_chs), 3), 1));

panel_lines   = { {pf_alpha_cp3, vf_alpha_cp3, pf_beta_cp3, vf_beta_cp3}, ...
                  {pf_alpha_cp4, vf_alpha_cp4, pf_beta_cp4, vf_beta_cp4} };
panel_titles6 = {'Contralateral', 'Ipsilateral'};

fig6 = figure('Units', 'centimeters', 'Position', [5 5 18 7], 'Color', 'w');
clf;
ax6 = gobjects(1, 2);
lh6 = gobjects(4, 2);

for col = 1:2
    ax6(col) = subplot(1, 2, col);
    hold on;

    patch([0 750 750 0],         [-4 -4 1 1], ...
        [0.75 0.75 0.75], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    patch([2000 3000 3000 2000], [-4 -4 1 1], ...
        [0.70 0.80 0.95], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    lines = panel_lines{col};
    lh6(1,col) = plot(tf_times_ms, lines{1}, '-',  'Color', col_pf6, 'LineWidth', 1.8);
    lh6(2,col) = plot(tf_times_ms, lines{2}, '-',  'Color', col_vf6, 'LineWidth', 1.8);
    lh6(3,col) = plot(tf_times_ms, lines{3}, '--', 'Color', col_pf6, 'LineWidth', 1.8);
    lh6(4,col) = plot(tf_times_ms, lines{4}, '--', 'Color', col_vf6, 'LineWidth', 1.8);

    yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    hold off;

    set(gca, 'XTick', [-1000 -500 0 1000 2000 3000], 'FontSize', 9);
    xlim([tf_times_ms(1) tf_times_ms(end)]);
    xlabel('Time (ms)', 'FontSize', 10);
    if col == 1, ylabel('Power [dB]', 'FontSize', 10); end
    title(panel_titles6{col}, 'FontSize', 10, 'FontWeight', 'bold');
    box off;
end

legend(lh6(:,2), {'\alpha PF', '\alpha VF', '\beta PF', '\beta VF'}, ...
    'Location', 'southwest', 'FontSize', 9, 'Box', 'off');
ylim(ax6(1), [-4 1]);
ylim(ax6(2), [-4 1]);

drawnow;
set(fig6, 'Color', 'none');
print(fig6, [figs_path 'figure4c.svg'], '-dsvg');
set(fig6, 'Color', 'w');


%% Part 7 — Figure 4d: ERD boxplots, PF vs VF, active & late windows, alpha & beta
% Self-contained: requires tf_db in workspace (run Parts 1-2 or load saved tf_db).
% Per-subject TF data = average of long and short articulation conditions.
close all;
stats_file_eeg = [figs_path 'statistics_eeg.txt'];
if exist(stats_file_eeg, 'file'), delete(stats_file_eeg); end

% --- Shared parameters (redeclare for self-contained running) ---
frex  = logspace(log10(2), log10(80), 50);
times = -3000:2:4999;
load EEG_chlocs.mat;
contra_labels = {'C1', 'C3', 'CP1', 'CP3'};
ipsi_labels   = {'C2', 'C4', 'CP2', 'CP4'};
contra_chs = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), contra_labels);
ipsi_chs   = cellfun(@(l) find(strcmpi({EEG_chlocs.labels}, l)), ipsi_labels);

alpha_idx7 = find(frex >= 8  & frex <= 13);
beta_idx7  = find(frex >= 13 & frex <= 30);

col_pf7 = [0.92, 0.46, 0.32];
col_vf7 = [0.28, 0.60, 0.42];

% Per-subject band×window×cluster average → [21×1]
% tf_db.*: [21 × 50 × 4000 × 60]
ep = @(data, fi, ti, ch) squeeze(mean(mean(mean(data(:, fi, ti, ch), 4), 2), 3));

active_idx7 = find(times >= 0    & times <= 750);
late_idx7   = find(times >= 2000 & times <= 3000);
pf_data     = tf_db.pf_exp;
vf_data     = tf_db.vf_exp;

row_labels7  = {'\alpha power [dB]', '\beta power [dB]'};
col_titles7  = {'Stim. window', 'Rebound window', 'Stim. window', 'Rebound window'};
p_idx_map7   = [1, 2, 5, 6; 3, 4, 7, 8];
comp_labels7 = {'\alpha active contra', '\alpha late contra', ...
                '\beta active contra',  '\beta late contra', ...
                '\alpha active ipsi',   '\alpha late ipsi', ...
                '\beta active ipsi',    '\beta late ipsi'};
n_comp7 = 8;

diary(stats_file_eeg);

    pf_aa_cp3 = ep(pf_data, alpha_idx7, active_idx7, contra_chs);
    pf_al_cp3 = ep(pf_data, alpha_idx7, late_idx7,   contra_chs);
    vf_aa_cp3 = ep(vf_data, alpha_idx7, active_idx7, contra_chs);
    vf_al_cp3 = ep(vf_data, alpha_idx7, late_idx7,   contra_chs);

    pf_ba_cp3 = ep(pf_data, beta_idx7, active_idx7, contra_chs);
    pf_bl_cp3 = ep(pf_data, beta_idx7, late_idx7,   contra_chs);
    vf_ba_cp3 = ep(vf_data, beta_idx7, active_idx7, contra_chs);
    vf_bl_cp3 = ep(vf_data, beta_idx7, late_idx7,   contra_chs);

    pf_aa_cp4 = ep(pf_data, alpha_idx7, active_idx7, ipsi_chs);
    pf_al_cp4 = ep(pf_data, alpha_idx7, late_idx7,   ipsi_chs);
    vf_aa_cp4 = ep(vf_data, alpha_idx7, active_idx7, ipsi_chs);
    vf_al_cp4 = ep(vf_data, alpha_idx7, late_idx7,   ipsi_chs);

    pf_ba_cp4 = ep(pf_data, beta_idx7, active_idx7, ipsi_chs);
    pf_bl_cp4 = ep(pf_data, beta_idx7, late_idx7,   ipsi_chs);
    vf_ba_cp4 = ep(vf_data, beta_idx7, active_idx7, ipsi_chs);
    vf_bl_cp4 = ep(vf_data, beta_idx7, late_idx7,   ipsi_chs);

    plot_pairs7 = { ...
        [pf_aa_cp3, vf_aa_cp3], [pf_al_cp3, vf_al_cp3], [pf_aa_cp4, vf_aa_cp4], [pf_al_cp4, vf_al_cp4]; ...
        [pf_ba_cp3, vf_ba_cp3], [pf_bl_cp3, vf_bl_cp3], [pf_ba_cp4, vf_ba_cp4], [pf_bl_cp4, vf_bl_cp4]  ...
    };

    pf_all7 = {pf_aa_cp3, pf_al_cp3, pf_ba_cp3, pf_bl_cp3, pf_aa_cp4, pf_al_cp4, pf_ba_cp4, pf_bl_cp4};
    vf_all7 = {vf_aa_cp3, vf_al_cp3, vf_ba_cp3, vf_bl_cp3, vf_aa_cp4, vf_al_cp4, vf_ba_cp4, vf_bl_cp4};

    p_raw7 = zeros(1, n_comp7);
    for k = 1:n_comp7
        p_raw7(k) = ranksum(pf_all7{k}, vf_all7{k}, 'tail', 'right');
    end
    p_bonf7 = min(p_raw7 * n_comp7, 1);

    % --- Figure ---
    fig7 = figure('Units', 'centimeters', 'Position', [5 5 26 12], 'Color', 'w');
    clf;

    for row = 1:2
        for col = 1:4
            subplot(2, 4, (row-1)*4 + col);

            data7  = plot_pairs7{row, col};
            p_idx7 = p_idx_map7(row, col);

            bh = boxplot(data7, {'PF', 'VF'}, 'Widths', 0.5);
            set(bh, 'LineWidth', 1.5);
            set(gca, 'FontSize', 8);
            hold on;

            h_box    = findobj(gca, 'Tag', 'Box');
            box_cols = [col_vf7; col_pf7];
            for j = 1:min(2, numel(h_box))
                patch(get(h_box(j), 'XData'), get(h_box(j), 'YData'), box_cols(j,:), ...
                    'FaceAlpha', 0.15, 'EdgeColor', box_cols(j,:), 'LineWidth', 1.5);
            end

            rng(0);
            n_pts = size(data7, 1);
            jit   = 0.08 * (rand(n_pts, 1) - 0.5);
            scatter(ones(n_pts,1)   + jit, data7(:,1), 20, col_pf7, 'filled', 'MarkerFaceAlpha', 0.5);
            scatter(2*ones(n_pts,1) + jit, data7(:,2), 20, col_vf7, 'filled', 'MarkerFaceAlpha', 0.5);

            yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            hold off;

            if row == 1, title(col_titles7{col}, 'FontSize', 8, 'FontWeight', 'bold'); end
            if col == 1, ylabel(row_labels7{row}, 'FontSize', 9, 'FontWeight', 'bold'); end
            box off;
        end
    end

    drawnow;
    pos_c1 = get(subplot(2,4,1), 'Position');
    pos_c2 = get(subplot(2,4,2), 'Position');
    pos_c3 = get(subplot(2,4,3), 'Position');
    pos_c4 = get(subplot(2,4,4), 'Position');
    contra_x = pos_c1(1);
    contra_w = (pos_c2(1) + pos_c2(3)) - pos_c1(1);
    ipsi_x   = pos_c3(1);
    ipsi_w   = (pos_c4(1) + pos_c4(3)) - pos_c3(1);
    grp_y    = pos_c1(2) + pos_c1(4) + 0.01;
    grp_h    = 0.05;
    annotation(fig7, 'textbox', [contra_x, grp_y, contra_w, grp_h], ...
        'String', 'Contralateral', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
    annotation(fig7, 'textbox', [ipsi_x, grp_y, ipsi_w, grp_h], ...
        'String', 'Ipsilateral', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

    drawnow;
    set(fig7, 'Color', 'none');
    print(fig7, [figs_path 'figure4d.svg'], '-dsvg');
    set(fig7, 'Color', 'w');

    % --- Statistics ---
    fprintf('\n========================================\n');
    fprintf('Figure 4d Statistics — PF vs VF ERD (alpha & beta)\n');
    fprintf('Active window: 0–750 ms  |  Late window: 2500–3000 ms\n');
    fprintf('Test: Wilcoxon rank-sum, one-tailed (H1: VF ERD more negative than PF, tail=right on PF vs VF)\n');
    fprintf('Multiple comparisons: Bonferroni correction, k = %d\n', n_comp7);
    fprintf('========================================\n');

    for k = 1:n_comp7
        x = pf_all7{k};  y = vf_all7{k};
        N = numel(x) + numel(y);
        [p_r, ~, s7] = ranksum(x, y, 'tail', 'right');
        z7    = s7.zval;
        r_eff = z7 / sqrt(N);
        se7   = 1 / sqrt(N - 3);
        ci_lo = tanh(atanh(r_eff) - 1.96*se7);
        ci_hi = tanh(atanh(r_eff) + 1.96*se7);
        fprintf('\n--- %s ---\n', comp_labels7{k});
        fprintf('  PF: median = %.3f dB, IQR = [%.3f, %.3f]\n', median(x), quantile(x,0.25), quantile(x,0.75));
        fprintf('  VF: median = %.3f dB, IQR = [%.3f, %.3f]\n', median(y), quantile(y,0.25), quantile(y,0.75));
        fprintf('  Z = %.4f\n', z7);
        fprintf('  p (raw)              = %.6f\n', p_r);
        fprintf('  p (Bonferroni k=%d)  = %.6f\n', n_comp7, p_bonf7(k));
        fprintf('  Effect size r        = %.4f,  95%% CI = [%.4f, %.4f]\n', r_eff, ci_lo, ci_hi);
    end

diary off;
