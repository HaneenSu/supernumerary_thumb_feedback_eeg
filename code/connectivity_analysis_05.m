% connectivity_analysis_05.m
%
% Weighted Phase Lag Index (wPLI) connectivity analysis for the sensory feedback study.
%
% Section 1 (PF) & Section 2 (VF): Computes time-resolved beta-band wPLI adjacency
%   matrices from epoched EEG using FieldTrip (mtmconvol + wpli), saves per-group
%   .mat files to ../derivatives/wPLI/.
%
% Small-Worldness sections: Loads wPLI adjacency,
%   computes small-worldness index (sigma) per subject using BCT, and compares
%   PF vs VF with a one-tailed Wilcoxon rank-sum test.
%
% Network Visualization: Plots beta wPLI brain networks for PF and VF,
%   with node size proportional to clustering coefficient.
%
% Expected folder layout (relative to code/):
%   ../derivatives/epoched/PF/  — preprocessed PF epochs
%   ../derivatives/epoched/VF/  — preprocessed VF epochs
%   ../derivatives/wPLI/        — wPLI adjacency .mat files (read & write)
%   ../figs/                    — figure and statistics output
%
% Dependencies: EEGLAB, FieldTrip, Brain Connectivity Toolbox (BCT)
%   Update the FieldTrip path in fieldtrip_candidates and the BCT path below.

clear; clc;


% ---- Toolbox setup ----
fieldtrip_candidates = { ...
    '../toolboxes/fieldtrip' ...  % <-- update to your FieldTrip installation path
    };

fieldtrip_root = '';
for iPath = 1:numel(fieldtrip_candidates)
    has_defaults = exist(fullfile(fieldtrip_candidates{iPath}, 'ft_defaults.m'), 'file');
    has_trackcfg = exist(fullfile(fieldtrip_candidates{iPath}, 'utilities', 'private', ...
        'ft_preamble_trackconfig.m'), 'file');
    if has_defaults && has_trackcfg
        fieldtrip_root = fieldtrip_candidates{iPath};
        break
    end
end

if isempty(fieldtrip_root)
    error(['A complete FieldTrip install was not found. ', ...
        'Update fieldtrip_candidates at the top of %s.'], mfilename('fullpath'));
end

% Remove conflicting FieldTrip copies that may already be on the MATLAB path.
all_paths = strsplit(path, pathsep);
for iPath = 1:numel(all_paths)
    this_path = all_paths{iPath};
    if ~isempty(strfind(lower(this_path), 'fieldtrip'))
        rmpath(this_path);
    end
end

% Prepend the selected FieldTrip directories explicitly so helper scripts in
% utilities/private come from the same install as ft_freqanalysis.
ft_subdirs = { ...
    fieldtrip_root, ...
    fullfile(fieldtrip_root, 'utilities'), ...
    fullfile(fieldtrip_root, 'specest'), ...
    fullfile(fieldtrip_root, 'connectivity'), ...
    fullfile(fieldtrip_root, 'preproc'), ...
    fullfile(fieldtrip_root, 'fileio'), ...
    fullfile(fieldtrip_root, 'forward'), ...
    fullfile(fieldtrip_root, 'inverse'), ...
    fullfile(fieldtrip_root, 'plotting'), ...
    fullfile(fieldtrip_root, 'trialfun') ...
    };

for iPath = 1:numel(ft_subdirs)
    if exist(ft_subdirs{iPath}, 'dir')
        addpath(ft_subdirs{iPath}, '-begin');
    end
end

ft_defaults;

assert(~isempty(which('ft_freqanalysis')), ...
    'ft_freqanalysis is not on the MATLAB path after FieldTrip setup.');
assert(~isempty(which('ft_specest_mtmconvol')), ...
    ['ft_specest_mtmconvol is not on the MATLAB path. ', ...
     'Remove conflicting FieldTrip copies and rerun ft_defaults.']);
assert(~isempty(which('ft_preamble')), ...
    'ft_preamble is not on the MATLAB path after FieldTrip setup.');
assert(exist(fullfile(fieldtrip_root, 'utilities', 'private', ...
    'ft_preamble_trackconfig.m'), 'file') == 2, ...
    'The selected FieldTrip install is incomplete: ft_preamble_trackconfig.m is missing.');
assert(~isempty(strfind(which('ft_preamble'), fieldtrip_root)), ...
    'MATLAB is still resolving ft_preamble from a different FieldTrip copy.');


figs_path = '../figs/';

%% ---- Pressure FB wPLI ----
load EEG_chlocs.mat
path_to_epoched = '../derivatives/epoched/PF/'; 

subjects = 1:21;
conds    = {'PF'};

% ---- Connectivity choices ----
conn_method = 'wpli';   % robust vs zero-lag volume conduction (recommended)

band = [13 30];                          % example band
foi  = 13:1:30   ;                  % frequencies of interest (adjust)

% ---- Time settings for time-resolved connectivity ----
toi = -3:0.1:5;                % 50 ms step (adjust)
twin = 1;                     % 500 ms window length (adjust)

A_all = []; labels_ref = [];

for si = 1:numel(subjects)
    sub = subjects(si);

    for ci = 1:numel(conds)
        cond = conds{ci};

        % ---- Load EEGLAB ----
        eeg_file   = ['SF' int2str(si) '_' cond '.set'];
        EEG = pop_loadset('filename', eeg_file, 'filepath', path_to_epoched);

        % Optional: resample to reduce compute (uncomment if desired)
        % EEG = pop_resample(EEG, 500);

        % ---- Convert to FieldTrip epoched format ----
        % 'preprocessing' keeps trials; 'none' avoids extra preprocessing here
        data_ft = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

        % ---- Spectral estimation over time (Fourier spectra per trial, freq, time) ----
        cfg = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'fourier';
        cfg.foi        = foi;
        cfg.toi        = toi;
        cfg.t_ftimwin  = twin * ones(size(foi));   % fixed window per frequency

        cfg.taper      = 'hanning';
        cfg.keeptrials = 'yes';
        freq = ft_freqanalysis(cfg, data_ft);

        % ---- Connectivity ----
        cfg = [];
        cfg.method = conn_method;

        conn = ft_connectivityanalysis(cfg, freq);

        % conn.<metric>spctrm shape: chan x chan x freq x time
        metric_field = [conn_method 'spctrm'];
        % X = conn.(metric_field); % Force absolute wPLI (0 to 1)
        X = abs(conn.(metric_field)); % Force absolute wPLI (0 to 1)
        X(isnan(X)) = 0;
        % band-average over frequency (dim 3)
        fidx  = conn.freq >= band(1) & conn.freq <= band(2);
        A = squeeze(mean(X(:,:,fidx,:), 3));       % 60x60x141
        % make sure diagonal is zero at each time
        nCh = size(A,1);
        nT  = size(A,3);
        for tt = 1:nT
            At = A(:,:,tt);
            At(1:nCh+1:end) = 0;
            A(:,:,tt) = At;
        end
        
        % store
        if isempty(A_all)
            A_all = zeros(numel(subjects), numel(conds), nCh, nCh, nT, 'single');
        end
        A_all(si, ci, :, :, :) = single(A);
        if isempty(labels_ref), labels_ref = conn.label; end
    end
end



save('../derivatives/wPLI/adj_wpli_beta_PF_ls.mat', ...
     'A_all','labels_ref','toi','foi','band','conn_method','twin','-v7.3');

%% ---- Vibration FB wPLI ----
load EEG_chlocs.mat
path_to_epoched = '../derivatives/epoched/VF/'; 

subjects = 1:21;
conds    = {'VF'};

% ---- Connectivity choices ----
conn_method = 'wpli';   % robust vs zero-lag volume conduction (recommended)

band = [13 30];                          % example band
foi  = 13:1:30   ;                  % frequencies of interest (adjust)

% ---- Time settings for time-resolved connectivity ----
toi = -3:0.1:5;                % 50 ms step (adjust)
twin = 1;                     % 500 ms window length (adjust)

A_all = []; labels_ref = [];

for si = 1:numel(subjects)
    sub = subjects(si);

    for ci = 1:numel(conds)
        cond = conds{ci};

        % ---- Load EEGLAB ----
        eeg_file   = ['SF' int2str(si) '_' cond '.set'];
        EEG = pop_loadset('filename', eeg_file, 'filepath', path_to_epoched);

        % Optional: resample to reduce compute (uncomment if desired)
        % EEG = pop_resample(EEG, 500);

        % ---- Convert to FieldTrip epoched format ----
        % 'preprocessing' keeps trials; 'none' avoids extra preprocessing here
        data_ft = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

        % ---- Spectral estimation over time (Fourier spectra per trial, freq, time) ----
        cfg = [];
        cfg.method     = 'mtmconvol';
        cfg.output     = 'fourier';
        cfg.foi        = foi;
        cfg.toi        = toi;
        cfg.t_ftimwin  = twin * ones(size(foi));   % fixed window per frequency

        cfg.taper      = 'hanning';
        cfg.keeptrials = 'yes';
        freq = ft_freqanalysis(cfg, data_ft);

        % ---- Connectivity ----
        cfg = [];
        cfg.method = conn_method;

        conn = ft_connectivityanalysis(cfg, freq);

        % conn.<metric>spctrm shape: chan x chan x freq x time
        metric_field = [conn_method 'spctrm'];
        % X = conn.(metric_field); % Force absolute wPLI (0 to 1)

        X = abs(conn.(metric_field)); % Force absolute wPLI (0 to 1)
        X(isnan(X)) = 0;
        % band-average over frequency (dim 3)
        fidx  = conn.freq >= band(1) & conn.freq <= band(2);
        A = squeeze(mean(X(:,:,fidx,:), 3));       % 60x60x141
        % make sure diagonal is zero at each time
        nCh = size(A,1);
        nT  = size(A,3);
        for tt = 1:nT
            At = A(:,:,tt);
            At(1:nCh+1:end) = 0;
            A(:,:,tt) = At;
        end
        
        % store
        if isempty(A_all)
            A_all = zeros(numel(subjects), numel(conds), nCh, nCh, nT, 'single');
        end
        A_all(si, ci, :, :, :) = single(A);
        if isempty(labels_ref), labels_ref = conn.label; end
    end
end



save('../derivatives/wPLI/adj_wpli_beta_VF_ls.mat', ...
     'A_all','labels_ref','toi','foi','band','conn_method','twin','-v7.3');


%% ========== Small-Worldness: PF vs VF ==========
% Uses BCT (Brain Connectivity Toolbox).
% This section is self-contained.
clc; close all

% ---- Add BCT to path ----
addpath('../toolboxes/BCT');  % <-- update to your BCT installation path

% ---- Load abs wPLI adjacency (unsigned) ----
data_PF_cpl = load('../derivatives/wPLI/adj_wpli_beta_PF_ls.mat');
data_VF_cpl = load('../derivatives/wPLI/adj_wpli_beta_VF_ls.mat');

% ---- Channel labels ----
if ~isempty(data_PF_cpl.labels_ref)
    labels_cpl = data_PF_cpl.labels_ref;
else
    load EEG_chlocs.mat
    labels_cpl = strtrim({EEG_chlocs.labels})';
end

% ---- Time setup ----
times_cpl   = data_PF_cpl.toi;
lateWin_cpl = [2 3];

[~,tL1c] = min(abs(times_cpl - lateWin_cpl(1)));
[~,tL2c] = min(abs(times_cpl - lateWin_cpl(2)));

% ---- Merged condition, late window only ----
A_PF_cS = squeeze(data_PF_cpl.A_all(:,1,:,:,:));  % single merged condition
A_VF_cS = squeeze(data_VF_cpl.A_all(:,1,:,:,:));
nSub_cpl = size(A_PF_cS, 1);

A_PF_cS_lat = squeeze(mean(A_PF_cS(:,:,:,tL1c:tL2c), 4));  % [nSub x 60 x 60]
A_VF_cS_lat = squeeze(mean(A_VF_cS(:,:,:,tL1c:tL2c), 4));

% ---- Compute small-worldness, clustering coefficient, and path length ----
sw_pf  = nan(nSub_cpl, 1);
sw_vf  = nan(nSub_cpl, 1);
cc_pf  = nan(nSub_cpl, 1);
cc_vf  = nan(nSub_cpl, 1);
cpl_pf = nan(nSub_cpl, 1);
cpl_vf = nan(nSub_cpl, 1);
nRand = 1000;

for s = 1:nSub_cpl
    for grp = 1:2
        if grp == 1
            W = double(squeeze(A_PF_cS_lat(s,:,:)));
        else
            W = double(squeeze(A_VF_cS_lat(s,:,:)));
        end
        W(isnan(W)) = 0;
        W = abs(W);
        W(1:size(W,1)+1:end) = 0;
        W = (W + W') / 2;

        % Threshold: keep top 20% edges
        thr = prctile(W(W>0), 80);
        W(W < thr) = 0;

        W_nrm = weight_conversion(W, 'normalize');
        C = mean(clustering_coef_wu(W_nrm));
        L_mat = zeros(size(W));
        L_mat(W > 0) = 1 ./ W(W > 0);
        Lpath = charpath(distance_wei(L_mat), 0, 0);

        C_rand_all = zeros(nRand, 1);
        L_rand_all = zeros(nRand, 1);
        for ri = 1:nRand
            R     = randmio_und(W, 5);
            R_nrm = weight_conversion(R, 'normalize');
            C_rand_all(ri) = mean(clustering_coef_wu(R_nrm));
            Lr = zeros(size(R));
            Lr(R > 0) = 1 ./ R(R > 0);
            L_rand_all(ri) = charpath(distance_wei(Lr), 0, 0);
        end
        C_rand = mean(C_rand_all);
        L_rand = mean(L_rand_all);
        sigma = NaN;
        if C_rand > 0 && L_rand > 0
            sigma = (C / C_rand) / (Lpath / L_rand);
        end
        if grp == 1
            sw_pf(s)  = sigma;
            cc_pf(s)  = C;
            cpl_pf(s) = Lpath;
        else
            sw_vf(s)  = sigma;
            cc_vf(s)  = C;
            cpl_vf(s) = Lpath;
        end
    end
    fprintf('Subject %d/%d done\n', s, nSub_cpl);
end

% ---- Stats: one-tailed (H1: PF sigma > VF sigma) ----
[p_sw, ~, st_sw] = ranksum(sw_pf, sw_vf, 'method', 'approximate', 'tail', 'right');
N_sw  = nSub_cpl * 2;
r_sw  = st_sw.zval / sqrt(N_sw);
se_sw = 1 / sqrt(N_sw - 3);
ci_lo = tanh(atanh(r_sw) - 1.96*se_sw);
ci_hi = tanh(atanh(r_sw) + 1.96*se_sw);

fid_stats = fopen([figs_path 'statistics.txt'], 'a');

fprintf(fid_stats, '\n========== SMALL-WORLDNESS INDEX: PF vs VF ( Late %.1f–%.1fs) ==========\n', lateWin_cpl);
fprintf(fid_stats, 'Test: Wilcoxon rank-sum, one-tailed (H1: PF sigma > VF sigma)\n');
fprintf(fid_stats, 'PF: median = %.3f, IQR = [%.3f, %.3f]\n', median(sw_pf,'omitnan'), quantile(sw_pf,0.25), quantile(sw_pf,0.75));
fprintf(fid_stats, 'VF: median = %.3f, IQR = [%.3f, %.3f]\n', median(sw_vf,'omitnan'), quantile(sw_vf,0.25), quantile(sw_vf,0.75));
fprintf(fid_stats, 'Z = %.4f,  p = %.6f\n', st_sw.zval, p_sw);
fprintf(fid_stats, 'Effect size r = %.4f,  95%% CI = [%.4f, %.4f]\n', r_sw, ci_lo, ci_hi);

fprintf('\n========== SMALL-WORLDNESS INDEX: PF vs VF ( Late %.1f–%.1fs) ==========\n', lateWin_cpl);
fprintf('Test: Wilcoxon rank-sum, one-tailed (H1: PF sigma > VF sigma)\n');
fprintf('PF: median = %.3f, IQR = [%.3f, %.3f]\n', ...
    median(sw_pf,'omitnan'), quantile(sw_pf,0.25), quantile(sw_pf,0.75));
fprintf('VF: median = %.3f, IQR = [%.3f, %.3f]\n', ...
    median(sw_vf,'omitnan'), quantile(sw_vf,0.25), quantile(sw_vf,0.75));
fprintf('Z = %.4f,  p = %.6f\n', st_sw.zval, p_sw);
fprintf('Effect size r = %.4f,  95%% CI = [%.4f, %.4f]\n', r_sw, ci_lo, ci_hi);

% ---- Plot ----
col_pf_sw = [0.92, 0.46, 0.32];   % coral
col_vf_sw = [0.28, 0.60, 0.42];   % forest green

fig_sw = figure('Units', 'centimeters', 'Position', [5 5 8 10], 'Color', 'w');

data_sw = [sw_pf, sw_vf];   % [nSub x 2]
bh = boxplot(data_sw, {'PF', 'VF'}, 'Widths', 0.5);
set(bh, 'LineWidth', 1.5);
set(gca, 'FontSize', 10);

hold on;

h_box    = findobj(gca, 'Tag', 'Box');
box_cols = [col_vf_sw; col_pf_sw];
for j = 1:min(2, numel(h_box))
    patch(get(h_box(j), 'XData'), get(h_box(j), 'YData'), box_cols(j,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', box_cols(j,:), 'LineWidth', 1.5);
end

rng(0);
n_pts = nSub_cpl;
jit   = 0.08 * (rand(n_pts, 1) - 0.5);
scatter(ones(n_pts,1)   + jit, sw_pf, 20, col_pf_sw, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(2*ones(n_pts,1) + jit, sw_vf, 20, col_vf_sw, 'filled', 'MarkerFaceAlpha', 0.5);

yline(1, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

if p_sw < 0.05
    y_sig = max(data_sw(:)) + 0.15 * (max(data_sw(:)) - min(data_sw(:)));
    text(1.5, y_sig, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end

hold off;

ylabel('\sigma (Small-Worldness)', 'FontSize', 10);
title(sprintf('Rebound window'), ...
    'FontSize', 10, 'FontWeight', 'bold');
box off;

drawnow;
set(fig_sw, 'Color', 'none');
print(fig_sw, [figs_path 'figure_sw_boxplot.svg'], '-dsvg');
set(fig_sw, 'Color', 'w');

%% ========== Network Visualization: Beta wPLI, PF vs VF (Short, Late) ==========
% Two side-by-side brain network plots.
% Node size = clustering coefficient; edge thickness/opacity = wPLI strength.
% Assumes adj_wpli_beta_PF.mat / adj_wpli_beta_VF.mat (already loaded as
% data_PF_cpl / data_VF_cpl from the Small-Worldness section above).
clc; close all;

addpath('../toolboxes/BCT');  % <-- update to your BCT installation path

if ~exist('data_PF_cpl', 'var')
    data_PF_cpl = load('../derivatives/wPLI/adj_wpli_beta_PF_ls.mat');
    data_VF_cpl = load('../derivatives/wPLI/adj_wpli_beta_VF_ls.mat');
end

% ---- Late window indices (merged condition, index 1) ----
times_net   = data_PF_cpl.toi;
lateWin_net = [2 3];
[~, tL1n] = min(abs(times_net - lateWin_net(1)));
[~, tL2n] = min(abs(times_net - lateWin_net(2)));

% Average across subjects (dim 1) and late time points (dim 5) → [60 x 60]
% A_all: [nSub x 1 x 60 x 60 x nTime]
W_pf_raw = squeeze(mean(mean(data_PF_cpl.A_all(:,1,:,:,tL1n:tL2n), 5), 1));
W_vf_raw = squeeze(mean(mean(data_VF_cpl.A_all(:,1,:,:,tL1n:tL2n), 5), 1));

nCh_net  = size(W_pf_raw, 1);   % 60

% Absolute value, symmetrise, zero diagonal
for W_tmp = {W_pf_raw, W_vf_raw}; end  % just to iterate cleanly below
W_pf_raw = abs(W_pf_raw); W_pf_raw = (W_pf_raw + W_pf_raw')/2; W_pf_raw(1:nCh_net+1:end) = 0;
W_vf_raw = abs(W_vf_raw); W_vf_raw = (W_vf_raw + W_vf_raw')/2; W_vf_raw(1:nCh_net+1:end) = 0;

% ---- Threshold: top 20% of connections per group ----
thr_pf = prctile(W_pf_raw(W_pf_raw > 0), 80);
thr_vf = prctile(W_vf_raw(W_vf_raw > 0), 80);
W_pf = W_pf_raw .* (W_pf_raw >= thr_pf);
W_vf = W_vf_raw .* (W_vf_raw >= thr_vf);

% ---- Clustering coefficient per node (BCT) ----
cc_pf_net = clustering_coef_wu(weight_conversion(W_pf, 'normalize'));
cc_vf_net = clustering_coef_wu(weight_conversion(W_vf, 'normalize'));

% ---- Shared node size scale across both groups ----
sz_min = 25; sz_max = 220;
cc_all_net = [cc_pf_net; cc_vf_net];
cc_lo = min(cc_all_net);  cc_hi = max(cc_all_net);
norm_cc = @(v) sz_min + (v - cc_lo) / max(cc_hi - cc_lo, eps) * (sz_max - sz_min);
sz_pf_net = norm_cc(cc_pf_net);
sz_vf_net = norm_cc(cc_vf_net);

% ---- Shared edge strength scale across both groups ----
e_vals = [W_pf(W_pf > 0); W_vf(W_vf > 0)];
e_min  = min(e_vals);  e_max = max(e_vals);
lw_min_net = 0.5;  lw_max_net = 5;
al_min_net = 0.08; al_max_net = 0.85;
cmap_net   = parula(256);

norm_e = @(v) (v - e_min) / max(e_max - e_min, eps);  % → [0,1]

% ---- Electrode positions (shared) ----
load('../code/systems.mat');
t_net = systems(18).layout;
t_net([5,10,21,27],:) = [];
a_head = 0.45; b_head = 0.5;
x0h = t_net(21,1); y0h = t_net(21,2);
tt_h = -pi:0.01:pi;

titles_net  = {'PF', 'VF'};
fnames_net  = {'figure_network_PF.svg', 'figure_network_VT.svg'};
W_grps      = {W_pf,      W_vf};
sz_grps     = {sz_pf_net, sz_vf_net};

for g = 1:2
    fig_g = figure('Units', 'centimeters', 'Position', [2 2 16 15], 'Color', 'w');
    ax_n  = axes(fig_g);

    % Gray scalp outline
    plot(ax_n, x0h + a_head*cos(tt_h), y0h + b_head*sin(tt_h), ...
        'Color', [0.80 0.80 0.80], 'LineWidth', 2);
    hold(ax_n, 'on');

    % Draw edges (below nodes)
    W_g  = W_grps{g};
    sz_g = sz_grps{g};
    [r_idx, c_idx] = find(triu(W_g > 0, 1));
    for ei = 1:numel(r_idx)
        i = r_idx(ei);  j = c_idx(ei);
        sn   = norm_e(W_g(i,j));
        lw_e = lw_min_net + sn * (lw_max_net - lw_min_net);
        al_e = al_min_net + sn * (al_max_net - al_min_net);
        ci_e = max(1, round(1 + sn * 255));
        edge_col = [cmap_net(ci_e,:), al_e];
        drawline_thick(ax_n, i, j, edge_col, lw_e);
    end

    % Draw nodes on top
    scatter(ax_n, t_net(:,1), t_net(:,2), sz_g, [0.15 0.15 0.15], 'filled');

    title(ax_n, titles_net{g}, 'FontSize', 12, 'FontWeight', 'bold');
    set(ax_n, 'YTickLabel',[], 'XTickLabel',[], 'YTick',[], 'XTick',[]);
    axis(ax_n, 'equal'); axis(ax_n, 'off');

    % Colorbar
    colormap(ax_n, parula);
    clim(ax_n, [e_min e_max]);
    cb = colorbar(ax_n, 'eastoutside');
    ylabel(cb, 'wPLI strength', 'FontSize', 10);
    set(cb, 'FontSize', 9);

    drawnow;
    set(fig_g, 'Color', 'w');
    print(fig_g, [figs_path fnames_net{g}], '-dpdf', '-bestfit');
end


