% SF_pre_bhv.m
%
% Behavioral analysis for the sensory feedback study.
% Computes explicit and implicit identification accuracy across for two groups:
%   - PF: pressure feedback
%   - VF: vibrotactile feedback
%
%
% Expected folder layout (relative to code/):
%   ../rawdata/SFP/   — raw PF EEG data (.vhdr files)
%   ../rawdata/SFV/   — raw VF EEG data (.vhdr files)
%   ../figs/          — output figures and statistics.txt
%
% Dependencies: EEGLAB

clearvars;
close all;
clc;


cfg = get_analysis_config();
if ~exist(cfg.figs_path, 'dir'), mkdir(cfg.figs_path); end

stats_file = fullfile(cfg.figs_path, 'statistics.txt');
if exist(stats_file, 'file'), delete(stats_file); end

vf = collect_group_metrics(cfg.vf_path, 'SFV', 'VF', cfg);
pf = collect_group_metrics(cfg.pf_path, 'SFP', 'PF', cfg);

print_descriptive_summary(vf, pf);

%% Figure 2: overall identification accuracy
close all;
plot_figure2(vf, pf, cfg.figs_path);
diary(stats_file); run_figure2_statistics(vf, pf); diary off;

%% Supplementary Figure 1: consistency matrices
close all; plot_supp1(vf, pf, cfg);

%% Supplementary Figure 2: implicit articulation length vs direction
close all; plot_supp2(vf, pf, cfg.figs_path);
diary(stats_file); run_supp2_statistics(vf, pf); diary off;

function cfg = get_analysis_config()
cfg = struct();
cfg.vf_path   = '../rawdata/SFV/';
cfg.pf_path   = '../rawdata/SFP/';
cfg.figs_path = '../figs';
cfg.n_subject = 21;
cfg.trials_per_session = 36;
cfg.sessions_per_condition = 2;
cfg.motion_codes = {'4', '5', '6', '7', '8', '9'};
cfg.motion_labels = {'P1 to P3', 'P3 to P1', 'P1 to P2', 'P2 to P3', 'P3 to P2', 'P2 to P1'};
cfg.reported_motion_codes = {'6', '9'};
end

function metrics = collect_group_metrics(path_to_rawdata, file_prefix, modality, cfg)
n_subject = cfg.n_subject;
n_motion_types = numel(cfg.motion_codes);
n_transitions = (cfg.trials_per_session - 1) * cfg.sessions_per_condition;

metrics = struct();
metrics.modality = modality;
metrics.explicit_accuracy = zeros(n_subject, 1);
metrics.implicit_accuracy = zeros(n_subject, 1);
metrics.explicit_continuity = false(n_subject, n_transitions);
metrics.implicit_continuity = false(n_subject, n_transitions);
metrics.implicit_length_accuracy = zeros(n_subject, 1);
metrics.implicit_direction_accuracy = zeros(n_subject, 1);

for sub = 1:n_subject
    subEEG = load_clean_subject_eeg(path_to_rawdata, file_prefix, sub, modality);
    trials = extract_trial_records(subEEG);

    explicit_trials = trials([trials.is_explicit]);
    implicit_trials = trials(~[trials.is_explicit]);

    metrics.explicit_accuracy(sub) = mean([explicit_trials.is_correct]);
    metrics.implicit_accuracy(sub) = mean([implicit_trials.is_correct]);
    metrics.explicit_continuity(sub, :) = compute_continuity_mask(explicit_trials, cfg.trials_per_session);
    metrics.implicit_continuity(sub, :) = compute_continuity_mask(implicit_trials, cfg.trials_per_session);

    [length_acc, direction_acc] = compute_length_direction_accuracy(implicit_trials);
    metrics.implicit_length_accuracy(sub) = length_acc;
    metrics.implicit_direction_accuracy(sub) = direction_acc;
end
end

function plot_figure2(vf, pf, figs_path)
col_imp = [0.60, 0.60, 0.60];  % gray         — implicit
col_pf  = [0.92, 0.46, 0.32];  % coral        — pressure
col_vf  = [0.28, 0.60, 0.42];  % forest green — vibrotactile
fig_size = [2, 2, 8, 9];       % [left bottom width height] in centimeters

% Figure 2A: PF explicit vs PF implicit (within-group)
figure('Name', 'Figure 2A - PF vs Implicit', 'Color', 'w', 'Units', 'centimeters', 'Position', fig_size);
styled_boxplot(100 * [pf.implicit_accuracy, pf.explicit_accuracy], ...
    {'Implicit', 'Pressure'}, 'Accuracy [%]', [], false, [col_imp; col_pf]);
ylim([0 110]); title('PF group');
drawnow;
exportgraphics(gcf, fullfile(figs_path, 'fig2a_pf_vs_implicit.pdf'), 'ContentType', 'vector');

% Figure 2B: VF explicit vs VF implicit (within-group)
figure('Name', 'Figure 2B - VF vs Implicit', 'Color', 'w', 'Units', 'centimeters', 'Position', fig_size);
styled_boxplot(100 * [vf.implicit_accuracy, vf.explicit_accuracy], ...
    {'Implicit', 'Vibrotactile'}, 'Accuracy [%]', [], false, [col_imp; col_vf]);
ylim([0 110]); title('VF group');
drawnow;
exportgraphics(gcf, fullfile(figs_path, 'fig2b_vf_vs_implicit.pdf'), 'ContentType', 'vector');

% Figure 2C: VF vs PF overall (between-group)
figure('Name', 'Figure 2C - VF vs PF Overall', 'Color', 'w', 'Units', 'centimeters', 'Position', fig_size);
styled_boxplot(100 * [vf.explicit_accuracy, pf.explicit_accuracy], ...
    {'Vibrotactile', 'Pressure'}, 'Accuracy [%]', [], false, [col_vf; col_pf]);
ylim([0 110]); title('Overall');
drawnow;
exportgraphics(gcf, fullfile(figs_path, 'fig2c_vf_vs_pf_overall.pdf'), 'ContentType', 'vector');
end

function run_figure2_statistics(vf, pf)
fprintf('\n========================================\n');
fprintf('Figure 2 Statistics\n');
fprintf('========================================\n');

fprintf('\n--- 2A: PF explicit vs PF implicit (within-group) ---\n');
report_wilcoxon_signed(pf.implicit_accuracy, pf.explicit_accuracy, ...
    'PF implicit', 'PF explicit', 'two-tailed', 'none');

fprintf('\n--- 2B: VF explicit vs VF implicit (within-group) ---\n');
report_wilcoxon_signed(vf.implicit_accuracy, vf.explicit_accuracy, ...
    'VF implicit', 'VF explicit', 'two-tailed', 'none');

p_overall = ranksum(vf.explicit_accuracy, pf.explicit_accuracy);

fprintf('\n--- 2C: VF vs PF overall accuracy (between-group) ---\n');
report_wilcoxon_ranksum(vf.explicit_accuracy, pf.explicit_accuracy, ...
    'VF overall', 'PF overall', 'two-tailed', 'none', p_overall);
end

function report_wilcoxon_signed(x, y, label_x, label_y, tails, correction)
n = numel(x);
[p_raw, ~, stats] = signrank(x, y);
z = stats.zval;
r = z / sqrt(n);
[ci_lo, ci_hi] = effect_size_ci(r, n);

fprintf('  Test:         Wilcoxon signed-rank, %s\n', tails);
fprintf('  Correction:   %s\n', correction);
fprintf('  n (pairs):    %d\n', n);
fprintf('  %-22s median = %.1f%%, IQR = [%.1f, %.1f]%%\n', ...
    label_x, 100*median(x), 100*quantile(x,0.25), 100*quantile(x,0.75));
fprintf('  %-22s median = %.1f%%, IQR = [%.1f, %.1f]%%\n', ...
    label_y, 100*median(y), 100*quantile(y,0.25), 100*quantile(y,0.75));
fprintf('  Z = %.4f\n', z);
fprintf('  p (raw) = %s\n', format_pvalue(p_raw));
fprintf('  Effect size r = Z/sqrt(N) = %.4f\n', r);
fprintf('  95%% CI on r = [%.4f, %.4f]\n', ci_lo, ci_hi);
end

function report_wilcoxon_ranksum(x, y, label_x, label_y, tails, correction, p_corr)
n1 = numel(x);
n2 = numel(y);
N  = n1 + n2;
[p_raw, ~, stats] = ranksum(x, y);
z = stats.zval;
r = z / sqrt(N);
[ci_lo, ci_hi] = effect_size_ci(r, N);

fprintf('  Test:         Wilcoxon rank-sum, %s\n', tails);
fprintf('  Correction:   %s\n', correction);
fprintf('  n per group:  %d, %d\n', n1, n2);
fprintf('  %-22s median = %.1f%%, IQR = [%.1f, %.1f]%%\n', ...
    label_x, 100*median(x), 100*quantile(x,0.25), 100*quantile(x,0.75));
fprintf('  %-22s median = %.1f%%, IQR = [%.1f, %.1f]%%\n', ...
    label_y, 100*median(y), 100*quantile(y,0.25), 100*quantile(y,0.75));
fprintf('  Z = %.4f\n', z);
fprintf('  p (raw)       = %s\n', format_pvalue(p_raw));
fprintf('  p (corrected) = %s\n', format_pvalue(p_corr));
fprintf('  Effect size r = Z/sqrt(N) = %.4f\n', r);
fprintf('  95%% CI on r = [%.4f, %.4f]\n', ci_lo, ci_hi);
end

function [ci_lo, ci_hi] = effect_size_ci(r, N)
% 95% CI on effect size r via Fisher z-transformation
z_r   = atanh(r);
se    = 1 / sqrt(N - 3);
ci_lo = tanh(z_r - 1.96 * se);
ci_hi = tanh(z_r + 1.96 * se);
end


function plot_supp1(vf, pf, cfg)
[vf_implicit_s1, vf_implicit_s2] = split_continuity_sessions(vf.implicit_continuity, cfg.trials_per_session);
[pf_implicit_s1, pf_implicit_s2] = split_continuity_sessions(pf.implicit_continuity, cfg.trials_per_session);
[vf_explicit_s1, vf_explicit_s2] = split_continuity_sessions(vf.explicit_continuity, cfg.trials_per_session);
[pf_explicit_s1, pf_explicit_s2] = split_continuity_sessions(pf.explicit_continuity, cfg.trials_per_session);

figure('Name', 'Supp Figure 1 - Implicit', 'Color', 'w');

subplot(2, 2, 1);
plot_consistency_matrix(vf_implicit_s1, {'Vibrotactile'; 'Session 1'}, false, true);

subplot(2, 2, 2);
plot_consistency_matrix(vf_implicit_s2, {'Vibrotactile'; 'Session 2'}, false, false);

subplot(2, 2, 3);
plot_consistency_matrix(pf_implicit_s1, {'Pressure'; 'Session 1'}, true, true);

subplot(2, 2, 4);
plot_consistency_matrix(pf_implicit_s2, {'Pressure'; 'Session 2'}, true, false);

set(gcf, 'Position', [100 100 800 600]);
exportgraphics(gcf, fullfile(cfg.figs_path, 'supp1a_consistency_implicit.pdf'), 'ContentType', 'vector');

figure('Name', 'Supp Figure 1 - Explicit', 'Color', 'w');

subplot(2, 2, 1);
plot_consistency_matrix(vf_explicit_s1, {'Vibrotactile'; 'Session 1'}, false, true);

subplot(2, 2, 2);
plot_consistency_matrix(vf_explicit_s2, {'Vibrotactile'; 'Session 2'}, false, false);

subplot(2, 2, 3);
plot_consistency_matrix(pf_explicit_s1, {'Pressure'; 'Session 1'}, true, true);

subplot(2, 2, 4);
plot_consistency_matrix(pf_explicit_s2, {'Pressure'; 'Session 2'}, true, false);

set(gcf, 'Position', [100 100 800 600]);
exportgraphics(gcf, fullfile(cfg.figs_path, 'supp1b_consistency_explicit.pdf'), 'ContentType', 'vector');
end

function plot_supp2(vf, pf, figs_path)
length_accuracy = [vf.implicit_length_accuracy; pf.implicit_length_accuracy];
direction_accuracy = [vf.implicit_direction_accuracy; pf.implicit_direction_accuracy];

figure('Name', 'Supp Figure 2');
styled_boxplot(100 * [length_accuracy, direction_accuracy], ...
    {'Articulation length', 'Articulation direction'}, 'Accuracy [%]', [], false);
set(gcf, 'Position', [100 100 550 420]);
ylim([0 110]);
exportgraphics(gcf, fullfile(figs_path, 'supp2_length_vs_direction.pdf'), 'ContentType', 'vector');
end

function run_supp2_statistics(vf, pf)
length_accuracy    = [vf.implicit_length_accuracy;    pf.implicit_length_accuracy];
direction_accuracy = [vf.implicit_direction_accuracy; pf.implicit_direction_accuracy];

fprintf('\n========================================\n');
fprintf('Supplementary Figure 2 Statistics\n');
fprintf('========================================\n');

fprintf('\n--- Supp2: Implicit length vs direction accuracy (within-subject, pooled groups) ---\n');
report_wilcoxon_signed(length_accuracy, direction_accuracy, ...
    'Implicit length', 'Implicit direction', 'two-tailed', 'none');
end

function print_descriptive_summary(vf, pf)
all_implicit = [vf.implicit_accuracy; pf.implicit_accuracy];

fprintf('Descriptive summary\n');
fprintf('Vibrotactile explicit: %.2f +/- %.2f %%\n', mean(100 * vf.explicit_accuracy), std(100 * vf.explicit_accuracy));
fprintf('Pressure explicit: %.2f +/- %.2f %%\n', mean(100 * pf.explicit_accuracy), std(100 * pf.explicit_accuracy));
fprintf('Implicit pooled: %.2f +/- %.2f %%\n', mean(100 * all_implicit), std(100 * all_implicit));
end

function print_swtest_if_available(label, values)
if exist('swtest', 'file') ~= 2
    return;
end

[h_value, p_value] = swtest(values);
fprintf('%s SW test: H=%d, p=%s\n', label, h_value, format_pvalue(p_value));
end

function formatted_p = format_pvalue(p_value)
formatted_p = sprintf('%.8f', p_value);
end


function plot_consistency_matrix(matrix_data, plot_title, is_bottom, is_left)
imagesc(double(matrix_data));
colormap([0 0 0; 1 1 1]);
caxis([0 1]);

set(gca, 'YDir', 'reverse');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', [], 'YTick', []);
title(plot_title, 'FontSize', 13, 'FontWeight', 'bold');

if nargin > 2 && is_bottom
    xlabel('Consecutive trial pair', 'FontSize', 14, 'FontWeight', 'bold');
end
if nargin > 3 && is_left
    ylabel('Participant', 'FontSize', 14, 'FontWeight', 'bold');
end
hold on;

n_rows = size(matrix_data, 1);
n_cols = size(matrix_data, 2);
grid_color = [0.7 0.7 0.7];

for x = 0.5:(n_cols + 0.5)
    plot([x x], [0.5 n_rows + 0.5], 'Color', grid_color, 'LineWidth', 0.5);
end

for y = 0.5:(n_rows + 0.5)
    plot([0.5 n_cols + 0.5], [y y], 'Color', grid_color, 'LineWidth', 0.5);
end

xlim([0.5 n_cols + 0.5]);
ylim([0.5 n_rows + 0.5]);
hold off;
end

function [session1, session2] = split_continuity_sessions(matrix_data, trials_per_session)
session_transitions = trials_per_session - 1;
session1 = matrix_data(:, 1:session_transitions);
session2 = matrix_data(:, session_transitions + 1:2 * session_transitions);
end

function subEEG = load_clean_subject_eeg(path_to_rawdata, file_prefix, sub, modality)
sub_str = sprintf('%04d', sub - 1);
set_file = fullfile(path_to_rawdata, [file_prefix sub_str '.vhdr']);
disp(['reading... ' set_file]);

EEG = pop_fileio(set_file);
subEEG = struct(EEG);
subEEG = remove_empty_markers(subEEG);
subEEG = apply_subject_specific_fixes(subEEG, sub, modality);
end

function subEEG = remove_empty_markers(subEEG)
del = 0;
for j = 1:size(subEEG.event, 2)
    if strcmp(subEEG.event(j - del).type, 'empty')
        subEEG.event(j - del) = [];
        del = del + 1;
    end
end
end

function subEEG = apply_subject_specific_fixes(subEEG, sub, modality)
if sub == 1
    subEEG = relabel_sub1_no_feedback_sessions(subEEG);
end

if strcmp(modality, 'PF')
    if sub == 7
        subEEG = remove_pf_extra_markers(subEEG);
    end

    if sub == 15
        subEEG = swap_pf_start_end_labels(subEEG);
    end
end
end

function subEEG = relabel_sub1_no_feedback_sessions(subEEG)
for j = 1:size(subEEG.event, 2)
    if j <= 180 || j > 540
        switch subEEG.event(j).type
            case 'S  4'
                subEEG.event(j).type = 'R  4';
            case 'S  5'
                subEEG.event(j).type = 'R  5';
            case 'S  6'
                subEEG.event(j).type = 'R  6';
            case 'S  7'
                subEEG.event(j).type = 'R  7';
            case 'S  8'
                subEEG.event(j).type = 'R  8';
            case 'S  9'
                subEEG.event(j).type = 'R  9';
        end
    end
end
end

function subEEG = remove_pf_extra_markers(subEEG)
del = 0;
for j = 1:size(subEEG.event, 2)
    if j > 180 && j < 190
        subEEG.event(j - del) = [];
        del = del + 1;
    end
end
end

function subEEG = swap_pf_start_end_labels(subEEG)
for j = 1:size(subEEG.event, 2)
    if strcmp(subEEG.event(j).type, 'S  1')
        subEEG.event(j).type = 'S  3';
    elseif strcmp(subEEG.event(j).type, 'S  3')
        subEEG.event(j).type = 'S  1';
    end
end
end

function trials = extract_trial_records(subEEG)
trials = struct('code', {}, 'start', {}, 'end', {}, 'reported_start', {}, ...
    'reported_end', {}, 'is_correct', {}, 'is_explicit', {}, 'is_long', {}, 'is_flex', {});

for j = 1:size(subEEG.event, 2)
    event_type = subEEG.event(j).type;
    motion = get_motion_definition(event_type);
    if isempty(motion)
        continue;
    end

    reported_start = answer_marker_to_state(subEEG.event(j + 2).type);
    reported_end = answer_marker_to_state(subEEG.event(j + 4).type);

    trials(end + 1) = struct( ...
        'code', event_type(end), ...
        'start', motion.start, ...
        'end', motion.end, ...
        'reported_start', reported_start, ...
        'reported_end', reported_end, ...
        'is_correct', double(reported_start == motion.start && reported_end == motion.end), ...
        'is_explicit', event_type(1) == 'S', ...
        'is_long', motion.isLong, ...
        'is_flex', motion.isFlex);
end
end

function accuracy = compute_motion_accuracy(trials, motion_codes)
accuracy = NaN(1, numel(motion_codes));
trial_codes = {trials.code};

for i = 1:numel(motion_codes)
    idx = strcmp(trial_codes, motion_codes{i});
    accuracy(i) = mean([trials(idx).is_correct]);
end
end

function continuity = compute_continuity_mask(trials, trials_per_session)
n_sessions = numel(trials) / trials_per_session;
continuity = false(1, (trials_per_session - 1) * n_sessions);

cursor = 1;
for session_idx = 1:n_sessions
    session_start = (session_idx - 1) * trials_per_session + 1;
    session_end = session_idx * trials_per_session;
    session_trials = trials(session_start:session_end);

    for trial_idx = 1:(trials_per_session - 1)
        continuity(cursor) = session_trials(trial_idx).reported_end == session_trials(trial_idx + 1).reported_start;
        cursor = cursor + 1;
    end
end
end

function [length_accuracy, direction_accuracy] = compute_length_direction_accuracy(trials)
true_length = [trials.is_long];
true_direction = [trials.is_flex];

reported_length = zeros(1, numel(trials));
reported_direction = zeros(1, numel(trials));

for i = 1:numel(trials)
    reported_length(i) = double(abs(trials(i).reported_end - trials(i).reported_start) == 2);
    reported_direction(i) = double(trials(i).reported_end > trials(i).reported_start);
end

length_accuracy = mean(reported_length == true_length);
direction_accuracy = mean(reported_direction == true_direction);
end

function corrected_p = benjamini_hochberg(raw_p)
[sorted_p, sort_idx] = sort(raw_p(:));
n_tests = numel(sorted_p);
sorted_adjusted = sorted_p .* n_tests ./ (1:n_tests)';
sorted_adjusted = flipud(cummin(flipud(sorted_adjusted)));
sorted_adjusted(sorted_adjusted > 1) = 1;

corrected_p = zeros(size(sorted_p));
corrected_p(sort_idx) = sorted_adjusted;
corrected_p = reshape(corrected_p, size(raw_p));
end

function state = answer_marker_to_state(marker_type)
state = NaN;

switch marker_type
    case 'S  1'
        state = 1;
    case 'S  2'
        state = 2;
    case 'S  3'
        state = 3;
end
end

function motion = get_motion_definition(event_type)
motion = [];

switch event_type(end)
    case '4'
        motion = struct('start', 1, 'end', 3, 'isLong', true, 'isFlex', true);
    case '5'
        motion = struct('start', 3, 'end', 1, 'isLong', true, 'isFlex', false);
    case '6'
        motion = struct('start', 1, 'end', 2, 'isLong', false, 'isFlex', true);
    case '7'
        motion = struct('start', 2, 'end', 3, 'isLong', false, 'isFlex', true);
    case '8'
        motion = struct('start', 3, 'end', 2, 'isLong', false, 'isFlex', false);
    case '9'
        motion = struct('start', 2, 'end', 1, 'isLong', false, 'isFlex', false);
end
end

function bh = styled_boxplot(data, labels, y_label, y_limits, add_grid, colors)
bh = boxplot(data, labels, 'Widths', 0.5);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
set(bh, 'LineWidth', 2);
ylabel(y_label, 'FontSize', 14, 'FontWeight', 'bold');

if ~isempty(y_limits)
    ylim(y_limits);
end

% Default palette — overridden per figure via optional colors argument
default_colors = [
    0.28, 0.60, 0.42;   % forest green
    0.92, 0.46, 0.32;   % coral
    0.60, 0.60, 0.60;   % gray
    0.15, 0.60, 0.68;   % teal
];
if nargin < 6 || isempty(colors)
    group_colors = default_colors;
else
    group_colors = colors;
end

n_groups = size(data, 2);
n_points = size(data, 1);

hold on;

% Color the boxes and overlay individual data points
rng(0);
h_box = findobj(gca, 'Tag', 'Box');
for j = 1:n_groups
    % Box index is reversed in MATLAB (last group = first handle)
    box_idx = n_groups - j + 1;
    c = group_colors(min(j, size(group_colors, 1)), :);
    
    % Colored transparent box fill
    if box_idx <= numel(h_box)
        patch(get(h_box(box_idx), 'XData'), get(h_box(box_idx), 'YData'), c, ...
            'FaceAlpha', 0.15, 'EdgeColor', c, 'LineWidth', 2);
    end
    
    % Overlay individual data points with jitter
    jitter = 0.08 * (rand(n_points, 1) - 0.5);
    scatter(ones(n_points, 1) * j + jitter, data(:, j), 40, ...
        c, 'filled', 'MarkerFaceAlpha', 0.5);
end

hold off;

if add_grid
    grid on;
end
end

