% SPF_pre_longshort.m
%
% EEG preprocessing pipeline for the Pressure Feedback (PF) group.
% For each subject: loads raw .vhdr, applies bandpass (0.5–85 Hz) and
% notch (50 Hz) filters, runs ASR artifact rejection with channel
% interpolation, re-references to CAR, removes eye components via
% ICA + ICLabel, and saves merged long+short articulation epochs.
%
% Output epochs (.set files) are saved to ../derivatives/epoched/PF/.
%
% Expected folder layout (relative to code/):
%   ../rawdata/SFP/              — raw PF EEG data (.vhdr files)
%   ../derivatives/epoched/PF/  — epoched output
%
% Dependencies: EEGLAB, clean_rawdata, ICLabel

clear all;
close all;
clc;

path_to_rawdata = '../rawdata/SFP/';
path_to_epoched = '../derivatives/epoched/PF/';
short_art = ['S  6' 'S  7' 'S  8' 'S  9' 'R  6' 'R  7' 'R  8' 'R  9'];
long_art = ['S  4' 'S  5' 'R  4' 'R  5'];

nSubject = 21;

events = {
    'S  4', 'PF';
    'S  5', 'PF';
    'S  6', 'PF';
    'S  7', 'PF';
    'S  8', 'PF';
    'S  9', 'PF';
    'R  4', 'NPF';
    'R  5', 'NPF';
    'R  6', 'NPF';
    'R  7', 'NPF';
    'R  8', 'NPF';
    'R  9', 'NPF';
};

baseline_period = [-1000 -200];
epoch_period = [-3 5];

load NewEasyCap63

for sub = 1:nSubject
    value_to_display = sub - 1;
    four_digit_text = sprintf('%04d', value_to_display);
    set_file = [path_to_rawdata 'SFP' four_digit_text '.vhdr'];
    disp(['reading...' set_file]);

    EEG = pop_fileio(set_file);
    EEG.chanlocs = struct(chanlocsEasyCapNoRef);
    EEG = eeg_checkset(EEG);

    EEG = pop_chanedit(EEG, 'append', 63, 'changefield', {64 'labels' 'FCz'}, ...
        'changefield', {64 'X' '0.383'}, 'changefield', {64 'Y' '0'}, ...
        'changefield', {64 'Z' '0.923'}, 'convert', {'cart2all'});
    EEG = eeg_checkset(EEG);

    ex_channels = [5 10 21 27];
    EEG = pop_select(EEG, 'nochannel', ex_channels);
    nCh = EEG.nbchan;

    % Band Pass Filter [0.5-85Hz]
    EEG = pop_eegfiltnew(EEG, 0.5, 85, 33000, 0, [], 1);
    EEG = eeg_checkset(EEG);
    close all;

    % Notch filter 50Hz
    EEG = pop_eegfiltnew(EEG, 49.5, 50.5, 8250, 1, [], 1);
    EEG = eeg_checkset(EEG);
    close all;

    % Artifact Subspace Reconstruction (ASR)
    infoCh = {EEG.chanlocs.labels};
    saveEEG = struct(EEG);
    EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 20, 'off');

    nChMiss = EEG.nbchan;
    [M, N] = size(EEG.data);

    % Find and Interpolate Missing Channels
    EEG.data = [EEG.data(:, :); zeros(nCh - nChMiss, N)];
    missChNum = [];
    for i = 1:nCh
        temp = max(strcmp(infoCh(i), {EEG.chanlocs.labels}));
        if temp == 0
            missChNum = [missChNum i];
            EEG.data(i + 1:end, :) = EEG.data(i:end - 1, :);
            EEG.data(i, :) = zeros(1, N);
        end
    end

    EEG.nbchan = saveEEG.nbchan;
    EEG.chanlocs = struct(saveEEG.chanlocs);
    EEG = eeg_checkset(EEG);

    EEG = pop_interp(EEG, missChNum, 'spherical');
    EEG = eeg_checkset(EEG);

    % CAR Re-referencing
    EEG = pop_reref(EEG, [], 'refloc', struct('theta', {0}, 'radius', {0.1252}, ...
        'labels', {'FCz'}, 'sph_theta', {0}, 'sph_phi', {67.4639}, 'X', {0.383}, ...
        'Y', {0}, 'Z', {0.923}, 'sph_radius', {0.99931}, 'type', {''}, 'ref', {''}, ...
        'urchan', {[]}, 'datachan', {0}));
    EEG = eeg_checkset(EEG);

    % ICLabel: Remove Eye Artifacts
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
    EEG = eeg_checkset(EEG);
    EEG = pop_iclabel(EEG, 'default');
    
    eye_idx = find(EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.9);
    fprintf('Subject %d: removing %d eye component(s)\n', sub, length(eye_idx));
    
    EEG = pop_subcomp(EEG, eye_idx, 0);
    EEG = eeg_checkset(EEG);

    % Remove empty markers
    subEEG = struct(EEG);
    del = 0;
    for j = 1:size(subEEG.event, 2)
        if strcmp(subEEG.event(j - del).type, 'empty')
            subEEG.event(j - del) = [];
            del = del + 1;
        end
    end

    % Subject 1 specific event renaming
    if sub == 1
        for j = 1:size(subEEG.event, 2)
            if j <= 180 || j > 540
                if strcmp(subEEG.event(j).type, 'S  4')
                    subEEG.event(j).type = 'R  4';
                elseif strcmp(subEEG.event(j).type, 'S  5')
                    subEEG.event(j).type = 'R  5';
                elseif strcmp(subEEG.event(j).type, 'S  6')
                    subEEG.event(j).type = 'R  6';
                elseif strcmp(subEEG.event(j).type, 'S  7')
                    subEEG.event(j).type = 'R  7';
                elseif strcmp(subEEG.event(j).type, 'S  8')
                    subEEG.event(j).type = 'R  8';
                elseif strcmp(subEEG.event(j).type, 'S  9')
                    subEEG.event(j).type = 'R  9';
                end
            end
        end
    end

    % Epoching
    for evt = 7:12
        if evt == 7 
            epoch_file = ['SF' int2str(sub) '_' events{evt, 2}];
            epoch = pop_epoch(subEEG, events(evt, 1), epoch_period, 'epochinfo', 'yes');
            epoch = pop_rmbase(epoch, baseline_period);
            epoch = eeg_checkset(epoch);
        else
            epochx = pop_epoch(subEEG, events(evt, 1), epoch_period, 'epochinfo', 'yes');
            epochx = pop_rmbase(epochx, baseline_period);
            epochx = eeg_checkset(epochx);
            epoch = pop_mergeset(epoch, epochx);
        end

        if evt == 12
            epoch = pop_saveset(epoch, 'filename', epoch_file, 'filepath', path_to_epoched, 'savemode', 'onefile');
        end
    end
end