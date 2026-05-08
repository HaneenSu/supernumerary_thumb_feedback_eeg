# Supernumerary Thumb EEG

Code repository for the EEG study accompanying the paper:

> **Pressure Feedback Enhances Supernumerary Robotic Thumb
Proprioception over Vibrotactile Feedback: Behavioral and EEG
Evidence**

This repository contains MATLAB scripts for behavioral analysis, EEG preprocessing, time-frequency analysis, and functional connectivity analysis comparing two sensory feedback modalities during supernumerary thumb use.

---

## Study Overview

Participants operated a supernumerary robotic thumb and received one of two types of sensory feedback:
- **PF** — Pressure feedback
- **VF** — Vibrotactile feedback

EEG was recorded while participants performed explicit and implicit thumb articulation tasks. Analyses focus on behavioral identification accuracy and EEG markers of sensorimotor learning: alpha/beta event-related desynchronization (ERD) and beta-band functional connectivity (small-worldness).

---

## Repository Structure

```
supernumerary_thumb_eeg/
├── code/                          — Analysis scripts (run in order)
│   ├── behavioral_analysis_01.m   — Behavioral accuracy (Figures 2, Supp1, Supp2)
│   ├── eeg_preprocessing_PF_02.m  — EEG preprocessing for PF group
│   ├── eeg_preprocessing_VF_03.m  — EEG preprocessing for VF group
│   ├── timefrequency_analysis_04.m — Morlet TF analysis & ERD figures
│   ├── connectivity_analysis_05.m  — wPLI connectivity & small-worldness
│   ├── EEG_chlocs.mat             — 60-channel electrode locations
│   ├── NewEasyCap63.mat           — EasyCap channel layout
│   ├── systems.mat                — Brain system layout for network plots
│   └── drawline_thick.m           — Helper for network edge rendering
│
├── rawdata/
│   ├── SFP/                       — Raw EEG data, PF group (.vhdr files)
│   └── SFV/                       — Raw EEG data, VF group (.vhdr files)
│
├── derivatives/
│   ├── epoched/
│   │   ├── PF/                    — Preprocessed & epoched EEG, PF group
│   │   └── VF/                    — Preprocessed & epoched EEG, VF group
│   ├── TF/                        — Time-frequency .mat files
│   └── wPLI/                      — Connectivity adjacency .mat files
│
└── figs/                          — Output figures and statistics.txt
```

---

## Usage

Scripts are designed to be run **in order** from the `code/` directory:

1. **`behavioral_analysis_01.m`** — Load raw EEG, extract trial records, compute explicit/implicit identification accuracy. Produces Figures 2, Supp1, Supp2 and writes statistics to `../figs/statistics.txt`.

2. **`eeg_preprocessing_PF_02.m`** — Preprocess raw PF EEG: bandpass filter (0.5–85 Hz), notch filter (50 Hz), ASR artifact rejection, channel interpolation, CAR re-referencing, ICA + ICLabel eye artifact removal, epoching. Saves `.set` files to `../derivatives/epoched/PF/`.

3. **`eeg_preprocessing_VF_03.m`** — Same pipeline as PF, with an additional +1000 ms latency correction for the VF hardware trigger delay. Saves to `../derivatives/epoched/VF/`.

4. **`timefrequency_analysis_04.m`** — Morlet wavelet convolution on epoched EEG, dB baseline normalisation, alpha/beta topoplot and band-power time-series figures comparing PF vs VF.

5. **`connectivity_analysis_05.m`** — Time-resolved beta-band wPLI adjacency matrices (FieldTrip), small-worldness index (BCT), PF vs VF Wilcoxon comparison, and network visualisation.

---

## Dependencies

| Toolbox | Used by | Notes |
|---|---|---|
| [EEGLAB](https://sccn.ucsd.edu/eeglab/) | Scripts 01–04 | Tested with EEGLAB 2021+ |
| [FieldTrip](https://www.fieldtriptoolbox.org/) | Script 05 | Update path in `fieldtrip_candidates` |
| [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) | Script 05 | Update `addpath` to your BCT location |

Toolboxes are **not** included in this repository. Download and install them separately, then update the paths at the top of the relevant scripts.

---

## Notes

- Raw data files follow the naming convention `SFP0000.vhdr` (PF) and `SFV0000.vhdr` (VF), with zero-indexed subject numbers.
- Scripts 02 and 03 are computationally intensive (ICA per subject). Run overnight or on a cluster.
- Script 05 small-worldness computation uses 1000 random network realisations per subject — also time-intensive.
