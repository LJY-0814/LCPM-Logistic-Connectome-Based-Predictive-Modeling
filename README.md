# LCPM: Weighted Logistic Connectome-Based Predictive Modeling

MATLAB code for predicting suicide risk in major depressive disorder (MDD)
from resting-state functional connectivity using weighted logistic
connectome-based predictive modeling (LCPM), together with permutation
testing, consensus-pattern characterization, and virtual lesion analysis.

## Overview

The pipeline implements the following analyses:

1. **Repeated k-fold cross-validation of the LCPM classifier**
   (`func/lcpm_repeated_cv.m`)
   - In each of `cfg.n_iterations` iterations, subjects are randomly split
     into `cfg.n_folds` folds.
   - **Feature selection** (training folds only): point-biserial
     correlation between each functional connection and the binary label;
     edges with `p < cfg.p_threshold` are kept and split into positively
     and negatively correlated sets.
   - **Model**: composite feature (sum of positive edges minus sum of
     negative edges) together with covariates (HAMD-17 total score, sex,
     age, years of education), standardized on the
     training folds, entered into a logistic regression weighted inversely
     proportional to class size.
   - **Evaluation**: AUC (`perfcurve`) and accuracy / sensitivity /
     specificity at the fixed probability cutoff `cfg.class_threshold`.
2. **Label-permutation test** (`func/lcpm_permutation.m`) — one
   permutation per original repeated-CV split (default 1,000). For each
   split the labels are shuffled once and a single k-fold CV is run; the
   resulting single-CV AUCs are compared to the observed mean AUC to obtain
   a one-tailed p value.
3. **Consensus-pattern characterization** (`func/lcpm_characterize_pattern.m`)
   - *High-frequency edges*: edges selected in >= `cfg.freq_threshold` of
     all iterations x folds.
   - *Node degree*: number of high-frequency edges incident to each ROI.
   - *Network-level distribution*: number of high-frequency edges within 
     and between functional networks.
4. **Virtual lesion analysis** (`func/lcpm_virtual_lesion.m`) — for each
   functional network, all edges touching that network are removed and the
   full pipeline is re-run on the same splits; a one-tailed Wilcoxon
   signed-rank test (Benjamini-Hochberg FDR across networks) identifies
   networks whose removal significantly lowers the AUC.

## Repository structure

```
LCPM/
├── main.m                              % entry point: runs the full pipeline
├── config.m                            % all paths and thresholds
├── func/
│   ├── load_data.m                     % dataset-specific data loader
│   ├── lcpm_check_errors.m             % input validation
│   ├── lcpm_prepare_edges.m            % upper-triangle vectorization + Fisher z
│   ├── lcpm_repeated_cv.m              % core repeated k-fold CV engine
│   ├── lcpm_permutation.m              % label-permutation test
│   ├── lcpm_characterize_pattern.m     % high-frequency edges / degrees / networks
│   └── lcpm_virtual_lesion.m           % leave-one-network-out analysis
├── README.md
└── .gitignore
```

## Usage

1. Edit `config.m`:
   - set `cfg.data_dir`, `cfg.roi_table_file`, `cfg.label_file` to your data;
   - adjust the thresholds if needed — the three thresholds are separate
     parameters (`p_threshold` for feature selection, `class_threshold`
     for classification, `freq_threshold` for the consensus pattern), so
     the pipeline is not tied to the manuscript settings.
2. Run `main` from the `LCPM/` directory.

### Input data

The provided `load_data.m` expects the layout used in the manuscript:

- one `ROISignals_S*.mat` file per subject, containing a `ROISignals`
  time-series matrix (DPARSF/DPABI output);
- an ROI Excel table whose column 1 lists the ROI column indices, column 2
  the region names, and column 8 the functional network labels;
- a label Excel spreadsheet (no header) whose columns are indexed by
  `cfg.label_column` (binary outcome) and `cfg.covariate_columns`.

To apply LCPM to a different dataset, replace `load_data.m` (or write your
own loader). The downstream functions only need a 3-D
`ROI x ROI x subjects` connectivity array, a binary label vector, a
covariate matrix, and per-ROI region/network labels.

### Outputs (written to `cfg.output_dir`)

| File | Content |
| --- | --- |
| `cv_summary.xlsx` | mean +/- SD of AUC, accuracy, sensitivity, specificity, edge count over valid iterations |
| `permutation_results.xlsx` | observed vs. permutation AUC and one-tailed p value |
| `selected_edges_high_frequency.xlsx` | consensus edges (ROIs, regions, selection count/frequency, sign) |
| `node_degree_high_frequency.xlsx` | degree of each ROI in the consensus network |
| `network_edge_count.xlsx` | number of high-frequency edges within/between networks |
| `network_edge_avg_strength.xlsx` | mean Fisher-z strength within/between networks |
| `virtual_lesion_results.xlsx` | per-network AUC after exclusion, Wilcoxon p, FDR significance |
| `lcpm_results.mat` | all result structs for downstream use |

## Citation