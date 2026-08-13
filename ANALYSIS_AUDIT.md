# Analysis Audit for the HRV–Metabolic Exercise Dataset

## Scope

This audit was performed by inspecting the code and data as they actually exist in the repository, rather than relying on README statements or earlier manuscript notes. The focus was to determine what the current pipeline really does, whether the Sample Entropy and DFA calculations are valid, and whether the reported inferential analysis is statistically defensible.

## 1) Exact dataset used

The primary analysis dataset is:

- `merged_data/hrv_metabolic_merged.csv`

This file is created by `merge_data.py` from:

- `python_analysis/hrv_python_new_advancement.csv` (HRV windows)
- `all_bxb_10secondEpochs.csv` (metabolic epochs)
- `Old-to-New-IDs.csv` (subject mapping)

The merged file contains 923 rows and 34 columns after filtering.

## 2) Number of subjects

Actual subject count in the working merged dataset:

- 22 subjects

This differs from the README narrative, which states 23 athletes. The subject count is determined by the actual merged data after filtering and ID mapping:

- `df_merged['subject_id_new'].nunique()` = 22

One subject is not represented in the final merged clean dataset after mapping and QC.

## 3) Observations per subject

From the merged data:

- Mean observations per subject: 41.95
- Median observations per subject: 39
- Min: 35
- Max: 61

Distribution summary:

- count = 22
- std = 7.54
- 25% = 37.25
- 75% = 45

This is a repeated-measures dataset with 22 clusters and 923 repeated windows, not 923 independent observations.

## 4) HRV window length and advancement

The HRV pipeline in `hrv_analysis.py` defines:

- `WINDOW_DURATION = 300.0` seconds
- `NEW_ADVANCEMENT = 10.0` seconds
- `MIN_WINDOW_RR_COUNT = 200`

The analysis script processes sliding windows of 300 seconds with 10-second advancement.

The code also retains the original 30-second advance pipeline for comparison, but the active manuscript analysis uses the 10-second advancement dataset.

This means consecutive windows are highly overlapping:

- each 300-second window overlaps the preceding one by 290 seconds
- adjacent windows are therefore not independent observations

## 5) Metabolic alignment method

The merge code uses a per-subject nearest-neighbor temporal alignment:

- `pd.merge_asof(df_h, df_m, on='timestamp', direction='nearest', tolerance=60)`
- performed separately for each subject
- HRV windows are aligned to the nearest metabolic epoch within 60 seconds

The code also does:

- `df_hrv['timestamp'] = df_hrv['window_start']`
- `df_met['timestamp'] = df_met['epoch'] * 10.0`

This is a nearest-neighbor alignment, not a true physiological synchronization of the entire exercise test.

## 6) Intensity-group definitions

The classification used by `Paper1/paper1_autonomic_metabolic_coupling_analysis.py` is explicitly based on VO2 quantiles:

```python
vo2_quantiles = self.df[vo2_col].quantile([0.1, 0.3, 0.7])

if vo2 < vo2_quantiles[0.1]:
    return 'Rest'
elif vo2 < vo2_quantiles[0.3]:
    return 'Moderate'
elif vo2 < vo2_quantiles[0.7]:
    return 'High'
else:
    return 'Maximal'
```

This is not a physiological threshold system (e.g., VT, RCP, %VO2peak, %HRpeak). It is a distribution-defined grouping based on observed VO2 values.

The code labels these groups as:

- Rest
- Moderate
- High
- Maximal

But the actual definitions are percentiles, not physiologically defined exercise zones.

## 7) Exact implementation of Sample Entropy

The Sample Entropy code is in `hrv_analysis.py`:

```python
@staticmethod
def sample_entropy(rr_intervals, m=2, r=None):
    if len(rr_intervals) < m + 1:
        return np.nan

    if r is None:
        r = 0.2 * np.std(rr_intervals, ddof=1)

    rr_norm = (rr_intervals - np.mean(rr_intervals)) / np.std(rr_intervals, ddof=1)

    def count_patterns(x, m, r):
        n = len(x)
        count = 0
        for i in range(n - m):
            for j in range(i + 1, n - m + 1):
                if np.max(np.abs(x[i:i+m] - x[j:j+m])) <= r:
                    count += 1
        return count

    phi_m = count_patterns(rr_norm, m, r)
    phi_m1 = count_patterns(rr_norm, m + 1, r)

    if phi_m1 == 0:
        return np.nan

    sampen = -np.log(phi_m1 / phi_m)
    return sampen
```

Important details:

- `m = 2` is the embedding dimension
- default tolerance is `r = 0.2 * SD(RR)`
- RR intervals are normalized by subtracting the mean and dividing by the standard deviation before counting matches
- self-matches are excluded because the inner loop starts at `j = i + 1`
- the implementation counts pairwise template matches among all subsequences
- if `phi_m1 == 0`, the value is returned as NaN

This is a standard algorithmic form, but the result is highly sensitive to how regular the RR sequence is, how the tolerance is scaled, and how many matching templates exist in a window.

## 8) Exact implementation of DFA

The DFA function is in `hrv_analysis.py`:

```python
@staticmethod
def detrended_fluctuation_analysis(rr_intervals, min_scale=10, max_scale=500):
    if len(rr_intervals) < max_scale * 2:
        return {'dfa1': np.nan, 'dfa2': np.nan}

    rr_norm = rr_intervals - np.mean(rr_intervals)
    y = np.cumsum(rr_norm)

    scales = np.logspace(np.log10(min_scale), np.log10(max_scale), 30).astype(int)
    scales = np.unique(scales)
    fluctuations = []

    for scale in scales:
        n_full = len(y) // scale
        y_forward = y[:n_full * scale].reshape(n_full, scale)

        fit_forward = np.polyfit(np.arange(scale), y_forward[0], 1)

        fluct = 0
        for i in range(n_full):
            trend = np.polyval(fit_forward, np.arange(scale))
            fluct += np.mean((y_forward[i] - trend) ** 2)

        fluctuations.append(np.sqrt(fluct / n_full))

    log_scales = np.log10(scales)
    log_fluct = np.log10(fluctuations)

    mid_idx = len(scales) // 2
    slope1 = np.polyfit(log_scales[:mid_idx], log_fluct[:mid_idx], 1)[0]
    slope2 = np.polyfit(log_scales[mid_idx:], log_fluct[mid_idx:], 1)[0]

    return {'dfa1': slope1, 'dfa2': slope2}
```

Problems in the implementation:

- the code requires `len(rr_intervals) >= max_scale * 2`, and `max_scale` is set to 500
- therefore a sequence shorter than 1000 beats is declared invalid
- most windows in this dataset have ~507–1460 RR intervals, with mean ~752 RR intervals
- this means a large fraction of the exercise windows fail the DFA precondition
- in addition, the trend is fitted only to `y_forward[0]` instead of each segment, which is not the standard DFA procedure

This is the direct reason the current analysis reports no valid DFA alpha1 observations.

## 9) Missing-data handling

The current pipeline does the following:

- `hrv_analysis.py` calculates HRV only if `len(rr_intervals) >= MIN_WINDOW_RR_COUNT` (200 RR intervals)
- `remove_ectopic_beats` removes RR values outside a ±30% median threshold
- `merge_data.py` drops rows with missing `meanrr` or `vo2_mlkgmin_median`
- `merge_data.py` also removes HR outliers using `meanhr` outside 40–200 bpm

Thus, the dataset is cleaned in multiple stages, but there is no formal subject-level missingness model or explicit handling of non-random missingness.

## 10) Artifact/outlier handling

Artifact handling is limited and simple:

- median-based ectopic filtering at ±30%
- HR outlier filter at `40 <= meanhr <= 200`
- no explicit validation of all HRV windows for stationarity, non-stationary noise, or ectopic beat correction quality

This approach is pragmatic but not comprehensive.

## 11) Current statistical models

The analysis script used for the paper currently relies on:

- pooled Pearson correlations
- one-way ANOVA across exercise-stage groups
- no subject-level random effects
- no correction for within-subject repeated windows
- no autocorrelation modeling

The relevant code is in `Paper1/paper1_autonomic_metabolic_coupling_analysis.py`:

- `correlation_analysis()` uses `scipy.stats.pearsonr`
- `intensity_effects()` uses `stats.f_oneway()`
- `create_intensity_stages()` defines groups by VO2 percentiles

This approach treats 923 repeated windows as independent observations, which is not statistically valid.

## 12) Discrepancies between code, README, tables, and manuscript claims

There are several important discrepancies.

### A. Subject count mismatch

- README claims 23 athletes
- actual merged analysis data contains 22 subjects

### B. Window count mismatch

- README says 930 windows per subject in some places
- the actual HRV output has 930 total windows across all subjects, with mean ~41 windows per subject

### C. Intensity groups are not physiological thresholds

- the code defines classes by VO2 percentiles, not by ventilatory threshold, RCP, or %VO2peak
- the manuscript labels therefore risk implying physiologic exercise zones that are not actually defined in the data

### D. DFA is absent by design

- the current code is not reporting valid DFA because it fails preconditions at the scale definition level
- the manuscript must not imply valid DFA measurements unless the calculation is repaired and validated

### E. Sample Entropy values are strongly parameter-sensitive

- the code uses a very small tolerance and normalizes the RR series
- the resulting values are near zero when the RR series becomes regular, but this must be checked against sensitivity analysis rather than assumed to be a biological result without validation

### F. Repeated windows are treated as independent

- pooled correlations and ANOVA are not defensible for this dataset structure
- the real inferential analysis must use a mixed-effects or clustered repeated-measures framework

## Immediate conclusion

The current analyses should not be treated as finalized inferential evidence. The critical issues are:

1. repeated windows are nested within subjects
2. pooled tests ignore clustering and serial dependence
3. Sample Entropy values are extremely small and require explicit validation
4. DFA alpha1 is not valid under the current implementation because the code blocks it before analysis
5. intensity group definitions are distribution-based rather than physiologic

The next stage is to replace pooled inferential analyses with subject-aware repeated-measures models, validate the entropy implementation across parameter choices, and only then decide what claims the manuscript can support.
