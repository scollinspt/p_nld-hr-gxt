# Final Analysis Summary

This is the finalized statistical handoff for Paper 1. No manuscript files or manuscript prose were changed.

## Dataset

- Primary source: `merged_data/hrv_metabolic_merged.csv`
- Validated analysis source: `analysis_validation/hrv_metabolic_merged_validated.csv`
- HRV windows generated: 930
- Usable merged observations: 923
- Subjects in the usable merged dataset: 22
- Observations per subject: mean 41.95, median 39, range 35–61
- HRV window: 300 seconds
- Window advancement: 10 seconds
- Metabolic alignment: per-subject nearest timestamp with a 60-second tolerance

The observations remain highly overlapping repeated measures. Adjacent windows overlap by 290 seconds.

## Serial-dependence solution

The primary serial-correlation-aware model was Gaussian GEE with:

```text
outcome ~ VO2_z + VO2_z^2
participant as the clustering variable
AR(1) working correlation within participant
```

VO2 was centered and scaled to one sample standard deviation for model fitting. The reported coefficients below are therefore changes in the outcome per one-SD increase in VO2, with a quadratic term included. AR(1) was estimated by grid search because direct optimization could not bracket the autoregressive parameter.

Estimated AR(1) working correlations were extremely high:

- mean HR: 0.971
- SDNN: 0.987
- RMSSD: 0.998
- LF power: 0.970
- HF power: 0.981
- LF/HF: 0.988
- Sample Entropy: 0.994

This confirms that the original pooled Pearson correlations and ANOVA substantially overstated the effective amount of independent information.

The detailed results are in:

- `analysis_validation/gee_serial_models.csv`
- `analysis_validation/thinning_sensitivity_models.csv`
- `analysis_validation/resid_lag1_summary.csv`

## Temporal-thinning sensitivity

Temporal thinning was performed by retaining the first available window within each subject and then selecting the next window at least 60, 150, or 300 seconds later. It was used only as a sensitivity analysis.

The direction of the VO2 association was stable across all thinning levels:

- mean HR: positive at 60, 150, and 300 seconds
- SDNN: negative at all thinning levels
- RMSSD: null at all thinning levels
- LF power: negative at all thinning levels
- HF power: negative at all thinning levels
- LF/HF: negative at all thinning levels
- corrected Sample Entropy: positive at all thinning levels

Thinning reduced the apparent precision and changed standardized coefficient magnitudes, as expected. The substantive direction and the important null finding for RMSSD did not change.

Sensitivity observation counts were:

- 60-second spacing: 163 observations
- 150-second spacing: 72 observations
- 300-second spacing: 45 observations

The complete 923-window AR(1)-GEE analysis remains primary because thinning discards data. The thinning results demonstrate that the principal directions are not created solely by the 10-second window spacing.

## Nonlinear model assessment

Linear and quadratic mixed models were compared using maximum likelihood, not REML, because their fixed-effect structures differ. AIC, BIC, likelihood-ratio comparisons, residual skewness, and interpretability were considered.

Quadratic VO2 models had lower AIC and BIC than their corresponding linear random-intercept models for all seven outcomes. Likelihood-ratio improvements were substantial for mean HR, SDNN, RMSSD, LF power, HF power, and Sample Entropy. The LF/HF improvement was smaller but still favored the quadratic specification.

The final functional form is therefore quadratic VO2 for all seven primary outcomes. This should be described as curvature in the observed VO2 range, not as evidence of a specific physiologic threshold.

A spline was not selected because the quadratic specification provided an interpretable low-complexity representation and the dataset contains only 22 participants. The quadratic model is already a sensitivity beyond the primary serial-correlation-aware framework, so additional flexible spline degrees of freedom were not justified for the primary paper analysis.

The ML comparison output is:

- `analysis_validation/model_comparisons_ml.csv`

## Random-slope assessment

Random-slope models were fit as:

```text
outcome ~ VO2_z + (VO2_z | participant)
```

and compared with random-intercept models using maximum likelihood.

All random-slope fits converged in the validation run and had lower AIC than their corresponding random-intercept models. The estimated random-slope variances were nonzero, indicating meaningful between-participant differences in VO2 response slopes. However, the random-slope model is not used as the primary inferential model because the required serial-correlation structure is handled more directly and transparently by the AR(1)-GEE.

Random-slope results are retained as a model-structure sensitivity analysis. The full intercept variance, slope variance, intercept-slope covariance, residual variance, and residual diagnostics are in:

- `analysis_validation/model_diagnostics.csv`

## Transformation decisions and diagnostics

The following outcomes were modeled on the raw scale:

- mean HR
- SDNN

The following strictly positive, right-skewed outcomes were modeled on the natural-log scale:

- RMSSD
- LF power
- HF power
- LF/HF
- Sample Entropy

The raw-scale models showed substantial residual skew for RMSSD, LF power, HF power, LF/HF, and Sample Entropy. Log transformation reduced residual skew and materially improved model fit. For example, the quadratic random-intercept residual skew was approximately:

- RMSSD: 0.66 on the log scale versus 2.11 on the raw scale
- LF power: -0.06 on the log scale versus 4.96 on the raw scale
- HF power: 0.09 on the log scale versus 2.33 on the raw scale
- LF/HF: -0.23 on the log scale versus 1.32 on the raw scale
- Sample Entropy: 0.36 on the log scale versus 0.58 on the raw scale

The log transformation was therefore selected for distributional reasons, not to obtain statistical significance.

The final GEE coefficients for the selected quadratic form are:

| Outcome | Scale | VO2 linear beta | VO2 quadratic beta | SE of linear beta | 95% CI for linear beta | p-value | N | Subjects |
|---|---|---:|---:|---:|---|---:|---:|---:|
| Mean HR | raw | 1.4678 | -0.0741 | 0.2176 | 1.0412 to 1.8944 | <0.001 | 923 | 22 |
| SDNN | raw | -0.6861 | -0.0607 | 0.1073 | -0.8964 to -0.4758 | <0.001 | 923 | 22 |
| RMSSD | log | -0.0093 | 0.0039 | 0.0052 | -0.0195 to 0.0010 | 0.076 | 923 | 22 |
| LF power | log | -0.1319 | -0.0003 | 0.0202 | -0.1714 to -0.0924 | <0.001 | 923 | 22 |
| HF power | log | -0.0626 | 0.0155 | 0.0159 | -0.0938 to -0.0315 | <0.001 | 923 | 22 |
| LF/HF | log | -0.0411 | -0.0139 | 0.0111 | -0.0628 to -0.0194 | <0.001 | 923 | 22 |
| Sample Entropy | log | 0.0289 | 0.0105 | 0.0099 | 0.0095 to 0.0483 | 0.009 | 923 | 22 |

The p-values above are the model-based AR(1)-GEE Wald p-values. The quadratic coefficient should also be reported when describing the fitted trajectories; the table's confidence interval is for the linear VO2 term only.

## Full Sample Entropy sensitivity analysis

The corrected implementation was run across all 930 generated HRV windows, not just three example windows or subjects.

The corrected implementation uses:

- embedding dimension `m = 2`
- standardized RR intervals
- tolerance `r = factor × SD(RR)` before standardization, equivalently `r = factor` on the standardized series
- pairwise template matching with Chebyshev distance
- self-matches excluded
- undefined results returned only when the required template counts are zero

All 930 windows produced defined estimates at every tested tolerance.

| r factor | N windows | N usable | Median | IQR | Mean | SD | Range | Undefined |
|---:|---:|---:|---:|---:|---:|---:|---|---:|
| 0.10 | 930 | 930 | 0.6706 | 0.5896 | 0.7954 | 0.4418 | 0.2387 to 2.0653 | 0.0% |
| 0.15 | 930 | 930 | 0.4160 | 0.4124 | 0.5494 | 0.3881 | 0.1361 to 2.0653 | 0.0% |
| 0.20 | 930 | 930 | 0.2843 | 0.3054 | 0.4103 | 0.3358 | 0.0923 to 1.5395 | 0.0% |
| 0.25 | 930 | 930 | 0.2116 | 0.2320 | 0.3261 | 0.3033 | 0.0708 to 1.3578 | 0.0% |

Pairwise correlations among the four SampEn parameterizations were high:

- r=0.10 versus 0.15: 0.975
- r=0.10 versus 0.20: 0.957
- r=0.10 versus 0.25: 0.918
- r=0.15 versus 0.20: 0.985
- r=0.15 versus 0.25: 0.965
- r=0.20 versus 0.25: 0.989

The corrected Sample Entropy association with VO2 remained positive for every parameterization. Random-intercept models using the corrected all-window values produced VO2 slopes of 0.0246, 0.0206, 0.0162, and 0.0133 for r factors 0.10, 0.15, 0.20, and 0.25, respectively; all were statistically positive with 923 observations and 22 subjects.

The corrected values differ materially from the old values near 0.003 because the original implementation calculated the tolerance in raw milliseconds and then applied it to a standardized series. That was a coding-scale error. The old entropy column should not be used for Paper 1.

Files:

- `analysis_validation/sample_entropy_sensitivity_all_windows.csv`
- `analysis_validation/sample_entropy_sensitivity_summary.csv`
- `analysis_validation/sample_entropy_sensitivity_correlations.csv`
- `analysis_validation/sample_entropy_sensitivity_vo2_models.csv`
- `analysis_validation/sample_entropy_sensitivity_all_windows.png`

## DFA decision

DFA was excluded because the original implementation was incompatible with the available window lengths and no independently validated replacement was established for this analysis.

DFA values must not be included in the final Paper 1 figures, tables, models, or claims.

## Final analysis files

- [analysis_validation/definitive_reanalysis.py](analysis_validation/definitive_reanalysis.py)
- [analysis_validation/hrv_metabolic_merged_validated.csv](analysis_validation/hrv_metabolic_merged_validated.csv)
- [analysis_validation/gee_serial_models.csv](analysis_validation/gee_serial_models.csv)
- [analysis_validation/model_comparisons_ml.csv](analysis_validation/model_comparisons_ml.csv)
- [analysis_validation/model_diagnostics.csv](analysis_validation/model_diagnostics.csv)
- [analysis_validation/thinning_sensitivity_models.csv](analysis_validation/thinning_sensitivity_models.csv)
- [analysis_validation/sample_entropy_sensitivity_summary.csv](analysis_validation/sample_entropy_sensitivity_summary.csv)
- [analysis_validation/sample_entropy_sensitivity_correlations.csv](analysis_validation/sample_entropy_sensitivity_correlations.csv)

## Findings supported for the manuscript

- Within this dataset, mean HR increases with VO2 after accounting for participant clustering and strong within-subject serial dependence.
- SDNN decreases with VO2 under the AR(1)-GEE analysis.
- LF power and HF power decrease with VO2 on the log scale.
- LF/HF decreases rather than increasing monotonically with VO2 in the fitted data.
- RMSSD does not show a statistically reliable VO2 association in the primary AR(1)-GEE model or in the temporal-thinning sensitivity analyses.
- Correctly implemented Sample Entropy is positive and robustly associated with VO2 across r values from 0.10 to 0.25 SD.
- The VO2 relationships are better represented by quadratic than linear functions within the observed VO2 range.
- Participant-specific response slopes vary substantially, supporting heterogeneity in individual trajectories.
- The final inferential sample is 923 observations from 22 subjects.

## Claims not supported by the analysis

- LF/HF as a simple monotonic measure of sympathetic dominance.
- Progressive Sample Entropy decline with exercise intensity.
- “Complexity collapse” as an established physiologic mechanism.
- DFA dynamical transitions or DFA alpha1 findings.
- Autonomic saturation at maximal exercise without a formal saturation analysis.
- Causal claims from these observational repeated-measures associations.
- Treating the 923 overlapping windows as independent cases.
- Calling the percentile-defined VO2 groups Rest, Moderate, High, and Maximal physiological zones.
