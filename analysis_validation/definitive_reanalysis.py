#!/usr/bin/env python3
"""Definitive Paper 1 statistical validation.

No manuscript files are touched. Outputs are written under analysis_validation/.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, skew
from statsmodels.formula.api import mixedlm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.cov_struct import Autoregressive, Independence
from statsmodels.genmod.families import Gaussian

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis_validation'
OUT.mkdir(exist_ok=True)
DATA = ROOT / 'merged_data' / 'hrv_metabolic_merged.csv'
HRV = ROOT / 'python_analysis' / 'hrv_python_new_advancement.csv'
RR_DIR = ROOT / 'rr'
ID_MAP = ROOT / 'Old-to-New-IDs.csv'
OUTCOMES = ['meanhr', 'sdnn', 'rmssd', 'lfpowfft', 'hfpowfft', 'lfhffft', 'sampen']
POSITIVE_OUTCOMES = ['rmssd', 'lfpowfft', 'hfpowfft', 'lfhffft', 'sampen']


def sampen(rr: np.ndarray, m: int = 2, r_factor: float = 0.2) -> float:
    """Exact Chebyshev template-count Sample Entropy, excluding self matches."""
    rr = np.asarray(rr, dtype=float)
    if len(rr) < m + 1:
        return np.nan
    sd = np.std(rr, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return np.nan
    z = (rr - np.mean(rr)) / sd
    radius = r_factor
    def count(embedding: int) -> int:
        templates = np.lib.stride_tricks.sliding_window_view(z, embedding)
        # query_pairs counts each unordered pair once, exactly matching i < j.
        return int(cKDTree(templates).query_pairs(radius, p=np.inf, output_type='ndarray').shape[0])
    b = count(m)
    a = count(m + 1)
    if b == 0 or a == 0:
        return np.nan
    return float(-np.log(a / b))


def rr_by_subject():
    mapping_df = pd.read_csv(ID_MAP)
    mapping = dict(zip(mapping_df['id-old'], mapping_df['id-new']))
    result = {}
    for old, new in mapping.items():
        path = RR_DIR / f'{old}.txt'
        if path.exists():
            result[int(new)] = np.loadtxt(path)
    return result


def full_entropy_sensitivity():
    hrv = pd.read_csv(HRV)
    rr_map = rr_by_subject()
    records = []
    for row in hrv.itertuples(index=False):
        times = rr_map.get(int(row.subject_id_new))
        if times is None:
            continue
        window = times[(times >= row.window_start) & (times <= row.window_end)]
        rr = np.diff(window) * 1000.0
        for rf in [0.10, 0.15, 0.20, 0.25]:
            value = sampen(rr, m=2, r_factor=rf)
            records.append({'subject_id_new': int(row.subject_id_new),
                            'window_start': row.window_start, 'rr_count': len(rr),
                            'r_factor': rf, 'sampen': value})
    result = pd.DataFrame(records)
    result.to_csv(OUT / 'sample_entropy_sensitivity_all_windows.csv', index=False)
    wide = result.pivot_table(index=['subject_id_new', 'window_start'], columns='r_factor', values='sampen')
    summary = []
    for rf in [0.10, 0.15, 0.20, 0.25]:
        x = result.loc[result.r_factor == rf, 'sampen']
        q = x.quantile([0.25, 0.75])
        summary.append({'r_factor': rf, 'n_windows': len(x), 'n_usable': x.notna().sum(),
                        'median': x.median(), 'iqr': q.loc[0.75] - q.loc[0.25],
                        'mean': x.mean(), 'sd': x.std(), 'min': x.min(), 'max': x.max(),
                        'undefined_proportion': x.isna().mean()})
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT / 'sample_entropy_sensitivity_summary.csv', index=False)
    corr_rows = []
    for a in [0.10, 0.15, 0.20, 0.25]:
        for b in [0.10, 0.15, 0.20, 0.25]:
            if a < b:
                pair = wide[[a, b]].dropna()
                corr_rows.append({'r_a': a, 'r_b': b, 'n_pairs': len(pair),
                                  'pearson_r': pearsonr(pair[a], pair[b]).statistic if len(pair) > 2 else np.nan,
                                  'p_value': pearsonr(pair[a], pair[b]).pvalue if len(pair) > 2 else np.nan})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(OUT / 'sample_entropy_sensitivity_correlations.csv', index=False)
    assoc = []
    merged = pd.read_csv(DATA)
    for rf in [0.10, 0.15, 0.20, 0.25]:
        e = result[result.r_factor == rf][['subject_id_new', 'window_start', 'sampen']]
        d = merged.merge(e, on=['subject_id_new', 'window_start'], how='left') if 'window_start' in merged else None
        if d is None or d['sampen_y'].notna().sum() == 0:
            # merged data stores the alignment timestamp as timestamp; HRV window_start is timestamp.
            d = merged.merge(e.rename(columns={'window_start': 'timestamp'}), on=['subject_id_new', 'timestamp'], how='left')
            entropy_col = 'sampen_y'
        else:
            entropy_col = 'sampen_y'
        d = d[['subject_id_new', 'vo2_mlkgmin_median', entropy_col]].dropna()
        if len(d):
            fit = mixedlm(f'{entropy_col} ~ vo2_mlkgmin_median', d, groups=d.subject_id_new).fit(method='lbfgs', maxiter=200, reml=False)
            assoc.append({'r_factor': rf, 'n_obs': len(d), 'n_subjects': d.subject_id_new.nunique(),
                          'beta_vo2': fit.params['vo2_mlkgmin_median'], 'se': fit.bse['vo2_mlkgmin_median'],
                          'ci_low': fit.conf_int().loc['vo2_mlkgmin_median', 0],
                          'ci_high': fit.conf_int().loc['vo2_mlkgmin_median', 1],
                          'p_value': fit.pvalues['vo2_mlkgmin_median']})
    pd.DataFrame(assoc).to_csv(OUT / 'sample_entropy_sensitivity_vo2_models.csv', index=False)
    fig, ax = plt.subplots(figsize=(8, 5))
    for rf, group in result.groupby('r_factor'):
        ax.scatter(group['rr_count'], group['sampen'], s=8, alpha=.25, label=f'r={rf:.2f}')
    ax.set(xlabel='RR intervals per window', ylabel='Sample Entropy', title='All-window Sample Entropy sensitivity')
    ax.legend()
    fig.tight_layout(); fig.savefig(OUT / 'sample_entropy_sensitivity_all_windows.png', dpi=200); plt.close(fig)
    return summary_df, corr_df, pd.DataFrame(assoc)


def prep(df, outcome, transform=False, thin=None):
    cols = ['subject_id_new', 'timestamp', 'vo2_mlkgmin_median', outcome]
    d = df[cols].dropna().sort_values(['subject_id_new', 'timestamp']).copy()
    if thin:
        kept = []
        for _, g in d.groupby('subject_id_new', sort=False):
            last = -np.inf
            for idx, row in g.iterrows():
                if row.timestamp - last >= thin:
                    kept.append(idx); last = row.timestamp
        d = d.loc[kept].sort_values(['subject_id_new', 'timestamp'])
    d['vo2_z'] = (d['vo2_mlkgmin_median'] - d['vo2_mlkgmin_median'].mean()) / d['vo2_mlkgmin_median'].std(ddof=1)
    d['y'] = np.log(d[outcome]) if transform else d[outcome]
    return d


def fit_ml(d, formula, re_formula='1'):
    model = mixedlm(formula, d, groups=d.subject_id_new, re_formula=re_formula)
    return model.fit(method='lbfgs', maxiter=500, reml=False, disp=False)


def model_comparisons(df):
    rows = []
    diagnostics = []
    for outcome in OUTCOMES:
        for transformed in [False, True] if outcome in POSITIVE_OUTCOMES else [False]:
            d = prep(df, outcome, transformed)
            label = 'log' if transformed else 'raw'
            try:
                linear = fit_ml(d, 'y ~ vo2_z')
                quad = fit_ml(d, 'y ~ vo2_z + I(vo2_z**2)')
                slope = fit_ml(d, 'y ~ vo2_z', re_formula='1 + vo2_z')
                lr = max(0.0, 2 * (quad.llf - linear.llf))
                rows += [
                    {'outcome': outcome, 'transform': label, 'model': 'linear_random_intercept', 'aic': linear.aic, 'bic': linear.bic, 'llf': linear.llf, 'converged': linear.converged, 'lr_vs_linear': np.nan},
                    {'outcome': outcome, 'transform': label, 'model': 'quadratic_random_intercept', 'aic': quad.aic, 'bic': quad.bic, 'llf': quad.llf, 'converged': quad.converged, 'lr_vs_linear': lr},
                    {'outcome': outcome, 'transform': label, 'model': 'linear_random_slope', 'aic': slope.aic, 'bic': slope.bic, 'llf': slope.llf, 'converged': slope.converged, 'lr_vs_linear': 2 * (slope.llf - linear.llf)},
                ]
                for fit, name in [(linear, 'linear_random_intercept'), (quad, 'quadratic_random_intercept'), (slope, 'linear_random_slope')]:
                    diagnostics.append({'outcome': outcome, 'transform': label, 'model': name,
                                         'resid_sd': np.std(fit.resid, ddof=1), 'resid_skew': skew(fit.resid),
                                         'random_intercept_var': fit.cov_re.iloc[0, 0] if fit.cov_re.shape else np.nan,
                                         'random_slope_var': fit.cov_re.iloc[1, 1] if fit.cov_re.shape[0] > 1 else np.nan,
                                         'intercept_slope_cov': fit.cov_re.iloc[0, 1] if fit.cov_re.shape[0] > 1 else np.nan,
                                         'residual_var': fit.scale})
            except Exception as exc:
                rows.append({'outcome': outcome, 'transform': label, 'model': 'failed', 'error': str(exc)})
    pd.DataFrame(rows).to_csv(OUT / 'model_comparisons_ml.csv', index=False)
    pd.DataFrame(diagnostics).to_csv(OUT / 'model_diagnostics.csv', index=False)


def gee_models(df):
    rows = []
    for outcome in OUTCOMES:
        transformed = outcome in POSITIVE_OUTCOMES
        d = prep(df, outcome, transformed)
        for form_name, formula in [('linear', 'y ~ vo2_z'), ('quadratic', 'y ~ vo2_z + I(vo2_z**2)')]:
          for corr_name, corr in [('AR1', Autoregressive(grid=True)), ('independence', Independence())]:
            try:
                fit = GEE.from_formula(formula, groups='subject_id_new', time='timestamp', cov_struct=corr, family=Gaussian(), data=d).fit()
                rows.append({'outcome': outcome, 'transform': 'log' if transformed else 'raw', 'form': form_name, 'correlation': corr_name,
                             'n_obs': len(d), 'n_subjects': d.subject_id_new.nunique(), 'beta_vo2_z': fit.params['vo2_z'],
                             'beta_vo2_z2': fit.params.get('I(vo2_z ** 2)', np.nan),
                             'se': fit.bse['vo2_z'], 'ci_low': fit.conf_int().loc['vo2_z', 0], 'ci_high': fit.conf_int().loc['vo2_z', 1],
                             'p_value': fit.pvalues['vo2_z'], 'working_alpha': getattr(fit.cov_struct, 'dep_params', np.nan),
                             'qic': fit.qic()[0] if corr_name == 'AR1' else np.nan})
            except Exception as exc:
                rows.append({'outcome': outcome, 'transform': 'log' if transformed else 'raw', 'form': form_name, 'correlation': corr_name, 'error': str(exc)})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / 'gee_serial_models.csv', index=False)
    return result


def thinning_models(df):
    rows = []
    for outcome in OUTCOMES:
        transformed = outcome in POSITIVE_OUTCOMES
        for spacing in [60, 150, 300]:
            d = prep(df, outcome, transformed, spacing)
            try:
                fit = fit_ml(d, 'y ~ vo2_z')
                rows.append({'outcome': outcome, 'transform': 'log' if transformed else 'raw', 'spacing_seconds': spacing,
                             'n_obs': len(d), 'n_subjects': d.subject_id_new.nunique(), 'beta_vo2_z': fit.params['vo2_z'],
                             'se': fit.bse['vo2_z'], 'ci_low': fit.conf_int().loc['vo2_z', 0], 'ci_high': fit.conf_int().loc['vo2_z', 1],
                             'p_value': fit.pvalues['vo2_z'], 'aic': fit.aic, 'converged': fit.converged})
            except Exception as exc:
                rows.append({'outcome': outcome, 'spacing_seconds': spacing, 'error': str(exc)})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / 'thinning_sensitivity_models.csv', index=False)
    return result


def main():
    df = pd.read_csv(DATA)
    print('Running full-window Sample Entropy sensitivity...')
    summary, correlations, associations = full_entropy_sensitivity()
    print(summary.round(6).to_string(index=False))
    print(associations.round(6).to_string(index=False))
    # Replace the legacy, incorrectly scaled entropy column with the validated
    # r=0.20 * SD implementation before fitting the final models.
    entropy = pd.read_csv(OUT / 'sample_entropy_sensitivity_all_windows.csv')
    entropy = entropy[entropy['r_factor'] == 0.20].rename(columns={'window_start': 'timestamp', 'sampen': 'sampen_validated'})
    entropy = entropy[['subject_id_new', 'timestamp', 'sampen_validated']]
    df = df.drop(columns=['sampen']).merge(entropy, on=['subject_id_new', 'timestamp'], how='left')
    df['sampen'] = df['sampen_validated']
    df = df.drop(columns=['sampen_validated'])
    df.to_csv(OUT / 'hrv_metabolic_merged_validated.csv', index=False)
    print('Validated merged dataset:', len(df), 'rows; entropy missing:', df['sampen'].isna().sum())
    print('Running ML model comparisons...')
    model_comparisons(df)
    print('Running GEE AR(1)/independence models...')
    gee = gee_models(df)
    print(gee[['outcome', 'correlation', 'beta_vo2_z', 'se', 'p_value', 'working_alpha']].to_string(index=False))
    print('Running temporal thinning sensitivity...')
    thin = thinning_models(df)
    print(thin[['outcome', 'spacing_seconds', 'n_obs', 'beta_vo2_z', 'se', 'p_value']].to_string(index=False))
    with open(OUT / 'dfa_decision.json', 'w') as f:
        json.dump({'decision': 'exclude', 'reason': 'DFA was excluded because the original implementation was incompatible with the available window lengths and no independently validated replacement was established for this analysis.'}, f, indent=2)


if __name__ == '__main__':
    main()
