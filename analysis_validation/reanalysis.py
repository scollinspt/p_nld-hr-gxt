#!/usr/bin/env python3
"""Subject-aware reanalysis for the HRV exercise dataset.

This script validates the prior Sample Entropy/DFA logic, checks model assumptions,
and fits repeated-measures mixed models using participant as the clustering variable.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.formula.api import mixedlm

PROJECT_DIR = Path(__file__).resolve().parents[1]
MERGED_PATH = PROJECT_DIR / 'merged_data' / 'hrv_metabolic_merged.csv'
HRV_PATH = PROJECT_DIR / 'python_analysis' / 'hrv_python_new_advancement.csv'
RR_DIR = PROJECT_DIR / 'rr'
ID_MAP_PATH = PROJECT_DIR / 'Old-to-New-IDs.csv'
OUT_DIR = PROJECT_DIR / 'analysis_validation'
OUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette('husl')


def sample_entropy(rr, m=2, r=0.2):
    rr = np.asarray(rr, dtype=float)
    if len(rr) < m + 1:
        return np.nan
    if rr.size < 3:
        return np.nan
    rr = rr - rr.mean()
    scale = rr.std(ddof=1)
    if not np.isfinite(scale) or scale <= 0:
        return np.nan
    rr = rr / scale

    def count_mixtures(seq, emb):
        n = len(seq)
        count = 0
        for i in range(n - emb):
            xi = seq[i:i + emb]
            for j in range(i + 1, n - emb + 1):
                if np.max(np.abs(xi - seq[j:j + emb])) <= r:
                    count += 1
        return count

    phi_m = count_mixtures(rr, m)
    phi_m1 = count_mixtures(rr, m + 1)
    if phi_m == 0 or phi_m1 == 0:
        return np.nan
    return -np.log(phi_m1 / phi_m)


def dfa_alpha1(rr, min_scale=10, max_scale=None):
    rr = np.asarray(rr, dtype=float)
    if len(rr) < 30:
        return np.nan
    if max_scale is None:
        max_scale = min(200, max(20, len(rr) // 4))
    rr = rr - np.mean(rr)
    y = np.cumsum(rr)
    scales = np.unique(np.round(np.logspace(np.log10(min_scale), np.log10(max_scale), 20))).astype(int)
    scales = scales[scales >= 2]
    scales = scales[scales < len(y) / 2]
    if len(scales) < 5:
        return np.nan

    fluctuations = []
    for scale in scales:
        n_segments = len(y) // scale
        if n_segments < 2:
            continue
        vals = []
        for seg_idx in range(n_segments):
            seg = y[seg_idx * scale:(seg_idx + 1) * scale]
            x = np.arange(len(seg), dtype=float)
            coef = np.polyfit(x, seg, 1)
            trend = np.polyval(coef, x)
            vals.append(np.sqrt(np.mean((seg - trend) ** 2)))
        fluctuations.append(np.nanmean(vals))
    if len(fluctuations) < 5:
        return np.nan
    xx = np.log10(scales[:len(fluctuations)])
    yy = np.log10(fluctuations)
    if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(yy)):
        return np.nan
    slope, _ = np.polyfit(xx, yy, 1)
    return float(slope)


def get_rr_times_by_subject():
    id_map = pd.read_csv(ID_MAP_PATH)
    mapping = dict(zip(id_map['id-old'], id_map['id-new']))
    rr_by_subj = {}
    for old_id in sorted(mapping):
        path = RR_DIR / f'{old_id}.txt'
        if path.exists():
            rr_times = np.loadtxt(path)
            rr_by_subj[int(mapping[old_id])] = rr_times
    return rr_by_subj


def build_window_entropy_table():
    df_hrv = pd.read_csv(HRV_PATH)
    rr_by_subj = get_rr_times_by_subject()
    rows = []
    for _, row in df_hrv.iterrows():
        subj = int(row['subject_id_new'])
        rr_times = rr_by_subj.get(subj)
        if rr_times is None:
            continue
        start = float(row['window_start'])
        end = float(row['window_end'])
        window = rr_times[(rr_times >= start) & (rr_times <= end)]
        rr = np.diff(window) * 1000.0
        if len(rr) < 20:
            continue
        for rfac in [0.10, 0.15, 0.20, 0.25]:
            tolerance = rfac * rr.std(ddof=1)
            rows.append({
                'subject_id_new': subj,
                'window_start': start,
                'window_end': end,
                'rr_count': len(rr),
                'm': 2,
                'rfactor': rfac,
                'tolerance': tolerance,
                'sampen': sample_entropy(rr, m=2, r=tolerance),
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / 'sample_entropy_sensitivity.csv', index=False)
    return out


def plot_sample_entropy_sensitivity():
    df = pd.read_csv(OUT_DIR / 'sample_entropy_sensitivity.csv')
    summary = df.groupby('rfactor')['sampen'].agg(['median', 'mean', 'std', 'min', 'max'])
    fig, ax = plt.subplots(figsize=(8, 5))
    summary['median'].plot(kind='bar', ax=ax, color='steelblue')
    ax.set_title('Sample Entropy sensitivity across r factors')
    ax.set_xlabel('r factor (× SD(RR))')
    ax.set_ylabel('Median SampEn across windows')
    ax.set_xticklabels([f'{x:g}' for x in summary.index], rotation=0)
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'sample_entropy_sensitivity.png', dpi=200)
    plt.close(fig)


def build_dfa_table():
    df_hrv = pd.read_csv(HRV_PATH)
    rr_by_subj = get_rr_times_by_subject()
    rows = []
    valid = 0
    for _, row in df_hrv.iterrows():
        subj = int(row['subject_id_new'])
        rr_times = rr_by_subj.get(subj)
        if rr_times is None:
            continue
        start = float(row['window_start'])
        end = float(row['window_end'])
        window = rr_times[(rr_times >= start) & (rr_times <= end)]
        rr = np.diff(window) * 1000.0
        if len(rr) < 30:
            continue
        alpha1 = dfa_alpha1(rr)
        rows.append({
            'subject_id_new': subj,
            'window_start': start,
            'window_end': end,
            'rr_count': len(rr),
            'dfa1': alpha1,
        })
        if np.isfinite(alpha1):
            valid += 1
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / 'dfa_alpha1_windows.csv', index=False)
    with open(OUT_DIR / 'dfa_validity_summary.json', 'w') as f:
        json.dump({'valid_alpha1': valid, 'total': len(out), 'fraction': valid / max(len(out), 1)}, f, indent=2)
    return out


def fit_mixed_models():
    df = pd.read_csv(MERGED_PATH)
    df = df[df['subject_id_new'].notna()].copy()
    outcomes = [
        'meanhr',
        'sdnn',
        'rmssd',
        'lfpowfft',
        'hfpowfft',
        'lfhffft',
        'sampen',
    ]
    results = []
    for outcome in outcomes:
        model_df = df[['subject_id_new', 'vo2_mlkgmin_median', outcome]].dropna().copy()
        if len(model_df) == 0:
            continue
        # Fit random-intercept model and quadratic check
        form1 = f'{outcome} ~ vo2_mlkgmin_median'
        form2 = f'{outcome} ~ vo2_mlkgmin_median + I(vo2_mlkgmin_median**2)'
        try:
            m1 = mixedlm(formula=form1, data=model_df, groups=model_df['subject_id_new'], re_formula='1')
            fit1 = m1.fit(method='lbfgs', maxiter=100)
            m2 = mixedlm(formula=form2, data=model_df, groups=model_df['subject_id_new'], re_formula='1')
            fit2 = m2.fit(method='lbfgs', maxiter=100)
            aic1 = fit1.aic
            aic2 = fit2.aic
            best_model = 'linear' if aic1 <= aic2 else 'quadratic'
            chosen = fit1 if aic1 <= aic2 else fit2
            coef = chosen.summary().tables[1]
            results.append({
                'outcome': outcome,
                'n_obs': len(model_df),
                'n_subjects': model_df['subject_id_new'].nunique(),
                'model_choice': best_model,
                'aic_linear': aic1,
                'aic_quadratic': aic2,
                'vo2_coef': float(chosen.params.get('vo2_mlkgmin_median', np.nan)),
                'vo2_p': float(chosen.pvalues.get('vo2_mlkgmin_median', np.nan)),
                'vo2_se': float(chosen.bse.get('vo2_mlkgmin_median', np.nan)),
                'vo2_ci_low': float(chosen.conf_int().loc['vo2_mlkgmin_median', 0]) if 'vo2_mlkgmin_median' in chosen.conf_int().index else np.nan,
                'vo2_ci_high': float(chosen.conf_int().loc['vo2_mlkgmin_median', 1]) if 'vo2_mlkgmin_median' in chosen.conf_int().index else np.nan,
            })
        except Exception as exc:  # pragma: no cover
            results.append({
                'outcome': outcome,
                'n_obs': len(model_df),
                'n_subjects': model_df['subject_id_new'].nunique(),
                'model_choice': 'failed',
                'aic_linear': np.nan,
                'aic_quadratic': np.nan,
                'vo2_coef': np.nan,
                'vo2_p': np.nan,
                'vo2_se': np.nan,
                'vo2_ci_low': np.nan,
                'vo2_ci_high': np.nan,
                'error': str(exc),
            })
    out = pd.DataFrame(results)
    out.to_csv(OUT_DIR / 'mixed_model_vo2_summary.csv', index=False)
    return out


def main():
    # Audit summary checks
    df = pd.read_csv(MERGED_PATH)
    print('Merged rows:', len(df))
    print('Subjects:', df['subject_id_new'].nunique())
    print('Obs per subject summary:')
    print(df.groupby('subject_id_new').size().describe())
    print('SampEn missing:', df['sampen'].isna().sum())
    print('DFA1 missing:', df['dfa1'].isna().sum())

    # Sample entropy sensitivity
    sample_df = build_window_entropy_table()
    print('\nSample entropy summary by r factor:')
    print(sample_df.groupby('rfactor')['sampen'].agg(['median', 'mean', 'std']).round(6))
    plot_sample_entropy_sensitivity()

    # DFA check
    dfa_df = build_dfa_table()
    print('\nDFA valid count:', int(np.isfinite(dfa_df['dfa1']).sum()), '/', len(dfa_df))
    print(dfa_df['dfa1'].describe())

    # Mixed model summary
    model_df = fit_mixed_models()
    print('\nMixed-model summary:')
    print(model_df[['outcome', 'n_obs', 'n_subjects', 'model_choice', 'vo2_coef', 'vo2_p']].to_string(index=False))


if __name__ == '__main__':
    main()
