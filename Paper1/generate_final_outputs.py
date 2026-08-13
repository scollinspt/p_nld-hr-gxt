#!/usr/bin/env python3
"""Generate final Paper 1 tables and figures from validated analysis inputs.

This script does not modify validated analysis files or manuscript prose.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from statsmodels.genmod.cov_struct import Autoregressive
from statsmodels.genmod.families import Gaussian
from statsmodels.genmod.generalized_estimating_equations import GEE

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'Paper1' / 'final_outputs'
OUT.mkdir(parents=True, exist_ok=True)
VALIDATED = ROOT / 'analysis_validation' / 'hrv_metabolic_merged_validated.csv'
TEST = ROOT / 'subject_TestData.csv'
ID_MAP = ROOT / 'Old-to-New-IDs.csv'
ENTROPY_ALL = ROOT / 'analysis_validation' / 'sample_entropy_sensitivity_all_windows.csv'
ENTROPY_SUMMARY = ROOT / 'analysis_validation' / 'sample_entropy_sensitivity_summary.csv'
ENTROPY_CORR = ROOT / 'analysis_validation' / 'sample_entropy_sensitivity_correlations.csv'
ENTROPY_VO2 = ROOT / 'analysis_validation' / 'sample_entropy_sensitivity_vo2_models.csv'
THINNING = ROOT / 'analysis_validation' / 'thinning_sensitivity_models.csv'
ML_MODELS = ROOT / 'analysis_validation' / 'model_comparisons_ml.csv'

OUTCOMES = [
    ('meanhr', 'Mean HR', 'bpm', 'raw'),
    ('sdnn', 'SDNN', 'ms', 'raw'),
    ('rmssd', 'RMSSD', 'ms', 'log'),
    ('lfpowfft', 'LF power', 'ms$^2$', 'log'),
    ('hfpowfft', 'HF power', 'ms$^2$', 'log'),
    ('lfhffft', 'LF/HF', 'ratio', 'log'),
    ('sampen', 'Sample Entropy', 'unitless', 'log'),
]
OUTCOME_MAP = {x[0]: x for x in OUTCOMES}


def pvalue(value):
    if not np.isfinite(value):
        return 'NA'
    return '<0.001' if value < 0.001 else f'{value:.3f}'


def fmt(value, digits=3):
    return 'NA' if not np.isfinite(value) else f'{value:.{digits}f}'


def latex_escape(text):
    return str(text).replace('&', r'\&').replace('%', r'\%').replace('_', r'\_')


def table_document(body: str, label: str, caption: str, note: str = '') -> str:
    note_block = f'\\begin{{tablenotes}}\\footnotesize\\item {note}\\end{{tablenotes}}\n' if note else ''
    notes = f'\\begin{{threeparttable}}\n{body}{note_block}\\end{{threeparttable}}'
    return '\\begin{table}[htbp]\n\\centering\n' + notes + f'\n\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{table}}\n'


def participant_table():
    validated = pd.read_csv(VALIDATED)
    ids = sorted(validated['subject_id_new'].dropna().astype(int).unique())
    test = pd.read_csv(TEST)
    participants = test[test['id'].isin(ids)].copy()
    if len(participants) != len(ids):
        raise ValueError(f'Participant reconciliation failed: expected {len(ids)}, found {len(participants)}')
    variables = [
        ('age', 'Age (yr)'), ('ht_m', 'Height (m)'), ('wt_kg', 'Body mass (kg)'),
        ('bmi', r'BMI (kg/m$^2$)'), ('bodyfat', 'Body fat (\%)'),
        ('vo2_mlkgmin_max', r'VO$_2$max (ml$\cdot$kg$^{-1}\cdot$min$^{-1}$)'),
        ('hr_MAX', 'Maximal HR (beats/min)'),
    ]
    lines = [r'\begin{tabular}{lcc}', r'\toprule', r'Characteristic & Mean $\pm$ SD & Range \\', r'\midrule']
    lines.append(f'N & {len(participants)} & -- \\\\')
    for col, label in variables:
        x = pd.to_numeric(participants[col], errors='coerce').dropna()
        lines.append(f'{label} & {x.mean():.2f} $\\pm$ {x.std(ddof=1):.2f} & {x.min():.2f}--{x.max():.2f} \\\\')
    # AT VO2 is sufficiently complete and explicitly named in source data.
    x = pd.to_numeric(participants['AnThmlkgmin'], errors='coerce').dropna()
    if len(x) >= 0.8 * len(participants):
            lines.append(r'Anaerobic-threshold VO$_2$ (ml$\cdot$kg$^{-1}\cdot$min$^{-1}$) & '
                     f'{x.mean():.2f} $\\pm$ {x.std(ddof=1):.2f} & {x.min():.2f}--{x.max():.2f} \\\\')
    lines += [r'\bottomrule', r'\end{tabular}']
    note = f'Values are mean $\\pm$ SD with observed range. N = {len(participants)} analyzed participants.'
    (OUT / 'table_participant_characteristics.tex').write_text(
        table_document('\n'.join(lines), 'tab:participant-characteristics', 'Participant characteristics', note), encoding='utf-8')
    return participants


def prepare(df, outcome, transform):
    d = df[['subject_id_new', 'timestamp', 'vo2_mlkgmin_median', outcome]].dropna().sort_values(['subject_id_new', 'timestamp']).copy()
    mean = d['vo2_mlkgmin_median'].mean()
    sd = d['vo2_mlkgmin_median'].std(ddof=1)
    d['vo2_z'] = (d['vo2_mlkgmin_median'] - mean) / sd
    d['y'] = np.log(d[outcome]) if transform == 'log' else d[outcome]
    return d, mean, sd


def fit_final_models(df):
    fits = {}
    rows = []
    for col, label, unit, transform in OUTCOMES:
        d, vo2_mean, vo2_sd = prepare(df, col, transform)
        formula = 'y ~ vo2_z + I(vo2_z**2)'
        gee = GEE.from_formula(formula, groups='subject_id_new', time='timestamp', cov_struct=Autoregressive(grid=True), family=Gaussian(), data=d).fit()
        # Formula naming is stable in statsmodels, but locate the quadratic term defensively.
        quad_name = next(name for name in gee.params.index if name.startswith('I(vo2_z'))
        ci = gee.conf_int()
        fits[col] = {'fit': gee, 'data': d, 'vo2_mean': vo2_mean, 'vo2_sd': vo2_sd, 'quad_name': quad_name}
        rows.append({
            'outcome': col, 'label': label, 'scale': transform, 'unit': unit,
            'beta1': gee.params['vo2_z'], 'se1': gee.bse['vo2_z'], 'ci1_low': ci.loc['vo2_z', 0], 'ci1_high': ci.loc['vo2_z', 1], 'p1': gee.pvalues['vo2_z'],
            'beta2': gee.params[quad_name], 'se2': gee.bse[quad_name], 'ci2_low': ci.loc[quad_name, 0], 'ci2_high': ci.loc[quad_name, 1], 'p2': gee.pvalues[quad_name],
            'ar1_alpha': float(gee.cov_struct.dep_params), 'n_obs': len(d), 'n_subjects': d.subject_id_new.nunique(),
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / 'table_gee_primary_results.csv', index=False)
    return fits, result


def gee_table(results):
    lines = [r'\begin{tabular}{lcccccccccc}', r'\toprule',
             r'Outcome & Scale & $\beta_1$ & SE$_1$ & 95\% CI$_1$ & $p_1$ & $\beta_2$ & SE$_2$ & 95\% CI$_2$ & $p_2$ & AR(1) $\alpha$ \\', r'\midrule']
    for row in results.itertuples(index=False):
        scale = 'Raw' if row.scale == 'raw' else 'Natural log'
        ci1 = f'[{fmt(row.ci1_low)}, {fmt(row.ci1_high)}]'
        ci2 = f'[{fmt(row.ci2_low)}, {fmt(row.ci2_high)}]'
        lines.append(f'{row.label} & {scale} & {fmt(row.beta1)} & {fmt(row.se1)} & {ci1} & {pvalue(row.p1)} & {fmt(row.beta2)} & {fmt(row.se2)} & {ci2} & {pvalue(row.p2)} & {fmt(row.ar1_alpha)} \\\\')
    lines += [r'\bottomrule', r'\end{tabular}']
    note = (r'VO$_2$ was standardized before fitting. $\beta_1$ is the first-order coefficient evaluated at centered VO$_2$; '
            r'$\beta_2$ represents curvature. RMSSD, LF power, HF power, LF/HF, and Sample Entropy were natural-log transformed. '
            r'N = 923 observations from 22 participants.')
    (OUT / 'table_gee_primary_results.tex').write_text(table_document('\n'.join(lines), 'tab:gee-results', 'Primary quadratic AR(1)-GEE results', note), encoding='utf-8')


def thinning_table():
    d = pd.read_csv(THINNING)
    lines = [r'\begin{tabular}{lccc}', r'\toprule', r'Outcome & 60 s (N=163) & 150 s (N=72) & 300 s (N=45) \\', r'\midrule']
    for col, label, unit, transform in OUTCOMES:
        x = d[d.outcome == col].set_index('spacing_seconds')
        cells = []
        for spacing in [60, 150, 300]:
            row = x.loc[spacing]
            cells.append(f'{fmt(row.beta_vo2_z)} [{fmt(row.ci_low)}, {fmt(row.ci_high)}]')
        lines.append(f'{label} & ' + ' & '.join(cells) + r' \\')
    lines += [r'\bottomrule', r'\end{tabular}']
    note = r'Entries are VO$_2$ coefficients per one SD increase with 95\% CI. Thinning was a sensitivity analysis; the complete AR(1)-GEE used all 923 observations.'
    (OUT / 'table_s1_thinning_sensitivity.tex').write_text(table_document('\n'.join(lines), 'tab:s1-thinning', 'Temporal-thinning sensitivity', note), encoding='utf-8')


def model_selection_table():
    d = pd.read_csv(ML_MODELS)
    selected = d[d['transform'].isin(['raw', 'log'])].copy()
    lines = [r'\begin{tabular}{lrrrrrrc}', r'\toprule', r'Outcome & Linear AIC & Quadratic AIC & $\Delta$AIC & Linear BIC & Quadratic BIC & LR $\chi^2$ & Converged \\', r'\midrule']
    for col, label, unit, transform in OUTCOMES:
        x = selected[(selected['outcome'] == col) & (selected['transform'] == transform)].set_index('model')
        lin = x.loc['linear_random_intercept']; quad = x.loc['quadratic_random_intercept']
        lines.append(f'{label} & {lin.aic:.1f} & {quad.aic:.1f} & {lin.aic - quad.aic:.1f} & {lin.bic:.1f} & {quad.bic:.1f} & {quad.lr_vs_linear:.1f} & {str(bool(lin.converged and quad.converged)).lower()} \\\\')
    lines += [r'\bottomrule', r'\end{tabular}']
    note = r'Models were fit by maximum likelihood. The likelihood-ratio statistic compares quadratic with linear random-intercept models.'
    (OUT / 'table_s2_model_selection.tex').write_text(table_document('\n'.join(lines), 'tab:s2-model-selection', 'Maximum-likelihood model selection', note), encoding='utf-8')


def entropy_table():
    summary = pd.read_csv(ENTROPY_SUMMARY).set_index('r_factor')
    assoc = pd.read_csv(ENTROPY_VO2).set_index('r_factor')
    lines = [r'\begin{tabular}{lrrrrrrr}', r'\toprule', r'$r$ & N usable & Median & IQR & Mean & SD & Range & VO$_2$ association (95\% CI; $p$) \\', r'\midrule']
    for rf, row in summary.iterrows():
        a = assoc.loc[rf]
        association = f'{a.beta_vo2:.3f} [{a.ci_low:.3f}, {a.ci_high:.3f}]; {pvalue(a.p_value)}'
        lines.append(f'{rf:.2f} & {int(row["n_usable"])} & {row["median"]:.3f} & {row["iqr"]:.3f} & {row["mean"]:.3f} & {row["sd"]:.3f} & {row["min"]:.3f}--{row["max"]:.3f} & {association} \\\\')
    lines += [r'\bottomrule', r'\end{tabular}']
    note = r'All 930 HRV windows produced defined estimates. VO$_2$ associations are from random-intercept sensitivity models using the corrected Sample Entropy values.'
    (OUT / 'table_s3_sampen_sensitivity.tex').write_text(table_document('\n'.join(lines), 'tab:s3-sampen', 'Sample Entropy sensitivity analysis', note), encoding='utf-8')


def prediction_grid(fit_info, outcome, x):
    fit = fit_info['fit']; d = fit_info['data']; mean = fit_info['vo2_mean']; sd = fit_info['vo2_sd']
    z = (x - mean) / sd
    exog = pd.DataFrame({'Intercept': 1.0, 'vo2_z': z, 'I(vo2_z ** 2)': z ** 2})
    # Align columns to statsmodels' parameter order.
    exog = exog[fit.params.index]
    pred = exog.to_numpy() @ fit.params.to_numpy()
    cov = fit.cov_params().loc[fit.params.index, fit.params.index].to_numpy()
    se = np.sqrt(np.einsum('ij,jk,ik->i', exog.to_numpy(), cov, exog.to_numpy()))
    lo = pred - 1.96 * se; hi = pred + 1.96 * se
    if OUTCOME_MAP[outcome][3] == 'log':
        return np.exp(pred), np.exp(lo), np.exp(hi)
    return pred, lo, hi


def styling():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})


def primary_figure(df, fits):
    styling(); fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.8), constrained_layout=True); axes = axes.ravel()
    xmin, xmax = df.vo2_mlkgmin_median.min(), df.vo2_mlkgmin_median.max(); grid = np.linspace(xmin, xmax, 250)
    for i, (col, label, unit, transform) in enumerate(OUTCOMES):
        ax = axes[i]; d = fits[col]['data']; original = d[col]
        ax.scatter(d.vo2_mlkgmin_median, original, s=9, alpha=.16, color='#34495e', linewidths=0, rasterized=True)
        pred, lo, hi = prediction_grid(fits[col], col, grid)
        ax.fill_between(grid, lo, hi, color='#d95f02', alpha=.20, linewidth=0)
        ax.plot(grid, pred, color='#b33f00', linewidth=2.0)
        ax.set_xlabel(r'VO$_2$ (ml kg$^{-1}$ min$^{-1}$)')
        ax.set_ylabel(f'{label} ({unit})' if unit != 'unitless' else label)
        ax.text(.03, .95, chr(65+i), transform=ax.transAxes, va='top', fontweight='bold')
        ax.spines[['top','right']].set_visible(False); ax.grid(alpha=.18)
    axes[-1].axis('off')
    fig.savefig(OUT / 'figure_primary_trajectories.png', dpi=400); fig.savefig(OUT / 'figure_primary_trajectories.pdf'); plt.close(fig)


def individual_figure(df, fits):
    styling(); selected = ['meanhr', 'sdnn', 'lfhffft', 'sampen']; fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(.08, .92, df.subject_id_new.nunique()))
    for ax, col in zip(axes.ravel(), selected):
        label = OUTCOME_MAP[col][1]; unit = OUTCOME_MAP[col][2]
        for color, (subject, g) in zip(colors, df.groupby('subject_id_new')):
            g = g.sort_values('vo2_mlkgmin_median')
            ax.plot(g.vo2_mlkgmin_median, g[col], color=color, alpha=.42, linewidth=.65)
        grid = np.linspace(df.vo2_mlkgmin_median.min(), df.vo2_mlkgmin_median.max(), 250)
        pred, lo, hi = prediction_grid(fits[col], col, grid)
        ax.fill_between(grid, lo, hi, color='#111111', alpha=.12)
        ax.plot(grid, pred, color='#111111', linewidth=2.2)
        ax.set_xlabel(r'VO$_2$ (ml kg$^{-1}$ min$^{-1}$)'); ax.set_ylabel(f'{label} ({unit})' if unit != 'unitless' else label)
        ax.spines[['top','right']].set_visible(False); ax.grid(alpha=.18)
    handles = [Line2D([0], [0], color='#111111', lw=2.2, label='Population-average fit'), Line2D([0], [0], color='#6a51a3', lw=.8, alpha=.6, label='Participant trajectories')]
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False)
    fig.savefig(OUT / 'figure_individual_trajectories.png', dpi=400); fig.savefig(OUT / 'figure_individual_trajectories.pdf'); plt.close(fig)


def entropy_figure(df):
    styling(); all_e = pd.read_csv(ENTROPY_ALL); summary = pd.read_csv(ENTROPY_SUMMARY); vo2 = df[['subject_id_new','timestamp','vo2_mlkgmin_median']]
    all_e = all_e.rename(columns={'window_start':'timestamp'}).merge(vo2, on=['subject_id_new','timestamp'], how='inner')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    labels = [f'r = {x:.2f}' for x in sorted(all_e.r_factor.unique())]
    axes[0].boxplot([all_e.loc[all_e.r_factor == x, 'sampen'] for x in sorted(all_e.r_factor.unique())], labels=labels, patch_artist=True,
                    boxprops={'facecolor':'#9ecae1','alpha':.75}, medianprops={'color':'#08519c'})
    axes[0].set_xlabel('Tolerance factor'); axes[0].set_ylabel('Sample Entropy'); axes[0].spines[['top','right']].set_visible(False); axes[0].grid(axis='y', alpha=.18)
    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a']
    for color, rf in zip(colors, sorted(all_e.r_factor.unique())):
        g = all_e[all_e.r_factor == rf]
        axes[1].scatter(g.vo2_mlkgmin_median, g.sampen, s=7, alpha=.16, color=color, rasterized=True)
        # Descriptive line is intentionally not the inferential GEE fit; this panel shows robustness.
        coef = np.polyfit(g.vo2_mlkgmin_median, g.sampen, 1)
        grid = np.linspace(g.vo2_mlkgmin_median.min(), g.vo2_mlkgmin_median.max(), 100)
        axes[1].plot(grid, np.polyval(coef, grid), color=color, lw=1.8, label=f'r = {rf:.2f}')
    axes[1].set_xlabel(r'VO$_2$ (ml kg$^{-1}$ min$^{-1}$)'); axes[1].set_ylabel('Sample Entropy'); axes[1].spines[['top','right']].set_visible(False); axes[1].grid(alpha=.18); axes[1].legend(frameon=False)
    fig.savefig(OUT / 'figure_sampen_sensitivity.png', dpi=400); fig.savefig(OUT / 'figure_sampen_sensitivity.pdf'); plt.close(fig)


def main():
    df = pd.read_csv(VALIDATED)
    if len(df) != 923 or df.subject_id_new.nunique() != 22:
        raise ValueError('Validated dataset does not match 923 observations and 22 subjects')
    if df['sampen'].isna().any() or df['dfa1'].notna().any() or df['dfa2'].notna().any():
        raise ValueError('Validated dataset failed Sample Entropy/DFA quality gate')
    participants = participant_table()
    fits, results = fit_final_models(df)
    gee_table(results); thinning_table(); model_selection_table(); entropy_table()
    primary_figure(df, fits); individual_figure(df, fits); entropy_figure(df)
    print(f'Generated outputs in {OUT}')
    print('Participants:', len(participants), 'observations:', len(df), 'subjects:', df.subject_id_new.nunique())
    print(results[['outcome','beta1','beta2','se1','se2','ar1_alpha']].round(6).to_string(index=False))


if __name__ == '__main__':
    main()
