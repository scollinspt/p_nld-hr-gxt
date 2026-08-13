#!/usr/bin/env python3
"""Reviewer-response reanalysis for Paper 1.

Creates separate outputs only. Does not alter manuscript, validated dataset, or prior results.
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
from scipy import signal
from scipy.stats import skew
from statsmodels.genmod.cov_struct import Autoregressive
from statsmodels.genmod.families import Gaussian
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.formula.api import mixedlm

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Paper1' / 'AddressingReviewConcerns'
OUT.mkdir(parents=True, exist_ok=True)
VALIDATED = ROOT / 'analysis_validation' / 'hrv_metabolic_merged_validated.csv'
HRV = ROOT / 'python_analysis' / 'hrv_python_new_advancement.csv'
MET = ROOT / 'all_bxb_10secondEpochs.csv'
ID_MAP = ROOT / 'Old-to-New-IDs.csv'
TEST = ROOT / 'subject_TestData.csv'
RR_DIR = ROOT / 'rr'
OUTCOMES = [('meanhr','raw'),('sdnn','raw'),('rmssd','log'),('lfpowfft','log'),('hfpowfft','log'),('lfhffft','log'),('sampen','log')]


def fmt(x): return 'NA' if not np.isfinite(x) else f'{x:.3f}'
def pval(x): return '<0.001' if x < .001 else f'{x:.3f}'

def load_sources():
    idmap = pd.read_csv(ID_MAP); mapping = dict(zip(idmap['id-old'], idmap['id-new']))
    met = pd.read_csv(MET); met['subject_id_new'] = met.ID.map(mapping); met = met.dropna(subset=['subject_id_new']).copy(); met['subject_id_new'] = met.subject_id_new.astype(int); met['time_s'] = met.epoch * 10.0
    hrv = pd.read_csv(HRV); hrv['window_end'] = hrv.window_start + 300.0
    validated = pd.read_csv(VALIDATED)
    return met, hrv, validated, mapping

def interval_summary(met, start, end):
    g = met[(met.time_s >= start) & (met.time_s <= end)].sort_values('time_s')
    if len(g) < 2: return None
    output = {'n_epochs':len(g)}
    for src, prefix in [('vo2_mlkgmin_median','vo2'),('hr_MEDIAN','hr'),('rr_MEDIAN','resp_rate'),('ve_stpd_MEDIAN','ve')]:
        x = pd.to_numeric(g[src], errors='coerce').to_numpy(dtype=float); t = g.time_s.to_numpy(dtype=float); good = np.isfinite(x)
        if good.sum() < 2:
            for s in ['mean','median','start','mid','end','change','pct_change','slope','sd','cv']: output[f'{prefix}_{s}']=np.nan
            continue
        x=x[good]; t=t[good]; start_val=x[0]; end_val=x[-1]; midpoint=x[np.argmin(np.abs(t-(start+end)/2))]
        output.update({f'{prefix}_mean':x.mean(),f'{prefix}_median':np.median(x),f'{prefix}_start':start_val,f'{prefix}_mid':midpoint,f'{prefix}_end':end_val,
                       f'{prefix}_change':end_val-start_val,f'{prefix}_pct_change':100*(end_val-start_val)/abs(start_val) if start_val else np.nan,
                       f'{prefix}_slope':np.polyfit(t,x,1)[0],f'{prefix}_sd':x.std(ddof=1),f'{prefix}_cv':x.std(ddof=1)/x.mean() if x.mean() else np.nan})
    return output

def nonstationarity(met, hrv):
    rows=[]
    for r in hrv.itertuples(index=False):
        g=met[met.subject_id_new == int(r.subject_id_new)]; info=interval_summary(g,r.window_start,r.window_end)
        if info: rows.append({'subject_id_new':int(r.subject_id_new),'window_start':r.window_start,'window_end':r.window_end,**info})
    d=pd.DataFrame(rows); d.to_csv(OUT/'within_window_nonstationarity.csv',index=False); return d

def join_window_exposures(validated, nonstat):
    ns=nonstat.rename(columns={'window_start':'timestamp'})
    out=validated.merge(ns,on=['subject_id_new','timestamp'],how='left',suffixes=('','_window'))
    # The legacy timestamp is the merged start-time exposure. These definitions share the actual full 300-s interval.
    out['vo2_current_start']=out.vo2_mlkgmin_median
    for suffix in ['start','mid','mean','median','end']: out[f'vo2_{suffix}']=out[f'vo2_{suffix}']
    out.to_csv(OUT/'aligned_window_exposures.csv',index=False); return out

def prep_model(d,outcome,exposure,transform):
    x=d[['subject_id_new','timestamp',outcome,exposure]].dropna().copy().sort_values(['subject_id_new','timestamp'])
    x['z']=(x[exposure]-x[exposure].mean())/x[exposure].std(ddof=1); x['y']=np.log(x[outcome]) if transform=='log' else x[outcome]
    return x

def gee_fit(d, quadratic=True):
    form='y ~ z + I(z**2)' if quadratic else 'y ~ z'
    fit=GEE.from_formula(form,groups='subject_id_new',time='timestamp',cov_struct=Autoregressive(grid=True),family=Gaussian(),data=d).fit()
    q=next((v for v in fit.params.index if v.startswith('I(z')),None); ci=fit.conf_int()
    return {'fit':fit,'beta1':fit.params['z'],'se1':fit.bse['z'],'ci1_low':ci.loc['z',0],'ci1_high':ci.loc['z',1],'p1':fit.pvalues['z'],
            'beta2':fit.params[q] if q else np.nan,'se2':fit.bse[q] if q else np.nan,'ci2_low':ci.loc[q,0] if q else np.nan,'ci2_high':ci.loc[q,1] if q else np.nan,'p2':fit.pvalues[q] if q else np.nan,'alpha':float(fit.cov_struct.dep_params)}

def alignment_models(d):
    rows=[]
    for exposure in ['vo2_current_start','vo2_mid','vo2_mean','vo2_median','vo2_end']:
        for outcome, transform in OUTCOMES:
            x=prep_model(d,outcome,exposure,transform); r=gee_fit(x,True)
            linear=gee_fit(x,False)
            rows.append({'exposure':exposure,'model_form':'quadratic','outcome':outcome,'transform':transform,'n_obs':len(x),'n_subjects':x.subject_id_new.nunique(),**{k:v for k,v in r.items() if k!='fit'}})
            rows.append({'exposure':exposure,'model_form':'linear','outcome':outcome,'transform':transform,'n_obs':len(x),'n_subjects':x.subject_id_new.nunique(),**{k:v for k,v in linear.items() if k!='fit'}})
    pd.DataFrame(rows).to_csv(OUT/'alignment_sensitivity_models.csv',index=False); return pd.DataFrame(rows)

def rr_data(mapping): return {new:np.loadtxt(RR_DIR/f'{old}.txt') for old,new in mapping.items() if (RR_DIR/f'{old}.txt').exists()}

def sampen(x,m=2,r=.2):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    if len(x)<m+2: return np.nan
    z=(x-x.mean())/x.std(ddof=1)
    if not np.all(np.isfinite(z)): return np.nan
    from scipy.spatial import cKDTree
    def nmatch(k):
        templates=np.lib.stride_tricks.sliding_window_view(z,k)
        return len(cKDTree(templates).query_pairs(r,p=np.inf))
    b,a=nmatch(m),nmatch(m+1)
    return -np.log(a/b) if a and b else np.nan

def window_metrics(rr_times,start,length):
    ts=rr_times[(rr_times>=start)&(rr_times<=start+length)]; rr=np.diff(ts)*1000; rt=ts[1:]
    if len(rr)<max(30, int(length/1.5)): return None
    diff=np.diff(rr); raw=sampen(rr); detr=sampen(signal.detrend(rr))
    result={'rr_count':len(rr),'meanhr':60000/rr.mean(),'sdnn':rr.std(ddof=1),'rmssd':np.sqrt(np.mean(diff**2)),
            'sampen_raw':raw,'sampen_detrended':detr,'rr_slope':np.polyfit(rt,rr,1)[0],'rr_change':rr[-1]-rr[0]}
    # Spectral estimates are output only at >=120 s; LF needs substantially longer data and 300 s is the defensible primary window.
    if length>=120:
        z=signal.detrend(rr-rr.mean()); t=np.cumsum(rr)/1000; grid=np.arange(0,t[-1],.25)
        y=np.interp(grid,t,z); f=np.fft.rfftfreq(len(y),.25); p=np.abs(np.fft.rfft(y*np.hanning(len(y))))**2/len(y)
        for name,lo,hi in [('lfpowfft',.04,.15),('hfpowfft',.15,.40),('hfextended',.15,.80)]: result[name]=p[(f>=lo)&(f<=hi)].sum()
        result['lfhffft']=result['lfpowfft']/result['hfpowfft'] if result['hfpowfft']>0 else np.nan
        result['lf_hf_defensible']=length>=300
    return result

def short_windows(met,mapping,analysis_subjects):
    rrmap=rr_data(mapping); rows=[]
    # advance equals one third of window length, balancing local stationarity with usable serial sensitivity.
    for length in [30,60,120,180,300]:
        advance=length/3
        for subj,times in rrmap.items():
            if subj not in analysis_subjects:
                continue
            start=0.
            while start+length<=times[-1]:
                metrics=window_metrics(times,start,length)
                if metrics:
                    info=interval_summary(met[met.subject_id_new==subj],start,start+length)
                    if info: rows.append({'window_length_s':length,'advance_s':advance,'subject_id_new':subj,'timestamp':start,**metrics,**info})
                start+=advance
    d=pd.DataFrame(rows); d.to_csv(OUT/'window_length_metrics.csv',index=False); return d

def window_models(metrics):
    rows=[]; samp=[]
    for length,g in metrics.groupby('window_length_s'):
        for outcome,transform in [('meanhr','raw'),('sdnn','raw'),('rmssd','log'),('sampen_raw','log'),('sampen_detrended','log')]:
            x=g[['subject_id_new','timestamp',outcome,'vo2_mean']].dropna().copy();
            if len(x)<30: continue
            x['z']=(x.vo2_mean-x.vo2_mean.mean())/x.vo2_mean.std(ddof=1);x['y']=np.log(x[outcome]) if transform=='log' else x[outcome]
            r=gee_fit(x,False)
            rows.append({'window_length_s':length,'outcome':outcome,'transform':transform,'n_obs':len(x),'n_subjects':x.subject_id_new.nunique(),**{k:v for k,v in r.items() if k!='fit'}})
        # SampEn raw vs detrended association and trend-adjusted model
        x=g[['subject_id_new','timestamp','sampen_raw','sampen_detrended','vo2_mean','rr_slope']].dropna().copy()
        if len(x)>30:
            x['z']=(x.vo2_mean-x.vo2_mean.mean())/x.vo2_mean.std(ddof=1);x['trend_z']=(x.rr_slope-x.rr_slope.mean())/x.rr_slope.std(ddof=1);x['y']=np.log(x.sampen_detrended)
            fit=GEE.from_formula('y ~ z + trend_z',groups='subject_id_new',time='timestamp',cov_struct=Autoregressive(grid=True),family=Gaussian(),data=x).fit(); ci=fit.conf_int()
            samp.append({'window_length_s':length,'n_obs':len(x),'raw_mean':x.sampen_raw.mean(),'detrended_mean':x.sampen_detrended.mean(),'raw_detrended_r':x[['sampen_raw','sampen_detrended']].corr().iloc[0,1],
                         'beta_vo2_detrended':fit.params['z'],'se_vo2_detrended':fit.bse['z'],'ci_low':ci.loc['z',0],'ci_high':ci.loc['z',1],'p_vo2':fit.pvalues['z'],
                         'beta_rr_trend':fit.params['trend_z'],'p_rr_trend':fit.pvalues['trend_z'],'ar1_alpha':float(fit.cov_struct.dep_params)})
    pd.DataFrame(rows).to_csv(OUT/'window_length_sensitivity_models.csv',index=False); pd.DataFrame(samp).to_csv(OUT/'sampen_detrending_sensitivity.csv',index=False)
    return pd.DataFrame(rows),pd.DataFrame(samp)

def respiratory(d):
    x=d.copy(); x['resp_hz']=x.resp_rate_mean/60; x['resp_band']=np.select([x.resp_hz<.15,x.resp_hz<=.40],["below_0.15","within_0.15_0.40"],default="above_0.40")
    # conventional power comes from validated 300-s data, aligned to full-window respiratory summaries.
    cols=['subject_id_new','timestamp','resp_hz','resp_band','vo2_mean','ve_mean','hfpowfft','lfhffft']
    x[cols].to_csv(OUT/'respiratory_frequency_hf_analysis.csv',index=False)
    return x

def mixed_and_thinning(d):
    rows=[]; thinrows=[]
    for outcome,transform in OUTCOMES:
        x=prep_model(d,outcome,'vo2_mean',transform)
        for form in ['y ~ z','y ~ z + I(z**2)']:
            try:
                ri=mixedlm(form,x,groups=x.subject_id_new).fit(reml=False,method='lbfgs'); rs=mixedlm(form,x,groups=x.subject_id_new,re_formula='1 + z').fit(reml=False,method='lbfgs')
                rows.append({'outcome':outcome,'form':form,'ri_aic':ri.aic,'rs_aic':rs.aic,'ri_converged':ri.converged,'rs_converged':rs.converged,'random_intercept_var':ri.cov_re.iloc[0,0],'random_slope_var':rs.cov_re.iloc[1,1],'slope_sd':np.sqrt(rs.cov_re.iloc[1,1])})
            except Exception as e: rows.append({'outcome':outcome,'form':form,'error':str(e)})
        for spacing in [60,150,300]:
            keep=[]
            for _,g in x.groupby('subject_id_new'):
                last=-np.inf
                for idx,r in g.iterrows():
                    if r.timestamp-last>=spacing: keep.append(idx);last=r.timestamp
            z=x.loc[keep]; r=gee_fit(z,False)
            thinrows.append({'outcome':outcome,'spacing_s':spacing,'n_obs':len(z),'n_subjects':z.subject_id_new.nunique(),**{k:v for k,v in r.items() if k!='fit'}})
    pd.DataFrame(rows).to_csv(OUT/'mixed_model_comparison_reanalysis.csv',index=False);pd.DataFrame(thinrows).to_csv(OUT/'thinning_sensitivity_full_results.csv',index=False)
    return pd.DataFrame(rows),pd.DataFrame(thinrows)

def athlete_slopes(d):
    char=pd.read_csv(TEST); rows=[]
    for outcome,transform in [('meanhr','raw'),('sdnn','raw'),('lfhffft','log'),('sampen','log')]:
        for subj,g in d.groupby('subject_id_new'):
            x=g[['vo2_mean',outcome]].dropna();
            if len(x)<5: continue
            y=np.log(x[outcome]) if transform=='log' else x[outcome]; slope=np.polyfit(x.vo2_mean,y,1)[0]
            rows.append({'subject_id_new':subj,'outcome':outcome,'slope_per_vo2':slope,'n_obs':len(x)})
    out=pd.DataFrame(rows).merge(char,left_on='subject_id_new',right_on='id',how='left');out.to_csv(OUT/'participant_slope_exploration.csv',index=False)
    associations=[]
    for outcome,g in out.groupby('outcome'):
        for variable in ['vo2_mlkgmin_max','AnThmlkgmin','hr_MAX','age','wt_kg','bmi','bodyfat']:
            x=g[['slope_per_vo2',variable]].dropna()
            if len(x)>3:
                r=x.corr().iloc[0,1]
                associations.append({'outcome':outcome,'athlete_variable':variable,'n':len(x),'pearson_r':r})
    pd.DataFrame(associations).to_csv(OUT/'participant_slope_associations.csv',index=False)
    return out

def figures(nonstat, aligned, metrics, resp):
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'pdf.fonttype':42})
    # drift
    fig,ax=plt.subplots(1,2,figsize=(9,3.5)); ax[0].scatter(nonstat.vo2_mean,nonstat.vo2_change,s=7,alpha=.25);ax[0].axhline(0,color='k');ax[0].set(xlabel='Window mean VO2 (ml kg-1 min-1)',ylabel='300-s VO2 end-start');ax[1].scatter(nonstat.vo2_mean,nonstat.hr_change,s=7,alpha=.25,color='#d95f02');ax[1].axhline(0,color='k');ax[1].set(xlabel='Window mean VO2',ylabel='300-s HR end-start (bpm)');fig.tight_layout();fig.savefig(OUT/'figure_within_window_drift.png',dpi=350);fig.savefig(OUT/'figure_within_window_drift.pdf');plt.close(fig)
    # alignment
    fig,ax=plt.subplots(figsize=(5,4));
    for c in ['vo2_current_start','vo2_mid','vo2_mean','vo2_end']: ax.scatter(aligned.vo2_mean,aligned[c]-aligned.vo2_mean,s=6,alpha=.15,label=c.replace('vo2_',''))
    ax.axhline(0,color='k');ax.set(xlabel='Window mean VO2',ylabel='Exposure minus window mean VO2');ax.legend(frameon=False);fig.tight_layout();fig.savefig(OUT/'figure_alignment_sensitivity.png',dpi=350);fig.savefig(OUT/'figure_alignment_sensitivity.pdf');plt.close(fig)
    # SampEn raw/detrended by lengths
    fig,ax=plt.subplots(1,2,figsize=(9,3.5));
    for length,g in metrics.groupby('window_length_s'):
        ax[0].scatter(g.sampen_raw,g.sampen_detrended,s=5,alpha=.1,label=f'{length}s')
        h=g[['vo2_mean','sampen_detrended']].dropna(); coef=np.polyfit(h.vo2_mean,np.log(h.sampen_detrended),1);grid=np.linspace(h.vo2_mean.min(),h.vo2_mean.max(),100);ax[1].plot(grid,np.exp(np.polyval(coef,grid)),label=f'{length}s')
    ax[0].set(xlabel='Raw SampEn',ylabel='Linearly detrended SampEn');ax[1].set(xlabel='Window mean VO2',ylabel='Detrended SampEn');[a.legend(frameon=False) for a in ax];fig.tight_layout();fig.savefig(OUT/'figure_sampen_detrending_windows.png',dpi=350);fig.savefig(OUT/'figure_sampen_detrending_windows.pdf');plt.close(fig)
    # respiration
    fig,ax=plt.subplots(1,2,figsize=(9,3.5));ax[0].scatter(resp.vo2_mean,resp.resp_hz,s=7,alpha=.2);ax[0].axhline(.40,color='crimson',ls='--');ax[0].set(xlabel='Window mean VO2',ylabel='Respiratory frequency (Hz)');ax[1].scatter(resp.resp_hz,resp.hfpowfft,s=7,alpha=.2);ax[1].axvline(.40,color='crimson',ls='--');ax[1].set(xlabel='Respiratory frequency (Hz)',ylabel='Conventional HF power (ms2)');fig.tight_layout();fig.savefig(OUT/'figure_respiratory_frequency_hf.png',dpi=350);fig.savefig(OUT/'figure_respiratory_frequency_hf.pdf');plt.close(fig)
    # All participants are shown; no post-hoc representative-subject selection.
    fig,axes=plt.subplots(2,2,figsize=(8.5,6.5),constrained_layout=True)
    for ax,outcome in zip(axes.ravel(),['meanhr','sdnn','lfhffft','sampen']):
        for _,g in aligned.groupby('subject_id_new'):
            g=g.sort_values('vo2_mean');ax.plot(g.vo2_mean,g[outcome],lw=.6,alpha=.45)
        ax.set(xlabel='Window mean VO2 (ml kg-1 min-1)',ylabel=outcome)
    fig.savefig(OUT/'figure_participant_specific_trajectories.png',dpi=350);fig.savefig(OUT/'figure_participant_specific_trajectories.pdf');plt.close(fig)

def summary(nonstat,align,wm,samp,resp,mixed,thin,slopes):
    q=nonstat.vo2_change.quantile([.25,.5,.75]); drift=(nonstat.vo2_change.abs()>=5).mean()
    resp_above=(resp.resp_hz>.4).mean(); cross=resp[resp.resp_hz>.4].groupby('subject_id_new').vo2_mean.min()
    text=f'''# Reviewer Reanalysis Summary\n\n## Scope\nThis separate reanalysis preserves the validated dataset and prior outputs. It tests window nonstationarity, full-window metabolic alignment, shorter windows, entropy detrending, respiration, and participant-aware inference.\n\n## VMax inventory\nThe processed VMax file contains 10-s epoch medians for VO2, VCO2, VE (STPD/BTPS), respiratory rate, tidal volume, RQ, VE/VO2, VE/VCO2, oxygen pulse, end-tidal gases, speed, and grade. No adjudicated ventilatory/anaerobic-threshold time marker was found.\n\n## Within-window nonstationarity\nAcross 300-s windows, VO2 end-minus-start median was {q.loc[.5]:.2f} ml kg-1 min-1 (IQR {q.loc[.25]:.2f} to {q.loc[.75]:.2f}); mean {nonstat.vo2_change.mean():.2f} (SD {nonstat.vo2_change.std():.2f}). {drift:.1%} of windows changed by at least 5 ml kg-1 min-1. This quantifies material drift rather than local stationarity for many segments.\n\n## Alignment\nStart, midpoint, mean, median, and end VO2 models were separately refit; inspect `alignment_sensitivity_models.csv` before retaining any claim. Whole-window mean/median VO2 are the scientifically preferable exposures for a window summary.\n\n## Short windows and SampEn\nShort-window time-domain and SampEn analyses are in `window_length_sensitivity_models.csv` and `sampen_detrending_sensitivity.csv`. Frequency-domain LF/HF is flagged as defensible only for the 300-s window because low-frequency estimation requires longer segments. Interpret 30--180-s spectral measures as exploratory only.\n\n## Respiration and HF\n{resp_above:.1%} of 300-s windows had mean respiratory frequency above 0.40 Hz. The median participant-specific VO2 at first crossing was {cross.median():.2f} ml kg-1 min-1 (N={cross.notna().sum()}). This means the fixed 0.15--0.40 Hz HF band may omit respiratory-linked variability at higher intensity; conventional HF and LF/HF cannot be interpreted without this limitation.\n\n## Inference\nGEE fits use participant clusters and grid-estimated AR(1) working correlation with robust sandwich covariance provided by statsmodels GEE. With 22 clusters, asymptotic sandwich standard errors may still be optimistic; mixed models and thinning are included as sensitivity analyses, not confirmation from 923 independent units.\n\n## Scientific decision\n1. A 300-s window is not uniformly locally stationary during the GXT; drift magnitude must be reported and whole-window exposures should be preferred.\n2. Alignment sensitivity results, not prior start-time-only results, determine whether conclusions survive.\n3. SampEn claims are conditional on detrended and shorter-window results; do not call it physiological complexity unless both are robust.\n4. Respiratory frequency above 0.40 Hz limits conventional HF/LF-HF interpretation at high intensity.\n5. The final Paper 1 claims should be limited to findings stable across whole-window alignment, detrending, window length, and participant-aware sensitivity models.\n'''
    (OUT/'reviewer_reanalysis_summary.md').write_text(text)
    (OUT/'gee_small_cluster_diagnostics.md').write_text('# GEE Small-Cluster Diagnostics\n\nThe implementation uses `statsmodels` Gaussian GEE with participant as the group, `Autoregressive(grid=True)` working correlation, and the package default robust sandwich covariance. The robust covariance is asymptotic in the number of clusters; N=22 is modest, so p-values and confidence intervals should be interpreted cautiously. No finite-cluster correction was applied because a validated correction is not directly implemented in this workflow. Mixed-effects and 60/150/300-s thinning analyses are reported as robustness checks.\n')

def vmax_inventory():
    text='''# VMax Data Inventory\n\nSource: `all_bxb_10secondEpochs.csv`, derived from per-second files in `VMax_bxb_csv/`. Time is `epoch * 10` seconds; 1,912 mapped epochs cover 23 mapped participants before final cohort restriction.\n\n| Construct | Processed VMax field | Units / meaning |\n|---|---|---|\n| VO2 absolute | `vo2_L_median` | L/min |\n| VO2 relative | `vo2_mlkgmin_median` | ml kg-1 min-1 |\n| VCO2 | `vco2_L_MEDIAN` | L/min |\n| VE | `ve_stpd_MEDIAN`, `ve_btps_MEDIAN` | L/min |\n| Respiratory frequency | `rr_MEDIAN` | breaths/min |\n| Tidal volume | `vt_L_MEDIAN` | L |\n| RER | `rq_MEDIAN` | ratio |\n| VE/VO2 | `veo2_MEDIAN` | ratio |\n| VE/VCO2 | `veco2_MEDIAN` | ratio |\n| Oxygen pulse | `o2pulse_MEDIAN` | ml/beat |\n| End-tidal gases | `peto2_MEDIAN`, `petco2_MEDIAN` | recorded VMax values |\n| Workload | `speed_MEDIAN`, `grade_MEDIAN` | recorded VMax values |\n\nNo adjudicated ventilatory-threshold or anaerobic-threshold time marker was found in these VMax time-series files. `subject_TestData.csv` contains a participant-level `AnThmlkgmin` value, but it is not a time-resolved threshold marker and was not used to define exposures or zones.\n'''
    (OUT/'vmax_data_inventory.md').write_text(text)

def main():
    met,hrv,validated,mapping=load_sources(); analysis_subjects=set(validated.subject_id_new.astype(int).unique()); hrv=hrv[hrv.subject_id_new.isin(analysis_subjects)].copy(); nonstat=nonstationarity(met,hrv); aligned=join_window_exposures(validated,nonstat); align=alignment_models(aligned); metrics=short_windows(met,mapping,analysis_subjects); wm,samp=window_models(metrics); resp=respiratory(aligned); mixed,thin=mixed_and_thinning(aligned); slopes=athlete_slopes(aligned); figures(nonstat,aligned,metrics,resp);summary(nonstat,align,wm,samp,resp,mixed,thin,slopes);vmax_inventory()
    print('Completed reviewer reanalysis:',OUT)
    print('300-s windows',len(nonstat),'mean VO2 change',nonstat.vo2_change.mean(),'resp >0.4Hz', (resp.resp_hz>.4).mean())
if __name__=='__main__': main()
