# VMax Data Inventory

Source: `all_bxb_10secondEpochs.csv`, derived from per-second files in `VMax_bxb_csv/`. Time is `epoch * 10` seconds; 1,912 mapped epochs cover 23 mapped participants before final cohort restriction.

| Construct | Processed VMax field | Units / meaning |
|---|---|---|
| VO2 absolute | `vo2_L_median` | L/min |
| VO2 relative | `vo2_mlkgmin_median` | ml kg-1 min-1 |
| VCO2 | `vco2_L_MEDIAN` | L/min |
| VE | `ve_stpd_MEDIAN`, `ve_btps_MEDIAN` | L/min |
| Respiratory frequency | `rr_MEDIAN` | breaths/min |
| Tidal volume | `vt_L_MEDIAN` | L |
| RER | `rq_MEDIAN` | ratio |
| VE/VO2 | `veo2_MEDIAN` | ratio |
| VE/VCO2 | `veco2_MEDIAN` | ratio |
| Oxygen pulse | `o2pulse_MEDIAN` | ml/beat |
| End-tidal gases | `peto2_MEDIAN`, `petco2_MEDIAN` | recorded VMax values |
| Workload | `speed_MEDIAN`, `grade_MEDIAN` | recorded VMax values |

No adjudicated ventilatory-threshold or anaerobic-threshold time marker was found in these VMax time-series files. `subject_TestData.csv` contains a participant-level `AnThmlkgmin` value, but it is not a time-resolved threshold marker and was not used to define exposures or zones.
