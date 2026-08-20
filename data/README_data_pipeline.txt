ECB SPF data for Chapter 4. Already cleaned; ready to use.

- "yhat.csv": forecasts. "forERR.csv": forecast errors. "ytrue.csv": true
  series. Read directly by empirical_spf_gdp.R - no further processing.
- Dimensions: T = 107 quarters (1999Q3-2026Q1), p = 35 forecasters.

To update or reproduce the data, follow these steps. This extends Lee & Seregina
(2025) (GitHub: ekat92/RD-FGL).

1. Download quarterly forecasts from the ECB SPF archive (one zip, all
   rounds, updated regularly):
   https://www.ecb.europa.eu/stats/prices/indic/forecast/shared/files/SPF_individual_forecasts.zip
   RD-FGL's repo already contains the older rounds - only add the csv
   files for newer rounds (e.g. from 2024Q2.csv
   onwards) into RD-FGL's DATA/ECB SPF data/Get Updated Data/New folder,
   alongside its existing files.
2. Run "Excel.py" then "Excel2.py" (Lee & Seregina's scripts, in the same
   folder) to merge the rounds.
3. Output: "second_edit.csv", all forecasts collected together
   (forecasters in rows, rounds in columns).
4. Download updated GDP actuals from Eurostat, dataset namq_10_gdp, with
   filters geo=EA20, na_item=B1GQ (GDP at market prices), unit=CLV_PCH_SM
   (chain-linked volume, % change on same quarter previous year),
   s_adj=SCA (seasonally and calendar adjusted), via their data browser.
   Save as "Actual_gdpgrowth_data_upd.xlsx".
5. Run "spf_clean.R" on both files. It applies a 40% missingness filter,
   imputes remaining gaps with missForest, and computes forecast error. 
   Produces "yhat.csv", "forERR.csv", "ytrue.csv".

Verification: checked against Lee & Seregina's own forERR.csv on the
p=35 forecasters common to both panels. Mean correlation 0.9965, mean
absolute difference 0.001. Differences are due to routine ECB data
revisions, not a processing error.
