"""Check monthly energy for the modeled plant IDs against the final EIA-923 plant survey."""
import argparse
from calendar import month_name
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from project_paths import INPUTS, TABLES


def extract(archive):
    modeled = pd.read_csv(TABLES/"fleet_audit.csv")[["eia_id","resource"]].drop_duplicates()
    with ZipFile(archive) as z:
        book = next(n for n in z.namelist() if "Schedules_2_3_4_5" in n)
        d = pd.read_excel(z.open(book),sheet_name="Page 1 Generation and Fuel Data",skiprows=5)
    d.columns = [" ".join(c.split()) for c in d.columns]
    d = d.rename(columns={"Plant Id":"eia_id"})
    d["resource"] = d["Reported Fuel Type Code"].map({"WND":"wind","SUN":"solar"})
    d = d.merge(modeled,on=["eia_id","resource"])
    d = d[d["Reported Prime Mover"].isin(["WT","WS","PV"])].copy()
    months = [f"Netgen {month_name[i]}" for i in range(1,13)]
    mwh = d[months].apply(pd.to_numeric,errors="coerce")
    annual = pd.to_numeric(d["Net Generation (Megawatthours)"],errors="raise")
    assert np.allclose(mwh.sum(axis=1),annual,atol=1.1)
    data = pd.concat([d[["eia_id","resource","Respondent Frequency"]],mwh],axis=1)
    data.to_csv(INPUTS/"eia923_pjm_2025.csv",index=False)


def compare():
    survey = pd.read_csv(INPUTS/"eia923_pjm_2025.csv")
    fleet = pd.read_csv(TABLES/"fleet_audit.csv")[["eia_id","resource","ac_mw"]]
    matches = survey[["eia_id","resource"]].drop_duplicates()
    print("Inventory without a matching survey record")
    print(fleet.merge(matches,on=["eia_id","resource"],how="left",indicator=True).query("_merge=='left_only'").groupby("resource").agg(plants=("eia_id","nunique"),mw=("ac_mw","sum")))
    results = pd.read_csv(TABLES/"hourly_validation.csv.gz",parse_dates=["utc"])
    results = results[results.region.eq("RTO")].copy()
    results["month"] = results.utc.dt.tz_convert("America/New_York").dt.month
    plants = pd.read_csv(TABLES/"modeled_plant_monthly.csv.gz")
    plants["resource"] = plants.model.str.split().str[-1]
    paired = plants.merge(matches, on=["eia_id","resource"]).groupby(["model","month"]).mwh.sum()
    rows = []
    for (model,month), g in results.groupby(["model","month"]):
        r = survey[survey.resource.eq(model.split()[-1])]
        rows.append({"model":model,"month":month,"paired_hours":len(g),
                     "eia923_mwh_full_month":r[f"Netgen {month_name[month]}"].sum(min_count=1),
                     "model_mwh_eia_matched_plants":paired.loc[(model,month)],
                     "pjm_mwh_paired":g.observed_mw.sum(),"model_mwh_paired":g.modeled_mw.sum()})
    df = pd.DataFrame(rows)
    df.to_csv(TABLES/"monthly_energy_check.csv",index=False)
    print("January through November (TWh)")
    print(df[df.month<12].groupby("model")[["model_mwh_paired","model_mwh_eia_matched_plants","eia923_mwh_full_month","pjm_mwh_paired"]].sum()/1e6)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eia-zip",type=Path)
    args=p.parse_args()
    if args.eia_zip:
        extract(args.eia_zip)
    compare()
