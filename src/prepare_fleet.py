"""Rebuild the small PJM generator inventory from EIA-860, then attach USGS locations."""
import argparse
from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import numpy as np
import pandas as pd

from project_paths import INPUTS, RAW_DATA
from wind_inventory import read_wind_turbines

MIDATL = set("BGE METED PENELEC PPL PECO PEPCO PSEG RECO DPL JCPL AECO".split())
WEST = set("DUQ APS OVEC COMED DAY DEOK ATSI EKPC AEP".split())
SOUTH = {"DOM", "MECK"}


def extract_eia(archive):
    with ZipFile(archive) as z:
        plants = pd.read_excel(z.open("2___Plant_Y2025.xlsx"), skiprows=1)
        plants = plants.loc[plants["Balancing Authority Code"].eq("PJM"),
                            ["Plant Code", "Latitude", "Longitude", "Plant Name", "State", "Transmission or Distribution System Owner"]]
        picked = []
        for resource, name in [("wind", "3_2_Wind"), ("solar", "3_3_Solar")]:
            book = pd.ExcelFile(z.open(f"{name}_Y2025.xlsx"))
            for sheet in book.sheet_names:
                df = pd.read_excel(book, sheet, skiprows=1).merge(plants, on=["Plant Code", "State", "Plant Name"])
                eligible = df.Status.eq("OP") if sheet == "Operable" else df.Status.eq("RE") & df["Retirement Year"].eq(2025)
                df = df[eligible & df.Technology.isin(["Solar Photovoltaic", "Onshore Wind Turbine", "Offshore Wind Turbine"])].copy()
                df["resource"] = resource
                picked.append(df)
        df = pd.concat(picked)
        fields = {"Transmission or Distribution System Owner":"td_owner", "Plant Code":"eia_id", "Generator ID":"generator", "Plant Name":"plant", "State":"state",
                  "Operating Year":"start_year", "Operating Month":"start_month", "Retirement Year":"end_year",
                  "Retirement Month":"end_month", "Nameplate Capacity (MW)":"ac_mw", "Latitude":"eia_lat",
                  "Longitude":"eia_lon", "DC Net Capacity (MW)":"dc_mw", "Turbine Hub Height (Feet)":"hub_ft",
                  "Fixed Tilt?":"fixed", "Single-Axis Tracking?":"single", "Dual-Axis Tracking?":"dual",
                  "Tilt Angle":"tilt", "Azimuth Angle":"azimuth", "Net Metering Agreement?":"net_metered"}
        df.rename(columns=fields)[["resource", *fields.values()]].to_csv(INPUTS / "eia_pjm_2025.csv", index=False)


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def prepare():
    f = pd.read_csv(INPUTS / "eia_pjm_2025.csv", dtype={"generator": str})
    assert not f.duplicated(["resource", "eia_id", "generator"]).any()
    for k in ["eia_lat","eia_lon","ac_mw","dc_mw","hub_ft","azimuth","tilt","start_year","start_month","end_year","end_month"]:
        f[k] = numeric(f[k])
    w = read_wind_turbines()
    w = w[w.eia_id.isin(f.eia_id) & w.p_year.le(2025)].copy()
    w["weight"] = w.t_cap.where(w.t_cap.gt(0), w.groupby("eia_id").t_cap.transform("median")).fillna(1)
    stats = []
    for code, g in w.groupby("eia_id"):
        hh = g.t_hh.between(10, 250)
        stats.append(dict(eia_id=code, usgs_lat=np.average(g.ylat,weights=g.weight),
                          usgs_lon=np.average(g.xlong,weights=g.weight),
                          usgs_hub=np.average(g.t_hh[hh],weights=g.weight[hh]) if hh.any() else np.nan))
    pv = pd.read_csv(next((RAW_DATA / "usgs/uspvdbCSV").glob("*.csv")))
    pv = pv[pv.eia_id.isin(f.eia_id) & pv.p_year.le(2025)].copy()
    assert not pv.eia_id.duplicated().any()
    pv = pv.rename(columns={"ylat":"usgs_lat","xlong":"usgs_lon"})
    pv["usgs_ratio"] = pv.p_cap_dc / pv.p_cap_ac
    props = pd.concat([pd.DataFrame(stats).assign(resource="wind"),
                       pv[["eia_id","usgs_lat","usgs_lon","p_axis","p_azimuth","p_tilt","usgs_ratio"]].assign(resource="solar")])
    f = f.merge(props, on=["eia_id","resource"], how="left", validate="many_to_one")
    f["coord_source"] = np.where(f.usgs_lat.notna(), "USGS", "EIA")
    f["lat"] = f.usgs_lat.fillna(f.eia_lat)
    f["lon"] = f.usgs_lon.fillna(f.eia_lon)
    # Prefer generator-specific geometry, then the matched USGS plant's attributes.
    wind = f.resource.eq("wind")
    height = (0.3048*f.hub_ft).where(lambda x:x.between(10,250)).fillna(f.usgs_hub)
    f["height_m"] = height.fillna(80).where(wind)
    solar = ~wind
    f["axis"] = np.select([f.single.eq("Y"), f.dual.eq("Y"), f.fixed.eq("Y"),
                           f.p_axis.str.contains("single",na=False), f.p_axis.eq("dual-axis")],
                          ["single","dual","fixed","single","dual"], default="fixed")
    f["tilt_deg"] = f.tilt.where(f.tilt.between(0,90)).fillna(f.p_tilt.where(f.p_tilt.between(0,90))).fillna(25).where(solar)
    f["azimuth_deg"] = f.azimuth.where(f.azimuth.between(0,360)).fillna(f.p_azimuth.where(f.p_azimuth.between(0,360))).fillna(180).where(solar)
    f["dc_model_mw"] = f.dc_mw.where(f.dc_mw.gt(0)).fillna(f.ac_mw*f.usgs_ratio.where(f.usgs_ratio.between(0.5,3))).fillna(f.ac_mw*1.25).where(solar)
    zones = gpd.read_file(INPUTS / "pjm_zone.geojson").to_crs(4326)
    assert set(zones.LOC_NAME) == MIDATL | WEST | SOUTH
    pts = gpd.GeoDataFrame(f, geometry=gpd.points_from_xy(f.lon,f.lat), crs=4326)
    joined = gpd.sjoin(pts, zones[["geometry","LOC_NAME"]], how="left", predicate="intersects")
    unique = joined.groupby(level=0).LOC_NAME.agg(lambda g: g.dropna().iloc[0] if g.nunique()==1 else "UNASSIGNED")
    f["polygon_zone"] = unique
    owner_zones = {"Commonwealth Edison Co": "COMED", "Indiana Michigan Power Co": "AEP",
                   "Virginia Electric & Power Co": "DOM", "The Potomac Edison Company": "APS",
                   "Potomac Electric Power Co": "PEPCO", "Delmarva Power": "DPL"}
    fallback = f.td_owner.map(owner_zones)
    f["zone"] = unique.where(unique.ne("UNASSIGNED"), fallback.fillna("UNASSIGNED"))
    f["region_source"] = np.select([unique.ne("UNASSIGNED"), fallback.notna()],
                                   ["ICARUS polygon", "EIA T&D owner"], default="unassigned")
    f["region"] = np.select([f.zone.isin(MIDATL), f.zone.isin(WEST), f.zone.isin(SOUTH)],["MIDATL","WEST","SOUTH"],default="UNASSIGNED")
    assert f[["eia_id","ac_mw","lat","lon","start_year","start_month"]].notna().all().all()
    assert f.ac_mw.gt(0).all() and f.start_month.between(1,12).all()
    f.to_csv(INPUTS/"fleet.csv",index=False)
    print(f.groupby(["resource","region"]).ac_mw.agg(["count","sum"]))
    print(f.groupby(["resource","coord_source"]).ac_mw.agg(["count","sum"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eia-zip",type=Path,help="Optional official 2025 EIA-860 archive to refresh the included PJM subset")
    args = p.parse_args()
    if args.eia_zip:
        extract_eia(args.eia_zip)
    prepare()
