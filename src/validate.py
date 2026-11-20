"""Convert 2025 weather to hourly power and compare on common PJM hour-beginning UTC intervals."""

import numpy as np
import pandas as pd
import pvlib
import xarray as xr

from project_paths import INPUTS, RAW_DATA, TABLES
from wind_inventory import read_wind_turbines

YEAR = pd.date_range("2025-01-01","2026-01-01",freq="h",inclusive="left",tz="UTC")
# Radiation is accumulated over the preceding hour. Wind uses both bounding snapshots.
HOURS = YEAR[:-1][YEAR[:-1].tz_convert("America/New_York").year == 2025]
EPT = HOURS.tz_convert("America/New_York")
REGIONS = ["RTO","WEST","MIDATL","SOUTH"]
GROUPS = REGIONS + ["UNASSIGNED"]
MODELS = ["ERA5 wind","HRRR-MET wind","ERA5 solar"]
DIMS = {"winds":{"u100":"m s**-1","v100":"m s**-1"},
        "temps":{"t2m":"K"},"solar":{"ssrd":"J m**-2","fdir":"J m**-2"}}


def observations():
    parts = []
    for resource in ("wind","solar"):
        d = pd.read_csv(RAW_DATA/"pjm"/f"{resource}_gen.csv")
        d["utc"] = pd.to_datetime(d.datetime_beginning_utc,format="%m/%d/%Y %I:%M:%S %p",utc=True)
        d["ept"] = pd.to_datetime(d.datetime_beginning_ept,format="%m/%d/%Y %I:%M:%S %p")
        assert (d.utc.dt.tz_convert("America/New_York").dt.tz_localize(None)==d.ept).all()
        wide = d.pivot(index="utc",columns="area",values=f"{resource}_generation_mw").sort_index()
        assert len(wide)==8760 and len(wide.columns)==6 and wide.notna().all().all()
        assert (wide.index.to_series().diff().dropna()==pd.Timedelta(hours=1)).all()
        audit = pd.DataFrame({"resource":resource,"area":wide.columns,"hours":len(wide),
                              "annual_mwh":wide.sum().values,"negative_hours":wide.lt(0).sum().values,
                              "min_mw":wide.min().values,"max_mw":wide.max().values})
        parts.append((wide,audit))
    pd.concat([a for _,a in parts]).to_csv(TABLES/"observed_area_audit.csv",index=False)
    return dict(zip(["wind","solar"],[w for w,_ in parts]))


def era_at_fleet(fleet):
    first = RAW_DATA/"era5/winds/2025/era5_winds_202501.nc"
    with xr.open_dataset(first) as d:
        lat, lon = d.latitude.values, d.longitude.values
    inside = fleet.lat.between(min(lat),max(lat)) & fleet.lon.between(min(lon),max(lon))
    fleet = fleet.loc[inside].copy()
    row = abs(fleet.lat.to_numpy()[:,None]-lat).argmin(axis=1)
    col = abs(fleet.lon.to_numpy()[:,None]-lon).argmin(axis=1)
    pairs, fleet["cell"] = np.unique(np.column_stack((row,col)),axis=0,return_inverse=True)
    latlon = np.column_stack((lat[pairs[:,0]],lon[pairs[:,1]]))
    weather = {}
    for kind, fields in DIMS.items():
        blocks = {key:[] for key in fields}
        for month in range(1,13):
            file = RAW_DATA/f"era5/{kind}/2025/era5_{kind}_2025{month:02}.nc"
            with xr.open_dataset(file) as d:
                assert np.array_equal(d.latitude,lat) and np.array_equal(d.longitude,lon)
                assert np.array_equal(d.valid_time.values,YEAR[YEAR.month==month].tz_localize(None).values)
                for key, unit in fields.items():
                    assert d[key].attrs["units"]==unit
                    blocks[key].append(d[key].transpose("valid_time","latitude","longitude").values[:,pairs[:,0],pairs[:,1]])
        for key in fields:
            weather[key] = np.concatenate(blocks[key])
            assert np.isfinite(weather[key]).all()
    return fleet, latlon, weather


def wind_curves():
    table = pd.read_csv(INPUTS/"wind_curves.csv")
    curves = {name:(g.speed_ms.to_numpy(),g.power_fraction.to_numpy()) for name,g in table.groupby("usgs_model")}
    w = read_wind_turbines()
    w = w[w.p_year.le(2025) & w.eia_id.notna() & w.t_cap.gt(0)].copy()
    w["direct"] = w.t_model.isin(curves)
    w["curve"] = w.t_model.where(w.direct,"GE1.5-77")
    # For unmatched machines use the GE reference at comparable rotor power loading.
    ratio = (1500 / 77**2) / (w.t_cap / w.t_rd**2)
    w["speed_scale"] = np.where(w.direct,1,ratio.where(w.t_rd.gt(0),1).fillna(1)**(1/3))
    grouped = w.groupby(["eia_id","curve","speed_scale","direct"],as_index=False).t_cap.sum()
    grouped["weight"] = grouped.t_cap/grouped.groupby("eia_id").t_cap.transform("sum")
    mix = {k:list(zip(g.curve,g.speed_scale,g.weight,g.direct)) for k,g in grouped.groupby("eia_id")}
    return curves,mix


def wind_power(speed, height, mix, curves):
    hub = speed*(height/100)**(1/7)
    cf = np.zeros(len(speed))
    for name, scale, weight, _ in mix:
        xs, ys = curves[name]
        cf += weight*np.interp(hub*scale,xs,ys,left=0,right=0)
    # Endpoints in UTC, then a trapezoidal hourly mean. Same losses for both products.
    return 0.85*(cf[:-1]+cf[1:])/2


def solar_geometry(lat,lon):
    mids = HOURS+pd.Timedelta(minutes=30)
    pos = pvlib.solarposition.get_solarposition(mids,lat,lon)
    zen, az = pos.apparent_zenith.to_numpy(),pos.azimuth.to_numpy()
    track = pvlib.tracking.singleaxis(zen,az,axis_tilt=0,axis_azimuth=0,max_angle=60,backtrack=True,gcr=0.35)
    return zen,az,track


def pv_power(g, geometry, weather):
    zen, az, track = geometry
    ghi, beam_horizontal, ambient = weather
    if g.axis=="single":
        tilt, face = np.nan_to_num(track["surface_tilt"]),np.nan_to_num(track["surface_azimuth"])
    elif g.axis=="dual":
        tilt, face = np.minimum(zen,90),az
    else:
        tilt, face = g.tilt_deg,g.azimuth_deg
    cosz = np.maximum(np.cos(np.radians(zen)),0)
    dni = np.clip(beam_horizontal/np.maximum(cosz,0.065),0,1367)
    dhi = np.maximum(ghi-dni*cosz,0)
    poa = pvlib.irradiance.get_total_irradiance(tilt,face,zen,az,dni,ghi,dhi,albedo=0.2)["poa_global"]
    poa = np.where(zen<90,np.maximum(poa,0),0)
    # Simple nominal operating cell temperature model with NOCT 45 C, gamma -0.4%/C.
    cell = ambient+25/800*poa
    dc = g.dc_model_mw*poa/1000*np.maximum(0,1-0.004*(cell-25))
    return np.clip(dc*0.90*0.96,0,g.ac_mw)


def run():
    TABLES.mkdir(parents=True, exist_ok=True)
    obs = observations()
    fleet = pd.read_csv(INPUTS/"fleet.csv")
    used, coords, w = era_at_fleet(fleet)
    curve, mixes = wind_curves()
    start = YEAR.get_loc(HOURS[0])
    hrrr = np.load(INPUTS/"hrrr_2025.npz")
    assert pd.DatetimeIndex(pd.to_datetime(hrrr["time_utc"],utc=True)).equals(YEAR)
    speed_100m = hrrr["speed_100m"][start:].astype(float)/hrrr["scale"]
    hwind = {k:speed_100m[:,i] for i,k in enumerate(hrrr["plant_id"])}
    era = np.hypot(w["u100"][start:], w["v100"][start:])
    sums = {m:pd.DataFrame(0.,index=HOURS,columns=GROUPS) for m in MODELS}
    sensitivity = pd.DataFrame(0.,index=HOURS,columns=["solar_net_metered_mw","solar_eia_location_mw"])
    capacities = pd.DataFrame(0.,index=HOURS,columns=["wind","solar"])
    audits = []
    plant_energy = []
    month = EPT.month.to_numpy()
    ymonth = EPT.year*12+EPT.month
    for cell, group in used.groupby("cell"):
        if group.resource.eq("solar").any():
            geo = solar_geometry(*coords[cell])
            irradiance = np.maximum(w["ssrd"][start+1:,cell]/3600,0)
            beam = np.clip(w["fdir"][start+1:,cell]/3600,0,irradiance)
            temp = (w["t2m"][start:-1,cell]+w["t2m"][start+1:,cell])/2-273.15
        for g in group.itertuples():
            # Monthly service dates: include all hours of both entry and retirement months.
            live = np.asarray(ymonth>=g.start_year*12+g.start_month)
            if pd.notna(g.end_year):
                live &= np.asarray(ymonth<=g.end_year*12+g.end_month)
            capacities[g.resource] += live*g.ac_mw
            match_share = np.nan
            if g.resource=="wind":
                mix = mixes.get(g.eia_id,[("GE1.5-77",1,1,False)])
                match_share = sum(weight for _,_,weight,direct in mix if direct)
                profiles = {m:g.ac_mw*wind_power(s,g.height_m,mix,curve)
                            for m,s in zip(MODELS[:2],(era[:,cell],hwind[g.eia_id]))}
            else:
                profiles = {"ERA5 solar":pv_power(g,geo,(irradiance,beam,temp))}
            for model, profile in profiles.items():
                power = np.asarray(profile)*live
                assert np.isfinite(power).all() and power.min()>=0 and power.max()<=g.ac_mw+1e-6
                monthly_energy = np.bincount(month, weights=power, minlength=13)[1:]
                plant_energy.extend((model,g.eia_id,i,mwh) for i,mwh in enumerate(monthly_energy,1))
                sums[model]["RTO"] += power
                sums[model][g.region] += power
                if model=="ERA5 solar":
                    if g.net_metered=="Y":
                        sensitivity["solar_net_metered_mw"] += power
                    if g.coord_source=="EIA":
                        sensitivity["solar_eia_location_mw"] += power
            audits.append([g.resource,g.eia_id,g.generator,g.region,g.ac_mw,live.sum(),match_share])
    pd.DataFrame(audits,columns=["resource","eia_id","generator","region","ac_mw","active_hours","wind_curve_direct_share"]).to_csv(TABLES/"fleet_audit.csv",index=False)
    plant_table = pd.DataFrame(plant_energy,columns=["model","eia_id","month","mwh"])
    plant_table.groupby(["model","eia_id","month"],as_index=False).mwh.sum().to_csv(TABLES/"modeled_plant_monthly.csv.gz",index=False,float_format="%.6f")
    output = []
    for model, totals in sums.items():
        assert np.allclose(totals.RTO,totals[GROUPS[1:]].sum(axis=1))
        resource = model.split()[-1]
        for region in REGIONS:
            frame = pd.DataFrame({"observed_mw":obs[resource][region].loc[HOURS],"modeled_mw":totals[region],
                                  "model":model,"region":region})
            output.append(frame)
    detail = pd.concat(output).rename_axis("utc").reset_index()
    detail.to_csv(TABLES/"hourly_validation.csv.gz",index=False,float_format="%.5f")
    for model in MODELS:
        sensitivity[model+"_unassigned_mw"] = sums[model].UNASSIGNED
    sensitivity.join(capacities.add_suffix("_active_mw")).rename_axis("utc").to_csv(TABLES/"hourly_fleet_diagnostics.csv.gz",float_format="%.5f")
    return detail


def metrics(g):
    d = g.modeled_mw-g.observed_mw
    return dict(hours=len(d),observed_mean_mw=g.observed_mw.mean(),modeled_mean_mw=g.modeled_mw.mean(),
                bias_mw=d.mean(),rmse_mw=np.sqrt(np.mean(d**2)),mae_mw=d.abs().mean(),
                bias_pct=100*d.sum()/g.observed_mw.sum() if g.observed_mw.sum()>0 else np.nan,
                pearson_r=g.modeled_mw.corr(g.observed_mw) if g.modeled_mw.nunique()>1 and g.observed_mw.nunique()>1 else np.nan)


def summarize(detail):
    times = detail.utc.dt.tz_convert("America/New_York")
    detail["season"] = np.array(["DJF","MAM","JJA","SON"])[(times.dt.month%12)//3]
    detail["hour_ept"] = times.dt.hour
    detail["month"] = times.dt.month
    zen,_,_ = solar_geometry(39.5,-79)
    sunlight = set(HOURS[zen<80])
    metrics_rows = []
    for (model,region),g in detail.groupby(["model","region"]):
        splits = [("all","all",g)]
        for key in ["season","hour_ept","month"]:
            splits.extend((key,k,sub) for k,sub in g.groupby(key))
        if model.endswith("solar"):
            splits.append(("sun_elevation","above_10_deg",g[g.utc.isin(sunlight)]))
        for split,label,sub in splits:
            metrics_rows.append({"model":model,"region":region,"split":split,"value":label,**metrics(sub)})
    table = pd.DataFrame(metrics_rows)
    table.to_csv(TABLES/"metrics.csv",index=False,float_format="%.6f")
    print(table[(table.split=="all")&(table.region=="RTO")].to_string(index=False))


if __name__ == "__main__":
    summarize(run())
