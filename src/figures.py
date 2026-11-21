"""Make figures and sensitivity tables from the validation output."""
from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"me637-matplotlib"))
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from project_paths import FIGURES, INPUTS, TABLES
from validate import MODELS, HOURS, REGIONS, metrics, solar_geometry

COLORS = {"PJM":"#202930","ERA5 wind":"#19768a","HRRR-MET wind":"#ab501a","ERA5 solar":"#19768a"}
STYLES = {"PJM": "-", "ERA5 wind": "--", "ERA5 solar": "--", "HRRR-MET wind": "-."}
WIND_LOSSES = (.05, .15, .25)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.labelcolor":"#202930","xtick.color":"#35424b","ytick.color":"#35424b",
                     "legend.frameon":False,"savefig.facecolor":"white","pdf.fonttype":42})


def save_figure(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=320, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", dpi=320, bbox_inches="tight", metadata={"CreationDate": None})
    plt.close(fig)


def figure_map():
    f = pd.read_csv(INPUTS/"fleet.csv")
    z = gpd.read_file(INPUTS/"pjm_zone.geojson").to_crs(4326)
    fig, axes = plt.subplots(1,2,figsize=(7.5,3.0),layout="constrained")
    assignments = [("ICARUS polygon", "Zone polygon", "#257182", "o"),
                   ("EIA T&D owner", "Transmission owner", "#b35b27", "^"),
                   ("unassigned", "Unassigned region", "#684780", "s")]
    for resource,ax in zip(["wind","solar"],axes):
        z.plot(ax=ax,color="#e8eeee",edgecolor="#9aa8ae",linewidth=0.35,rasterized=True)
        group = f[(f.resource==resource)&f.lat.between(34.25,43)&f.lon.between(-91,-71.25)].groupby("eia_id").agg(
            lat=("lat","first"),lon=("lon","first"),ac_mw=("ac_mw","sum"),method=("region_source","first"))
        for method, _, color, marker in assignments:
            rows=group[group.method==method]
            ax.scatter(rows.lon,rows.lat,s=3+0.15*rows.ac_mw,c=color,marker=marker,alpha=.80,linewidths=0)
        ax.set_title("(a) Wind plants" if resource=="wind" else "(b) Solar PV plants",loc="left",fontweight="bold")
        ax.set(xlim=(-91,-73),ylim=(34.5,43))
        ax.set_xticks([-90,-85,-80,-75], ["90° W","85° W","80° W","75° W"])
        ax.set_yticks([36,39,42], ["36° N","39° N","42° N"])
        ax.tick_params(labelsize=8, length=2)
        ax.set_aspect(1/np.cos(np.radians(39)))
        for spine in ax.spines.values():spine.set_visible(False)
    handles = [Line2D([], [], marker=m, ls="", color=c, markersize=5, label=l)
               for _, l, c, m in assignments]
    sizes = [Line2D([], [], marker="o", ls="", color="#58646c", markersize=np.sqrt(3+.15*mw),
                    markeredgewidth=0, label=f"{mw} MW") for mw in (10,100,300)]
    fig.legend(handles=handles+sizes, loc="outside lower center", ncol=6,
               fontsize=8, handletextpad=.3, columnspacing=.8)
    save_figure(fig, "plant_locations")


def figure_seasonal(m):
    fig, axes=plt.subplots(2,2,figsize=(7.5,5.3),layout="constrained")
    for i,(tech,mods) in enumerate([("Wind",MODELS[:2]),("Solar PV",MODELS[2:])]):
        for j,split in enumerate(["month","hour_ept"]):
            ax=axes[i,j]
            d=m[(m.region=="RTO")&(m.split==split)&m.model.isin(mods)].copy()
            d["x"]=d.value.astype(int)
            for key in ["PJM",*mods]:
                batch=d[d.model==(mods[0] if key=="PJM" else key)].sort_values("x")
                y=batch.observed_mean_mw if key=="PJM" else batch.modeled_mean_mw
                ax.plot(batch.x,y/1000,label=key.replace(" wind","").replace(" solar",""),lw=1.8,
                        ls=STYLES[key],color=COLORS[key])
            panel = "abcd"[2*i+j]
            ax.set_title(f"({panel}) {tech} by "+("month" if j==0 else "hour"),loc="left",fontweight="bold",fontsize=10)
            ax.set_ylabel("Mean generation (GW)");ax.set_ylim(bottom=-.15);ax.grid(axis="y",alpha=.2)
            if j==0:ax.set_xticks([1,3,5,7,9,11],["Jan","Mar","May","Jul","Sep","Nov"])
            else:ax.set_xticks([0,4,8,12,16,20,23]);ax.set_xlabel("Hour beginning (Eastern prevailing time)")
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="outside lower center")
    save_figure(fig, "seasonal_diurnal")



def example_week(data):
    import matplotlib.dates as md
    from zoneinfo import ZoneInfo
    tz=ZoneInfo("America/New_York")
    start=pd.Timestamp("2025-03-03",tz=tz)
    end=pd.Timestamp("2025-03-10",tz=tz)
    week=data[data.region.eq("RTO") & data.utc.between(start,end,inclusive="left")].copy()
    fig,axes=plt.subplots(2,1,figsize=(7.4,3.6),sharex=True,layout="constrained")
    for ax,mods in zip(axes,[MODELS[:2],MODELS[2:]]):
        for key in ["PJM",*mods]:
            frame=week[week.model.eq(mods[0] if key=="PJM" else key)].sort_values("utc")
            power=frame.observed_mw if key=="PJM" else frame.modeled_mw
            ax.plot(frame.utc,power/1000,label=key.split()[0],color=COLORS[key],linewidth=1.25,
                    linestyle=STYLES[key])
        ax.set_title("(a) Wind" if len(mods)==2 else "(b) Solar PV",loc="left",fontweight="bold",fontsize=10)
        ax.set_ylabel("Wind (GW)" if len(mods)==2 else "Solar (GW)")
        ax.set_ylim(bottom=-.1);ax.grid(axis="y",alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3)
    axes[1].xaxis.set_major_locator(md.DayLocator(tz=tz))
    axes[1].xaxis.set_major_formatter(md.DateFormatter("%b %d",tz=tz))
    axes[1].set_xlim(start,end)
    axes[1].set_xlabel("Date in 2025 (Eastern prevailing time)")
    save_figure(fig, "hourly_example")


def figure_wind_losses(data):
    wind = data[data.model.isin(MODELS[:2])]
    difference = pd.DataFrame(index=REGIONS, columns=WIND_LOSSES, dtype=float)
    for loss in WIND_LOSSES:
        residual = wind.modeled_mw * (1 - loss) / .85 - wind.observed_mw
        rmse = residual.pow(2).groupby([wind.region, wind.model]).mean().pow(.5).unstack()
        difference[loss] = rmse[MODELS[0]] - rmse[MODELS[1]]

    fig, ax = plt.subplots(figsize=(6.2, 2.8), layout="constrained")
    limit = np.abs(difference.to_numpy()).max()
    cells = ax.imshow(difference, cmap="PuOr", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(WIND_LOSSES)), [f"{loss*100:.0f} %" for loss in WIND_LOSSES])
    ax.set_yticks(range(len(REGIONS)), REGIONS)
    ax.set_xlabel("Assumed wind plant loss")
    ax.tick_params(length=0)
    for row, region in enumerate(REGIONS):
        for col, loss in enumerate(WIND_LOSSES):
            value = difference.loc[region, loss]
            rgb = np.array(cells.cmap(cells.norm(value))[:3])
            linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
            color = "white" if linear @ [.2126, .7152, .0722] < .18 else "black"
            ax.text(col, row, f"{value:+.0f}".replace("-", "−"), ha="center", va="center", color=color)
    bar = fig.colorbar(cells, ax=ax, pad=.04)
    bar.set_label("RMSE difference (MW)")
    ax.set_title("ERA5 RMSE minus HRRR-MET RMSE", loc="left", fontsize=10, fontweight="bold")
    save_figure(fig, "wind_loss_comparison")

def diagnostics(d):
    d=d[d.region=="RTO"].copy()
    solar=d[d.model==MODELS[2]].set_index("utc")
    extra=pd.read_csv(TABLES/"hourly_fleet_diagnostics.csv.gz",parse_dates=["utc"]).set_index("utc")
    sensitivity=[]
    for label,field in [("exclude net-metered generators","solar_net_metered_mw"),("USGS-located generators only","solar_eia_location_mw")]:
        alt=solar.copy()
        alt.modeled_mw-=extra[field]
        sensitivity.append({"model":"ERA5 solar","scenario":label,**metrics(alt)})
    paired=[]
    for model in MODELS[:2]:
        g=d[d.model==model].sort_values("utc").copy()
        paired.append((g.modeled_mw-g.observed_mw).to_numpy()**2)
        baseline=g.modeled_mw.copy()
        for loss in WIND_LOSSES:
            g.modeled_mw=baseline*(1-loss)/.85
            sensitivity.append({"model":model,"scenario":f"{loss:.0%} wind loss",**metrics(g)})
    pd.DataFrame(sensitivity).to_csv(TABLES/"sensitivity.csv",index=False)
    # Paired circular moving-block bootstrap preserves synoptic-scale serial dependence.
    n=len(paired[0]); block=168; rng=np.random.default_rng(637)
    def segments(errors):
        v=np.concatenate([errors,errors[:block]])
        cum=np.r_[0,np.cumsum(v)]
        return cum[block:]-cum[:-block]
    sums=np.stack([segments(x)[:n] for x in paired])
    nblock,tail=divmod(n,block)
    starts=rng.integers(0,n,size=(2000,nblock))
    sq=sums[:,starts].sum(axis=2)
    if tail:
        starts=rng.integers(0,n,size=2000)
        for i,e in enumerate(paired):
            p=np.r_[0,np.cumsum(np.r_[e,e[:tail]])]
            sq[i]+=p[starts+tail]-p[starts]
    delta=np.sqrt(sq[0]/n)-np.sqrt(sq[1]/n)
    delta_point=np.sqrt(np.mean(paired[0]))-np.sqrt(np.mean(paired[1]))
    pd.DataFrame([{"era5_minus_hrrr_rmse_mw":delta_point,"p025_mw":np.quantile(delta,.025),
                   "p975_mw":np.quantile(delta,.975),"replicates":len(delta),"block_hours":block}]).to_csv(TABLES/"paired_wind_rmse.csv",index=False)
    rows=[]
    for model,g in d.groupby("model"):
        g=g.sort_values("utc")
        if model.endswith("solar"):
            zen,_,_=solar_geometry(39.5,-79)
            g=g[g.utc.isin(HOURS[zen<80])]
        q=g.observed_mw.quantile(.1)
        low=g[g.observed_mw<=q]
        rows.append(dict(model=model,observed_p10_mw=q,**metrics(low)))
    pd.DataFrame(rows).to_csv(TABLES/"low_output.csv",index=False)


if __name__=="__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    df=pd.read_csv(TABLES/"hourly_validation.csv.gz",parse_dates=["utc"])
    diagnostics(df)
    example_week(df)
    figure_wind_losses(df)
    figure_map()
    figure_seasonal(pd.read_csv(TABLES/"metrics.csv"))
