import calendar
import datetime
import pathlib

import pandas as pd
import numpy as np

from data2dd_funcs import (euro_demand2dd, scen2dd, temporal2dd, trans_links,
                           co2target2dd, add_vre_connection_costs)
                           

root = pathlib.Path(snakemake.output[0]).parent
data_root = root
out = pathlib.Path(".")


co2target2dd(
    snakemake.input.co2_targets_db,
    snakemake.output.co2_target,
    snakemake.params.co2_target_scenario,
    snakemake.params.co2_target_type,
    snakemake.params.co2_target_extent)

psys_scen = snakemake.wildcards.psys_scenario
esys_scen = "BASE"

scen_db = snakemake.input[1]
f_techno = snakemake.input[2]

params_to_write = {}

params_to_write["gen"] = {}
params_to_write["gen"]["set"] = [
    "g",
    "non_vre",
    "vre",
    "hydro_res",
    "uc_int",
    "uc_lin",
    "quick"
]
params_to_write["gen"]["parameter"] = {}
params_to_write["gen"]["parameter"]["all"] = [
    "emis fac",
    "max ramp",
    "min down",
    "min up",
    "startup cost",
    "inertia",
    "unit size",
    "capex2050",
    "varom",
    "fom",
    "fuelcost2050",
    "cap2area",
    "af",
    "min gen",
    "peak af",
]

params_to_write["store"] = {}
params_to_write["store"]["set"] = ["s", "uc_lin","quick"]
params_to_write["store"]["parameter"] = {}
params_to_write["store"]["parameter"]["all"] = [
    "max_freq",
    "max_res",
    "eff_in",
    "eff_out",
    "loss_per_hr",
    "p_to_e",
    "varom",
    "e capex2050",
    "p capex2050",
    "fom",
    "af",
    "min down",
    "min up",
    "startup cost",
    "inertia",
    "unit size",
    "min gen",
    "max ramp",
]


if snakemake.wildcards.spatial == "focus":
    
    zones=np.array(snakemake.params.aggregated_regions)
    for (key,z) in snakemake.params.focus_countries.items():
        zones=np.append(zones[zones!=key],z)
        
    if snakemake.params.aggregated_countries is not None:
        for key,z in snakemake.params.aggregated_countries.items():
            zones=np.append(zones[np.isin(zones,z,invert=True)],key)
            
    agg_countries=snakemake.params.aggregated_countries
    
    zones=zones.tolist()
    
else:
    zones=snakemake.params.aggregated_regions
    agg_countries=None
    
    
scen2dd(
    snakemake.output[1],
    root,
    f_techno,
    params_to_write,
    scen_db,
    psys_scen,
    esys_scen,
    zones,
    out=out,
    esys_cap=False,
    exist_cap=True,
    exist_agg=snakemake.wildcards.spatial,
    agg_countries=agg_countries
)

trans_links(
    root,
    f_techno,
    aggregated_regions=zones,
    out=out,
    agg_countries=agg_countries
)

add_vre_connection_costs(
    root,
    out,
    f_techno,
    psys_scen,
    snakemake.input.vre_connection_dists
    )
    

years = [snakemake.wildcards.year]
date_range = snakemake.params.date_range

if snakemake.wildcards.spatial == "focus":
    focus_dem_shares=snakemake.input.focus_dem_shares
else:
    focus_dem_shares=None

for yr in years:
    date_range = [yr + "-" + date for date in date_range]
    yr = int(yr)

    dstart = datetime.datetime.fromisoformat(date_range[0])
    dstop = datetime.datetime.fromisoformat(date_range[1]) + datetime.timedelta(
        hours=23
    )

    if calendar.isleap(yr):
        rleap = False
    else:
        rleap = True

    temporal2dd(dstart, dstop, root / out, snakemake.output[0])

    euro_demand2dd(
        snakemake.input["europedemandcsvlocation"],
        snakemake.params.aggregated_regions,
        root,
        root / out,
        dstart,
        dstop,
        scen_db,
        esys_scen,
        yr,
        focus_dem_shares=focus_dem_shares,
        focus_countries=snakemake.params.focus_countries,
        agg_countries=agg_countries
    )
