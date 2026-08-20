import calendar
import datetime
import pathlib

import pandas as pd

from data2dd_funcs import (euro_demand2dd, scen2dd, temporal2dd, trans_links,
                           co2target2dd, add_vre_connection_costs)
                           

root = pathlib.Path(snakemake.output[0]).parent
data_root = root
out = pathlib.Path(".")

model_year = int(snakemake.wildcards.model_year)
co2target2dd(
    snakemake.input.co2_targets_db,
    snakemake.output.co2_target,
    snakemake.wildcards.esys_scenario, # snakemake.params.co2_target_scenario
    snakemake.params.co2_target_type,
    snakemake.params.co2_target_extent,
    model_year)

psys_scen = snakemake.wildcards.psys_scenario
esys_scen = snakemake.wildcards.esys_scenario

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
    f"capex{model_year}",
    "varom",
    "fom",
    f"fuelcost{model_year}",
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
    f"e capex{model_year}",
    f"p capex{model_year}",
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

zones = pd.read_csv(snakemake.input[0]).loc[:, "zone"]

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
    model_year=model_year,
    model_years=snakemake.params.model_years,
)

trans_links(
    root,
    f_techno,
    aggregated_regions=snakemake.params.aggregated_regions,
    out=out,
    model_year=model_year,
)

add_vre_connection_costs(
    root,
    out,
    f_techno,
    psys_scen,
    snakemake.input.vre_connection_dists,
    model_year
    )
    

years = [str(model_year)]
date_range = snakemake.params.date_range

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

    print(rleap)
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
        model_year,
    )
