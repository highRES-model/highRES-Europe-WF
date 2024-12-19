import pandas as pd
import pathlib
import sqlite3
import atlite

onshore_turbine = snakemake.params.windturbines.get('onshore')
offshore_bottom_turbine = snakemake.params.windturbines.get('offshore_bottom')
onshore_size = atlite.resource.windturbine_rated_capacity_per_unit(onshore_turbine)
offshore_size = atlite.resource.windturbine_rated_capacity_per_unit(offshore_bottom_turbine)

turbinedict = {
    'onshore' : onshore_size,
    'offshore' : offshore_size,
}

df_ef = pd.read_csv(snakemake.input.employmentfactors)

path2 = pathlib.Path(snakemake.input.resultsdb)

con = sqlite3.connect(path2)
pcap = pd.read_sql_query("SELECT z, g, level FROM var_tot_pcap_z", con)

df_pcap = (
    pcap
    .assign(level = lambda x : x.level.mul(1000))
    .query("g == ['Windonshore','Windoffshore'] & level > 1")
    .rename(columns={'z' : 'zone','g' : 'turbine_type'})
    .replace({'Windonshore' : 'onshore','Windoffshore' : 'offshore'})
)

df_merged = (
    df_pcap
    .merge(df_ef,on = 'turbine_type',how='inner')
    .assign(
        number_turbines = lambda x : x.level.div(x.turbine_type.map(turbinedict)),
        job_years = lambda x : x.number_turbines*x.value
    )
)

df_merged.loc[:,['zone','turbine_type','project_stage','job_type','job_years']].to_csv(snakemake.output[0],index=False)