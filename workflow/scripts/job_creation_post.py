import pandas as pd
import pathlib
import sqlite3

turbine_size_dict = snakemake.params.turbinedict

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
        number_turbines = lambda x : x.level.div(x.turbine_type.map(turbine_size_dict)),
        job_years = lambda x : x.number_turbines*x.value
    )
)

df_merged.loc[:,['zone','turbine_type','project_stage','job_type','job_years']].to_csv(snakemake.output[0],index=False)