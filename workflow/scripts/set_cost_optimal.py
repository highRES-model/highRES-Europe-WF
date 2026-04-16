import sqlite3
import pandas as pd

costs=[]
con = sqlite3.connect(snakemake.input.resultspath)
(
    costs
    .append(
        (
            round(pd
            .read_sql_query("SELECT * from scalarvariables", con),3)
        )
    )
)
con.close()
df_costs = pd.concat(costs)

df_costs.to_csv(snakemake.output.cost_opt,sep='\t')