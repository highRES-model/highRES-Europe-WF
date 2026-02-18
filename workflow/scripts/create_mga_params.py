import pandas as pd
from pathlib import Path

ALL_ZONES = snakemake.params.all_zones # Only necessary to declare the focus zones. 
cost_opt_file = snakemake.input.cost_optimal
mga_params_file = Path(snakemake.output.mga_params)
mga_mode_file = Path(snakemake.output.mgaMode)

slack_level = snakemake.wildcards.slack
objective = snakemake.wildcards.objective
focus_techs = snakemake.params.techs 
focus_zones = snakemake.params.zones

if focus_zones == ["*"]:
    focus_zones = ALL_ZONES

df_costOpt = pd.read_csv(cost_opt_file, sep='\t')
c_opt = df_costOpt.query("name == 'costs'")["level"].iloc[0] # Value of objective function

tech_lines = [f"  {t} 1" for t in focus_techs]

par_focus_g_block = "parameter par_focus_g(g) /\n"
par_focus_g_block += ",\n".join(tech_lines) + "\n/\n" if tech_lines else "/\n"

zone_lines = [f"  {z} 1" for z in focus_zones]

par_focus_z_block = "parameter par_focus_z(z) /\n"
par_focus_z_block += ",\n".join(zone_lines) + "\n/\n" if zone_lines else "/\n"

content_dd = f"""scalar
par_optimal_cost /
{c_opt}
/

parameter
par_slack /
{slack_level}
/

{par_focus_g_block}

{par_focus_z_block}
"""

mga_params_file.parent.mkdir(parents=True, exist_ok=True)
mga_params_file.write_text(content_dd, encoding="utf-8", newline="\n")

# 2) mgaMode.gms: compile-time macro for minimizing/maximizing
content_mode = f"$setglobal mode {objective}\n"
mga_mode_file.write_text(content_mode, encoding="utf-8", newline="\n")