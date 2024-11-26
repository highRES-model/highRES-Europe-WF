import pandas as pd
import atlite

onshore_turbine = snakemake.params.windturbines.get('onshore')
offshore_bottom_turbine = snakemake.params.windturbines.get('offshore_bottom')
onshore_size = atlite.resource.windturbine_rated_capacity_per_unit(onshore_turbine)
offshore_size = atlite.resource.windturbine_rated_capacity_per_unit(offshore_bottom_turbine)

turbinedict = {
    'onshore' : onshore_size,
    'offshore' : offshore_size,
}

new_project_stage = ['manufacturing','construction','OM','decomissioning']

job_type_new = 'direct'

OS_value = {'onshore' : 1, 'offshore' : 0} 

# Fixed employment factors
ef_dict = {
    'development_direct' : 0.69,
    'development_indirect' : 0.59,
    'development_induced' : 0.59,
    'manufacturing_indirect' : 5.68,
    'manufacturing_induced' : 3.65,
    'construction_indirect' : 3.08,
    'construction_induced' : 4.16,
    'OM_indirect' : 9.36,
    'OM_induced' : 12.53,
    'decomissioning_indirect' : 0.9,
    'decomissioning_induced' : 1.01,
}

df = pd.DataFrame(list(ef_dict.items()), columns=['combined', 'value'])
df[['project_stage', 'job_type']] = df['combined'].str.split('_', expand=True)

#Drop the 'combined' column
df = df.drop(columns=['combined'])

df_onshore = df.copy()
df_onshore['turbine_type'] = 'onshore'

df_offshore = df.copy()
df_offshore['turbine_type'] = 'offshore'

# Concatenate onshore and offshore DataFrames
df_final = pd.concat([df_onshore, df_offshore], ignore_index=True)

# Reorder the columns for better clarity
df_final = df_final[['turbine_type', 'project_stage', 'job_type', 'value']]

for turbine in ['onshore','offshore']:
    turbine_size = turbinedict.get(turbine)
    for stage in new_project_stage:
        if stage == 'manufacturing':
            new_value = 13.37 - (9.5 * OS_value.get(turbine))
        if stage == 'construction':
            new_value = 9.46-(3.65*OS_value.get(turbine))-4.01
        if stage == 'OM':
            new_value = 19.38 + (0.92*turbine_size) - 15.74
        if stage == 'decomissioning':
            new_value = 2.82-(2.11*OS_value.get(turbine))
     
        df_final.loc[len(df_final)] = {'turbine_type' : turbine, 'project_stage' : stage, 'job_type' : job_type_new, 'value' : round(new_value,3)}
        del(new_value)

df_final.to_csv(snakemake.output[0], index=False)        
