import numpy as np
from data2dd_funcs import wrapdd


zones=np.array(snakemake.params.aggregated_regions)

if snakemake.params.focus_countries is not None:
 
    for key,z in snakemake.params.focus_countries.items():
    
        zones=np.append(zones[zones!=key],z)
        
if snakemake.params.aggregated_countries is not None:
    
    for key,z in snakemake.params.aggregated_countries.items():
        
        zones=np.append(zones[np.isin(zones,z,invert=True)],key)
           

zones_out=wrapdd(zones.reshape(-1,1),"z","set")

if snakemake.params.unit_commitment == "ON":
    uc_zones=snakemake.params.unit_commitment_zones
    if type(uc_zones)==list:
        uc_zones=wrapdd(
            np.array(uc_zones).reshape(-1,1)
            ,"uc_z",
            "set"
            )
    if type(uc_zones)==str and uc_zones=="ALL":
        uc_zones=wrapdd(zones.copy().reshape(-1,1)
                        ,"uc_z","set") 
   
    np.savetxt(snakemake.output[0],
           np.concatenate((zones_out,uc_zones),axis=0)
           ,fmt="%s")
else:
    np.savetxt(snakemake.output[0],
               zones_out,
               fmt="%s")
