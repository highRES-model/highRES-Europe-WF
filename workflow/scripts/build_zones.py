import numpy as np


zones=np.array(snakemake.params.aggregated_regions)

if snakemake.params.focus_countries is not None:
 
    for key,z in snakemake.params.focus_countries.items():
    
        zones=np.append(zones[zones!=key],z)
        
if snakemake.params.aggregated_countries is not None:
    
    for key,z in snakemake.params.aggregated_countries.items():
        
        zones=np.append(zones[np.isin(zones,z,invert=True)],key)
        
np.savetxt(snakemake.output[0],zones,fmt="%s")
