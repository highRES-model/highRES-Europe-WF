import numpy as np


zones=np.array(snakemake.params.aggregated_regions)

if snakemake.params.focus_countries:
 
    for (key,z) in snakemake.params.focus_countries.items():
    
        zones=np.append(zones[zones!=key],z)
        
np.savetxt(snakemake.output[0],zones,fmt="%s")
