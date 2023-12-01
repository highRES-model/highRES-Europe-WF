# highRES-model

To run the full workflow, two datapackages are needed they can be downloaded from:

1. (~80MB compressed, ~300MB uncompressed) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1>.
2. (~10GB) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1>

## Windows
1. Clone the repository
2. Install snakemake
    - Download miniforge windows exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>
    - Install Minforge
    - Run the minimal install of the snakemake environment `mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas zstd`
3. Activate the snakemake environment
4. Navigate to the repository in your snakemake conda environment shell
4. Get the required input files
    ```
   curl -o shared_input.tar -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1" -o resources.tar.zst -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1"
   ```
5. Extract the required input files
    ```
    zstd -d resources.tar.zst
    mkdir resources
    tar xf resources.tar -C resources
    mkdir shared_input
    tar xf shared_input.tar -C shared_input
    ```
5. Run snakemake -c --use-conda

## Notes (to be removed)

This file collects all documentation from this repository and all documentation from "documentation-in-progress" in the files:

- README.md
- doc/Investigation.md
- doc/README.md
- doc/Unit commitment notes.md
- units_highRES.txt

## What is highRES?

Welcome to the repository for the European version of the high temporal and spatial resolution electricity system model (highRES-Europe). The model is used to plan least-cost electricity systems for Europe and specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage.

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data.

## How to run the model

This repository contains all GAMS code and necessary text files/GDX format input data for a 8760 hour model run. To execute the code:

1. GAMS must be installed and licensed. This version was tested/developed with GAMS version 27.2.0.
2. All files must be in the same directory and then open highres.gms, the main driving script for the model, in the GAMS IDE and hit run.
3. Full model outputs are written into the file "hR_dev.gdx" which is written into the same directory as the code/data and can be viewed using the GAMS IDE. Outputs include: the capacity of generation, storage and transmission by node, the hourly operation of these assets (including flows into and out of storage plus the storage level and total system costs.)
4. The GDX output file can be converted to SQLite using the command line utility gdx2sqlite which is distributed with GAMS. From the command line do "gdx2sqlite -i hR_dev.gdx -o hR_dev.db -fast". This SQLite database can then be easily read by Python using, e.g., Pandas.

## Data

This model uses the following data sources:

- ERA5 for on and offshore wind capacity factors and runoff data (<https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>)

- CMSAF-SARAH2 for solar PV capacity factors (<https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>)

- Demand data is from ENTSO-E for the year 2013

- Cost and technical data is taken from UKTM (<https://www.ucl.ac.uk/energy-models/models/uk-times>) with data from the JRC report "Cost development of low carbon energy technologies: Scenario-based cost trajectories to 2050" (<https://publications.jrc.ec.europa.eu/repository/handle/JRC109894>) used to update some areas.

- Data on run-of-river, reservoir and pumped hydro power capacities is taken from <https://transparency.entsoe.eu/>, <https://www.entsoe.eu/data/power-stats/> and <https://github.com/energy-modelling-toolkit/hydro-power-database>

- Energy storage capacities for reservoir and pumped storage are taken from Schlachtberger et al. (2017) and Geth et al. (2015) respectively (see <https://doi.org/10.1016/j.energy.2017.06.004> and <https://doi.org/10.1016/j.rser.2015.07.145>).

- Exclusion data

## highRES-model publications

This repository contains the power system expansion and dispatch model "highRES". A version of it has been used in the following paper:

- Zeyringer, Marianne, James Price, Birgit Fais, Pei-Hao Li, and Ed Sharp. ‘Designing Low-Carbon Power Systems for Great Britain in 2050 That Are Robust to the Spatiotemporal and Inter-Annual Variability of Weather’. Nature Energy 3, no. 5 (May 2018): 395–403. <https://doi.org/10.1038/s41560-018-0128-x>.

      Note it is currently set to run with a carbon price and budget. Only one weather/demand year (2013) - the demand is not rescaled, i.e. it is historic.

      Another version of it is available on [GitHub](https://github.com/highRES-model/highRES-Europe) under the MIT license and is mirrored here in the branch "james".

      Documentation that has been found on the internet is collected in the folder "doc".
      An incomplete documentation which stems from trying to understand the code can be found there as well in a file called ["Investigation.md"](https://github.com/highRES-model/documentation-in-progress/blob/master/doc/Investigation.md). "Investigation.md" is now incorporated into this document.

Furthermore, the following article provides a high level overview:

- Price, James, and Marianne Zeyringer. 2022. ‘HighRES-Europe: The High Spatial and Temporal Resolution Electricity System Model for Europe’. SoftwareX 17 (January): 101003. <https://doi.org/10.1016/j.softx.2022.101003>.

### Publications using highRES

- Zeyringer, Marianne, Birgit Fais, and James Price. 2016. ‘“New” or “Old” Technologies to Decarbonize UK’s Electricity System? A Long-Term High Spatial and Temporal Resolution Assessment for Marine and Wind Energy’. In 2016 13th International Conference on the European Energy Market (Eem). Vol. 2016-July. <https://doi.org/10.1109/EEM.2016.7521318>.

- Zeyringer, Marianne, Worrell, Ernst, Schmid, Erwin, and University Utrecht. 2017. ‘Temporal and Spatial Explicit Modelling of Renewable Energy Systems : Modelling Variable Renewable Energy Systems to Address Climate Change Mitigation and Universal Electricity Access’. Utrecht University. <https://dspace.library.uu.nl/handle/1874/350879>.

- Moore, Andy, James Price, and Marianne Zeyringer. 2018. ‘The Role of Floating Offshore Wind in a Renewable Focused Electricity System for Great Britain in 2050’. Energy Strategy Reviews 22 (November): 270–78. <https://doi.org/10.1016/j.esr.2018.10.002>.

- Price, James, Marianne Zeyringer, Dennis Konadu, Zenaida Sobral Mourao, Andy Moore, and Ed Sharp. 2018. ‘Low Carbon Electricity Systems for Great Britain in 2050: An Energy-Land-Water Perspective’. Applied Energy 228 (October): 928–41. <https://doi.org/10.1016/j.apenergy.2018.06.127>.

- Zeyringer, Marianne, Birgit Fais, Ilkka Keppo, and James Price. 2018. ‘The Potential of Marine Energy Technologies in the UK - Evaluation from a Systems Perspective’. Renewable Energy 115 (January): 1281–93. <https://doi.org/10.1016/j.renene.2017.07.092>.

- Zeyringer, Marianne, James Price, Birgit Fais, Pei-Hao Li, and Ed Sharp. 2018. ‘Designing Low-Carbon Power Systems for Great Britain in 2050 That Are Robust to the Spatiotemporal and Inter-Annual Variability of Weather’. Nature Energy 3 (5): 395–403. <https://doi.org/10.1038/s41560-018-0128-x>.

## Investigation

This section contains the findings from learning how to use the model.

### Filetypes

The repository contains the following file extensions:

| Filetype           | Description                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **.dd**            | Seems to contain mainly input data for the model (sets, parameters, variables?).                                                                                                                                                                                                                                                                                                                                                        |
| **.gms**           | Seems to contain model code.                                                                                                                                                                                                                                                                                                                                                                                                            |
| **.gpr**           | Seems to contain metadata for the model. File containing the configuration of a project (it is automatically generated when using GAMS IDLE).                                                                                                                                                                                                                                                                                           |
| **.gdx**           | "A GDX file is a file that stores the values of one or more GAMS symbols such as sets, parameters variables and equations. GDX files can be used to prepare data for a GAMS model, present results of a GAMS model, store results of the same model using different parameters etc. A GDX file does not store a model formulation or executable statements." Source: [GAMS Documentation](https://www.gams.com/latest/docs/UG_GDX.html) |
| **.md**            | [Markdown](https://en.wikipedia.org/wiki/Markdown) files used to document the project.                                                                                                                                                                                                                                                                                                                                                  |
| **.yml** **.yaml** |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| **.py**            |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| **.ipynb**         |                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| **.opt**           |                                                                                                                                                                                                                                                                                                                                                                                                                                         |

### GAMS components

The follow components are the basic building blocks of a GAMS model.
| Component | Description |
|-|-|
| sets | Sets may be seen as categories or groups of elements and are assigned names. |
| data | Definition of the input data through parameter declaration and definition. |
| variables | Elements that define a solution to the optimization problem?. They should be declared in order to be used. |
| Equations | Objective function(s), constraints of the problem, and additional functions that define the optimization problem. Each equation must be declared with a name and subsequently defined. |
| Model and solver definition | Definition of which equations are considered part of the model, type of problem, and options for the solver. |

For more information see the GAMS documentation <https://www.gams.com/latest/docs/>.

### Documentation of individual **.dd** files

- **2013_temporal.dd** contains the

  - set h, which contains 8760 hours of a full year.
    This set starts counting at zero.
  - The set yr contains only one value, the digit 0, representing the considered years (0 is the base case?).
  - The set hr2yr_map is equal to the set h but has a 0. in front of every value. The notation is year.hour.

  Since highRES is able to be run with multiple consecutive years in one go, it needs to understand which hour belongs to which year. This is done by the set hr2yr_map, which contains one dimension for the year and one dimension for the hours. Both dimensions have to be declared beforehand.

- **dev_regions.dd** contains the

  - 31 regions of the model (countries in this case), abbreviated by two letters.
    These regions are the higher resolution spatial layer and can be aggregated to zones. This aggregation needs a mapping, which is given in the file

- **NewPl_min80_co2_budget.dd** contains the

  - CO~2~ Budget for the model run.
    Presumably it is assumed constant for the whole year.
    TODO check and clarify unit.

- **NewPl_min80_gen_parameters.dd** specifies

  - gen_lim_pcap_z describes the non VRE capacity limit for each zone
  - gen_exist_pcap_z describes the existing zonal power capacity
  - gen_emisfac describes the generational emission factor
  - gen_lim_ecap_z describes the generation capacity limit for each zone (HydroRes)
  - gen_exist_ecap_z describes the generation capacity for each existing zonal generation technology (HydroRes)
  - gen_maxramp describes the ramp-up time (?) of NG and Nuclear (MW/min)
  - gen_capex2050 describes the CAPEX of the different technologies in 2050 (?)
  - gen_varom describes the variable cost (O\&M) of the different technologies
  - gen_fom describes the fixed costs (O\&M) of the different technologies
  - gen_fuelcost2050 describes the expected fuel costs in 2050 (?)
  - gen_cap2area contains the convertion factor from capacity to area. Debugging the code, I think this factor is a convertion from area to capacity. TODO: check the equations where it is used.
  - gen_af describes the fractional availability factor for each technology
  - gen_mingen
  - gen_peakaf describes the fractional peak availability factor (?)

- **NewPl_min80_gen_set_all.dd** specifies

  - the 8 sets (technologies) of the model

- **NewPl_min80_gen_set_nonvre.dd** check

  - Subset containing the non VRE technologies

- **NewPl_min80_gen_set_vre.dd** check

  - Subset containing the VRE technologies

- **NewPl_min80_gen_set_hydrores.dd** specifies

- **NewPl_min80_nonrescale_demand_2013.dd** check

  - Electricity demand for each country/zone in 2013

- **NewPl_min80_store_parameters.dd** specifies

  - store_lim_ecap_z contains the energy? storage limit capacity
  - store_exist_pcap_z defines the existing storage capacity
  - store_exist_ecap_z defines the existing energy? storage capacity
  - store_max_freq describes the fractional maximum constribution to frequency response (definition obtained from the highres_storage_setup.gms file)
  - store_max_res describes the fractional maximum contribution to operating reserve
  - store_eff_in defines the fractional charge efficiency
  - store_eff_out defines the fractional discharge efficiency
  - store_loss_per_hr contains the fractional energy loss per hour (self discharge)
  - store_p_tp_e describes the power to energy ratio
  - store_varom describes the variable O&M (GBP k per MWh)
  - store_ecapex2050 describes the annualised energy capex (GBP k per MWh)
  - store_pcapex2050 describes the annualised power capex (GBP k per MW)
  - store_fom defines the fixed O&M (GBP k per MW per yr)
  - store_af defines the fractional availability factor

- **NewPl_min80_store_set_all.dd** check

  - Set containing the storage technologies (?)

- **trans_links.dd**

  - contains information on the transmission capacity between the regions of the model (see dev*regions.dd). Each row appears to be structured ~~\_FROM.TO.TYPE*~~ _TO.FROM.TYPE_ (Debugging the model, I realized that the notation is different. Check the highres.gms file, where balance is carry out: import and export) . For example, AT.HU.HVAC400KV, describes that there is an 400kV HVAC transmission line from ~~Austria~~Hungary to ~~Hungary~~Austria. This is a subset depending on the regions and type. TODO: Check this notation

- **trans_parameters.dd**

  - contains several parameters for the transmission lines (see trans_links.dd).
    - trans_links_dist describe the length (km) of the transmission cable.
    - trans*links_cap \_probably* describes the capacity of the transmission cable (currently set to 0).
    - trans*loss \_probably* describes the transmission losses (\%) for 400kV HVAC cables and subsea HVDC cables.
    - trans_capex describes the capital expenditures (CAPEX) of 400kV HVAC cables and subsea HVDC cables.
    - trans_varom describes the variable operation and maintenance cost for 400kV HVAC cables and subsea HVDC cables (currently set to 0).

- **trans_set_all.dd**

  - set containing two elements: HVAC400KV and HVDCSubsea

- **vre_areas_2013_dev.dd**

  - contains the area available in 2013 for additional (?) VRE generation in the different regions. TODO: Check if it is _addition_

- **vre_areas_2017_dev.dd** contains
  - the area available in 2017 for additional (?) VRE generation in the different regions
- **zones.dd** contains the zones the model runs on, which are the lower resolution spatial layer, on which transmission, storage and demand/supply balancing happens. For runs where zones and regions differ, the regions get aggregated to zonal level. The mapping between zones and regions is done in

### Documentation of **.gms** files

- **highres.gms** \
  Contains the main model code.

  - It starts with [$ontext](https://www.gams.com/latest/docs/UG_DollarControlOptions.html#DOLLARonofftext).
  This is a [comment](https://www.gams.com/latest/docs/UG_GAMSPrograms.html#UG_GAMSPrograms_CommentsBlock), which hides the option called profile which would be set to one from GAMS.

  - The [profile](https://www.gams.com/latest/docs/UG_GamsCall.html#GAMSAOprofile) option is useful for performance analysis, as it prints time and memory statistics for different parts of the program.
      It appears as if the comment declaration $ontext is commented out itself by the star in front of it, meaning profiling is enabled.

  - Next is the option [limrow](https://www.gams.com/latest/docs/UG_GamsCall.html#GAMSAOlimrow) which changes the GAMS default behaviour to display the first three equations of each block, to supressing the equation output completely as explained [here](https://www.gams.com/latest/docs/UG_GAMSOutput.html).

  - The option [limcol](https://www.gams.com/latest/docs/UG_GamsCall.html#GAMSAOlimcol) is very similar, as it supresses the coeficcient as documented [here](https://www.gams.com/latest/docs/UG_GAMSOutput.html#UG_GAMSOutput_TheColumnListing). TODO: Understand the column listing

  - The option [solPrint](https://www.gams.com/latest/docs/UG_GamsCall.html#GAMSAOsolprint) disables the [solution listing](https://www.gams.com/latest/docs/UG_GAMSOutput.html#UG_GAMSOutput_TheSolutionListing). TODO: Understand the solution listing

  - The option [decimals](https://www.gams.com/latest/docs/UG_GamsCall.html#GAMSAOdecimals) shows how many decimal places should be displayed.

  - The option [$offlisting](https://www.gams.com/34/docs/UG_DollarControlOptions.html#DOLLARonofflisting) tells GAMS to not replicate the model code after this statement in the \*.lst file that is generated when running the model.

  - The option [$ONMULTI](https://www.gams.com/34/docs/UG_DollarControlOptions.html#DOLLARonoffmulti) controls the program flow and enables that parameters (called data statements in the documentation) are redefined. All the new entries are merged with the old values.
      This is not allowed by default.
      It maybe messes with the program flow .
      TODO: find out if the warning about two pass processing is relevant to us!

  - The option [$onEps](https://www.gams.com/35/docs/UG_DollarControlOptions.html#DOLLARonoffeps) is used to interpret zero values as EPS. It is useful when binary variables are used to express existence.

  - The option [$offDigit](https://www.gams.com/35/docs/UG_DollarControlOptions.html#DOLLARonoffdigit) controls the internal precision of numbers. This tells GAMS to use as much precision as possible and ignore the rest of the numbers.

  - The option [$setGlobal Varname text](https://www.gams.com/35/docs/UG_DollarControlOptions.html#DOLLARsetglobal) defines a global compile-time variable.
      TODO: understand the difference between Global, Local, and scoped compile-time variables [Here](https://www.gams.com/35/docs/UG_DollarControlOptions.html#DOLLARset) there is an example.

- **highres_data_input.gms** \
  Contains the declaration and definition of the sets and parameters used in the main model _highres.gms_. This file uses the information contained in the \*.dd files. The demand and vre capacity information is loaded in the model.

- **highres_storage_setup.gms** \
  Contains the declaration and definition of the sets, parameters and variables for the storage module. This file is used in main model code if the module _storage_ is active (_$setglobal storage "ON"_)

### Documentation of **.gdx** files

- **hr_DEV.gdx** \
  Contains the outputs produced by the model run.
  It can be converted to an SQLITE file with a tool supplied with GAMS.
  In this conversion, parameters with 0 dimensions are summarized in a table called scalars in the resulting *.sqlite file.

- **vre2013_dev.gdx** \
  Contains the information of the VRE capacities for each technology, zone and hour. This file contains the vre_gen parameter which is structured as _HR.VRE.ZONE cap_fac.

### Documentation of **.ipynb** files

- **highRES-build_weather.py.ipynb** \
  Exclusions:
  - Exclusions from WDPA.
  - Explain the reasoning behind which solar codes are included and exluded.
  - Areas above 2000 m above sea level and 15 degree inclination are excluded.
  - Offshore regions choosen by exclusive economic area and depth. Explain reasoning behind the different depths: \
  Bottom: 0-70 m \
  Floating: 70 - 1000 m \
  Future: more than 1000 m (ignore)


### Documentation of **.py** files

### Documentation of **Snakemake**

- Scenario definition.

## highRES versions

- highRES-UK \
  This repo contains the UK version of highRES. The first paper using it is under review at Energy and the preprint of that paper is available at <https://arxiv.org/abs/2109.15173>.

- highRES-Europe

- highRES-AtLAST

- highRES-GB

- highRES-Norway

## Units of parameters in highRES

- gen_exist_pcap_z: GW (most likely)

## Unit commitment notes

Conceptually the unit commitment implementation in this model can be seen in the
following way:  
To simplify the modelling of unit commitment, which we have to do to keep the
model computationally feasible, there are 4 options:

1. **"binary/boolean"** unit commitment: Here every powerplant will have its own
   on/off switch. We currently do not use this in our implementation at all.

2. **"integer unit commitment"**: Plants of same design are clustered. You can
   have 0,1 or 2 and so forth on/off but not 1.5 plants. The difference in this
   approach from the one above is that here all the plants are uniform (they have
   the same installed generation capacity, etc.)

3. **"linear unit commitment"**: Similar to 2) but all the integers are relaxed
   to be continuous (floats). This means that you can have 1.5 plants on/off but
   still we look at startup cost. The main reason for having this is to model are
   the following 2 things: - frequency response (50 HZ) which has a 10 second response window - operating reserve, which has a 20 minute response window and is then
   required to run for a ?full hour?

4. **"linear no unit commitment"**: Here we can have 1.5 plants on/off and it
   does not cost anything to switch them on/off.

The availibility factor also plays into this.
