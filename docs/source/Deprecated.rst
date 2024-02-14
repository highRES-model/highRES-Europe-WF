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
