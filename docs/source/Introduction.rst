Introduction
==============
The high spatial and temporal resolution electricity system model, highRES, is used to design cost effective, flexible and weather resilient electricity systems for Great Britain and Europe. The model is specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. As the proportion of renewables in electricity generation increases, there will be increasing imbalances between electricity demand and supply. highRES is a high resolution electricity system model that simultaneously considers infrastructure planning (investment) and operational (dispatch) decisions to identify the most cost-effective strategies to cope with growing shares of intermittent renewables. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage. 

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of unit and system constraints. It can model a variety of technical characteristics of thermal generators (e.g. ramping restrictions, minimum stable generation, startup costs, minimum up and down times) depending on the requirements of the research question, their CO2 emissions, and the technical characteristics of a variety of energy storage options. The transmission grid is represented using a linear transport model.

.. image:: figures/highRES_framework.jpg
   :alt: The highRES modelling framework.


highRES-Europe is further described in `Price & Zeyringer (2022) <https://doi.org/10.1016/j.softx.2022.101003>`_.


Data
-------------

Although highRES-Europe is customisable and any data can be used, as long as it is formatted in the correct way, there is some default data which is commonly used. 

Demand 
~~~~~~~~~~~~~~

Demand time series are originally based on historical data from the `European Network of Transmission System Operators for Electricity (ENTSO-E) Transparency Platform <https://transparency.entsoe.eu/dashboard/show>`_ , but need to be adjusted to account for inconsistencies and missing data. `van der Most <https://doi.org/10.1016/j.rser.2022.112987>`_ uses climate data and applies a logistic smooth transmission regression (LSTR) model to the ENTSO-E dataset to correlate historical electricity demand to temperature and generate daily electricity demand for a set of European countries. Subsequently, `Frysztacki, van der Most and Neumann <https://zenodo.org/records/7070438#.Y2OfViYo9hE>`_ uses hourly profiles from the `Open Power Systems Database <https://data.open-power-system-data.org/time_series/>`_ to disaggregate the daily electricity demand to an hourly resolution, on a country level.

- ERA5 for on and offshore wind capacity factors and runoff data (<https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>)
- CMSAF-SARAH2 for solar PV capacity factors (<https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>)
- Cost and technical data is taken from UKTM (<https://www.ucl.ac.uk/energy-models/models/uk-times>) with data from the JRC report "Cost development of low carbon energy technologies: Scenario-based cost trajectories to 2050" (<https://publications.jrc.ec.europa.eu/repository/handle/JRC109894>) used to update some areas.
- Data on run-of-river, reservoir and pumped hydro power capacities is taken from <https://transparency.entsoe.eu/>, <https://www.entsoe.eu/data/power-stats/> and <https://github.com/energy-modelling-toolkit/hydro-power-database>
- Energy storage capacities for reservoir and pumped storage are taken from Schlachtberger et al. (2017) and Geth et al. (2015) respectively (see <https://doi.org/10.1016/j.energy.2017.06.004> and <https://doi.org/10.1016/j.rser.2015.07.145>).
- Exclusion data

