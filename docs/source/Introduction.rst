Introduction
==============
The high spatial and temporal resolution electricity system model, highRES, is used to design cost effective, flexible and weather resilient electricity systems for Great Britain and Europe. The model is specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. As the proportion of renewables in electricity generation increases, there will be increasing imbalances between electricity demand and supply. highRES is a high resolution electricity system model that simultaneously considers infrastructure planning (investment) and operational (dispatch) decisions to identify the most cost-effective strategies to cope with growing shares of intermittent renewables. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage. 

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of unit and system constraints. It can model a variety of technical characteristics of thermal generators (e.g. ramping restrictions, minimum stable generation, startup costs, minimum up and down times) depending on the requirements of the research question, their CO2 emissions, and the technical characteristics of a variety of energy storage options. The transmission grid is represented using a linear transport model.

.. image:: figures/softwareX_fig1.jpg
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

highRES-model publications
-------------
highRES has been used in a number of publications. An occasionally updated list is provided below, sorted by newest first. 

Viole, I., Shen, L., Camargo, L. R., Zeyringer, M., & Sartori, S. (2023). Sustainable Astronomy: A comparative Life Cycle Assessment of Off-grid Hybrid Energy Systems to supply large Telescopes [Preprint]. In Review. https://doi.org/10.21203/rs.3.rs-3281965/v2

Valenzuela-Venegas, G., Lode, M. L., Viole, I., Felice, A., Alonso, A. M., Camargo, L. R., Sartori, S., & Zeyringer, M. (2023). Designing renewable and socially accepted energy systems for astronomical telescopes: A move towards energy justice [Preprint]. In Review. https://doi.org/10.21203/rs.3.rs-3181969/v1

Viole, I., Valenzuela-Venegas, G., Zeyringer, M., & Sartori, S. (2023). A renewable power system for an off-grid sustainable telescope fueled by solar power, batteries and green hydrogen. Energy, 128570. https://doi.org/10.1016/j.energy.2023.128570

Price, J., Keppo, I., & Dodds, P. E. (2023). The role of new nuclear power in the UK’s net-zero emissions energy system. Energy, 262, 125450. https://doi.org/10.1016/j.energy.2022.125450

Price, J., & Zeyringer, M. (2022). highRES-Europe: The high spatial and temporal Resolution Electricity System model for Europe. SoftwareX, 17, 101003. https://doi.org/10.1016/j.softx.2022.101003

Price, J., Mainzer, K., Petrović, S., Zeyringer, M., & McKenna, R. (2022). The Implications of Landscape Visual Impact on Future Highly Renewable Power Systems: A Case Study for Great Britain. IEEE Transactions on Power Systems, 37(4), 3311–3320. https://doi.org/10.1109/TPWRS.2020.2992061

Price, J., Zeyringer, M., Konadu, D., Sobral Mourão, Z., Moore, A., & Sharp, E. (2018). Low carbon electricity systems for Great Britain in 2050: An energy-land-water perspective. Applied Energy, 228, 928–941. https://doi.org/10.1016/j.apenergy.2018.06.127

Zeyringer, M., Fais, B., Keppo, I., & Price, J. (2018). The potential of marine energy technologies in the UK – Evaluation from a systems perspective. Renewable Energy, 115, 1281–1293. https://doi.org/10.1016/j.renene.2017.07.092

Moore, A., Price, J., & Zeyringer, M. (2018). The role of floating offshore wind in a renewable focused electricity system for Great Britain in 2050. Energy Strategy Reviews, 22, 270–278. https://doi.org/10.1016/j.esr.2018.10.002

Zeyringer, M., Price, J., Fais, B., Li, P.-H., & Sharp, E. (2018). Designing low-carbon power systems for Great Britain in 2050 that are robust to the spatiotemporal and inter-annual variability of weather. Nature Energy, 3(5), Article 5. https://doi.org/10.1038/s41560-018-0128-x

Price, J., Zeyringer, M., Konadu, D., Sobral Mourão, Z., Moore, A., & Sharp, E. (2018). Low carbon electricity systems for Great Britain in 2050: An energy-land-water perspective. Applied Energy, 228, 928–941. https://doi.org/10.1016/j.apenergy.2018.06.127

Zeyringer, M., Fais, B., & Price, J. (2016). “New” or “old” technologies to decarbonize UK’s electricity system? 2016 13th International Conference on the European Energy Market (EEM), 1–5. https://doi.org/10.1109/EEM.2016.7521318
