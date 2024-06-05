Basic examples
=================

.. note::

   This page is work-in-progress.

highRES comes shipped with a default setup, with enough data to run the model. 



Data
-------------
The default data is based on the following sources:

Demand 
~~~~~~~~~

Demand time series are originally based on historical data from the `European Network of Transmission System Operators for Electricity (ENTSO-E) Transparency Platform <https://transparency.entsoe.eu/dashboard/show>`_ , but need to be adjusted to account for inconsistencies and missing data. `van der Most <https://doi.org/10.1016/j.rser.2022.112987>`_ uses climate data and applies a logistic smooth transmission regression (LSTR) model to the ENTSO-E dataset to correlate historical electricity demand to temperature and generate daily electricity demand for a set of European countries. Subsequently, `Frysztacki, van der Most and Neumann <https://zenodo.org/records/7070438#.Y2OfViYo9hE>`_ uses hourly profiles from the `Open Power Systems Database <https://data.open-power-system-data.org/time_series/>`_ to disaggregate the daily electricity demand to an hourly resolution, on a country level.

Weather data
~~~~~~~~~~~~~~

Weather data is based on the `ERA5 reanalysis dataset <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ from the European Centre for Medium-Range Weather Forecasts (ECMWF). The ERA5 dataset provides hourly data on a 0.25° grid for the period 1979-01-01 to 2021-12-31. The data is used to calculate capacity factors for wind and solar power, as well as runoff data for hydro power. However, the ERA5 dataset performs poorly in some regions, particularly for solar irradiance. As such, there is an option to use the `SARAH dataset <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`_ for solar PV capacity factors instead, but this is covered in the :doc:`advanced examples <advanced_examples>`.

- Data on run-of-river, reservoir and pumped hydro power capacities is taken from <https://transparency.entsoe.eu/>, <https://www.entsoe.eu/data/power-stats/> and <https://github.com/energy-modelling-toolkit/hydro-power-database>
- Energy storage capacities for reservoir and pumped storage are taken from Schlachtberger et al. (2017) and Geth et al. (2015) respectively (see <https://doi.org/10.1016/j.energy.2017.06.004> and <https://doi.org/10.1016/j.rser.2015.07.145>).
- Exclusion data




- Cost and technical data is taken from UKTM (<https://www.ucl.ac.uk/energy-models/models/uk-times>) with data from the JRC report "Cost development of low carbon energy technologies: Scenario-based cost trajectories to 2050" (<https://publications.jrc.ec.europa.eu/repository/handle/JRC109894>) used to update some areas.

