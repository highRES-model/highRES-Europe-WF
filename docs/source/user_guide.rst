Detailed user guide
====================

Nomenclature
-------------
* Zone | A higher spatial level, typically countries
* Region | A lower spatial level, for example NUTS2 or NUTS3 regions. 
* .gms files | .gms files are the primary file format used in GAMS. They contain the modeling code written in the GAMS language. A .gms file typically includes the model definition, sets, parameters, variables, equations, objective functions, and any additional GAMS-specific syntax required to structure the model.
* .dd files | .dd files contain data input which is used in the model. 

Abbreviations
--------------
* Res | Reservoir 
* RoR | Run-off-River
* UC | Unit commitment
* WF | Workflow

highRES-Europe consists of two modules, a GAMS module and a workflow (WF) module.  

.. _workflow-label:

Workflow
------------
Work-in-progress

GAMS
------------

The general algebraic modeling system (GAMS) is the modelling system for optimisation which highRES is written in. The main GAMS file of highRES is ``highres.gms``. Here, the essential variables and equations are declared and defined. 

* The **objective equation** details the total system cost of the model, which is to be minimised. This includes capital expenditures, fixed operation and maintenance and variable operation and maintenance for generation, storage and transmission infrastructure. 
* The **demand balance equation** ensures that the supply ≥ demand for every hour in every zone. 
* The **transmission equations** allows for electricity to flow between zones for every hour. 
* Additional **miscellaneous equations** 

For more descriptions of GAMS components and syntax, see the `GAMS documentation <https://www.gams.com/latest/docs/>`_.

| **Objective equation**
| The objective equation of the model governs the central objective of the model. By default this is to minimise the total system cost, but it can be changed, as in the case of Modelling to Generate Alternatives (MGA).

The objective equation (``eq_obj``) and the total system cost is composed of generation, storage and transmission costs. Depending on the setup, start up costs (from UC) as well as penalty generation (value of lost load) may be included. Cost are divided into capital expenditure (Capex), fixed operation and maintenance costs (FOM) and variable operation and maintenance (VOM). There are no VOM costs included for transmission. 

.. math::
   :nowrap:

   \begin{gather*}
   \text{generation costs} = \sum_{g,z}(gen\_capex_{g} \times gen\_capacity_{g,z}) + \sum_{g,z,h}(gen\_VOM_{g,h} \times gen_{g,z,h}) + \\ \sum_{g,z}(gen\_FOM_{g} \times gen\_capacity_{g,z}) \\

   \text{storage costs} = \sum_{g,z}(store\_capex_{g} \times store\_capacity) + \sum_{g,z,h}(store\_gen_{g,z,h} \times store\_VOM_{g,h}) + \\ \sum_{g,z}(store\_FOM_{g} \times store\_capacity_{g,z}) \\

   \text{transmission costs} = \sum_{g,z}(trans\_capex_{g} \times trans\_cap_{g}) \\ + \sum_{g,z}(trans\_FOM_{g} \times trans\_cap_{g}) \\

   \text{penalty generation costs} = \sum_{z,h}(pgen\_cost \times pgen_{z,h}) \\
   \end{gather*}


The total system cost is then the sum of these different components, which, typically, are to be minimised. 

.. math::
   \min \text{total system cost} = \text{generation costs} + \text{storage costs} + \\ \text{transmission costs} + \text{penalty generation costs}

| **Demand balance equation**
| The demand balance equation (``eq_elc_balance(h,z)``) ensures that the demand is met in each of the zones (*z*) and for every hour (*h*) of the model. The demand can be met by in-region electricity generation, imported electricity from neighbouring regions through transmission infrastructure or discharging either of the storage technologies. At a high cost, the model can, if penalty generation is turned on, shed load. 

| **Transmission equations**
| The electricity transmission of highRES is represented using a computationally efficient linear transshipment formulation, where electricity flows similarly to fuel transport in pipelines. The benefit with a transshipment formulation compared to e.g. an *direct current optimal flow model* is that it is simpler `(Matar and Elshurafa, 2019) <https://doi.org/10.1016/j.egyr.2018.04.004>`_. 

The flow of electricity is constrained to not exceed the transmission capacity (``eq_trans_flow``) and bidirectionality is required (``eq_trans_bidirect``).  

| **Miscellaneous equations**
| One important miscellaneous equation is the CO₂ constraint (``eq_co2_budget``). It limits the total CO₂ emissions to be lower than a given value. The constraints scale with demand and as such indicate a maximum average emission intensity. By default, the intensity is 2gCO₂/kWh. 

Additionally, the model includes a set of submodules, containing various features. In general, these can be controlled by an IF statement. 

Module for data input 
~~~~~~~~~~~~~~~~~~~~~~

Whereas **highres.gms** contains the essential variables and equations, the data input submodule (``highres_data_input.gms``) contains the data input. This includes, among other things, the demand, the generation, the storage and the transmission data.

.. code-block:: gams

   $INCLUDE highres_data_input.gms

Within ``highres_data_input.gms`` numerous data files are loaded, such as for the defined spatial levels (regions and zones) as well as the temporal extent, technoeconomic generation and transmission data as well as the demand data. These are generated in the :ref:`workflow <workflow-label>`. 

The files are loaded through the following code:

::

       r regions /
       $BATINCLUDE %datafolderpath%/%vre_restrict%_regions.dd
       /

       z zones /
       $BATINCLUDE %datafolderpath%/zones.dd
       /

       $INCLUDE %datafolderpath%/%weather_yr%_temporal.dd

       $INCLUDE %datafolderpath%/%psys_scen%_gen.dd

       $INCLUDE %datafolderpath%/trans.dd

       $INCLUDE %datafolderpath%/%esys_scen%_demand_%dem_yr%.dd


Note that ``%datafolderpath%``, and other % enclosed variables are defined through Snakemake (see :ref:`workflow <workflow-label>` for further details). 

Before we go through the contents of those files, we need to introduce an important set, namely *lt*. 

.. code-block:: gams
   Sets

   lt / UP, LO, FX /

*lt* defines three types of limits that are loaded together with the technoeconomic input data. These are the upper limit (UP), the lower limit (LO) and the fixed limit (FX). These are used, for example in ``parameter gen_lim_pcap_z(z,g,lt);``. For example, in the line ``DK.HydroRoR.UP 0.009`` in ``gen.dd``, the upper limit for the generation capacity of run-off-river hydropower in Denmark is set to 0.009. This means that the model is allowed to build up to 0.009 GW of run-off-river hydropower in Denmark. If on the contrary, UP would be replaced by FX, the model would be forced to build exactly 0.009 GW of run-off-river hydropower in Denmark. 

Now, to the input data files.

.. code-block:: gams

   r regions /
   $BATINCLUDE %datafolderpath%/%vre_restrict%_regions.dd
   /

The regions.dd file contains the regions, which are the lower spatial level. 

.. code-block:: gams

   z zones /
   $BATINCLUDE %datafolderpath%/zones.dd
   /
   ;

The zones.dd file contains the zones, which are the higher spatial level.

.. code-block:: gams

   $INCLUDE %datafolderpath%/%weather_yr%_temporal.dd

The temporal.dd file contains the set h, for the temporal dimension in the model. Typically, this is a range between 0 and 8759, representing the hours of the year. 

.. code-block:: gams

   $INCLUDE %datafolderpath%/%psys_scen%_gen.dd

The gen.dd file contain information on generation technologies and their characteristics. It includes the ``set g``, with the different generation technologies, as well as subsets for, among other things, which technologies are variable (``set_vre(g)``) or not (``set_nonvre(g)``). Additionally, there are power capacity limits and existing infrastructure through the parameter ``gen_lim_pcap_z`` and ``gen_exist_pcap_z``, respectively. Similarly, there are energy capacity limits (storage) and existing infrastructure for reservoir hydro through the parameter ``gen_lim_ecap_z`` and ``gen_exist_ecap_z``, respectively. 

There are a few additional parameters, such as emission factors (``gen_emisfac``), cost parameters (``gen_capex``, ``gen_varom``, ``gen_fom``, ``gen_fuelcost``) and features related to unit commitment, if that is turned on. 

.. code-block:: gams

   $INCLUDE %datafolderpath%/trans.dd

The trans.dd file contains the ``set trans`` which includes the types of transmission technologies (typically HVAC400KV and HVDCSubsea) as well as the transmission links available to the model ``set trans_links`` and their associated distance ``parameter trans_links_dist`` and capacity limit ``parameter trans_links_cap``. 

.. code-block:: gams

      $INCLUDE %datafolderpath%/%esys_scen%_demand_%dem_yr%.dd

This file contains the demand, stored in the parameter ``demand(z,h)``. The demand is given in MWh for every hour and zone.

Module for storage
~~~~~~~~~~~~~~~~~~~~~~~~

The option of modelling storage in highRES is controlled in the $setglobal statement, whereas the IF statement loads the external storage submodule.

.. code-block:: gams

   $setglobal storage "ON"

   $IF "%storage%" == ON $INCLUDE highres_storage_setup.gms

By default, storage is turned on. 

A few important equations is the storage balance equation, the maximum storage level constraint and the storage end constraint.

The storage balance equation (``eq_store_balance(h,s_lim(z,s))``) models the storage level of each storage technology (*s*) for every hour (*h*) and zone (*z*). Essentially, the storage level (``var_store_level(h,z,s)``) is based on the electricity of the previous hour, with additionally stored electricity going into the storage level and electricity used for consumption subtracted from it. Additionally, there are efficiency losses and self-discharge. 

The storage level is constrained (``eq_store_level(s_lim(z,s),h)``) to always be lower or equal to the maximum storage capacity. Furthermore, the storage technologies are set to be cyclical (``eq_store_end_level``), meaning that they are not necessarily empty in the first hour of the model, but that they need to end at the same level as they started. 

Module for reservoir hydropower
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whereas run-off-river hydropower functions the same as other VREs, reservoir hydropower functions differently. Again, the $setglobal controls whether it is activated or not, and the IF statement loads the submodule (``highres_hydro.gms``).

.. code-block:: gams
    
   $setglobal hydrores "ON"

    $IF "%hydrores%" == ON $INCLUDE highres_hydro.gms

Reservoir hydropower functions similar to a storage technology, but with a natural inflow of energy (electricity) ``parameter hydro_inflow(h,z,hydro_res)``, as opposed to charging electricity from the grid. The storage level ``var_hydro_level`` at any given hour is the storage level in the previous hour, plus the inflow of water (in energy units), minus the electricity generated and water which is "spilled" if it is necessary to e.g. not overflow the reservoir. The inflow is loaded as an input, and generated in the :ref:`workflow <workflow-label>`. 

Additional equations ensure that the level of the reservoir does not exceed the maximum storage level ``eq_hydro_level(h,gen_lim(z,hydro_res))`` and not generate more electricity than the maximum power capacity ``eq_hydro_gen_max(h,gen_lim(z,hydro_res))``.

highRES does not include any cascading effects, meaning that the outflow of one reservoir is not the inflow of another. Rather, the model sees one large reservoir at the zonal or regional level, depending on the setup. However, the hydro power inflow is normalised, based on historical production data, to ensure that the total electricity available corresponds with reality. See the :ref:`workflow <workflow-label>` for more details. 

Module for EV flexibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Upcoming work.

References
-------------

Matar, W., & Elshurafa, A. M. (2018). Electricity transmission formulations in multi-sector national planning models: An illustration using the KAPSARC energy model. Energy Reports, 4, 328–340. https://doi.org/10.1016/j.egyr.2018.04.004
