Detailed user guide
====================

Nomenclature
-------------
* Zone | A higher spatial level, typically countries
* Region | A lower spatial level, for example NUTS2 or NUTS3 regions. 

Abbreviations
--------------
* UC | Unit commitment

highRES-Europe consists of two modules, a GAMS module and a workflow (WF) module.  

Workflow
------------


GAMS
------------

The general algebraic modeling system (GAMS) is the modelling system which highRES-Europe is written in. The main GAMS file of highRES-Europe is **highres.gms**. Here, the essential variables and equations are declared and defined. 

* The **objective equation** details the total system cost of the model, which is to be minimised. This includes capital expenditures, fixed operation and maintenance and variable operation and maintenance for generation, storage and transmission infrastructure. 
* The **demand balance equation** ensures that the supply ≥ demand for every hour in every zone. 
* The **transmission equations** allows for electricity to flow between zones for every hour. 
* Additional **miscellaneous equations** 

**Objective equation**

The objective equation of the model governs the central objective of the model. By default this is to minimise the total system cost, but it can be changed, as in the case of Modelling to Generate Alternatives (MGA).

.. math::

   \text{min} \sum_z{placeholder}


The objective equation (eq_obj) and the total system cost is composed of generation, storage and transmission costs. Depending on the setup, start up costs (from UC) as well as penalty generation (value of lost load) may be turned on. Cost are divided into capital expenditure (Capex), fixed operation and maintenance costs (FOM) and variable operation and maintenance (VOM). There are no VOM included for transmission. 

**Demand balance equation**

The demand balance equation (**eq_elc_balance(h,z)**) ensures that the demand is met in each of the zones (*z*) and for every hour (*h*) of the model. The demand can be met by in-zone electricity generation, imported electricity from neighbouring zones through transmission infrastructure or discharging either of the storage technologies. At a high cost, the model can, if penalty generation is turned on, shed load. 

**Transmission equations**

The flow of electricity is constrained to not exceed the transmission capacity (**eq_trans_flow**) and bidirectionality is required (**eq_trans_bidirect**).  

**Miscellaneous equations**

One important miscellanneous equation is the CO2 constraint (**eq_co2_budget**). It limits the total CO2 emissions to be lower than a given value. The constraints scale with demand and as such indicate a maximum average emission intensity. By default, the intensity is 2gCO2/kWh. 

Additionally, the model includes a set of submodules, containing various features. In general, these can be controlled by an IF statement. 

Module for storage
~~~~~~~~~~~~~~~~~~~~~~~~

The option of modelling storage in highRES is controlled in the $setglobal statement, whereas the IF statement loads the external storage submodule.

.. code-block:: gams

   $setglobal storage "ON"

   $IF "%storage%" == ON $INCLUDE highres_storage_setup.gms

By default, storage is turned on. 

A few important equations is the storage balance equation, the maximum storage level constraint and the storage end constraint.

The storage balance equation (**eq_store_balance(h,s_lim(z,s))**) models the storage level of each storage technology (*s*) for every hour (*h*) and zone (*z*). Essentially, the storage level (**var_store_level(h,z,s)**) is based on the electricity of the previous hour, with additionally stored electricity going into the storage level and electricity used for consumption subtracted from it. Additionally, there are efficiency losses and self-discharge. 

The storage level is constrained (**eq_store_level(s_lim(z,s),h)**) to always be lower or equal to the maximum storage capacity. Furthermore, the storage technologies are set to be cyclical (**eq_store_end_level**), meaning that they are not necessarily empty in the first hour of the model, but that they need to end at the same level as they started. 

Module for reservoir hydropower
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whereas run-off-river hydropower functions the same as other VREs, reservoir hydropower functions differently. Again, the $setglobal controls whether it is activated or not, and the IF statement loads the submodule (**highres_hydro.gms**).

.. code-block:: gams
    
   $setglobal hydrores "ON"

    $IF "%hydrores%" == ON $INCLUDE highres_hydro.gms

Module for EV flexibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~