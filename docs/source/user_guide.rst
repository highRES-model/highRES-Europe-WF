Detailed user guide
====================

Nomenclature
-------------
* Zone
* Region

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

The objective equation of the model governs the central objective of the model. By default this is to minimise the total system cost, but in the case of Modelling to Generate Alternatives (MGA), this might be different. 

.. math::

   \text{min} \sum_z{placeholder}


The objective equation (eq_obj) is composed of generation, storage and transmission costs and, depending on the setup, start up costs (if UC is on) as well as penalty generation (value of lost load). The generation costs 

**Demand balance equation**


**Transmission equations**


**Miscellaneous equations**

One important miscellanneous equation is the CO:subscript:`2` constraint **eq_co2_budget**. It limits the total CO:subscript:`2` emissions to be lower than a given value. The constraints scale with demand and as such indicate a maximum average emission intensity. By default, the intensity is 2gCO:subscript:`2`/kWh.   

Module for hydropower
~~~~~~~~~~~~~~~~~~~~~~~~

Additionally, the main model file loads submodules, such as **highres_hydro.gms** for hydropower and **highres_storage_setup.gms** for storage, if those features should be included in the model. This is typically controlled through switches

.. code-block:: gams
    
    $IF "%hydrores%" == ON $INCLUDE highres_hydro.gms

