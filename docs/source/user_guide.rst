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

**Objective equation**
The objective equation of the model governs the central objective of the model. By default this is to minimise the total system cost, but in the case of Modelling to Generate Alternatives (MGA), this might be different. 

.. math::

   \text{eq\_obj : } \text{costs} = \sum_{z} ( \text{costs\_gen\_capex}(z) + \text{costs\_gen\_fom}(z) + \text{costs\_gen\_varom}(z) \\
   \quad \quad \quad \quad + \text{costs\_gen\_start}(z) \text{ if } \text{\%UC\%} == \text{ON} + \text{costs\_store\_capex}(z) + \text{costs\_store\_fom}(z) + \text{costs\_store\_varom}(z) \text{ if } \text{\%store\_uc\%} == \text{ON} \\
   \quad \quad \quad \quad + \text{costs\_store\_start}(z) \text{ if } \text{\%pen\_gen\%} == \text{ON} + \text{costs\_pgen}(z) + \text{costs\_trans\_capex}(z) + \text{costs\_trans\_fom}(z) )


The objective equation (eq_obj) is composed of generation, storage and transmission costs and occasionally by start up costs (if UC is on) and penalty generation (value of lost load). The generation costs  

Module for hydropower
~~~~~~~~~~~~~~~~~~~~~~~~

Additionally, the main model file loads submodules, such as **highres_hydro.gms** for hydropower and **highres_storage_setup.gms** for storage, if those features should be included in the model. This is typically controlled through switches

.. code-block:: gams
    
    $IF "%hydrores%" == ON $INCLUDE highres_hydro.gms

