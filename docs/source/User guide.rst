User guide
============

highRES-Europe consists of two modules, a GAMS module and a workflow (WF) module.  

GAMS
------------

The general algebraic modeling system (GAMS) is the modelling system which highRES-Europe is written in. The main GAMS file of highRES-Europe is **highres.gms**. Here, the essential variables and equations are declared and defined. 

* The **objective equation** details the total system cost of the model, which is to be minimised. This includes capital expenditures, fixed operation and maintenance and variable operation and maintenance for generation, storage and transmission infrastructure. 
* The **demand balance equation** ensures that the supply ≥ demand for every hour in every zone. 

Additionally, the main model file loads submodules, such as **highres_hydro.gms** for hydropower and **highres_storage_setup.gms** for storage, if those features should be included in the model. This is typically controlled through switches

.. code-block:: gams
    
    $IF "%hydrores%" == ON $INCLUDE highres_hydro.gms

