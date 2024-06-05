Welcome to highRES documentation!
===========================================

Welcome to the repository for the European version of the high temporal and spatial resolution electricity system model (highRES). The model is used to plan least-cost electricity systems for Europe and specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage.

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data.

Check out the :doc:`Getting started <getting_started>` section for further information on how to get started using the model. 

.. note::

   This project is under active development.

To inspect the results files produced (.db or .parquet), we recommend the
program "tadviewer" https://www.tadviewer.com/ .

Contents
---------

.. toctree::
   :maxdepth: 2

   Introduction <introduction>
   Getting started <getting_started>
   User guide <user_guide>
   Basic examples <basic_examples.rst>
   Advanced examples <advanced_examples.rst>
   Publications <publications>


.. toctree::
   :hidden:

   Deprecated
   
   
