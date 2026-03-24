Advanced examples
=================
.. note::

   This page is work-in-progress.


Replacing ERA5 with SARAH
~~~~~~~~~~~~~~~~~~~~~~~~~~


Modelling to generate alternatives (MGA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Transmission grid network and net transfer capacities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the components of the technoeconomic input is the transmission infrastructure that connects zones and allows for the flow of electricity between them. Different models represent transmission differently. Models like PyPSA-Eur employ a power flow model, whereas highRES models transmission through a transhipment formulation. 

The currently most comprehensive publicly available dataset of the European high-voltage grid is a `representation of the European high-voltage grid based on OpenStreetMap <https://doi.org/10.5281/zenodo.12799201>`_ processed and described by `Xiong et al. (2025) <https://www.nature.com/articles/s41597-025-04550-7>`_. The dataset covers AC lines from 220 to 750 kV and all DC lines. Power capacity data is in nominal apparent power capacity (S_nom) [MV A], and represents the maximum technical capacity of a line. However, when representing the transmission grid through a transhipment formulation, we need to consider how much net power can be safely transferred accross borders, taking into account security constraints and internal bottlenecks. This is often referred to as Net Transfer Capacities (NTCs), which can be acquired from ENTSO-E's Transparency platform between bidding zones at different time horizons. However, there is no European-wide high-resolution dataset on NTCs at higher spatial resolution. 

To go from the nominal apparent power capacity data in the European high-voltage grid, we need to make approximations. 

First, we can acquire data from ENTSO-E on available NTC for evaluation. Data can be downloaded directly from the `Transparency portal <https://transparency.entsoe.eu/>`_, or more systematically with the ENTSO-E API (e.g. via `entsoe-py <https://github.com/EnergieID/entsoe-py>`_). 

Table 1. Exemplary NTC data from ENTSO-E
==========  ==========  ==========
zone_from     zone_to      NTC_MW  
==========  ==========  ==========
AT             CH          1200
AT             CZ          900
AT             DE-LU       4900
AT             HU          300
==========  ==========  ==========

The data acquired on NTCs from ENTSO-E could a bit like Table 1. 

Table 2. Raw data snippet from Xiong et al. 
========  ========  ========  ========
line_id   voltage    s_nom    geometry  
========  ========  ========  ========
merged_relation/10264161-225+1      225      502.728     LINESTRING (5.97276 43.15786, 5.97909 43.15992...
merged_relation/11762422-400-b+2    400 	   1787.476 	LINESTRING (6.05553 46.25095, 6.06287 46.25738...
========  ========  ========  ========

Table 2 show a few of the column in the high-voltage grid data from Xiong et al. that we are interested in. The geometry column can be used to identify which lines connect different nodes. 


.. image:: /_static/figures/voltage_grid_bins.jpg
   :alt: Overview of the European high-voltage grid.
