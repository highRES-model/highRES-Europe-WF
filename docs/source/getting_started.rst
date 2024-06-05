Getting started
================

Requirements
-------------

To use highRES, a GAMS license is required. See the `GAMS website <https://www.gams.com/latest/docs/UG_License.html>`_ for more information.

Data bundle
------------

To run the full workflow, two datapackages are needed they can be downloaded from:

1. (~80MB compressed, ~300MB uncompressed) https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1.
2. (~10GB) https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1

Windows
----------------
1. Clone the repository
2. Install snakemake

   - Download miniforge windows exe https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe
   - Install Miniforge
   - Run the minimal install of the snakemake environment mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas zstd
3. Activate the snakemake environment
4. Navigate to the repository in your snakemake conda environment shell
5. Get the required input files

.. code-block:: console

   curl -o shared_input.tar -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-   OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1" -o resources.tar.zst -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/   tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1"

6. Extract the required input files

.. code-block:: console

    zstd -d resources.tar.zst
    mkdir resources
    tar xf resources.tar -C resources
    mkdir shared_input
    tar xf shared_input.tar -C shared_input

5. Run snakemake -c --use-conda
