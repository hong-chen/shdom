========
Case LES
========

Please make sure `data/` contains the following files before run any scripts

* ``data/les_lwc.txt``
* ``data/sfc_lambertian.txt``
* ``data/sfc_ocean.txt``
* ``data/shdom-inp_les.nml``


1. **Generate property file**: run the following command

   .. code-block:: bash

      ./01_gen.sh

   Once run successfully, the following files will be generated

   * ``tmp-data/water_w065c_sca.txt``
   * ``tmp-data/les_prp.txt``

2. **Run SHDOM**: run the following command

   .. code-block:: bash

      ./02_run.sh

   Once run successfully, the following file will be generated

   * ``tmp-data/shdom-out_les.txt``

3. **Plot the results** from SHDOM output: run the following command

   .. code-block:: bash

      python 03_plot.py
