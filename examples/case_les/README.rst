========
Case LES
========

Please make sure `data/` contains the following files before run any scripts

* ``data/les.lwc``
* ``data/ocean.sfc``


1. **Generate property file**: run the following command

   .. code-block:: bash

      ./01_gen.sh

   Once run successfully, the following files will be generated

   * ``data/water_w065c.scat``
   * ``data/les.prp``

2. **Run SHDOM**: run the following command

   .. code-block:: bash

      ./02_run.sh

   Once run successfully, the following file will be generated

   * ``shdom-out.txt``

3. **Plot the results** from SHDOM output: run the following command

   .. code-block:: bash

      python 03_plot.py
