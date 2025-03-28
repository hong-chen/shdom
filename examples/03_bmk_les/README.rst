========
Case LES
========

Please make sure ``<shdom repo>/data/`` and ``data/`` contains the following files before run any scripts

* ``<shdom repo>/data/shdom-lwc_les.txt``
* ``<shdom repo>/data/shdom-sfc_land-lsrt.txt``
* ``data/shdom-inp_les.nml``


1. **Generate property file**: run the following command

   .. code-block:: bash

      ./01_gen.sh

   Once run successfully, the following files will be generated

   * ``tmp-data/shdom-sca-065_les.txt``
   * ``tmp-data/shdom-prp_les.txt``

2. **Run SHDOM**: run the following command

   .. code-block:: bash

      ./02_run.sh

   Once run successfully, the following file will be generated

   * ``tmp-data/shdom-out_les.txt``

3. **Plot the results** from SHDOM output: run the following command

   .. code-block:: bash

      python 03_plot.py

   Once run successfully, the following image sequence will be generated

   .. image:: https://github.com/hong-chen/shdom/blob/master/examples/03_bmk_les/data/les.gif
       :width: 400
