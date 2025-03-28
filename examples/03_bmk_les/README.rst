========
Case LES
========

Please make sure ``<shdom repo>/data/`` and ``data/`` contains the following files before run any scripts

* ``<shdom repo>/data/shdom-lwc_les-fast.txt``
* ``<shdom repo>/data/shdom-sfc_land-lsrt-fast.txt``
* ``data/shdom-inp_les-fast_save.nml``
* ``data/shdom-inp_les-fast_rad.nml``
* ``data/shdom-inp_les-fast_flux.nml``


1. **Generate property file**: run the following command

   .. code-block:: bash

      ./01_gen.sh

   Once run successfully, the following files will be generated

   * ``tmp-data/shdom-sca-065_les-fast.txt``
   * ``tmp-data/shdom-prp_les-fast.txt``

2. **Run SHDOM**: run the following command

   .. code-block:: bash

      ./02_run.sh

   Once run successfully, the following file will be generated

   * ``tmp-data/shdom-out_les-fast_r-all.txt``
   * ``tmp-data/shdom-out_les-fast_r-all.txt.sHdOm-out``
   * ``tmp-data/shdom-out_les-fast_r-nad.txt``
   * ``tmp-data/shdom-out_les-fast_r-nad.txt.sHdOm-out``
   * ``tmp-data/shdom-out_les-fast_f1.txt``
   * ``tmp-data/shdom-out_les-fast_f1.txt.sHdOm-out``
   * ``tmp-data/shdom-out_les-fast_f4.txt``
   * ``tmp-data/shdom-out_les-fast_f4.txt.sHdOm-out``

3. **Plot the results** from SHDOM output: run the following command

   .. code-block:: bash

      python 03_plot.py

   Once run successfully, the following image sequence will be generated

   .. image:: https://github.com/hong-chen/shdom/blob/master/examples/03_bmk_les/data/les-fast.gif
       :width: 400
