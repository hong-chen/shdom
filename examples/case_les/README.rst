========
Case LES
========

Gnerate Property File
---------------------

Please make sure `data/` contains the following files before run any scripts

* `data/les.lwc`

* `data/ocean.sfc`


1. Run the following command

   .. code-block:: bash

      ./gen.sh

   Once run successfully, the following files will be generated

   * ``data/water_w065c.scat``
   * ``data/les.prp``

2. Run the following command

   .. code-block:: bash

      ./run.sh

   Once run successfully, the following files will be generated

   * ``shdom-out.txt``

2. Run the following command

   .. code-block:: bash

      python plot.py
