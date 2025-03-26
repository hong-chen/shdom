Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://github.com/hong-chen/shdom/actions/workflows/build_gfortran.yml/badge.svg
   :target: https://github.com/hong-chen/shdom/actions/workflows/build_gfortran.yml
   :alt: Build Status

.. image:: https://img.shields.io/badge/dynamic/json?color=blue&label=unique%20clone&query=uniques&url=https://gist.githubusercontent.com/hong-chen/8d66aae71ef31274c9cebbe74c548a82/raw/clone.json&logo=github
    :target: https://github.com/hong-chen/shdom
    :alt: Git Clone Counts

|

SHDOM was initially developed by `Frank Evans <evans@nit.colorado.edu>`_ and available `here <https://coloradolinux.com/shdom/>`_.

**Documents**:

* `shdom <https://github.com/hong-chen/shdom/blob/master/docs/shdom.rst>`_

* `propgen <https://github.com/hong-chen/shdom/blob/master/docs/propgen.rst>`_

|

**Major changes:**

* Added ``CMake`` support (contributed by `Dr. Cheng Da <https://www.cdalab.org>`_);

* Fixed a segmentation error when saving "state" from calculation;

* Changed ASCII output to binary format to reduce file size and accelerate processing speed;

* Added support for LSRT BRDF model;

* Fixed domain shift issue in 3D radiance calculations;

* Added examples (including legacy cases) for benchmarking purpose;

* Changed the legacy code from ``F77`` to ``F90``;

* Formatted the code.

|

How to Install
--------------

If you encountered any issue, please feel free to reach out at `Discord SUPPORT/installation <https://discord.com/channels/681619528945500252/1123343093417119754>`__
for community support.

**Dependencies**

* Require a ``fortran`` compiler, we recommend `GFortran <https://fortran-lang.org/learn/os_setup/install_gfortran/>`_.

* We recommend ``CMake``, but it is optional.

|

**Clone the Repository**

#. Open a terminal, type in the following

    .. code-block:: bash

       git clone https://github.com/hong-chen/shdom.git

|

**Build SHDOM (based on CMake)** (contributed by `Dr. Cheng Da <https://www.cdalab.org>`_)

#. Under newly cloned ``<shdom repo>``, where it contains ``LICENSE``, type in the following

    .. code-block:: bash

       make

    * By default, neither ``netCDF`` nor ``MPI`` is chosen;

    * If using the ``MPI``, run ``make mpi`` instead;

    * If using the ``netCDF``, run ``make nc`` instead. However, you will need
      to modify ``-DNC_INC=<NETCDF_INCLUDE_DIR>`` and ``-DNC_LIBS=<NETCDF_LIBRARY_DIR> -lnetcdff``
      within ``makefile``, for example (current settings)

        .. code-block:: bash

           -DNC_INC="/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/include"
           -DNC_LIBS="-L/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/lib -lnetcdff"

OR **Build SHDOM (based on GNU Make)**

#. Under ``<shdom repo>/src``, where it contains ``_legacy/``, type in the following

    .. code-block:: bash

       make

|

**Post-Installation**

#. Add the following to the system source file, e.g. ``~/.bashrc`` for ``bash`` (replace ``<shdom repo>`` with a valid path),

    .. code-block:: bash

       export PATH="<shdom repo>/build/bin:$PATH"

#. Run the case under ``examples/03_bmk_les`` (see detailed instructions `there <https://github.com/hong-chen/shdom/blob/master/examples/03_bmk_les>`_).

|

Tips
----

* Warning messages are expected during the compilation process, and will not affect the installation of ``shdom``.
  Once the ``shdom`` is successfully compiled, the following executables will be generated

   * ``bin/shdom`` (namelist input)

   * ``bin/shdom_`` (user input)

   * ``bin/propgen``

   * ``bin/make_mie_table``

   * ``bin/make_tmatrix_table``

   * ``bin/plotscattab``

   * ``bin/put``


* For Mac users, add the following to ``~/.bashrc``,

  .. code-block:: bash

     export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

|

File-to-File Reference
----------------------

* ``src/makefile``: was ``makefile``

* ``src/shdom.f90`` and ``src/shdom_.f90``: was ``shdom.f90``

* ``src/shdom_mpi.f90``: was ``shdom_mpi.f90``

* ``src/shdom_nompi.f90``: was ``shdom_nompi.f``

* ``src/shdSolv.f90``: was ``shdomsub1.f``

* ``src/shdComp.f90``: was ``shdomsub2.f``

* ``src/shdIO_inp.f90``: partially from ``shdomsub3.f``

* ``src/shdIO_out.f90``: partially from ``shdomsub3.f``

* ``src/shdSurf_brdf.f90``: partially from ``shdomsub2.f`` and combined with ``ocean_brdf.f``

* ``src/fftpack.f90``: was ``fftpack.f``

|

* ``src/tools/shdTool_propgen.f90``: was ``propgen.f90``

* ``src/tools/shdTool_make_mie_tab.f90``: was ``make_mie_table.f90``

* ``src/tools/shdScat_ind_wat.f90``: was ``indexwatice.f``

* ``src/tools/shdScat_mie_wig.f90``: was ``miewig.f``

* ``src/tools/shdTool_make_tmx_tab.f90``: was ``make_tmatrix_table.f90``

* ``src/tools/shdScat_tmx_wig.f90``: was ``tmatrixwig.f``

* ``src/tools/shdScat_tmx_par.f90``: was ``tmd.par.f``

* ``src/tools/_lapack.f``: was ``lpd.f``

* ``src/tools/shdTool_plot_sca_tab.f90``: was ``plotscattab.f90``

* ``src/tools/shdTool_put.f90``: was ``put.c``

|

* Planning to remove support for ``netCDF``, thus ``netCDF`` related codes are not presented.

|

* ``docs/README.rst``: was ``README``

* ``docs/shdom.rst``: was ``shdom.txt``

* ``docs/propgen.rst``: was ``propgen.txt``

* ``HISTORY.rst``: was ``Updatelist.txt``

|

Other Resources
---------------

* `SHDOM <https://coloradolinux.com/shdom/>`_: legacy SHDOM

* `AT3D <https://github.com/CloudTomography/AT3D>`_: a Python-interfaced SHDOM tailored towards cloud tomography applications
