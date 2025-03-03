======================================================================================
Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
======================================================================================

All of the code here was initially developed by `Frank Evans <evans@nit.colorado.edu>`_ and available `here <https://coloradolinux.com/shdom/>`_.

Documents:

* `shdom <https://github.com/hong-chen/shdom/blob/master/docs/shdom.rst>`_

* `propgen <https://github.com/hong-chen/shdom/blob/master/docs/propgen.rst>`_

|

How to Install
--------------

**Dependencies**

1. Require ``fortran`` and ``c`` compilers (``c`` is usually built-in, pick **one** from the following for ``fortran``)

    * `GFortran <https://fortran-lang.org/learn/os_setup/install_gfortran/>`_ (free, more popular);

    * `Intel Fortran <https://docs.conda.io/en/latest/miniconda.html>`_ (costly, faster).

   We recommend ``gfortran`` (`GFortran <https://fortran-lang.org/learn/os_setup/install_gfortran/>`_) and current ``Makefile`` is developed around ``gfortran``.

**SHDOM**

1. Open a terminal, type in the following

    .. code-block:: bash

       git clone https://github.com/hong-chen/shdom.git

2. Under newly cloned ``<shdom repo>``, where it contains ``README.rst``, type in the following

    .. code-block:: bash

       cd src

3. Under ``<shdom repo>/src``, where it contains ``Makefile``, type in the following

    .. code-block:: bash

       make

   Warning messages are expected during the compilation process, and will not affect the installation of ``shdom``.

   Once the ``shdom`` is successfully compiled, the following executables will be generated

   * ``bin/shdom``

   * ``bin/propgen``

   * ``bin/make_mie_table``

   * ``bin/make_tmatrix_table``

   * ``bin/plotscattab``

   * ``bin/put``

4. Add the following to the system source file, e.g. ``~/.bashrc`` for ``bash`` (replace ``<shdom repo>`` with a valid path),

    .. code-block:: bash

       export PATH="<shdom repo>/bin:$PATH"

5. Run the case under ``examples/case_les`` (see detailed instructions there).

|

Tips
----

* For Mac users, add the following to ``~/.bashrc``,

  .. code-block:: bash

     export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

If you encountered any error, please feel free to reach out at `Discord SUPPORT/installation <https://discord.com/channels/681619528945500252/1123343093417119754>`__
for community support.

|

Document File Change
--------------------

* ``docs/README.rst``: was ``README``

  * changed to rST format

* ``docs/shdom.rst``: was ``shdom.txt``

  * changed to rST format

* ``docs/propgen.rst``: was ``propgen.txt``

  * changed to rST format

* ``HISTORY.rst``: was ``Updatelist.txt``

  * changed to rST format

|

Source Code File Change
-----------------------

* ``src/Makefile``: was ``makefile``

  * added ``-fallow-argument-mismatc`` to ``FFLAGS``

  * added ``clean`` option for ``make``

  * changed the module names accordingly

* ``src/shdom.f90``: was ``shdom.f90``

  * formatted the code

  * changed the default from ``user_input`` to ``namelist_input``

* ``src/shdom_mpi.f90``: was ``shdom_mpi.f90``

* ``src/shdom_netcdf.f90``: was ``shdom_netcdf.f90``

* ``src/shdom_nompi.f90``: was ``shdom_nompi.f``

  * converted from F77 to F90
  * formatted the code

* ``src/shdom_nonetcdf.f90``: was ``shdom_nonetcdf.f90``

  * formatted the code

* ``src/shdSca_indexwatice.f``: was ``indexwatice.f``

* ``src/shdSca_miewig.f``: was ``miewig.f``

* ``src/shdSca_tmatrixwig.f``: was ``tmatrixwig.f``

* ``src/shdUtil_rte1.f``: was ``shdomsub1.f``

* ``src/shdUtil_rte2.f``: was ``shdomsub2.f``

* ``src/shdUtil_io.f90``: was ``shdomsub3.f``

  * converted from F77 to F90
  * formatted the code

* ``src/shdUtil_ocean.f``: was ``ocean_brdf.f``

* ``src/shdTool_propgen.f90``: was ``propgen.f90``

  * formatted the code

* ``src/shdTool_make_mie_table.f90``: was ``make_mie_table.f90``

  * formatted the code

* ``src/shdTool_make_tmatrix_table.f90``: was ``make_tmatrix_table.f90``

* ``src/shdTool_put.c``: was ``put.c``

  * modified the code to resolve a warning from compilation

* ``src/shdTool_plotscattab.f90``: was ``plotscattab.f90``

* ``src/shdTool_txt2nc.c``: was ``AsciiPropertyFileToNetcdf.f90``

* ``src/_tmd.par.f``: was ``tmd.par.f``

* ``src/_fftpack.f``: was ``fftpack.f``

* ``src/_lapack.f``: was ``lpd.f``

|

Other Resources
---------------

* `SHDOM <https://coloradolinux.com/shdom/>`_: legacy SHDOM

* `AT3D <https://github.com/CloudTomography/AT3D>`_: a Python-based SHDOM tailored towards cloud tomography applications
