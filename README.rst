======================================================================================
Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
======================================================================================

All of the code here was initially developed by `Frank Evans <evans@nit.colorado.edu>`_ and available `here <https://coloradolinux.com/shdom/>`_.

Installation Tweaks
-------------------

* For Mac users, add the following to ``~/.bashrc``,

  .. code-block:: bash

    export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

* Add the following to ``~/.bashrc`` (replace ``<shdom repo>`` with a valid path),

  .. code-block:: bash

    export PATH="<shdom repo>/bin:$PATH"

Document File Change
--------------------
* ``docs/README.rst``: was ``README``

  * changed to reStructuredText format

* ``docs/shdom.rst``: was ``shdom.txt``

  * changed to reStructuredText format

* ``docs/propgen.rst``: was ``propgen.txt``

  * changed to reStructuredText format

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

* ``src/shdom_nompi.f``: was ``shdom_nompi.f``

* ``src/shdom_nonetcdf.f90``: was ``shdom_nonetcdf.f90``

  * formatted the code

* ``src/shdSca_indexwatice.f``: was ``indexwatice.f``

* ``src/shdSca_miewig.f``: was ``miewig.f``

* ``src/shdSca_tmatrixwig.f``: was ``tmatrixwig.f``

* ``src/shdUtil_rte1.f``: was ``shdomsub1.f``

* ``src/shdUtil_rte2.f``: was ``shdomsub2.f``

* ``src/shdUtil_io.f``: was ``shdomsub3.f``

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


Other Resources
---------------

* `SHDOM <https://coloradolinux.com/shdom/>`_: legacy SHDOM

* `AT3D <https://github.com/CloudTomography/AT3D>`_: a Python-based SHDOM tailored towards cloud tomography applications
