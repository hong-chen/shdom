=================
SHDOM Source Code
=================

All of the code here was initially developed by `Frank Evans <evans@nit.colorado.edu>`_ and available `here <https://coloradolinux.com/shdom/>`_.

===================
Installation Tweaks
===================
* For Mac users, add the following to ``~/.bashrc``,

  .. code-block:: bash

     export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

* Add the following to ``~/.bashrc`` (replace ``<shdom repo>`` with a valid path),

  .. code-block:: bash

     export PATH="<shdom repo>/bin:$PATH"


===========
File Change
===========

* ``Makefile``: was ``makefile``

  * added ``-fallow-argument-mismatc`` to ``FFLAGS``

  * added ``clean`` option for ``make``

  * changed the module names accordingly

* ``shdom.f90``: was ``shdom.f90``

  * formatted the code

  * changed the default from ``user_input`` to ``namelist_input``

* ``shdom_mpi.f90``: was ``shdom_mpi.f90``

* ``shdom_netcdf.f90``: was ``shdom_netcdf.f90``

* ``shdom_nompi.f``: was ``shdom_nompi.f``

* ``shdom_nonetcdf.f90``: was ``shdom_nonetcdf.f90``

  * formatted the code

* ``shdSca_indexwatice.f``: was ``indexwatice.f``

* ``shdSca_miewig.f``: was ``miewig.f``

* ``shdSca_tmatrixwig.f``: was ``tmatrixwig.f``

* ``shdUtil_rte1.f``: was ``shdomsub1.f``

* ``shdUtil_rte2.f``: was ``shdomsub2.f``

* ``shdUtil_io.f``: was ``shdomsub3.f``

* ``shdUtil_ocean.f``: was ``ocean_brdf.f``

* ``shdTool_propgen.f90``: was ``propgen.f90``

  * formatted the code

* ``shdTool_make_mie_table.f90``: was ``make_mie_table.f90``

  * formatted the code

* ``shdTool_make_tmatrix_table.f90``: was ``make_tmatrix_table.f90``

* ``shdTool_put.c``: was ``put.c``

  * modified the code to resolve a warning from compilation

* ``shdTool_plotscattab.f90``: was ``plotscattab.f90``

* ``shdTool_txt2nc.c``: was ``AsciiPropertyFileToNetcdf.f90``

* ``_tmd.par.f``: was ``tmd.par.f``

* ``_fftpack.f``: was ``fftpack.f``

* ``_lapack.f``: was ``lpd.f``
