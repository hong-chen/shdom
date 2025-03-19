======================================================================================
Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
======================================================================================

SHDOM was initially developed by `Frank Evans <evans@nit.colorado.edu>`_ and available `here <https://coloradolinux.com/shdom/>`_.

Documents:

* `shdom <https://github.com/hong-chen/shdom/blob/master/docs/shdom.rst>`_

* `propgen <https://github.com/hong-chen/shdom/blob/master/docs/propgen.rst>`_

|

How to Install
--------------

**Dependencies**

Require a ``fortran`` compiler, we recommend `GFortran <https://fortran-lang.org/learn/os_setup/install_gfortran/>`_.

**SHDOM**

1. Open a terminal, type in the following

    .. code-block:: bash

       git clone https://github.com/hong-chen/shdom.git

2. Under newly cloned ``<shdom repo>``, where it contains ``README.rst``, type in the following

    .. code-block:: bash

       cd src

3. Under ``<shdom repo>/src``, where it contains ``makefile``, type in the following

    .. code-block:: bash

       make

   Warning messages are expected during the compilation process, and will not affect the installation of ``shdom``.

   Once the ``shdom`` is successfully compiled, the following executables will be generated

   * ``bin/shdom`` (namelist input)

   * ``bin/shdom_`` (user input)

   * ``bin/propgen``

   * ``bin/make_mie_table``

   * ``bin/make_tmatrix_table``

   * ``bin/plotscattab``

   * ``bin/put``

4. Add the following to the system source file, e.g. ``~/.bashrc`` for ``bash`` (replace ``<shdom repo>`` with a valid path),

    .. code-block:: bash

       export PATH="<shdom repo>/bin:$PATH"

5. Run the case under ``examples/03_bmk_les`` (see detailed instructions there).

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

**Major changes:**

#. changed the format from `txt` to `rst`.

The following is a list of file-to-file reference:

* ``docs/README.rst``: was ``README``

* ``docs/shdom.rst``: was ``shdom.txt``

* ``docs/propgen.rst``: was ``propgen.txt``

* ``HISTORY.rst``: was ``Updatelist.txt``

|

Source Code File Change
-----------------------

**Major changes:**

#. added support for LSRT BRDF model;

#. fixed domain shift issue in 3D radiance calculations;

#. added examples for benchmarking purpose;

#. changed the legacy code from `F77` to `F90`;

#. formatted the code.

|

The following is a list of file-to-file reference:

* ``src/makefile``: was ``makefile``

* ``src/shdom.f90`` and ``src/shdom_.f90``: was ``shdom.f90``

* ``src/shdom_mpi.f90``: was ``shdom_mpi.f90``

* ``src/shdom_nompi.f90``: was ``shdom_nompi.f``

* ``src/shdom_nc.f90``: was ``shdom_netcdf.f90``

* ``src/shdom_nonc.f90``: was ``shdom_nonetcdf.f90``

* ``src/shdSolv_1.f90``: was ``shdomsub1.f``

* ``src/shdSolv_2.f90``: was ``shdomsub2.f``

* ``src/shdUtil_io.f90``: was ``shdomsub3.f``

* ``src/shdSurf_brdf.f90``: separated from ``shdomsub2.f``

* ``src/shdSurf_ocean.f90``: was ``ocean_brdf.f``

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

* ``src/tools/shdTool_txt2nc.f90``: was ``AsciiPropertyFileToNetcdf.f90``

* ``src/tools/shdTool_put.f90``: was ``put.c``

|

Other Resources
---------------

* `SHDOM <https://coloradolinux.com/shdom/>`_: legacy SHDOM

* `AT3D <https://github.com/CloudTomography/AT3D>`_: a Python-interfaced SHDOM tailored towards cloud tomography applications
