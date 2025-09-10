Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://github.com/hong-chen/shdom/actions/workflows/build_gfortran.yml/badge.svg?branch=master
   :target: https://github.com/hong-chen/shdom/actions/workflows/build_gfortran.yml
   :alt: Build Status

.. image:: https://img.shields.io/badge/dynamic/json?color=blue&label=unique%20clone&query=uniques&url=https://gist.githubusercontent.com/hong-chen/8d66aae71ef31274c9cebbe74c548a82/raw/clone.json&logo=github
    :target: https://github.com/hong-chen/shdom
    :alt: Git Clone Counts

|

SHDOM was originally developed by `Frank Evans <evans@nit.colorado.edu>`_ and is available `here <https://coloradolinux.com/shdom/>`_.

**Documents**:

* `shdom <https://github.com/hong-chen/shdom/blob/master/docs>`_

* `propgen <https://github.com/hong-chen/shdom/blob/master/docs/tools/propgen.rst>`_

|

**Major changes:**

* Added support for snow BRDF model (by `Vikas Nataraja <https://github.com/vikasnataraja>`_) and for mixed BRDF surface;

* Conducted benchmarking of radiance (by `Yu-Wen Chen <https://lasp.colorado.edu/people/yu-wen-chen/>`_) and irradiance (by `Ken Hirata <https://kenhira.gitlab.io>`_) calculations (SHDOM vs MCARaTS vs libRadtran [DISORT]);

* Removed layer averaging in the z-direction to ensure consistency in model-to-model comparisons;

* Modified ``propgen`` to accept direct inputs for atmospheric extinction coefficients;

* Fixed the "single pixel shift issue" in surface calculations;

* Modularized and restructured the code;

* Fixed the "dead process issue" when running MPI version with 1D property field input;

* Added ``CMake`` support and auto-build workflow on Github (contributed by `Dr. Cheng Da <https://www.cdalab.org>`_);

* Fixed a segmentation error when saving "state" from calculation;

* Replaced floating-point number representation with scientific notation;

* Converted most I/O from ASCII to binary format;

* Added support for LSRT BRDF model;

* Fixed the "domain shift issue" in 3D radiance calculations (also for MPI version);

* Added examples (including legacy cases) for benchmarking purpose;

* Fixed most of the warnings during code compilation;

* Modernized ``F77`` legacy code to ``F90``;

* Formatted the code.

|

How to Install
--------------

If you encounter any issue, please feel free to reach out at `Discord SUPPORT/installation <https://discord.com/channels/681619528945500252/1123343093417119754>`__
for community support.

**Dependencies**

* Require a ``fortran`` compiler, we recommend `GFortran <https://gcc.gnu.org/fortran>`_ (``gfortran``).

* For MPI use, we recommend `OpenMPI <https://www.open-mpi.org>`_ (``mpif90``, optional).

* For code compilation, we recommend `CMake <https://cmake.org>`_ (``cmake``, optional).

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

    * For ``MPI`` use, run ``make mpi`` instead of ``make``.

or **Build SHDOM (based on GNU Make)**

#. Under ``<shdom repo>/src``, where it contains ``_legacy/``, type in the following

    .. code-block:: bash

       make

    * For ``MPI`` use, comment out ``FC = gfortran`` and uncomment ``FC = mpif90`` in ``makefile``, and then run ``make``.

|

**Post-Installation**

#. Add the following to system source file, e.g. ``~/.bashrc`` for ``bash`` (replace ``<shdom repo>`` with the exact path),

    .. code-block:: bash

       export PATH="<shdom repo>/build/bin:$PATH"

#. Run the case under ``examples/03_bmk_les`` (see detailed instructions `there <https://github.com/hong-chen/shdom/blob/master/examples/03_bmk_les>`_).

|

Tips
----

* For Mac users, add the following to ``~/.bashrc``,

  .. code-block:: bash

     export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

* Warning messages are expected during the compilation process, and will not affect the installation of ``shdom``.
  Once the ``shdom`` is successfully compiled, the following executables will be generated under ``<shdom repo>/build/bin/``

   * ``shdom`` (namelist input)

   * ``shdom_`` (user input)

   * ``propgen``

   * ``make_mie_table``

   * ``make_tmatrix_table``

   * ``plotscattab``

   * ``put``

|

File-to-File Reference
----------------------

* ``src/makefile``: Makefile for code compilation and installation

  * was ``makefile``

* ``src/shdom.f90`` and ``src/shdom_.f90``: main programs

  * was ``shdom.f90``

* ``src/shdUtil.f90``: relates to utility functions

  * partially from ``shdomsub2.f``

* ``src/shdSfc.f90``: relates to surface (e.g., BRDF)

  * partially from ``shdomsub2.f`` and combined with ``ocean_brdf.f``

* ``src/shdPha.f90``: relates to phase function

  * partially from ``shdomsub2.f``

* ``src/shdInit.f90``: relates to initialization

  * partially from ``shdomsub2.f``

* ``src/shdSHDO.f90``: relates to SH (spherical harmonic) and DO (discrete ordinate) operations

  * partially from ``shdomsub1.f``

* ``src/shdSrc.f90``: relates to source function

  * partially from ``shdomsub2.f``

* ``src/shdRad.f90``: relates to radiation quantities (radiance/flux)

  * partially from ``shdomsub2.f``

* ``src/shdGrid.f90``: relates to grid/cell/point operations

  * partially from ``shdomsub1.f``

* ``src/shdIntg.f90``: relates to source function integration for solver

  * partially from ``shdomsub1.f``

* ``src/shdSolv.f90``: relates to solver

  * partially from ``shdomsub1.f``

* ``src/shdInp.f90``: relates to I/O input

  * partially from ``shdomsub3.f``

* ``src/shdOut.f90``: relates to I/O output

  * partially from ``shdomsub3.f``

* ``src/shdMPI.f90``: relates to MPI

  * was ``shdom_mpi.f90`` and ``shdom_nompi.f``

* ``src/globals.f90``: relates to global variables

  * partially from ``shdom.f90`` and mostly from ``hparx`` (``MCARaTS`` related)

* ``src/fftpack.f90``: relates to fast Fourier transform

  * was ``fftpack.f``

|

* ``src/tools/shdTool_propgen.f90``:

  * was ``propgen.f90``

* ``src/tools/shdTool_make_mie_tab.f90``:

  * was ``make_mie_table.f90``

* ``src/tools/shdScat_ind_wat.f90``:

  * was ``indexwatice.f``

* ``src/tools/shdScat_mie_wig.f90``:

  * was ``miewig.f``

* ``src/tools/shdTool_make_tmx_tab.f90``:

  * was ``make_tmatrix_table.f90``

* ``src/tools/shdScat_tmx_wig.f90``:

  * was ``tmatrixwig.f``

* ``src/tools/shdScat_tmx_par.f90``:

  * was ``tmd.par.f``

* ``src/tools/_lapack.f``:

  * was ``lpd.f``

* ``src/tools/shdTool_plot_sca_tab.f90``:

  * was ``plotscattab.f90``

* ``src/tools/shdTool_put.f90``:

  * was ``put.c``

|

* Planning to remove support for ``netCDF``, thus ``netCDF`` related codes are not presented.

|

* ``docs/README.rst``:

  * was ``shdom.txt``

* ``docs/propgen.rst``:

  * was ``propgen.txt``

* ``HISTORY.rst``:

  * was ``Updatelist.txt``

|

Other Resources
---------------

* `SHDOM <https://coloradolinux.com/shdom/>`_: legacy SHDOM

* `AT3D <https://github.com/CloudTomography/AT3D>`_: a Python-interfaced SHDOM tailored towards cloud tomography applications
