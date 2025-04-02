* 2024-04-02 v0.1.1

    * Fixed "dead process issue" when running MPI version with 1D property field input;

    * Fixed domain shift issue in 3D radiance calculations (also for MPI version);

    * Added ``CMake`` support and auto-build workflow on Github (contributed by Dr. Cheng Da);

    * Fixed a segmentation error when saving "state" from calculation;

    * Changed ASCII output to binary format to reduce file size and accelerate processing speed;

    * Added support for LSRT BRDF model;

    * Added examples (including legacy cases) for benchmarking purpose.

* 2024-03-05 v0.1.0

    * Converted the majority of the Fortran 77 code to Fortran 90.

    * Resolved the "domain shifting" issue for 3D radiance calculations (see `PR #1 <https://github.com/hong-chen/shdom/pull/1>`_).

    * Fixed bugs in input/output (I/O) operations.

    * Addressed warnings during code compilation.

    * Benchmarked with legacy test cases:

        - ``run_brdf`` (passed)
        - ``run_mono_les`` (passed)
        - ``run_polarized`` (passed)
        - ``run_visualize`` (passed)
        - ``run_multproc`` (not yet verified)

* 2015-05-13 by Frank Evans

    Fixed bugs concerning CALC_SOURCE_PNT and CALC_SOURCE_PNT_UNPOL.

* 2014-07-10 by Frank Evans

    Fixed a few bugs.
    There are still known bugs that cause crashes when using MPI.

* 2014-03-07 by Frank Evans

    Major upgrade release: beta distribution of polarized SHDOM.
    Parallelized (with MPI) visualization output option.
    Fortran 90 compiler is now required (main program is shdom.f90).
    Fixed bugs in visualization output (camera images flipped in
    sample direction, cross-track images flipped in line).
    Fixed bug in Flux output format 2 that made garbage output.
    Dropped make_ice_table.f90, cloudprp.f, ckdfu.f, ckdswrrtm.f,
    ckdlwrrtm.f, and several old scripts from distribution.

* 2010-05-03 by Frank Evans

    Separated large sw_ice_scatter.db from the distribution tar file.
    Added warning about ocean surface reflection model to ``shdom.txt``.

* 2009-11-15 by Frank Evans

    Added memory parameter input file for multiple processors
    on shared memory computers.
    Fixed several rare MPI bugs.

* 2009-03-09 by Frank Evans

    Major upgrade release (SHDOM input parameter list changed).
    Added ability to run SHDOM on multiple processors with MPI.
    Added option to read property inputs from and write most
    output formats to netcdf files (thanks to Robert Pincus).

* 2008-03-03 by Frank Evans

    Fixed open boundary condition bug in DIRECT_BEAM_PROP.
    (Thanks to Richard Siddans)

* 2005-06-06 by Frank Evans

    Limited SHDOM acceleration parameter to be less than 10; having
    no limit might have caused divergence for optically thick IPA.
    Fixed make_ice_table.f90 database reading allocation bug.
    Fixed lowest level indexing error in ckdlwrrtm.f (27-JAN-2004).

* 2003-07-10 by Frank Evans

    Corrected problem with ice scattering database (sw_ice_scatter.db)
    by rescaling extinction efficiency, single scattering albedo,
    and asymmetry parameter for "delta transmission" of rays
    through parallel ice crystal faces.

* 2003-06-20 by Frank Evans

    Major upgrade release.
    Added new optical property file generation system consisting
    of make_mie_table.f90 and make_ice_table.f90 for creating
    single scattering table files for spherical particles
    and 8 shapes of ice crystals (shortwave only), and
    propgen.f90 for producing property files with arbitrary
    mixtures of particles at each grid point.
    Added modified version of 6S ocean surface reflectance BRDF
    model to SHDOM (wavelength is now a required input parameter).
    Modified RPV BRDF surface reflection function according to
    Alexei Lyapustin for added robustness.
    Added camera and cross-track scanning image "visualization"
    output to SHDOM and example scripts and IDL programs.
    Added longwave and shortwave RRTM based k-distribution programs.

* 1999-07-27 by Frank Evans

    Modified cloudprp to change aerosol input file format (added
    distribution width and index of refraction dependence on wavelength)

* 1999-06-01 by Frank Evans

    Major upgrade release.
    Added Fortran 90 main program (shdom90) with allocatable arrays.
    New version of cloudprp.f.  Different modes for Mie table
    generation. Faster Mie table generation.  Added aerosol
    properties input from file.  Added lognormal size distributions
    (courtesy of Lin Chambers).  Fixed bug for effective radius
    off end of table.   (2/99)
    Included program to plot phase functions in Mie tables.
    Included documentation file for cloudprp and ckdfu.
    Included more example scripts in distribution.
    Included IDL SHDOM file input and visualization procedures.

    Fixed horizontal range of radiance output in open boundary conditions
    Note: this changed the input parameters for radiance output.
    Added input parameter checking subroutine.
    Fixed bugs for large optical depth per grid cell in    (5/98)
    k-distributions by initializing for each k (INRADFLAG=.FALSE.).
    Fixed bug for isotropic source function with NPHI>32.  (1/99)
    Added checks for negative absorption in CKD file, and
    negative phase function from severe Legendre truncation.


* 1997-06-20 by Frank Evans

    Added general surface bidirectional reflectance functions.
    Added open horizontal boundary conditions (new input - BCFLAG).
    Improved adaptive grid cell splitting method (affects SPLITACC).
    Improved response when adaptive grid generation runs out of memory.
    Added code option for outputting rms of higher order radiance.
    Added double gaussian discrete ordinate set option.

* 1996-12-20 by Frank Evans

    Put in MAXRADOUT constant for radiance output array size.
    Fixed heating rate bug for thermal emission with SH truncation

* 1996-08-02 by Frank Evans

    Modified ckdfu.f to output CKD info at atmosphere file levels.

* 1996-07-30 by Frank Evans

    Fixed bug that happened when NPHI was not a good FFT size (44,56,).
    Fixed negative source bug: "Grid point has no contribution".

* 1996-07-11 by Frank Evans

    Put in reading Z base grid from zgrid.inp file (GRIDTYPE='F').
    Allowed delta-M mode for thermal source.
    Referred to online journal article in documentation.

* 1996-05-16 by Frank Evans

    Put in transmission cutoff to speed up optically thick radiances
    Added J (source function) output

* 1996-05-07 by Frank Evans

    Fixed cell integration formula for large cell optical path
    Fixed initialization of l=0 m=0 radiance term
    cloudprp.f: fixed number of size integration steps

* 1996-04-23 by Frank Evans

    Azimuth transform coefficient Nphi>32 bug fixed.
    Echo parameter values in input routine.
    Solar direct beam in mean radiance, net flux output

* 1996-04-10 by Frank Evans

    Added variable Lambertian surface, input surface file,
    changed GNDEMIS to GNDALBEDO.
    Changed output file header.
    Stored grid point sweeping order for speed up.

* 1996-03-29 by Frank Evans

    First beta version
