!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!
!
!   Copyright (C) 2025-
!   Hong Chen, hong.chen-1@colorado.edu
!
!   SHDOMNG is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   SHDOMNG is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SHDOMNG. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
!
! Acknowledgement:
! SHDOM was originally developed by Frank Evans and is available at https://coloradolinux.com/shdom.
!╰────────────────────────────────────────────────────────────────────────────╯!

! Computes the light scattering properties of a single randomly
! oriented spheroid or cylinder using the T-matrix code of M. Mishc
! Spheroids are calculated with cylflag=.false. and cylinders with
! cylflag=.true.  The axis ratio (axratio) is > 1 for an oblate
! spheroid and < 1 for a prolate spheroid.   For cylinders axratio
! is the ratio of the diameter to the length.   The radius is
! the equal volume radius. The imaginary part of the complex index
! of refraction specified in mindex can have either sign, as it is
! set to positive when used.
! The Wigner d-function coefficients of the six unique elements of
! the scattering phase matrix are returned in wigcoef(6,0:maxleg) i
! order alpha1, alpha2, alpha3, alpha4, beta1, beta2.
! The inputs and outputs to this subroutine are in REAL*4, but
! the T-matrix calculations are done in REAL*8.
!
! This interface was written by Frank Evans  (May 2011)
!   Modified for Wigner d-functions  (January 2014)
!
! Compile: pgf77 -c tmatrixdwig.f lpd.f  (shdScat_tmx_par.f90 include file
!
! The main code was modified to put in the subroutine interface at
! the beginning, set parameters for a single particle, and calculat
! the Legendre coefficients of the six unique elements of the phase
! matrix (first make Gaussian quadrature angles, call MATR to get t
! scattering matrix at these angles, and then do the integration ov
! Wigner functions).  The MATR subroutine was modified to output th
! scattering matrix at specified angles.


!   Comments below are from the original tmd.lp.f code.

!   New release including the LAPACK matrix inversion procedure.
!   We thank Cory Davis (University of Edinburgh) for pointing
!   out the possibility of replacing the proprietary NAG matrix
!   inversion routine by the public-domain LAPACK equivalent.

!   CALCULATION OF LIGHT SCATTERING BY POLYDISPERSE, RANDOMLY
!   ORIENTED PARTICLES OF IDENTICAL AXIALLY SYMMETRIC SHAPE

!   This version of the code uses DOUBLE PRECISION variables
!   and must be used along with the accompanying files shdScat_tmx_par.f90
!   and lpd.f.

!   Last update 08/06/2005

!   The code has been developed by Michael Mishchenko at the NASA
!   Goddard Institute for Space Studies, New York. This research
!   was funded by the NASA Radiation Sciences Program.

!   The code can be used without limitations in any not-for-
!   profit scientific research.  We only request that in any
!   publication using the code the source of the code be acknowledged
!   and relevant references (see below) be made.

!   This version of the code is applicable to spheroids,
!   Chebyshev particles, and finite circular cylinders.

!   The computational method is based on the Watermsn's T-matrix
!   approach and is described in detail in the following papers:
!
!   1.  M. I. Mishchenko, Light scattering by randomly oriented
!       axially symmetric particles, J. Opt. Soc. Am. A,
!       vol. 8, 871-882 (1991).
!
!   2.  M. I. Mishchenko, Light scattering by size-shape
!       distributions of randomly oriented axially symmetric
!       particles of a size comparable to a wavelength,
!       Appl. Opt., vol. 32, 4652-4666 (1993).
!
!   3.  M. I. Mishchenko and L. D. Travis, T-matrix computations
!       of light scattering by large spheroidal particles,
!       Opt. Commun., vol. 109, 16-21 (1994).
!
!   4.  M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering
!       of light by polydisperse, randomly oriented, finite
!       circular cylinders, Appl. Opt., vol. 35, 4927-4940 (1996).
!
!   5.  D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson,
!       Improved T-matrix computations for large, nonabsorbing and
!       weakly absorbing nonspherical particles and comparison
!       with geometrical optics approximation, Appl. Opt., vol. 36,
!       4305-4313 (1997).
!
!   A general review of the T-matrix approach can be found in
!
!   6.  M. I. Mishchenko, L. D. Travis, and D. W. Mackowski,
!       T-matrix computations of light scattering by nonspherical
!       particles: a review, J. Quant. Spectrosc. Radiat.
!       Transfer, vol. 55, 535-575 (1996).
!
!   The following paper provides a detailed user guide to the
!   T-matrix code:
!
!   7.  M. I. Mishchenko and L. D. Travis, Capabilities and
!       limitations of a current FORTRAN implementation of the
!       T-matrix method for randomly oriented, rotationally
!       symmetric scatterers, J. Quant. Spectrosc. Radiat. Transfer,
!       vol. 60, 309-324 (1998).
!
!   These papers are available in the .pdf format at the web site
!
!   http://www.giss.nasa.gov/~crmim/publications/
!
!   or in hardcopy upon request from Michael Mishchenko
!   Please e-mail your request to crmim@giss.nasa.gov.
!
!   A comprehensive book "Scattering, Absorption, and Emission of
!   Light by Small Particles" (Cambridge University Press, Cambridge,
!   2002) is also available in the .pdf format at the web site
!
!   http://www.giss.nasa.gov/~crmim/books.html

!   Analytical averaging over particle orientations (Ref. 1) makes
!   this method the fastest exact technique currently available.
!   The use of an automatic convergence procedure
!   (Ref. 2) makes the code convenient in massive computations.
!   Ref. 4 describes features specific for finite cylinders as
!   particles with sharp rectangular edges.  Ref. 5 describes further
!   numerical improvements.

!   The use of extended precision variables (Ref. 3) can
!   significantly increase the maximal convergent equivalent-sphere
!   size parameter and make it greater than 200 (depending on
!   refractive index and aspect ratio). The extended-precision code
!   is also available. However, the use of extended precision varibales
!   results in a greater consumption of CPU time.
!   On IBM RISC workstations, that code is approximately
!   five times slower than this double-precision code.  The
!   CPU time difference between the double-precision and extended-
!   precision codes can be larger on supercomputers.

!   This is the first part of the full T-matrix code.  The second part,
!   lpd.f, is completely independent of the firsti part. It contains no
!   T-matrix-specific subroutines and can be compiled separately.
!   The second part of the code replaces the previously implemented
!   standard matrix inversion scheme based on Gaussian elimination
!   by a scheme based on the LU factorization technique.
!   As described in Ref. 5 above, the use of the LU factorization is
!   especially beneficial for nonabsorbing or weakly absorbing particles
!   In this code we use the LAPACK implementation of the LU factorizatio
!   scheme. LAPACK stands for Linear Algebra PACKage. The latter is
!   publicly available at the following internet site:
!
!   http://www.netlib.org/lapack/

!   INPUT PARAMETERS:
!
!      RAT = 1 - particle size is specified in terms of the
!                equal-volume-sphere radius
!      RAT.NE.1 - particle size is specified in terms of the
!                equal-surface-area-sphere radius
!      NDISTR specifies the distribution of equivalent-sphere radii
!      NDISTR = 1 - modified gamma distribution
!           [Eq. (40) of Ref. 7]
!               AXI=alpha
!               B=r_c
!               GAM=gamma
!      NDISTR = 2 - log-normal distribution
!           [Eq. 41) of Ref. 7]
!               AXI=r_g
!               B=[ln(sigma_g)]**2
!      NDISTR = 3 - power law distribution
!           [Eq. (42) of Ref. 7]
!                AXI=r_eff (effective radius)
!                B=v_eff (effective variance)
!                Parameters R1 and R2 (see below) are calculated
!                automatically for given AXI and B
!      NDISTR = 4 - gamma distribution
!           [Eq. (39) of Ref. 7]
!                AXI=a
!                B=b
!      NDISTR = 5 - modified power law distribution
!         [Eq. (24) in M. I. Mishchenko et al.,
!         Bidirectional reflectance of flat,
!         optically thick particulate laters: an efficient radiative
!         transfer solution and applications to snow and soil surfaces,
!         J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999)
!                B=alpha
!
!      The code computes NPNAX size distributions of the same type
!      and with the same values of B and GAM in one run.
!      The parameter AXI varies from AXMAX to AXMAX/NPNAX in steps of
!      AXMAX/NPNAX.  To compute a single size distribution, use
!      NPNAX=1 and AXMAX equal to AXI of this size distribution.
!
!      R1 and R2 - minimum and maximum equivalent-sphere radii
!           in the size distribution. They are calculated automatically
!           for the power law distribution with given AXI and B
!           but must be specified for other distributions
!           after the lines
!
!             DO 600 IAX=1,NPNAX
!                AXI=AXMAX-DAX*DFLOAT(IAX-1)
!
!           in the main program.
!           For the modified power law distribution (NDISTR=5), the
!           minimum radius is 0, R2 is the maximum radius,
!           and R1 is the intermediate radius at which the
!           n(r)=const dependence is replaced by the power law
!           dependence.
!
!      NKMAX.LE.988 is such that NKMAX+2 is the
!           number of Gaussian quadrature points used in
!           integrating over the size distribution for particles with
!           AXI=AXMAX.  For particles with AXI=AXMAX-AXMAX/NPNAX,
!           AXMAX-2*AXMAX/NPNAX, etc. the number of Gaussian points
!           linearly decreases.
!           For the modified power law distribution, the number
!           of integration points on the interval [0,R1] is also
!           equal to NKMAX.
!
!      LAM - wavelength of light
!      MRR and MRI - real and imaginary parts of the refractive
!                  index (MRI.GE.0)
!      EPS and NP - specify the shape of the particles.
!             For spheroids NP=-1 and EPS is the ratio of the
!                 horizontal to rotational axes.  EPS is larger than
!                 1 for oblate spheroids and smaller than 1 for
!                 prolate spheroids.
!             For cylinders NP=-2 and EPS is the ratio of the
!                 diameter to the length.
!             For Chebyshev particles NP must be positive and
!                 is the degree of the Chebyshev polynomial, while
!                 EPS is the deformation parameter
!                 [Eq. (33) of Ref. 7].
!      DDELT - accuracy of the computations
!      NPNA - number of equidistant scattering angles (from 0
!             to 180 deg) for which the scattering matrix is
!             calculated.
!      NDGS - parameter controlling the number of division points
!             in computing integrals over the particle surface.
!             For compact particles, the recommended value is 2.
!             For highly aspherical particles larger values (3, 4,...)
!             may be necessary to obtain convergence.
!             The code does not check convergence over this parameter.
!             Therefore, control comparisons of results obtained with
!             different NDGS-values are recommended.

!   OUTPUT PARAMETERS:
!
!      REFF and VEFF - effective radius and effective variance of
!          the size distribution as defined by Eqs. (43)-(45) of
!          Ref. 7.
!      CEXT - extinction cross section per particle
!      CSCA - scattering cross section per particle
!      W - single scattering albedo
!      <cos> - asymmetry parameter of the phase function
!      ALPHA1,...,BETA2 - coefficients appearing in the expansions
!          of the elements of the scattering matrix in
!          generalized spherical functions
!          [Eqs. (11)-(16) of Ref. 7].
!      F11,...,F44 - elements of the normalized scattering matrix [as
!          defined by Eqs. (1)-(3) of Ref. 7] versus scattering angle

!   Note that LAM, r_c, r_g, r_eff, a, R1, and R2 must
!   be given in the same units of length, and that
!   the dimension of CEXT and CSCA is that of LAM squared (e.g., if
!   LAM and AXI are given in microns, then CEXT and CSCA are
!   calculated in square microns).

!   The physical correctness of the computed results is tested using
!   the general inequalities derived by van der Mee and Hovenier,
!   Astron. Astrophys., vol. 228, 559-568 (1990).  Although
!   the message that the test of van der Mee and Hovenier is satisfied
!   does not guarantee that the results are absolutely correct,
!   the message that the test is not satisfied can mean that something
!   is wrong.

!   The convergence of the T-matrix method for particles with
!   different sizes, refractive indices, and aspect ratios can be
!   dramatically different.  Usually, large sizes and large aspect
!   ratios cause problems.  The user of this code
!   should first experiment with different input parameters in
!   order to get an idea of the range of applicability of this
!   technique.  Sometimes decreasing the aspect ratio
!   from 3 to 2 can increase the maximum convergent equivalent-
!   sphere size parameter by a factor of several (Ref. 7).
!   The CPU time required rapidly increases with increasing ratio
!   radius/wavelength and/or with increasing particle asphericity.
!   This should be taken into account in planning massive computations.
!   Using an optimizing compiler on IBM RISC workstations saves
!   about 70% of CPU time.

!   Execution can be automatically terminated if dimensions of certain
!   arrays are not big enough or if the convergence procedure decides
!   that the accuracy of double precision variables is insufficient
!   to obtain a converged T-matrix solution for given particles.
!   In all cases, a message appears explaining the cause of termination.

!   The message
!        "WARNING:  W IS GREATER THAN 1"
!   means that the single-scattering albedo exceeds the maximum
!   possible value 1.  If W is greater than 1 by more than
!   DDELT, this message can be an indication of numerical
!   instability caused by extreme values of particle parameters.

!   The message "WARNING: NGAUSS=NPNG1" means that convergence over
!   the parameter NG (see Ref. 2) cannot be obtained for the NPNG1
!   value specified in the PARAMETER statement in the file shdScat_tmx_par.f90
!   Often this is not a serious problem, especially for compact
!   particles.

!   Larger and/or more aspherical particles may require larger
!   values of the parameters NPN1, NPN4, and NPNG1 in the file
!   shdScat_tmx_par.f90.  It is recommended to keep NPN1=NPN4+25 and
!   NPNG1=3*NPN1.  Note that the memory requirement increases
!   as the third power of NPN4. If the memory of
!   a computer is too small to accomodate the code in its current
!   setting, the parameters NPN1, NPN4, and NPNG1 should be
!   decreased. However, this will decrease the maximum size parameter
!   that can be handled by the code.

!   In some cases any increases of NPN1 will not make the T-matrix
!   computations convergent.  This means that the particle is just
!   too "bad" (extreme size parameter and/or extreme aspect ratio
!   and/or extreme refractive index; see Ref. 7).
!   The main program contains several PRINT statements which are
!   currently commentd out.  If uncommented, these statements will
!   produce numbers which show the convergence rate and can be
!   used to determine whether T-matrix computations for given particle
!   parameters will converge at all.

!   Some of the common blocks are used to save memory rather than
!   to transfer data.  Therefore, if a compiler produces a warning
!   message that the lengths of a common block are different in
!   different subroutines, this is not a real problem.

!   The recommended value of DDELT is 0.001.  For bigger values,
!   false convergence can be obtained.

!   In computations for spheres use EPS=1.000001 instead of EPS=1.
!   The use of EPS=1 can cause overflows in some rare cases.

!   To calculate a monodisperse particle, use the options
!        NPNAX=1
!        AXMAX=R
!        B=1D-1
!        NKMAX=-1
!        NDISTR=4
!        ...
!        DO 600 IAX=1,NPNAX
!           AXI=AXMAX-DAX*DFLOAT(IAX-1)
!           R1=0.9999999*AXI
!           R2=1.0000001*AXI
!        ...
!   where R is the equivalent-sphere radius.

!   It is recommended to use the power law rather than the
!   gamma size distribution, because in this case convergent solution
!   can be obtained for larger REFF and VEFF assuming the same
!   maximal R2 (Mishchenko and Travis, Appl. Opt., vol. 33, 7206-7225,
!   1994).

!   For some compilers, DACOS must be raplaced by DARCOS and DASIN
!   by DARSIN.

!   If many different size distributions are computed and the
!   refractive index is fixed, then another approach can be more
!   efficient than running this code many times.  Specifically,
!   scattering results should be computed for monodisperse particles
!   with sizes ranging from essentially zero to some maximum value
!   with a small step size (say, 0.02 microns).  These results
!   should be stored on disk and can be used along with spline
!   interpolation to compute scattering by particles with intermediate
!   sizes.  Scattering patterns for monodisperse nonspherical
!   particles in random orientation are (much) smoother than for
!   monodisperse spheres, and spline interpolation usually gives good
!   results. In this way, averaging over any size distribution is a
!   matter of seconds.  For more on size averaging, see Refs. 2 and 4.

!   We would highly appreciate informing me of any problems encountered
!   with this code.  Please send your message to the following
!   e-mail address:  CRMIM@GISS.NASA.GOV.

!   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
!   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
!   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
!   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
!   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
!   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM.

subroutine TMATRIXWIG(wavelength, mindex, radius,&
      & cylflag, axratio, maxrank,&
      & extcross, scatcross, ssalbedo, asymmetry,&
      & nrank, wigcoef)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      integer maxrank, nrank
      logical cylflag
      real wavelength, radius, axratio
      complex*8 mindex
      real extcross, scatcross, ssalbedo, asymmetry

      real wigcoef(6, 0:maxrank)

      real*8 LAM, MRR, MRI, X(NPNG2), W(NPNG2), S(NPNG2), SS(NPNG2),&
          & AN(NPN1), R(NPNG2), DR(NPNG2),&
          & DDR(NPNG2), DRR(NPNG2), DRI(NPNG2), ANN(NPN1, NPN1)
      real*8 XG(1000), WG(1000), TR1(NPN2, NPN2), TI1(NPN2, NPN2),&
          & ALPH1(NPL), ALPH2(NPL), ALPH3(NPL), ALPH4(NPL), BET1(NPL),&
          & BET2(NPL), XG1(2000), WG1(2000),&
          & AL1(NPL), AL2(NPL), AL3(NPL), AL4(NPL), BE1(NPL), BE2(NPL)
      real*8&
          & RT11(NPN6, NPN4, NPN4), RT12(NPN6, NPN4, NPN4),&
          & RT21(NPN6, NPN4, NPN4), RT22(NPN6, NPN4, NPN4),&
          & IT11(NPN6, NPN4, NPN4), IT12(NPN6, NPN4, NPN4),&
          & IT21(NPN6, NPN4, NPN4), IT22(NPN6, NPN4, NPN4)
      parameter(maxnang=500)
      real*8 muang(maxnang), wts(maxnang), scatmat(6, maxnang)
      real*8 coef(6, 0:NPL), backscatvec(6)
      real*8 D00(0:maxnang), D20(0:maxnang)
      real*8 D22P(0:maxnang), D22M(0:maxnang)
      real*8 A1, A2, A3, A4, B1, B2, F

      common /CT/ TR1, TI1
      common /TMAT/ RT11, RT12, RT21, RT22, IT11, IT12, IT21, IT22
      P = DACOS(-1d0)

      !  OPEN FILES *******************************************************

      ! open (6,FILE='test')
      ! open (10,FILE='tmatr.write')

      extcross = 0.0
      scatcross = 0.0
      ssalbedo = 0.0
      asymmetry = 0.0

      !  INPUT DATA ********************************************************
      LAM = wavelength
      MRR = real(mindex)
      MRI = abs(aimag(mindex))
      ! NP=-1 for spheroids, NP=-2 for cylinders
      if (cylflag) then
          NP = -2
      else
          NP = -1
      end if
      EPS = axratio

      RAT = 1.0d0
      NDISTR = 4
      AXMAX = radius
      NPNAX = 1
      B = 0.1
      GAM = 0.5d0
      NKMAX = -1
      DDELT = 0.001d0
      NPNA = 19
      NDGS = 2

      NCHECK = 0
      if (NP .eq. -1 .or. NP .eq. -2) NCHECK = 1
      if (NP .gt. 0 .and. (-1)**NP .eq. 1) NCHECK = 1
      ! write (6,5454) NCHECK
5454  format('NCHECK=', I1)
      DAX = AXMAX/NPNAX
      if (DABS(RAT - 1d0) .gt. 1d-8 .and. NP .eq. -1) call SAREA(EPS, RAT)
      if (DABS(RAT - 1d0) .gt. 1d-8 .and. NP .ge. 0) call SURFCH(NP, EPS, RAT)
      if (DABS(RAT - 1d0) .gt. 1d-8 .and. NP .eq. -2) call SAREAC(EPS, RAT)
      ! print 8000, RAT
8000  format('RAT=', F8.6)
      ! if(NP.eq.-1.and.EPS.ge.1D0) print 7000,EPS
      ! if(NP.eq.-1.and.EPS.lt.1D0) print 7001,EPS
      if (NP .ge. 0) print 7100, NP, EPS
      ! if(NP.eq.-2.and.EPS.ge.1D0) print 7150,EPS
      ! if(NP.eq.-2.and.EPS.lt.1D0) print 7151,EPS
      ! print 7400, LAM,MRR,MRI
      ! print 7200, DDELT
7000  format('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=', F11.7)
7001  format('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=', F11.7)
7100  format('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',&
          & I1, '(', F5.2, ')')
7150  format('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=', F11.7)
7151  format('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=', F11.7)
7200  format('ACCURACY OF COMPUTATIONS DDELT = ', D8.2)
7400  format('LAM=', F10.6, 3x, 'MRR=', D10.4, 3x, 'MRI=', D10.4)
      DDELT = 0.1d0*DDELT

      do 600 IAX = 1, NPNAX
          AXI = AXMAX - DAX*DFLOAT(IAX - 1)
          R1 = 0.9999999*AXI
          R2 = 1.0000001*AXI
          NK = IDINT(AXI*NKMAX/AXMAX + 2)
          if (NK .gt. 1000) print 8001, NK
          if (NK .gt. 1000) stop
          if (NDISTR .eq. 3) call POWER(AXI, B, R1, R2)
8001      format('NK=', I4, ' I.E., IS GREATER THAN 1000. ',&
              & 'EXECUTION TERMINATED.')
          call GAUSS(NK, 0, 0, XG, WG)
          Z1 = (R2 - R1)*0.5d0
          Z2 = (R1 + R2)*0.5d0
          Z3 = R1*0.5d0
          if (NDISTR .eq. 5) GO TO 3
          do I = 1, NK
              XG1(I) = Z1*XG(I) + Z2
              WG1(I) = WG(I)*Z1
          end do
          GO TO 4
3         do I = 1, NK
              XG1(I) = Z3*XG(I) + Z3
              WG1(I) = WG(I)*Z3
          end do
          do I = NK + 1, 2*NK
              II = I - NK
              XG1(I) = Z1*XG(II) + Z2
              WG1(I) = WG(II)*Z1
          end do
          NK = NK*2
4         call DISTRB(NK, XG1, WG1, NDISTR, AXI, B, GAM, R1, R2,&
              & REFF, VEFF, P)
          ! print 8002,R1,R2
8002      format('R1=', F10.6, '   R2=', F10.6)
          ! if (dabs(RAT-1D0).le.1D-6) print 8003, REFF,VEFF
          ! if (dabs(RAT-1D0).gt.1D-6) print 8004, REFF,VEFF
8003      format('EQUAL-VOLUME-SPHERE REFF=', F8.4, '   VEFF=', F7.4)
8004      format('EQUAL-SURFACE-AREA-SPHERE REFF=', F8.4,&
              & '   VEFF=', F7.4)
          ! print 7250,NK
7250      format('NUMBER OF GAUSSIAN QUADRATURE POINTS ',&
              & 'IN SIZE AVERAGING =', I4)

          do I = 1, NPL
              ALPH1(I) = 0d0
              ALPH2(I) = 0d0
              ALPH3(I) = 0d0
              ALPH4(I) = 0d0
              BET1(I) = 0d0
              BET2(I) = 0d0
          end do
          CSCAT = 0d0
          CEXTIN = 0d0
          L1MAX = 0
          do 500 INK = 1, NK
              I = NK - INK + 1
              A = RAT*XG1(I)
              XEV = 2d0*P*A/LAM
              IXXX = XEV + 4.05d0*XEV**0.333333d0
              INM1 = MAX0(4, IXXX)
              if (INM1 .ge. NPN1) then
                  print 7333, NPN1
                  print '(A,F9.4,A,F5.2,A,2F8.4)',&
                      & 'radius=', radius, '  axratio=', axratio, '  m=', mindex
                  return
              end if
7333          format('CONVERGENCE IS NOT OBTAINED FOR NPN1=', I3,&
                  & '.  EXECUTION TERMINATED')
              QEXT1 = 0d0
              QSCA1 = 0d0
              do 50 NMA = INM1, NPN1
                  NMAX = NMA
                  MMAX = 1
                  NGAUSS = NMAX*NDGS
                  if (NGAUSS .gt. NPNG1) print 7340, NGAUSS
                  if (NGAUSS .gt. NPNG1) stop
7340              format('NGAUSS =', I3, ' I.E. IS GREATER THAN NPNG1.',&
                      & '  EXECUTION TERMINATED')
7334              format(' NMAX =', I3, '  DC2=', D8.2, '   DC1=', D8.2)
7335              format('                              NMAX1 =', I3, '  DC2=', D8.2,&
                      & '  DC1=', D8.2)
                  call CONST(NGAUSS, NMAX, MMAX, P, X, W, AN, ANN, S, SS, NP, EPS)
                  call VARY(LAM, MRR, MRI, A, EPS, NP, NGAUSS, X, P, PPI, PIR, PII, R,&
                      & DR, DDR, DRR, DRI, NMAX)
                  call TMATR0(NGAUSS, X, W, AN, ANN, S, SS, PPI, PIR, PII, R, DR,&
                      & DDR, DRR, DRI, NMAX, NCHECK)
                  QEXT = 0d0
                  QSCA = 0d0
                  do N = 1, NMAX
                      N1 = N + NMAX
                      TR1NN = TR1(N, N)
                      TI1NN = TI1(N, N)
                      TR1NN1 = TR1(N1, N1)
                      TI1NN1 = TI1(N1, N1)
                      DN1 = DFLOAT(2*N + 1)
                      QSCA = QSCA + DN1*(TR1NN*TR1NN + TI1NN*TI1NN&
                          & + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
                      QEXT = QEXT + (TR1NN + TR1NN1)*DN1
                  end do
                  DSCA = DABS((QSCA1 - QSCA)/QSCA)
                  DEXT = DABS((QEXT1 - QEXT)/QEXT)
                  ! print 7334, NMAX,DSCA,DEXT
                  QEXT1 = QEXT
                  QSCA1 = QSCA
                  NMIN = DFLOAT(NMAX)/2d0 + 1d0
                  do 10 N = NMIN, NMAX
                      N1 = N + NMAX
                      TR1NN = TR1(N, N)
                      TI1NN = TI1(N, N)
                      TR1NN1 = TR1(N1, N1)
                      TI1NN1 = TI1(N1, N1)
                      DN1 = DFLOAT(2*N + 1)
                      DQSCA = DN1*(TR1NN*TR1NN + TI1NN*TI1NN&
                          & + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
                      DQEXT = (TR1NN + TR1NN1)*DN1
                      DQSCA = DABS(DQSCA/QSCA)
                      DQEXT = DABS(DQEXT/QEXT)
                      NMAX1 = N
                      if (DQSCA .le. DDELT .and. DQEXT .le. DDELT) GO TO 12
10                continue
12                continue
                  ! PRINT 7335, NMAX1,DQSCA,DQEXT
                  if (DSCA .le. DDELT .and. DEXT .le. DDELT) GO TO 55
                  if (NMA .eq. NPN1) then
                      print 7333, NPN1
                      print '(A,F9.4,A,F5.2,A,2F8.4)',&
                          & 'radius=', radius, '  axratio=', axratio, '  m=', mindex
                      return
                  end if
50            continue
55            NNNGGG = NGAUSS + 1
              if (NGAUSS .eq. NPNG1) print 7336
              MMAX = NMAX1
              do 150 NGAUS = NNNGGG, NPNG1
                  NGAUSS = NGAUS
                  NGGG = 2*NGAUSS
7336              format('WARNING: NGAUSS=NPNG1')
7337              format(' NG=', I3, '  DC2=', D8.2, '   DC1=', D8.2)
                  call CONST(NGAUSS, NMAX, MMAX, P, X, W, AN, ANN, S, SS, NP, EPS)
                  call VARY(LAM, MRR, MRI, A, EPS, NP, NGAUSS, X, P, PPI, PIR, PII, R,&
                      & DR, DDR, DRR, DRI, NMAX)
                  call TMATR0(NGAUSS, X, W, AN, ANN, S, SS, PPI, PIR, PII, R, DR,&
                      & DDR, DRR, DRI, NMAX, NCHECK)
                  QEXT = 0d0
                  QSCA = 0d0
                  do 104 N = 1, NMAX
                      N1 = N + NMAX
                      TR1NN = TR1(N, N)
                      TI1NN = TI1(N, N)
                      TR1NN1 = TR1(N1, N1)
                      TI1NN1 = TI1(N1, N1)
                      DN1 = DFLOAT(2*N + 1)
                      QSCA = QSCA + DN1*(TR1NN*TR1NN + TI1NN*TI1NN&
                          & + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
                      QEXT = QEXT + (TR1NN + TR1NN1)*DN1
104               continue
                  DSCA = DABS((QSCA1 - QSCA)/QSCA)
                  DEXT = DABS((QEXT1 - QEXT)/QEXT)
                  ! PRINT 7337, NGGG,DSCA,DEXT
                  QEXT1 = QEXT
                  QSCA1 = QSCA
                  if (DSCA .le. DDELT .and. DEXT .le. DDELT) GO TO 155
                  if (NGAUS .eq. NPNG1) print 7336
150           continue
155           continue
              QSCA = 0d0
              QEXT = 0d0
              NNM = NMAX*2
              do 204 N = 1, NNM
                  QEXT = QEXT + TR1(N, N)
204           continue
              if (NMAX1 .gt. NPN4) print 7550, NMAX1
7550          format('NMAX1 = ', I3, ', i.e. greater than NPN4.',&
                  & ' Execution terminated')
              if (NMAX1 .gt. NPN4) stop
              do 213 N2 = 1, NMAX1
                  NN2 = N2 + NMAX
                  do 313 N1 = 1, NMAX1
                      NN1 = N1 + NMAX
                      ZZ1 = TR1(N1, N2)
                      RT11(1, N1, N2) = ZZ1
                      ZZ2 = TI1(N1, N2)
                      IT11(1, N1, N2) = ZZ2
                      ZZ3 = TR1(N1, NN2)
                      RT12(1, N1, N2) = ZZ3
                      ZZ4 = TI1(N1, NN2)
                      IT12(1, N1, N2) = ZZ4
                      ZZ5 = TR1(NN1, N2)
                      RT21(1, N1, N2) = ZZ5
                      ZZ6 = TI1(NN1, N2)
                      IT21(1, N1, N2) = ZZ6
                      ZZ7 = TR1(NN1, NN2)
                      RT22(1, N1, N2) = ZZ7
                      ZZ8 = TI1(NN1, NN2)
                      IT22(1, N1, N2) = ZZ8
                      QSCA = QSCA + ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4&
                          & + ZZ5*ZZ5 + ZZ6*ZZ6 + ZZ7*ZZ7 + ZZ8*ZZ8
313               continue
213           continue
              ! PRINT 7800,0,DABS(QEXT),QSCA,NMAX
              do 220 M = 1, NMAX1
                  call TMATR(M, NGAUSS, X, W, AN, ANN, S, SS, PPI, PIR, PII, R, DR,&
                      & DDR, DRR, DRI, NMAX, NCHECK)
                  NM = NMAX - M + 1
                  NM1 = NMAX1 - M + 1
                  M1 = M + 1
                  QSC = 0d0
                  do 214 N2 = 1, NM1
                      NN2 = N2 + M - 1
                      N22 = N2 + NM
                      do 314 N1 = 1, NM1
                          NN1 = N1 + M - 1
                          N11 = N1 + NM
                          ZZ1 = TR1(N1, N2)
                          RT11(M1, NN1, NN2) = ZZ1
                          ZZ2 = TI1(N1, N2)
                          IT11(M1, NN1, NN2) = ZZ2
                          ZZ3 = TR1(N1, N22)
                          RT12(M1, NN1, NN2) = ZZ3
                          ZZ4 = TI1(N1, N22)
                          IT12(M1, NN1, NN2) = ZZ4
                          ZZ5 = TR1(N11, N2)
                          RT21(M1, NN1, NN2) = ZZ5
                          ZZ6 = TI1(N11, N2)
                          IT21(M1, NN1, NN2) = ZZ6
                          ZZ7 = TR1(N11, N22)
                          RT22(M1, NN1, NN2) = ZZ7
                          ZZ8 = TI1(N11, N22)
                          IT22(M1, NN1, NN2) = ZZ8
                          QSC = QSC + (ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4&
                              & + ZZ5*ZZ5 + ZZ6*ZZ6 + ZZ7*ZZ7 + ZZ8*ZZ8)*2d0
314                   continue
214               continue
                  NNM = 2*NM
                  QXT = 0d0
                  do 215 N = 1, NNM
                      QXT = QXT + TR1(N, N)*2d0
215               continue
                  QSCA = QSCA + QSC
                  QEXT = QEXT + QXT
                  ! PRINT 7800,M,DABS(QXT),QSC,NMAX
7800              format(' m=', I3, '  qxt=', d12.6, '  qsc=', d12.6,&
                      & '  nmax=', I3)
220           continue
              COEFF1 = LAM*LAM*0.5d0/P
              CSCA = QSCA*COEFF1
              CEXT = -QEXT*COEFF1
              ! PRINT 7880, NMAX,NMAX1
7880          format('nmax=', I3, '   nmax1=', I3)
              call GSP(NMAX1, CSCA, LAM, AL1, AL2, AL3, AL4, BE1, BE2, LMAX)
              L1M = LMAX + 1
              L1MAX = max(L1MAX, L1M)
              WGII = WG1(I)
              WGI = WGII*CSCA
          do 250 L1 = 1, L1M
              ALPH1(L1) = ALPH1(L1) + AL1(L1)*WGI
              ALPH2(L1) = ALPH2(L1) + AL2(L1)*WGI
              ALPH3(L1) = ALPH3(L1) + AL3(L1)*WGI
              ALPH4(L1) = ALPH4(L1) + AL4(L1)*WGI
              BET1(L1) = BET1(L1) + BE1(L1)*WGI
              BET2(L1) = BET2(L1) + BE2(L1)*WGI
250       continue
          CSCAT = CSCAT + WGI
          CEXTIN = CEXTIN + CEXT*WGII
          ! PRINT 6070, I,NMAX,NMAX1,NGAUSS
6070      format(4i6)
500       continue
          do 510 L1 = 1, L1MAX
              ALPH1(L1) = ALPH1(L1)/CSCAT
              ALPH2(L1) = ALPH2(L1)/CSCAT
              ALPH3(L1) = ALPH3(L1)/CSCAT
              ALPH4(L1) = ALPH4(L1)/CSCAT
              BET1(L1) = BET1(L1)/CSCAT
              BET2(L1) = BET2(L1)/CSCAT
510       continue
          WALB = CSCAT/CEXTIN
          call HOVENR(L1MAX, ALPH1, ALPH2, ALPH3, ALPH4, BET1, BET2)
          ASYMM = ALPH1(2)/3d0
          ! PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
9100      format('CEXT=', D12.6, 2x, 'CSCA=', D12.6, 2x,&
              & 2x, 'W=', D12.6, 2x, '<COS>=', D12.6)
          if (WALB .gt. 1d0) print 9111
9111      format('WARNING: W IS GREATER THAN 1')

          extcross = CEXTIN
          scatcross = CSCAT
          ssalbedo = WALB
          asymmetry = ASYMM

          ! WRITE (10,580) WALB,L1MAX
          ! DO L=1,L1MAX
          !     WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),
          !     & BET1(L),BET2(L)
          ! ENDDO
575       format(6d14.7)
580       format(D14.8, I8)
          LMAX = L1MAX - 1

          ! Calculate the backscattering cross sectioun
          ! CALL MATR (ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,LMAX,&
              ! & 1, -1.0D0, backscatvec)
          ! backscat = scatcross*backscatvec(1)

          ! Code to calculate Wigner d coefficients of the phase matrix (E
          Nrank = min(L1MAX, maxrank)
          Nquad = 2*L1MAX
          if (Nquad .gt. maxnang)&
              & stop 'tmatrixwig: Number of angles exceeded'
          ! Make the Gaussian quadrature angles and get scattering matrix
          call GAUSS(Nquad, 0, 0, muang, wts)
          call MATR(ALPH1, ALPH2, ALPH3, ALPH4, BET1, BET2, LMAX,&
              & Nquad, muang, scatmat)
          do I = 1, Nquad
              do J = 1, 6
                  scatmat(J, I) = scatmat(J, I)*wts(I)
              end do
          end do
          do L = 0, Nrank
              do J = 1, 6
                  coef(J, L) = 0.0d0
              end do
          end do
          do I = 1, Nquad
              ! Calculate the Wigner functions D20, D22, D2-2 and D00 for
              ! l=0 to NRANK for muang(I)
              call WIGNERFCT(muang(I), Nrank, 0, 0, D00)
              call WIGNERFCT(muang(I), Nrank, 2, 2, D22P)
              call WIGNERFCT(-muang(I), Nrank, 2, 2, D22M)
              call WIGNERFCT(muang(I), Nrank, 2, 0, D20)
              ! scatmat elements in old Evans order (JQSRT 1991)
              A1 = scatmat(1, I)
              B1 = scatmat(2, I)
              A3 = scatmat(3, I)
              B2 = scatmat(4, I)
              A2 = scatmat(5, I)
              A4 = scatmat(6, I)
              ! Calculate the Wigner d-function series expansion coefficient
              do L = 0, Nrank
                  F = (-1.0d0)**L
                  coef(1, L) = coef(1, L) + (L + 0.5d0)*A1*D00(L)
                  coef(2, L) = coef(2, L) + (L + 0.5d0)*(A2 + A3)*D22P(L)
                  coef(3, L) = coef(3, L) + (L + 0.5d0)*(A2 - A3)*F*D22M(L)
                  coef(4, L) = coef(4, L) + (L + 0.5d0)*A4*D00(L)
                  coef(5, L) = coef(5, L) - (L + 0.5d0)*B1*D20(L)
                  ! Flip the sign of the B2 (F34) coefficient to match my con
                  coef(6, L) = coef(6, L) + (L + 0.5d0)*B2*D20(L)
              end do
          end do
          do L = 0, Nrank
              A2 = 0.5d0*(coef(2, L) + coef(3, L))
              A3 = 0.5d0*(coef(2, L) - coef(3, L))
              coef(2, L) = A2
              coef(3, L) = A3
          end do
          ! Normalize the coefficients (return Wigner d coefs times scatt
          do L = 0, Nrank
              do J = 1, 6
                  wigcoef(J, L) = SNGL(scatcross*coef(J, L)/coef(1, 0))
              end do
          end do

600   end do
      ! ITIME=MCLOCK()
      ! TIME=DFLOAT(ITIME)/6000D0
      ! PRINT 1001,TIME
1001  format(' time =', F8.2, ' min')
      return
end subroutine

!**********************************************************************

subroutine CONST(NGAUSS, NMAX, MMAX, P, X, W, AN, ANN, S, SS, NP, EPS)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)

      real*8 X(NPNG2), W(NPNG2), X1(NPNG1), W1(NPNG1),&
          & X2(NPNG1), W2(NPNG1),&
          & S(NPNG2), SS(NPNG2),&
          & AN(NPN1), ANN(NPN1, NPN1), DD(NPN1)

      do 10 N = 1, NMAX
          NN = N*(N + 1)
          AN(N) = DFLOAT(NN)
          D = DSQRT(DFLOAT(2*N + 1)/DFLOAT(NN))
          DD(N) = D
          do 101 N1 = 1, N
              DDD = D*DD(N1)*0.5d0
              ANN(N, N1) = DDD
              ANN(N1, N) = DDD
101       continue
10    continue
      NG = 2*NGAUSS
      if (NP .eq. -2) GO TO 11
      call GAUSS(NG, 0, 0, X, W)
      GO TO 19
11    NG1 = DFLOAT(NGAUSS)/2d0
      NG2 = NGAUSS - NG1
      XX = -DCOS(DATAN(EPS))
      call GAUSS(NG1, 0, 0, X1, W1)
      call GAUSS(NG2, 0, 0, X2, W2)
      do 12 I = 1, NG1
          W(I) = 0.5d0*(XX + 1d0)*W1(I)
          X(I) = 0.5d0*(XX + 1d0)*X1(I) + 0.5d0*(XX - 1d0)
12    end do
      do 14 I = 1, NG2
          W(I + NG1) = -0.5d0*XX*W2(I)
          X(I + NG1) = -0.5d0*XX*X2(I) + 0.5d0*XX
14    end do
      do 16 I = 1, NGAUSS
          W(NG - I + 1) = W(I)
          X(NG - I + 1) = -X(I)
16    end do
19    do 20 I = 1, NGAUSS
          Y = X(I)
          Y = 1d0/(1d0 - Y*Y)
          SS(I) = Y
          SS(NG - I + 1) = Y
          Y = DSQRT(Y)
          S(I) = Y
          S(NG - I + 1) = Y
20    end do
      return
end subroutine

!**********************************************************************

subroutine VARY(LAM, MRR, MRI, A, EPS, NP, NGAUSS, X, P, PPI, PIR, PII,&
      & R, DR, DDR, DRR, DRI, NMAX)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 X(NPNG2), R(NPNG2), DR(NPNG2), MRR, MRI, LAM,&
          & Z(NPNG2), ZR(NPNG2), ZI(NPNG2),&
          & J(NPNG2, NPN1), Y(NPNG2, NPN1), JR(NPNG2, NPN1),&
          & JI(NPNG2, NPN1), DJ(NPNG2, NPN1),&
          & DJR(NPNG2, NPN1), DJI(NPNG2, NPN1), DDR(NPNG2),&
          & DRR(NPNG2), DRI(NPNG2),&
          & DY(NPNG2, NPN1)
      common /CBESS/ J, Y, JR, JI, DJ, DY, DJR, DJI
      NG = NGAUSS*2
      if (NP .eq. -1) call RSP1(X, NG, NGAUSS, A, EPS, NP, R, DR)
      if (NP .ge. 0) call RSP2(X, NG, A, EPS, NP, R, DR)
      if (NP .eq. -2) call RSP3(X, NG, NGAUSS, A, EPS, R, DR)
      PI = P*2d0/LAM
      PPI = PI*PI
      PIR = PPI*MRR
      PII = PPI*MRI
      V = 1d0/(MRR*MRR + MRI*MRI)
      PRR = MRR*V
      PRI = -MRI*V
      TA = 0d0
      do 10 I = 1, NG
           VV = DSQRT(R(I))
           V = VV*PI
           TA = max(TA, V)
           VV = 1d0/V
           DDR(I) = VV
           DRR(I) = PRR*VV
           DRI(I) = PRI*VV
           V1 = V*MRR
           V2 = V*MRI
           Z(I) = V
           ZR(I) = V1
           ZI(I) = V2
10    end do
      if (NMAX .gt. NPN1) print 9000, NMAX, NPN1
      if (NMAX .gt. NPN1) stop
9000  format(' NMAX = ', I2, ', i.e., greater than ', I3)
      TB = TA*DSQRT(MRR*MRR + MRI*MRI)
      TB = DMAX1(TB, DFLOAT(NMAX))
      NNMAX1 = 1.2d0*DSQRT(DMAX1(TA, DFLOAT(NMAX))) + 3d0
      NNMAX2 = (TB + 4d0*(TB**0.33333d0) + 1.2d0*DSQRT(TB))
      NNMAX2 = NNMAX2 - NMAX + 5
      call BESS(Z, ZR, ZI, NG, NMAX, NNMAX1, NNMAX2)
      return
end subroutine

!**********************************************************************

subroutine RSP1(X, NG, NGAUSS, REV, EPS, NP, R, DR)
      implicit real*8(A - H, O - Z)
      real*8 X(NG), R(NG), DR(NG)
      A = REV*EPS**(1d0/3d0)
      AA = A*A
      EE = EPS*EPS
      EE1 = EE - 1d0
      do 50 I = 1, NGAUSS
          C = X(I)
          CC = C*C
          SS = 1d0 - CC
          S = DSQRT(SS)
          RR = 1d0/(SS + EE*CC)
          R(I) = AA*RR
          R(NG - I + 1) = R(I)
          DR(I) = RR*C*S*EE1
          DR(NG - I + 1) = -DR(I)
50    end do
      return
end subroutine

!**********************************************************************

subroutine RSP2(X, NG, REV, EPS, N, R, DR)
      implicit real*8(A - H, O - Z)
      real*8 X(NG), R(NG), DR(NG)
      DNP = DFLOAT(N)
      DN = DNP*DNP
      DN4 = DN*4d0
      EP = EPS*EPS
      A = 1d0 + 1.5d0*EP*(DN4 - 2d0)/(DN4 - 1d0)
      I = (DNP + 0.1d0)*0.5d0
      I = 2*I
      if (I .eq. N) A = A - 3d0*EPS*(1d0 + 0.25d0*EP)/&
          & (DN - 1d0) - 0.25d0*EP*EPS/(9d0*DN - 1d0)
      R0 = REV*A**(-1d0/3d0)
      do 50 I = 1, NG
         XI = DACOS(X(I))*DNP
         RI = R0*(1d0 + EPS*DCOS(XI))
         R(I) = RI*RI
         DR(I) = -R0*EPS*DNP*DSIN(XI)/RI
50    end do
      return
end subroutine

!**********************************************************************

subroutine RSP3(X, NG, NGAUSS, REV, EPS, R, DR)
      implicit real*8(A - H, O - Z)
      real*8 X(NG), R(NG), DR(NG)
      H = REV*((2d0/(3d0*EPS*EPS))**(1d0/3d0))
      A = H*EPS
      do 50 I = 1, NGAUSS
          CO = -X(I)
          SI = DSQRT(1d0 - CO*CO)
          if (SI/CO .gt. A/H) GO TO 20
          RAD = H/CO
          RTHET = H*SI/(CO*CO)
          GO TO 30
20        RAD = A/SI
          RTHET = -A*CO/(SI*SI)
30        R(I) = RAD*RAD
          R(NG - I + 1) = R(I)
          DR(I) = -RTHET/RAD
          DR(NG - I + 1) = -DR(I)
50    end do
      return
end subroutine

!***********************************************************************

subroutine BESS(X, XR, XI, NG, NMAX, NNMAX1, NNMAX2)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 X(NG), XR(NG), XI(NG),&
          & J(NPNG2, NPN1), Y(NPNG2, NPN1), JR(NPNG2, NPN1),&
          & JI(NPNG2, NPN1), DJ(NPNG2, NPN1), DY(NPNG2, NPN1),&
          & DJR(NPNG2, NPN1), DJI(NPNG2, NPN1),&
          & AJ(NPN1), AY(NPN1), AJR(NPN1), AJI(NPN1),&
          & ADJ(NPN1), ADY(NPN1), ADJR(NPN1),&
          & ADJI(NPN1)
      common /CBESS/ J, Y, JR, JI, DJ, DY, DJR, DJI

      do 10 I = 1, NG
          XX = X(I)
          call RJB(XX, AJ, ADJ, NMAX, NNMAX1)
          call RYB(XX, AY, ADY, NMAX)
          YR = XR(I)
          YI = XI(I)
          call CJB(YR, YI, AJR, AJI, ADJR, ADJI, NMAX, NNMAX2)
          do 103 N = 1, NMAX
              J(I, N) = AJ(N)
              Y(I, N) = AY(N)
              JR(I, N) = AJR(N)
              JI(I, N) = AJI(N)
              DJ(I, N) = ADJ(N)
              DY(I, N) = ADY(N)
              DJR(I, N) = ADJR(N)
              DJI(I, N) = ADJI(N)
103       continue
10    continue
      return
end subroutine

!**********************************************************************

subroutine RJB(X, Y, U, NMAX, NNMAX)
      implicit real*8(A - H, O - Z)
      real*8 Y(NMAX), U(NMAX), Z(800)
      L = NMAX + NNMAX
      XX = 1d0/X
      Z(L) = 1d0/(DFLOAT(2*L + 1)*XX)
      L1 = L - 1
      do 5 I = 1, L1
          I1 = L - I
          Z(I1) = 1d0/(DFLOAT(2*I1 + 1)*XX - Z(I1 + 1))
5     end do
      Z0 = 1d0/(XX - Z(1))
      Y0 = Z0*DCOS(X)*XX
      Y1 = Y0*Z(1)
      U(1) = Y0 - Y1*XX
      Y(1) = Y1
      do 10 I = 2, NMAX
          YI1 = Y(I - 1)
          YI = YI1*Z(I)
          U(I) = YI1 - DFLOAT(I)*YI*XX
          Y(I) = YI
10    end do
      return
end subroutine

!**********************************************************************

subroutine RYB(X, Y, V, NMAX)
      implicit real*8(A - H, O - Z)
      real*8 Y(NMAX), V(NMAX)
      C = DCOS(X)
      S = DSIN(X)
      X1 = 1d0/X
      X2 = X1*X1
      X3 = X2*X1
      Y1 = -C*X2 - S*X1
      Y(1) = Y1
      Y(2) = (-3d0*X3 + X1)*C - 3d0*X2*S
      NMAX1 = NMAX - 1
      do 5 I = 2, NMAX1
          Y(I + 1) = DFLOAT(2*I + 1)*X1*Y(I) - Y(I - 1)
5     continue
      V(1) = -X1*(C + Y1)
      do 10 I = 2, NMAX
          V(I) = Y(I - 1) - DFLOAT(I)*X1*Y(I)
10    continue
      return
end subroutine

!**********************************************************************
!                                                                     *
!   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!                                                                     *
!**********************************************************************

subroutine CJB(XR, XI, YR, YI, UR, UI, NMAX, NNMAX)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 YR(NMAX), YI(NMAX), UR(NMAX), UI(NMAX)
      real*8 CYR(NPN1), CYI(NPN1), CZR(1200), CZI(1200),&
          & CUR(NPN1), CUI(NPN1)
      L = NMAX + NNMAX
      XRXI = 1d0/(XR*XR + XI*XI)
      CXXR = XR*XRXI
      CXXI = -XI*XRXI
      QF = 1d0/DFLOAT(2*L + 1)
      CZR(L) = XR*QF
      CZI(L) = XI*QF
      L1 = L - 1
      do I = 1, L1
          I1 = L - I
          QF = DFLOAT(2*I1 + 1)
          AR = QF*CXXR - CZR(I1 + 1)
          AI = QF*CXXI - CZI(I1 + 1)
          ARI = 1d0/(AR*AR + AI*AI)
          CZR(I1) = AR*ARI
          CZI(I1) = -AI*ARI
      end do
      AR = CXXR - CZR(1)
      AI = CXXI - CZI(1)
      ARI = 1d0/(AR*AR + AI*AI)
      CZ0R = AR*ARI
      CZ0I = -AI*ARI
      CR = DCOS(XR)*DCOSH(XI)
      CI = -DSIN(XR)*DSINH(XI)
      AR = CZ0R*CR - CZ0I*CI
      AI = CZ0I*CR + CZ0R*CI
      CY0R = AR*CXXR - AI*CXXI
      CY0I = AI*CXXR + AR*CXXI
      CY1R = CY0R*CZR(1) - CY0I*CZI(1)
      CY1I = CY0I*CZR(1) + CY0R*CZI(1)
      AR = CY1R*CXXR - CY1I*CXXI
      AI = CY1I*CXXR + CY1R*CXXI
      CU1R = CY0R - AR
      CU1I = CY0I - AI
      CYR(1) = CY1R
      CYI(1) = CY1I
      CUR(1) = CU1R
      CUI(1) = CU1I
      YR(1) = CY1R
      YI(1) = CY1I
      UR(1) = CU1R
      UI(1) = CU1I
      do I = 2, NMAX
          QI = DFLOAT(I)
          CYI1R = CYR(I - 1)
          CYI1I = CYI(I - 1)
          CYIR = CYI1R*CZR(I) - CYI1I*CZI(I)
          CYII = CYI1I*CZR(I) + CYI1R*CZI(I)
          AR = CYIR*CXXR - CYII*CXXI
          AI = CYII*CXXR + CYIR*CXXI
          CUIR = CYI1R - QI*AR
          CUII = CYI1I - QI*AI
          CYR(I) = CYIR
          CYI(I) = CYII
          CUR(I) = CUIR
          CUI(I) = CUII
          YR(I) = CYIR
          YI(I) = CYII
          UR(I) = CUIR
          UI(I) = CUII
      end do
      return
end subroutine

!**********************************************************************

subroutine TMATR0(NGAUSS, X, W, AN, ANN, S, SS, PPI, PIR, PII, R, DR, DDR,&
      & DRR, DRI, NMAX, NCHECK)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)

      real*8 X(NPNG2), W(NPNG2), AN(NPN1), S(NPNG2), SS(NPNG2),&
          & R(NPNG2), DR(NPNG2), SIG(NPN2),&
          & J(NPNG2, NPN1), Y(NPNG2, NPN1),&
          & JR(NPNG2, NPN1), JI(NPNG2, NPN1), DJ(NPNG2, NPN1),&
          & DY(NPNG2, NPN1), DJR(NPNG2, NPN1),&
          & DJI(NPNG2, NPN1), DDR(NPNG2), DRR(NPNG2),&
          & D1(NPNG2, NPN1), D2(NPNG2, NPN1),&
          & DRI(NPNG2), DS(NPNG2), DSS(NPNG2), RR(NPNG2),&
          & DV1(NPN1), DV2(NPN1)

      real*8 R11(NPN1, NPN1), R12(NPN1, NPN1),&
          & R21(NPN1, NPN1), R22(NPN1, NPN1),&
          & I11(NPN1, NPN1), I12(NPN1, NPN1),&
          & I21(NPN1, NPN1), I22(NPN1, NPN1),&
          & RG11(NPN1, NPN1), RG12(NPN1, NPN1),&
          & RG21(NPN1, NPN1), RG22(NPN1, NPN1),&
          & IG11(NPN1, NPN1), IG12(NPN1, NPN1),&
          & IG21(NPN1, NPN1), IG22(NPN1, NPN1),&
          & ANN(NPN1, NPN1),&
          & QR(NPN2, NPN2), QI(NPN2, NPN2),&
          & RGQR(NPN2, NPN2), RGQI(NPN2, NPN2),&
          & TQR(NPN2, NPN2), TQI(NPN2, NPN2),&
          & TRGQR(NPN2, NPN2), TRGQI(NPN2, NPN2)
      real*8 TR1(NPN2, NPN2), TI1(NPN2, NPN2)
      real*4 PLUS(NPN6*NPN4*NPN4*8)
      common /TMAT0/ PLUS,&
          & R11, R12, R21, R22, I11, I12, I21, I22, RG11, RG12, RG21, RG22,&
          & IG11, IG12, IG21, IG22
      common /CBESS/ J, Y, JR, JI, DJ, DY, DJR, DJI
      common /CT/ TR1, TI1
      common /CTT/ QR, QI, RGQR, RGQI
      MM1 = 1
      NNMAX = NMAX + NMAX
      NG = 2*NGAUSS
      NGSS = NG
      FACTOR = 1d0
      if (NCHECK .eq. 1) then
          NGSS = NGAUSS
          FACTOR = 2d0
      else
          continue
      end if
      SI = 1d0
      do 5 N = 1, NNMAX
          SI = -SI
          SIG(N) = SI
5     end do
20    do 25 I = 1, NGAUSS
          I1 = NGAUSS + I
          I2 = NGAUSS - I + 1
          call VIG(X(I1), NMAX, 0, DV1, DV2)
          do 125 N = 1, NMAX
              SI = SIG(N)
              DD1 = DV1(N)
              DD2 = DV2(N)
              D1(I1, N) = DD1
              D2(I1, N) = DD2
              D1(I2, N) = DD1*SI
              D2(I2, N) = -DD2*SI
125       continue
25    continue
30    do 40 I = 1, NGSS
           RR(I) = W(I)*R(I)
40    end do

      do 300 N1 = MM1, NMAX
          AN1 = AN(N1)
          do 333 N2 = MM1, NMAX
              AN2 = AN(N2)
              AR12 = 0d0
              AR21 = 0d0
              AI12 = 0d0
              AI21 = 0d0
              GR12 = 0d0
              GR21 = 0d0
              GI12 = 0d0
              GI21 = 0d0
              if (NCHECK .eq. 1 .and. SIG(N1 + N2) .lt. 0d0) GO TO 205

              do 200 I = 1, NGSS
                  D1N1 = D1(I, N1)
                  D2N1 = D2(I, N1)
                  D1N2 = D1(I, N2)
                  D2N2 = D2(I, N2)
                  A12 = D1N1*D2N2
                  A21 = D2N1*D1N2
                  A22 = D2N1*D2N2
                  AA1 = A12 + A21

                  QJ1 = J(I, N1)
                  QY1 = Y(I, N1)
                  QJR2 = JR(I, N2)
                  QJI2 = JI(I, N2)
                  QDJR2 = DJR(I, N2)
                  QDJI2 = DJI(I, N2)
                  QDJ1 = DJ(I, N1)
                  QDY1 = DY(I, N1)

                  C1R = QJR2*QJ1
                  C1I = QJI2*QJ1
                  B1R = C1R - QJI2*QY1
                  B1I = C1I + QJR2*QY1

                  C2R = QJR2*QDJ1
                  C2I = QJI2*QDJ1
                  B2R = C2R - QJI2*QDY1
                  B2I = C2I + QJR2*QDY1

                  DDRI = DDR(I)
                  C3R = DDRI*C1R
                  C3I = DDRI*C1I
                  B3R = DDRI*B1R
                  B3I = DDRI*B1I

                  C4R = QDJR2*QJ1
                  C4I = QDJI2*QJ1
                  B4R = C4R - QDJI2*QY1
                  B4I = C4I + QDJR2*QY1

                  DRRI = DRR(I)
                  DRII = DRI(I)
                  C5R = C1R*DRRI - C1I*DRII
                  C5I = C1I*DRRI + C1R*DRII
                  B5R = B1R*DRRI - B1I*DRII
                  B5I = B1I*DRRI + B1R*DRII

                  URI = DR(I)
                  RRI = RR(I)

                  F1 = RRI*A22
                  F2 = RRI*URI*AN1*A12
                  AR12 = AR12 + F1*B2R + F2*B3R
                  AI12 = AI12 + F1*B2I + F2*B3I
                  GR12 = GR12 + F1*C2R + F2*C3R
                  GI12 = GI12 + F1*C2I + F2*C3I

                  F2 = RRI*URI*AN2*A21
                  AR21 = AR21 + F1*B4R + F2*B5R
                  AI21 = AI21 + F1*B4I + F2*B5I
                  GR21 = GR21 + F1*C4R + F2*C5R
                  GI21 = GI21 + F1*C4I + F2*C5I
200           continue

205           AN12 = ANN(N1, N2)*FACTOR
              R12(N1, N2) = AR12*AN12
              R21(N1, N2) = AR21*AN12
              I12(N1, N2) = AI12*AN12
              I21(N1, N2) = AI21*AN12
              RG12(N1, N2) = GR12*AN12
              RG21(N1, N2) = GR21*AN12
              IG12(N1, N2) = GI12*AN12
              IG21(N1, N2) = GI21*AN12
333       continue
300   continue

      TPIR = PIR
      TPII = PII
      TPPI = PPI

      NM = NMAX
      do 310 N1 = MM1, NMAX
          K1 = N1 - MM1 + 1
          KK1 = K1 + NM
          do 311 N2 = MM1, NMAX
              K2 = N2 - MM1 + 1
              KK2 = K2 + NM

              TAR12 = I12(N1, N2)
              TAI12 = -R12(N1, N2)
              TGR12 = IG12(N1, N2)
              TGI12 = -RG12(N1, N2)

              TAR21 = -I21(N1, N2)
              TAI21 = R21(N1, N2)
              TGR21 = -IG21(N1, N2)
              TGI21 = RG21(N1, N2)

              TQR(K1, K2) = TPIR*TAR21 - TPII*TAI21 + TPPI*TAR12
              TQI(K1, K2) = TPIR*TAI21 + TPII*TAR21 + TPPI*TAI12
              TRGQR(K1, K2) = TPIR*TGR21 - TPII*TGI21 + TPPI*TGR12
              TRGQI(K1, K2) = TPIR*TGI21 + TPII*TGR21 + TPPI*TGI12

              TQR(K1, KK2) = 0d0
              TQI(K1, KK2) = 0d0
              TRGQR(K1, KK2) = 0d0
              TRGQI(K1, KK2) = 0d0

              TQR(KK1, K2) = 0d0
              TQI(KK1, K2) = 0d0
              TRGQR(KK1, K2) = 0d0
              TRGQI(KK1, K2) = 0d0

              TQR(KK1, KK2) = TPIR*TAR12 - TPII*TAI12 + TPPI*TAR21
              TQI(KK1, KK2) = TPIR*TAI12 + TPII*TAR12 + TPPI*TAI21
              TRGQR(KK1, KK2) = TPIR*TGR12 - TPII*TGI12 + TPPI*TGR21
              TRGQI(KK1, KK2) = TPIR*TGI12 + TPII*TGR12 + TPPI*TGI21
311       continue
310   continue

      NNMAX = 2*NM
      do 320 N1 = 1, NNMAX
          do 322 N2 = 1, NNMAX
              QR(N1, N2) = TQR(N1, N2)
              QI(N1, N2) = TQI(N1, N2)
              RGQR(N1, N2) = TRGQR(N1, N2)
              RGQI(N1, N2) = TRGQI(N1, N2)
322       continue
320   continue
      call TT(NMAX, NCHECK)
      return
end subroutine

!**********************************************************************

subroutine TMATR(M, NGAUSS, X, W, AN, ANN, S, SS, PPI, PIR, PII, R, DR, DDR,&
      & DRR, DRI, NMAX, NCHECK)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 X(NPNG2), W(NPNG2), AN(NPN1), S(NPNG2), SS(NPNG2),&
          & R(NPNG2), DR(NPNG2), SIG(NPN2),&
          & J(NPNG2, NPN1), Y(NPNG2, NPN1),&
          & JR(NPNG2, NPN1), JI(NPNG2, NPN1), DJ(NPNG2, NPN1),&
          & DY(NPNG2, NPN1), DJR(NPNG2, NPN1),&
          & DJI(NPNG2, NPN1), DDR(NPNG2), DRR(NPNG2),&
          & D1(NPNG2, NPN1), D2(NPNG2, NPN1),&
          & DRI(NPNG2), DS(NPNG2), DSS(NPNG2), RR(NPNG2),&
          & DV1(NPN1), DV2(NPN1)

      real*8 R11(NPN1, NPN1), R12(NPN1, NPN1),&
          & R21(NPN1, NPN1), R22(NPN1, NPN1),&
          & I11(NPN1, NPN1), I12(NPN1, NPN1),&
          & I21(NPN1, NPN1), I22(NPN1, NPN1),&
          & RG11(NPN1, NPN1), RG12(NPN1, NPN1),&
          & RG21(NPN1, NPN1), RG22(NPN1, NPN1),&
          & IG11(NPN1, NPN1), IG12(NPN1, NPN1),&
          & IG21(NPN1, NPN1), IG22(NPN1, NPN1),&
          & ANN(NPN1, NPN1),&
          & QR(NPN2, NPN2), QI(NPN2, NPN2),&
          & RGQR(NPN2, NPN2), RGQI(NPN2, NPN2),&
          & TQR(NPN2, NPN2), TQI(NPN2, NPN2),&
          & TRGQR(NPN2, NPN2), TRGQI(NPN2, NPN2)
      real*8 TR1(NPN2, NPN2), TI1(NPN2, NPN2)
      real*4 PLUS(NPN6*NPN4*NPN4*8)
      common /TMAT0/ PLUS,&
          & R11, R12, R21, R22, I11, I12, I21, I22, RG11, RG12, RG21, RG22,&
          & IG11, IG12, IG21, IG22
      common /CBESS/ J, Y, JR, JI, DJ, DY, DJR, DJI
      common /CT/ TR1, TI1
      common /CTT/ QR, QI, RGQR, RGQI
      MM1 = M
      QM = DFLOAT(M)
      QMM = QM*QM
      NG = 2*NGAUSS
      NGSS = NG
      FACTOR = 1d0
      if (NCHECK .eq. 1) then
          NGSS = NGAUSS
          FACTOR = 2d0
      else
          continue
      end if
      SI = 1d0
      NM = NMAX + NMAX
      do 5 N = 1, NM
          SI = -SI
          SIG(N) = SI
5     end do
20    do 25 I = 1, NGAUSS
          I1 = NGAUSS + I
          I2 = NGAUSS - I + 1
          call VIG(X(I1), NMAX, M, DV1, DV2)
          do 225 N = 1, NMAX
              SI = SIG(N)
              DD1 = DV1(N)
              DD2 = DV2(N)
              D1(I1, N) = DD1
              D2(I1, N) = DD2
              D1(I2, N) = DD1*SI
              D2(I2, N) = -DD2*SI
225       continue
25    continue
30    do 40 I = 1, NGSS
          WR = W(I)*R(I)
          DS(I) = S(I)*QM*WR
          DSS(I) = SS(I)*QMM
          RR(I) = WR
40    end do

      do 300 N1 = MM1, NMAX
          AN1 = AN(N1)
          do 330 N2 = MM1, NMAX
              AN2 = AN(N2)
              AR11 = 0d0
              AR12 = 0d0
              AR21 = 0d0
              AR22 = 0d0
              AI11 = 0d0
              AI12 = 0d0
              AI21 = 0d0
              AI22 = 0d0
              GR11 = 0d0
              GR12 = 0d0
              GR21 = 0d0
              GR22 = 0d0
              GI11 = 0d0
              GI12 = 0d0
              GI21 = 0d0
              GI22 = 0d0
              SI = SIG(N1 + N2)

              do 200 I = 1, NGSS
                  D1N1 = D1(I, N1)
                  D2N1 = D2(I, N1)
                  D1N2 = D1(I, N2)
                  D2N2 = D2(I, N2)
                  A11 = D1N1*D1N2
                  A12 = D1N1*D2N2
                  A21 = D2N1*D1N2
                  A22 = D2N1*D2N2
                  AA1 = A12 + A21
                  AA2 = A11*DSS(I) + A22
                  QJ1 = J(I, N1)
                  QY1 = Y(I, N1)
                  QJR2 = JR(I, N2)
                  QJI2 = JI(I, N2)
                  QDJR2 = DJR(I, N2)
                  QDJI2 = DJI(I, N2)
                  QDJ1 = DJ(I, N1)
                  QDY1 = DY(I, N1)

                  C1R = QJR2*QJ1
                  C1I = QJI2*QJ1
                  B1R = C1R - QJI2*QY1
                  B1I = C1I + QJR2*QY1

                  C2R = QJR2*QDJ1
                  C2I = QJI2*QDJ1
                  B2R = C2R - QJI2*QDY1
                  B2I = C2I + QJR2*QDY1

                  DDRI = DDR(I)
                  C3R = DDRI*C1R
                  C3I = DDRI*C1I
                  B3R = DDRI*B1R
                  B3I = DDRI*B1I

                  C4R = QDJR2*QJ1
                  C4I = QDJI2*QJ1
                  B4R = C4R - QDJI2*QY1
                  B4I = C4I + QDJR2*QY1

                  DRRI = DRR(I)
                  DRII = DRI(I)
                  C5R = C1R*DRRI - C1I*DRII
                  C5I = C1I*DRRI + C1R*DRII
                  B5R = B1R*DRRI - B1I*DRII
                  B5I = B1I*DRRI + B1R*DRII

                  C6R = QDJR2*QDJ1
                  C6I = QDJI2*QDJ1
                  B6R = C6R - QDJI2*QDY1
                  B6I = C6I + QDJR2*QDY1

                  C7R = C4R*DDRI
                  C7I = C4I*DDRI
                  B7R = B4R*DDRI
                  B7I = B4I*DDRI

                  C8R = C2R*DRRI - C2I*DRII
                  C8I = C2I*DRRI + C2R*DRII
                  B8R = B2R*DRRI - B2I*DRII
                  B8I = B2I*DRRI + B2R*DRII

                  URI = DR(I)
                  DSI = DS(I)
                  DSSI = DSS(I)
                  RRI = RR(I)

                  if (NCHECK .eq. 1 .and. SI .gt. 0d0) GO TO 150

                  E1 = DSI*AA1
                  AR11 = AR11 + E1*B1R
                  AI11 = AI11 + E1*B1I
                  GR11 = GR11 + E1*C1R
                  GI11 = GI11 + E1*C1I
                  if (NCHECK .eq. 1) GO TO 160

150               F1 = RRI*AA2
                  F2 = RRI*URI*AN1*A12
                  AR12 = AR12 + F1*B2R + F2*B3R
                  AI12 = AI12 + F1*B2I + F2*B3I
                  GR12 = GR12 + F1*C2R + F2*C3R
                  GI12 = GI12 + F1*C2I + F2*C3I

                  F2 = RRI*URI*AN2*A21
                  AR21 = AR21 + F1*B4R + F2*B5R
                  AI21 = AI21 + F1*B4I + F2*B5I
                  GR21 = GR21 + F1*C4R + F2*C5R
                  GI21 = GI21 + F1*C4I + F2*C5I
                  if (NCHECK .eq. 1) GO TO 200

160               E2 = DSI*URI*A11
                  E3 = E2*AN2
                  E2 = E2*AN1
                  AR22 = AR22 + E1*B6R + E2*B7R + E3*B8R
                  AI22 = AI22 + E1*B6I + E2*B7I + E3*B8I
                  GR22 = GR22 + E1*C6R + E2*C7R + E3*C8R
                  GI22 = GI22 + E1*C6I + E2*C7I + E3*C8I
200           continue
              AN12 = ANN(N1, N2)*FACTOR
              R11(N1, N2) = AR11*AN12
              R12(N1, N2) = AR12*AN12
              R21(N1, N2) = AR21*AN12
              R22(N1, N2) = AR22*AN12
              I11(N1, N2) = AI11*AN12
              I12(N1, N2) = AI12*AN12
              I21(N1, N2) = AI21*AN12
              I22(N1, N2) = AI22*AN12
              RG11(N1, N2) = GR11*AN12
              RG12(N1, N2) = GR12*AN12
              RG21(N1, N2) = GR21*AN12
              RG22(N1, N2) = GR22*AN12
              IG11(N1, N2) = GI11*AN12
              IG12(N1, N2) = GI12*AN12
              IG21(N1, N2) = GI21*AN12
              IG22(N1, N2) = GI22*AN12
330       continue
300   continue
      TPIR = PIR
      TPII = PII
      TPPI = PPI
      NM = NMAX - MM1 + 1
      do 310 N1 = MM1, NMAX
          K1 = N1 - MM1 + 1
          KK1 = K1 + NM
          do 340 N2 = MM1, NMAX
              K2 = N2 - MM1 + 1
              KK2 = K2 + NM

              TAR11 = -R11(N1, N2)
              TAI11 = -I11(N1, N2)
              TGR11 = -RG11(N1, N2)
              TGI11 = -IG11(N1, N2)

              TAR12 = I12(N1, N2)
              TAI12 = -R12(N1, N2)
              TGR12 = IG12(N1, N2)
              TGI12 = -RG12(N1, N2)

              TAR21 = -I21(N1, N2)
              TAI21 = R21(N1, N2)
              TGR21 = -IG21(N1, N2)
              TGI21 = RG21(N1, N2)

              TAR22 = -R22(N1, N2)
              TAI22 = -I22(N1, N2)
              TGR22 = -RG22(N1, N2)
              TGI22 = -IG22(N1, N2)

              TQR(K1, K2) = TPIR*TAR21 - TPII*TAI21 + TPPI*TAR12
              TQI(K1, K2) = TPIR*TAI21 + TPII*TAR21 + TPPI*TAI12
              TRGQR(K1, K2) = TPIR*TGR21 - TPII*TGI21 + TPPI*TGR12
              TRGQI(K1, K2) = TPIR*TGI21 + TPII*TGR21 + TPPI*TGI12

              TQR(K1, KK2) = TPIR*TAR11 - TPII*TAI11 + TPPI*TAR22
              TQI(K1, KK2) = TPIR*TAI11 + TPII*TAR11 + TPPI*TAI22
              TRGQR(K1, KK2) = TPIR*TGR11 - TPII*TGI11 + TPPI*TGR22
              TRGQI(K1, KK2) = TPIR*TGI11 + TPII*TGR11 + TPPI*TGI22

              TQR(KK1, K2) = TPIR*TAR22 - TPII*TAI22 + TPPI*TAR11
              TQI(KK1, K2) = TPIR*TAI22 + TPII*TAR22 + TPPI*TAI11
              TRGQR(KK1, K2) = TPIR*TGR22 - TPII*TGI22 + TPPI*TGR11
              TRGQI(KK1, K2) = TPIR*TGI22 + TPII*TGR22 + TPPI*TGI11

              TQR(KK1, KK2) = TPIR*TAR12 - TPII*TAI12 + TPPI*TAR21
              TQI(KK1, KK2) = TPIR*TAI12 + TPII*TAR12 + TPPI*TAI21
              TRGQR(KK1, KK2) = TPIR*TGR12 - TPII*TGI12 + TPPI*TGR21
              TRGQI(KK1, KK2) = TPIR*TGI12 + TPII*TGR12 + TPPI*TGI21
340       continue
310   continue

      NNMAX = 2*NM
      do 320 N1 = 1, NNMAX
          do 323 N2 = 1, NNMAX
              QR(N1, N2) = TQR(N1, N2)
              QI(N1, N2) = TQI(N1, N2)
              RGQR(N1, N2) = TRGQR(N1, N2)
              RGQI(N1, N2) = TRGQI(N1, N2)
323       continue
320   continue

      call TT(NM, NCHECK)

      return
end subroutine

!*****************************************************************

subroutine VIG(X, NMAX, M, DV1, DV2)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 DV1(NPN1), DV2(NPN1)

      A = 1d0
      QS = DSQRT(1d0 - X*X)
      QS1 = 1d0/QS
      do N = 1, NMAX
          DV1(N) = 0d0
          DV2(N) = 0d0
      end do
      if (M .ne. 0) GO TO 20
      D1 = 1d0
      D2 = X
      do N = 1, NMAX
          QN = DFLOAT(N)
          QN1 = DFLOAT(N + 1)
          QN2 = DFLOAT(2*N + 1)
          D3 = (QN2*X*D2 - QN*D1)/QN1
          DER = QS1*(QN1*QN/QN2)*(-D1 + D3)
          DV1(N) = D2
          DV2(N) = DER
          D1 = D2
          D2 = D3
      end do
      return
20    QMM = DFLOAT(M*M)
      do I = 1, M
          I2 = I*2
          A = A*DSQRT(DFLOAT(I2 - 1)/DFLOAT(I2))*QS
      end do
      D1 = 0d0
      D2 = A
      do N = M, NMAX
          QN = DFLOAT(N)
          QN2 = DFLOAT(2*N + 1)
          QN1 = DFLOAT(N + 1)
          QNM = DSQRT(QN*QN - QMM)
          QNM1 = DSQRT(QN1*QN1 - QMM)
          D3 = (QN2*X*D2 - QNM*D1)/QNM1
          DER = QS1*(-QN1*QNM*D1 + QN*QNM1*D3)/QN2
          DV1(N) = D2
          DV2(N) = DER
          D1 = D2
          D2 = D3
      end do
      return
end subroutine

!**********************************************************************
!                                                                     *
!   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
!                                                                     *
!   INPUT INFORTMATION IS IN COMMON /CTT/                             *
!   OUTPUT INFORMATION IS IN COMMON /CT/                              *
!                                                                     *
!**********************************************************************

subroutine TT(NMAX, NCHECK)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 F(NPN2, NPN2), B(NPN2), WORK(NPN2),&
          & QR(NPN2, NPN2), QI(NPN2, NPN2),&
          & RGQR(NPN2, NPN2), RGQI(NPN2, NPN2),&
          & A(NPN2, NPN2), C(NPN2, NPN2), D(NPN2, NPN2), E(NPN2, NPN2)
      real*8 TR1(NPN2, NPN2), TI1(NPN2, NPN2)
      complex*16 ZQ(NPN2, NPN2), ZW(NPN2)
      integer IPIV(NPN2), IPVT(NPN2)
      common /CT/ TR1, TI1
      common /CTT/ QR, QI, RGQR, RGQI
      NDIM = NPN2
      NNMAX = 2*NMAX

      ! Matrix inversion from LAPACK

      do I = 1, NNMAX
      do J = 1, NNMAX
          ZQ(I, J) = DCMPLX(QR(I, J), QI(I, J))
      end do
      end do
      INFO = 0
      call ZGETRF(NNMAX, NNMAX, ZQ, NPN2, IPIV, INFO)
      if (INFO .ne. 0) write (6, 1100) INFO
      call ZGETRI(NNMAX, ZQ, NPN2, IPIV, ZW, NPN2, INFO)
      if (INFO .ne. 0) write (6, 1100) INFO

1100  format('WARNING: matrix inversion error info=', i6)
      do I = 1, NNMAX
          do J = 1, NNMAX
              TR = 0d0
              TI = 0d0
              do K = 1, NNMAX
                  ARR = RGQR(I, K)
                  ARI = RGQI(I, K)
                  AR = ZQ(K, J)
                  AI = DIMAG(ZQ(K, J))
                  TR = TR - ARR*AR + ARI*AI
                  TI = TI - ARR*AI - ARI*AR
              end do
              TR1(I, J) = TR
              TI1(I, J) = TI
          end do
      end do
      return
end subroutine

!********************************************************************
!                                                                   *
!   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
!   REPRESENTATION.                                                 *
!                                                                   *
!   INPUT PARAMETERS:                                               *
!                                                                   *
!      LAM - WAVELENGTH OF LIGHT                                    *
!      CSCA - SCATTERING CROSS SECTION                              *
!      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
!                  COMMON /CTM/                                     *
!      NMAX - DIMENSION OF T(M)-MATRICES                            *
!                                                                   *
!   OUTPUT INFORTMATION:                                            *
!                                                                   *
!      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
!      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
!                                                                   *
!********************************************************************

subroutine GSP(NMAX, CSCA, LAM, ALF1, ALF2, ALF3, ALF4, BET1, BET2, LMAX)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - B, D - H, O - Z), complex*16(C)

      real*8 LAM, SSIGN(900)

      real*8 CSCA, SSI(NPL), SSJ(NPN1),&
          & ALF1(NPL), ALF2(NPL), ALF3(NPL),&
          & ALF4(NPL), BET1(NPL), BET2(NPL),&
          & TR1(NPL1, NPN4), TR2(NPL1, NPN4),&
          & TI1(NPL1, NPN4), TI2(NPL1, NPN4),&
          & G1(NPL1, NPN6), G2(NPL1, NPN6),&
          & AR1(NPN4), AR2(NPN4), AI1(NPN4), AI2(NPN4),&
          & FR(NPN4, NPN4), FI(NPN4, NPN4), FF(NPN4, NPN4)

      real*4 B1R(NPL1, NPL1, NPN4), B1I(NPL1, NPL1, NPN4),&
          & B2R(NPL1, NPL1, NPN4), B2I(NPL1, NPL1, NPN4),&
          & D1(NPL1, NPN4, NPN4), D2(NPL1, NPN4, NPN4),&
          & D3(NPL1, NPN4, NPN4), D4(NPL1, NPN4, NPN4),&
          & D5R(NPL1, NPN4, NPN4), D5I(NPL1, NPN4, NPN4),&
          & PLUS1(NPN6*NPN4*NPN4*8)

      real*8&
          & TR11(NPN6, NPN4, NPN4), TR12(NPN6, NPN4, NPN4),&
          & TR21(NPN6, NPN4, NPN4), TR22(NPN6, NPN4, NPN4),&
          & TI11(NPN6, NPN4, NPN4), TI12(NPN6, NPN4, NPN4),&
          & TI21(NPN6, NPN4, NPN4), TI22(NPN6, NPN4, NPN4)

      !real*4&
      !    & RT11(NPN6, NPN4, NPN4), RT12(NPN6, NPN4, NPN4),&
      !    & RT21(NPN6, NPN4, NPN4), RT22(NPN6, NPN4, NPN4),&
      !    & IT11(NPN6, NPN4, NPN4), IT12(NPN6, NPN4, NPN4),&
      !    & IT21(NPN6, NPN4, NPN4), IT22(NPN6, NPN4, NPN4)

      complex*16 CIM(NPN1)

      common /TMAT/ TR11, TR12, TR21, TR22, TI11, TI12, TI21, TI22

      common /CBESS_G/ B1R, B1I, B2R, B2I
      common /SS/ SSIGN

      equivalence(PLUS1(1), TR11(1, 1, 1))

      equivalence(D1(1, 1, 1), PLUS1(1)),&
          & (D2(1, 1, 1), PLUS1(NPL1*NPN4*NPN4 + 1)),&
          & (D3(1, 1, 1), PLUS1(NPL1*NPN4*NPN4*2 + 1)),&
          & (D4(1, 1, 1), PLUS1(NPL1*NPN4*NPN4*3 + 1)),&
          & (D5R(1, 1, 1), PLUS1(NPL1*NPN4*NPN4*4 + 1))

      call FACT
      call SIGNUM
      LMAX = 2*NMAX
      L1MAX = LMAX + 1
      CI = (0d0, 1d0)
      CIM(1) = CI
      do 2 I = 2, NMAX
          CIM(I) = CIM(I - 1)*CI
2     end do
      SSI(1) = 1d0
      do 3 I = 1, LMAX
          I1 = I + 1
          SI = DFLOAT(2*I + 1)
          SSI(I1) = SI
          if (I .le. NMAX) SSJ(I) = DSQRT(SI)
3     end do
      CI = -CI
      do 5 I = 1, NMAX
          SI = SSJ(I)
          CCI = CIM(I)
          do 4 J = 1, NMAX
              SJ = 1d0/SSJ(J)
              CCJ = CIM(J)*SJ/CCI
              FR(J, I) = CCJ
              FI(J, I) = CCJ*CI
              FF(J, I) = SI*SJ
4         continue
5     end do
      NMAX1 = NMAX + 1

      ! *****  CALCULATION OF THE ARRAYS B1 AND B2  *****

      K1 = 1
      K2 = 0
      K3 = 0
      K4 = 1
      K5 = 1
      K6 = 2

      ! print 3300, B1,B2
3300  format(' B1 AND B2')
      do 100 N = 1, NMAX

      ! *****  CALCULATION OF THE ARRAYS T1 AND T2  *****

          do 10 NN = 1, NMAX
              M1MAX = MIN0(N, NN) + 1
              do 6 M1 = 1, M1MAX
                  M = M1 - 1
                  L1 = NPN6 + M

                  TT1 = TR11(M1, N, NN)
                  TT2 = TR12(M1, N, NN)
                  TT3 = TR21(M1, N, NN)
                  TT4 = TR22(M1, N, NN)
                  TT5 = TI11(M1, N, NN)
                  TT6 = TI12(M1, N, NN)
                  TT7 = TI21(M1, N, NN)
                  TT8 = TI22(M1, N, NN)

                  T1 = TT1 + TT2
                  T2 = TT3 + TT4
                  T3 = TT5 + TT6
                  T4 = TT7 + TT8
                  TR1(L1, NN) = T1 + T2
                  TR2(L1, NN) = T1 - T2
                  TI1(L1, NN) = T3 + T4
                  TI2(L1, NN) = T3 - T4
                  if (M .eq. 0) GO TO 6
                  L1 = NPN6 - M
                  T1 = TT1 - TT2
                  T2 = TT3 - TT4
                  T3 = TT5 - TT6
                  T4 = TT7 - TT8
                  TR1(L1, NN) = T1 - T2
                  TR2(L1, NN) = T1 + T2
                  TI1(L1, NN) = T3 - T4
                  TI2(L1, NN) = T3 + T4
6             continue
10        continue

          !  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****

          NN1MAX = NMAX1 + N
          do 40 NN1 = 1, NN1MAX
              N1 = NN1 - 1

              !  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****

              call CCG(N, N1, NMAX, K1, K2, G1)
              NNMAX = MIN0(NMAX, N1 + N)
              NNMIN = MAX0(1, IABS(N - N1))
              KN = N + NN1
              do 15 NN = NNMIN, NNMAX
                  NNN = NN + 1
                  SIG = SSIGN(KN + NN)
                  M1MAX = MIN0(N, NN) + NPN6
                  AAR1 = 0d0
                  AAR2 = 0d0
                  AAI1 = 0d0
                  AAI2 = 0d0
                  do 13 M1 = NPN6, M1MAX
                      M = M1 - NPN6
                      SSS = G1(M1, NNN)
                      RR1 = TR1(M1, NN)
                      RI1 = TI1(M1, NN)
                      RR2 = TR2(M1, NN)
                      RI2 = TI2(M1, NN)
                      if (M .eq. 0) GO TO 12
                      M2 = NPN6 - M
                      RR1 = RR1 + TR1(M2, NN)*SIG
                      RI1 = RI1 + TI1(M2, NN)*SIG
                      RR2 = RR2 + TR2(M2, NN)*SIG
                      RI2 = RI2 + TI2(M2, NN)*SIG
12                    AAR1 = AAR1 + SSS*RR1
                      AAI1 = AAI1 + SSS*RI1
                      AAR2 = AAR2 + SSS*RR2
                      AAI2 = AAI2 + SSS*RI2
13                continue
                  XR = FR(NN, N)
                  XI = FI(NN, N)
                  AR1(NN) = AAR1*XR - AAI1*XI
                  AI1(NN) = AAR1*XI + AAI1*XR
                  AR2(NN) = AAR2*XR - AAI2*XI
                  AI2(NN) = AAR2*XI + AAI2*XR
15            continue

              !  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****

              call CCG(N, N1, NMAX, K3, K4, G2)
              M1 = MAX0(-N1 + 1, -N)
              M2 = MIN0(N1 + 1, N)
              M1MAX = M2 + NPN6
              M1MIN = M1 + NPN6
              do 30 M1 = M1MIN, M1MAX
                  BBR1 = 0d0
                  BBI1 = 0d0
                  BBR2 = 0d0
                  BBI2 = 0d0
                  do 25 NN = NNMIN, NNMAX
                      NNN = NN + 1
                      SSS = G2(M1, NNN)
                      BBR1 = BBR1 + SSS*AR1(NN)
                      BBI1 = BBI1 + SSS*AI1(NN)
                      BBR2 = BBR2 + SSS*AR2(NN)
                      BBI2 = BBI2 + SSS*AI2(NN)
25                continue
                  B1R(NN1, M1, N) = BBR1
                  B1I(NN1, M1, N) = BBI1
                  B2R(NN1, M1, N) = BBR2
                  B2I(NN1, M1, N) = BBI2
30            continue
40        continue
100   end do

      !  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****

      !  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****

      ! print 3301
3301  format(' D1, D2, ...')
      do 200 N = 1, NMAX
          do 190 NN = 1, NMAX
              M1 = MIN0(N, NN)
              M1MAX = NPN6 + M1
              M1MIN = NPN6 - M1
              NN1MAX = NMAX1 + MIN0(N, NN)
              do 180 M1 = M1MIN, M1MAX
                  M = M1 - NPN6
                  NN1MIN = IABS(M - 1) + 1
                  DD1 = 0d0
                  DD2 = 0d0
                  do 150 NN1 = NN1MIN, NN1MAX
                      XX = SSI(NN1)
                      X1 = B1R(NN1, M1, N)
                      X2 = B1I(NN1, M1, N)
                      X3 = B1R(NN1, M1, NN)
                      X4 = B1I(NN1, M1, NN)
                      X5 = B2R(NN1, M1, N)
                      X6 = B2I(NN1, M1, N)
                      X7 = B2R(NN1, M1, NN)
                      X8 = B2I(NN1, M1, NN)
                      DD1 = DD1 + XX*(X1*X3 + X2*X4)
                      DD2 = DD2 + XX*(X5*X7 + X6*X8)
150               continue
                  D1(M1, NN, N) = DD1
                  D2(M1, NN, N) = DD2
180           continue
              MMAX = MIN0(N, NN + 2)
              MMIN = MAX0(-N, -NN + 2)
              M1MAX = NPN6 + MMAX
              M1MIN = NPN6 + MMIN
              do 186 M1 = M1MIN, M1MAX
                  M = M1 - NPN6
                  NN1MIN = IABS(M - 1) + 1
                  DD3 = 0d0
                  DD4 = 0d0
                  DD5R = 0d0
                  DD5I = 0d0
                  M2 = -M + 2 + NPN6
                  do 183 NN1 = NN1MIN, NN1MAX
                      XX = SSI(NN1)
                      X1 = B1R(NN1, M1, N)
                      X2 = B1I(NN1, M1, N)
                      X3 = B2R(NN1, M1, N)
                      X4 = B2I(NN1, M1, N)
                      X5 = B1R(NN1, M2, NN)
                      X6 = B1I(NN1, M2, NN)
                      X7 = B2R(NN1, M2, NN)
                      X8 = B2I(NN1, M2, NN)
                      DD3 = DD3 + XX*(X1*X5 + X2*X6)
                      DD4 = DD4 + XX*(X3*X7 + X4*X8)
                      DD5R = DD5R + XX*(X3*X5 + X4*X6)
                      DD5I = DD5I + XX*(X4*X5 - X3*X6)
183               continue
                  D3(M1, NN, N) = DD3
                  D4(M1, NN, N) = DD4
                  D5R(M1, NN, N) = DD5R
                  D5I(M1, NN, N) = DD5I
186           continue
190       continue
200   end do

      !  *****  END OF THE CALCULATION OF THE D-ARRAYS *****

      !  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****

      ! print 3303
3303  format(' G1, G2, ...')

      DK = LAM*LAM/(4d0*CSCA*DACOS(-1d0))
      do 300 L1 = 1, L1MAX
          G1L = 0d0
          G2L = 0d0
          G3L = 0d0
          G4L = 0d0
          G5LR = 0d0
          G5LI = 0d0
          L = L1 - 1
          SL = SSI(L1)*DK
          do 290 N = 1, NMAX
              NNMIN = MAX0(1, IABS(N - L))
              NNMAX = MIN0(NMAX, N + L)
              if (NNMAX .lt. NNMIN) GO TO 290
              call CCG(N, L, NMAX, K1, K2, G1)
              if (L .ge. 2) call CCG(N, L, NMAX, K5, K6, G2)
              NL = N + L
              do 280 NN = NNMIN, NNMAX
                  NNN = NN + 1
                  MMAX = MIN0(N, NN)
                  M1MIN = NPN6 - MMAX
                  M1MAX = NPN6 + MMAX
                  SI = SSIGN(NL + NNN)
                  DM1 = 0d0
                  DM2 = 0d0
                  do 270 M1 = M1MIN, M1MAX
                      M = M1 - NPN6
                      if (M .ge. 0) SSS1 = G1(M1, NNN)
                      if (M .lt. 0) SSS1 = G1(NPN6 - M, NNN)*SI
                      DM1 = DM1 + SSS1*D1(M1, NN, N)
                      DM2 = DM2 + SSS1*D2(M1, NN, N)
270               continue
                  FFN = FF(NN, N)
                  SSS = G1(NPN6 + 1, NNN)*FFN
                  G1L = G1L + SSS*DM1
                  G2L = G2L + SSS*DM2*SI
                  if (L .lt. 2) GO TO 280
                  DM3 = 0d0
                  DM4 = 0d0
                  DM5R = 0d0
                  DM5I = 0d0
                  MMAX = MIN0(N, NN + 2)
                  MMIN = MAX0(-N, -NN + 2)
                  M1MAX = NPN6 + MMAX
                  M1MIN = NPN6 + MMIN
                  do 275 M1 = M1MIN, M1MAX
                      M = M1 - NPN6
                      SSS1 = G2(NPN6 - M, NNN)
                      DM3 = DM3 + SSS1*D3(M1, NN, N)
                      DM4 = DM4 + SSS1*D4(M1, NN, N)
                      DM5R = DM5R + SSS1*D5R(M1, NN, N)
                      DM5I = DM5I + SSS1*D5I(M1, NN, N)
275               continue
                  G5LR = G5LR - SSS*DM5R
                  G5LI = G5LI - SSS*DM5I
                  SSS = G2(NPN4, NNN)*FFN
                  G3L = G3L + SSS*DM3
                  G4L = G4L + SSS*DM4*SI
280           continue
290       continue
          G1L = G1L*SL
          G2L = G2L*SL
          G3L = G3L*SL
          G4L = G4L*SL
          G5LR = G5LR*SL
          G5LI = G5LI*SL
          ALF1(L1) = G1L + G2L
          ALF2(L1) = G3L + G4L
          ALF3(L1) = G3L - G4L
          ALF4(L1) = G1L - G2L
          BET1(L1) = G5LR*2d0
          BET2(L1) = G5LI*2d0
          LMAX = L
          if (DABS(G1L) .lt. 1d-6) GO TO 500
300   end do
500   continue
      return
end subroutine

!****************************************************************

!   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
!   0.LE.N.LE.899

subroutine FACT
      real*8 F(900)
      common /FAC/ F
      F(1) = 0d0
      F(2) = 0d0
      do 2 I = 3, 900
          I1 = I - 1
          F(I) = F(I1) + 0.5d0*DLOG(DFLOAT(I1))
2     end do
      return
end subroutine

!************************************************************

!   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
!   0.LE.N.LE.899

subroutine SIGNUM
      real*8 SSIGN(900)
      common /SS/ SSIGN
      SSIGN(1) = 1d0
      do 2 N = 2, 899
          SSIGN(N) = -SSIGN(N - 1)
2     end do
      return
end subroutine

!******************************************************************
!
!   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
!   (N,M:N1,M1/NN,MM)
!   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
!   MM=M*K1+K2
!
!   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
!                               N.LE.NMAX
!                               N.GE.1
!                               N1.GE.0
!                               N1.LE.N+NMAX
!   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
!                                       COEFFICIENTS
!                               /M/.LE.N
!                               /M1/=/M*(K1-1)+K2/.LE.N1
!                               NN.LE.MIN(N+N1,NMAX)
!                               NN.GE.MAX(/MM/,/N-N1/)
!   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N

subroutine CCG(N, N1, NMAX, K1, K2, GG)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      real*8 GG(NPL1, NPN6), CD(0:NPN5), CU(0:NPN5)
      if (NMAX .le. NPN4 .AND.&
          & 0 .le. N1 .AND. &
          & N1 .le. NMAX + N .AND.&
          & N .ge. 1 .AND. N .le. NMAX) GO TO 1
      print 5001
      stop
5001  format(' ERROR IN SUBROUTINE CCG')
1     NNF = MIN0(N + N1, NMAX)
      MIN = NPN6 - N
      MF = NPN6 + N
      if (K1 .eq. 1 .and. K2 .eq. 0) MIN = NPN6
      do 100 MIND = MIN, MF
          M = MIND - NPN6
          MM = M*K1 + K2
          M1 = MM - M
          if (IABS(M1) .gt. N1) GO TO 90
          NNL = MAX0(IABS(MM), IABS(N - N1))
          if (NNL .gt. NNF) GO TO 90
          NNU = N + N1
          NNM = (NNU + NNL)*0.5d0
          if (NNU .eq. NNL) NNM = NNL
          call CCGIN(N, N1, M, MM, C)
          CU(NNL) = C
          if (NNL .eq. NNF) GO TO 50
          C2 = 0d0
          C1 = C
          do 7 NN = NNL + 1, MIN0(NNM, NNF)
              A = DFLOAT((NN + MM)*(NN - MM)*(N1 - N + NN))
              A = A*DFLOAT((N - N1 + NN)*(N + N1 - NN + 1)*(N + N1 + NN + 1))
              A = DFLOAT(4*NN*NN)/A
              A = A*DFLOAT((2*NN + 1)*(2*NN - 1))
              A = DSQRT(A)
              B = 0.5d0*DFLOAT(M - M1)
              D = 0d0
              if (NN .eq. 1) GO TO 5
              B = DFLOAT(2*NN*(NN - 1))
              B = DFLOAT((2*M - MM)*NN*(NN - 1) - MM*N*(N + 1) +&
                  & MM*N1*(N1 + 1))/B
              D = DFLOAT(4*(NN - 1)*(NN - 1))
              D = D*DFLOAT((2*NN - 3)*(2*NN - 1))
              D = DFLOAT((NN - MM - 1)*(NN + MM - 1)*(N1 - N + NN - 1))/D
              D = D*DFLOAT((N - N1 + NN - 1)*(N + N1 - NN + 2)*(N + N1 + NN))
              D = DSQRT(D)
5             C = A*(B*C1 - D*C2)
              C2 = C1
              C1 = C
              CU(NN) = C
7         continue
          if (NNF .le. NNM) GO TO 50
          call DIRECT(N, M, N1, M1, NNU, MM, C)
          CD(NNU) = C
          if (NNU .eq. NNM + 1) GO TO 50
          C2 = 0d0
          C1 = C
          do 12 NN = NNU - 1, NNM + 1, -1
              A = DFLOAT((NN - MM + 1)*(NN + MM + 1)*(N1 - N + NN + 1))
              A = A*DFLOAT((N - N1 + NN + 1)*(N + N1 - NN)*(N + N1 + NN + 2))
              A = DFLOAT(4*(NN + 1)*(NN + 1))/A
              A = A*DFLOAT((2*NN + 1)*(2*NN + 3))
              A = DSQRT(A)
              B = DFLOAT(2*(NN + 2)*(NN + 1))
              B = DFLOAT((2*M - MM)*(NN + 2)*(NN + 1) - MM*N*(N + 1)&
                  & + MM*N1*(N1 + 1))/B
              D = DFLOAT(4*(NN + 2)*(NN + 2))
              D = D*DFLOAT((2*NN + 5)*(2*NN + 3))
              D = DFLOAT((NN + MM + 2)*(NN - MM + 2)*(N1 - N + NN + 2))/D
              D = D*DFLOAT((N - N1 + NN + 2)*(N + N1 - NN - 1)*(N + N1 + NN + 3))
              D = DSQRT(D)
              C = A*(B*C1 - D*C2)
              C2 = C1
              C1 = C
              CD(NN) = C
12        continue
50        do 9 NN = NNL, NNF
              if (NN .le. NNM) GG(MIND, NN + 1) = CU(NN)
              if (NN .gt. NNM) GG(MIND, NN + 1) = CD(NN)
              ! write (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
9         continue
90        continue
100   end do
      return
end subroutine

!*********************************************************************

subroutine DIRECT(N, M, N1, M1, NN, MM, C)
      implicit real*8(A - H, O - Z)
      real*8 F(900)
      common /FAC/ F
      C = F(2*N + 1) + F(2*N1 + 1) + F(N + N1 + M + M1 + 1) + F(N + N1 - M - M1 + 1)
      C = C - F(2*(N + N1) + 1) - F(N + M + 1) - F(N - M + 1) - F(N1 + M1 + 1) - F(N1 - M1 + 1)
      C = DEXP(C)
      return
end subroutine

!*********************************************************************
!
!   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
!   G=(N,M:N1,MM-M/NN,MM)
!   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
!                               /M/.LE.N
!                               /MM-M/.LE.N1
!                               /MM/.LE.N+N1

subroutine CCGIN(N, N1, M, MM, G)
      implicit real*8(A - H, O - Z)
      real*8 F(900), SSIGN(900)
      common /SS/ SSIGN
      common /FAC/ F
      M1 = MM - M
      if (N .ge. IABS(M) .AND.&
          & N1 .ge. IABS(M1) .AND.&
          & IABS(MM) .le. (N + N1)) GO TO 1
      print 5001
      stop
5001  format(' ERROR IN SUBROUTINE CCGIN')
1     if (IABS(MM) .gt. IABS(N - N1)) GO TO 100
      L1 = N
      L2 = N1
      L3 = M
      if (N1 .le. N) GO TO 50
      K = N
      N = N1
      N1 = K
      K = M
      M = M1
      M1 = K
50    N2 = N*2
      M2 = M*2
      N12 = N1*2
      M12 = M1*2
      G = SSIGN(N1 + M1 + 1)&
          & *DEXP(F(N + M + 1) + F(N - M + 1) + F(N12 + 1) + F(N2 - N12 + 2) - F(N2 + 2)&
          & - F(N1 + M1 + 1) - F(N1 - M1 + 1) - F(N - N1 + MM + 1) - F(N - N1 - MM + 1))
      N = L1
      N1 = L2
      M = L3
      return
100   A = 1d0
      L1 = M
      L2 = MM
      if (MM .ge. 0) GO TO 150
      MM = -MM
      M = -M
      M1 = -M1
      A = SSIGN(MM + N + N1 + 1)
150   G = A*SSIGN(N + M + 1)&
          & *DEXP(F(2*MM + 2) + F(N + N1 - MM + 1) + F(N + M + 1) + F(N1 + M1 + 1)&
          & - F(N + N1 + MM + 2) - F(N - N1 + MM + 1) - F(-N + N1 + MM + 1) - F(N - M + 1)&
          & - F(N1 - M1 + 1))
      M = L1
      MM = L2
      return
end subroutine

!*****************************************************************

subroutine SAREA(D, RAT)
      implicit real*8(A - H, O - Z)
      if (D .ge. 1) GO TO 10
      E = DSQRT(1d0 - D*D)
      R = 0.5d0*(D**(2d0/3d0) + D**(-1d0/3d0)*DASIN(E)/E)
      R = DSQRT(R)
      RAT = 1d0/R
      return
10    E = DSQRT(1d0 - 1d0/(D*D))
      R = 0.25d0*(2d0*D**(2d0/3d0) + D**(-4d0/3d0)*DLOG((1d0 + E)/(1d0 - E))&
          & /E)
      R = DSQRT(R)
      RAT = 1d0/R
      return
end subroutine

!****************************************************************

subroutine SURFCH(N, E, RAT)
      implicit real*8(A - H, O - Z)
      real*8 X(60), W(60)
      DN = DFLOAT(N)
      E2 = E*E
      EN = E*DN
      NG = 60
      call GAUSS(NG, 0, 0, X, W)
      S = 0d0
      V = 0d0
      do 10 I = 1, NG
          XI = X(I)
          DX = DACOS(XI)
          DXN = DN*DX
          DS = DSIN(DX)
          DSN = DSIN(DXN)
          DCN = DCOS(DXN)
          A = 1d0 + E*DCN
          A2 = A*A
          ENS = EN*DSN
          S = S + W(I)*A*DSQRT(A2 + ENS*ENS)
          V = V + W(I)*(DS*A + XI*ENS)*DS*A2
10    end do
      RS = DSQRT(S*0.5d0)
      RV = (V*3d0/4d0)**(1d0/3d0)
      RAT = RV/RS
      return
end subroutine

!********************************************************************

subroutine SAREAC(EPS, RAT)
      implicit real*8(A - H, O - Z)
      RAT = (1.5d0/EPS)**(1d0/3d0)
      RAT = RAT/DSQRT((EPS + 2d0)/(2d0*EPS))
      return
end subroutine

!********************************************************************

!  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
!  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B

subroutine POWER(A, B, R1, R2)
      implicit real*8(A - H, O - Z)
      external F
      common AA, BB
      AA = A
      BB = B
      AX = 1d-5
      BX = A - 1d-5
      R1 = ZEROIN(AX, BX, F, 0d0)
      R2 = (1d0 + B)*2d0*A - R1
      return
end subroutine

!***********************************************************************

double precision function ZEROIN(AX, BX, F, TOL)
      implicit real*8(A - H, O - Z)
      EPS = 1d0
10    EPS = 0.5d0*EPS
      TOL1 = 1d0 + EPS
      if (TOL1 .gt. 1d0) GO TO 10
15    A = AX
      B = BX
      FA = F(A)
      FB = F(B)
20    C = A
      FC = FA
      D = B - A
      E = D
30    if (DABS(FC) .ge. DABS(FB)) GO TO 40
35    A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
40    TOL1 = 2d0*EPS*DABS(B) + 0.5d0*TOL
      XM = 0.5d0*(C - B)
      if (DABS(XM) .le. TOL1) GO TO 90
44    if (FB .eq. 0d0) GO TO 90
45    if (DABS(E) .lt. TOL1) GO TO 70
46    if (DABS(FA) .le. DABS(FB)) GO TO 70
47    if (A .ne. C) GO TO 50
48    S = FB/FA
      P = 2d0*XM*S
      Q = 1d0 - S
      GO TO 60
50    Q = FA/FC
      R = FB/FC
      S = FB/FA
      P = S*(2d0*XM*Q*(Q - R) - (B - A)*(R - 1d0))
      Q = (Q - 1d0)*(R - 1d0)*(S - 1d0)
60    if (P .gt. 0d0) Q = -Q
      P = DABS(P)
      if ((2d0*P) .ge. (3d0*XM*Q - DABS(TOL1*Q))) GO TO 70
64    if (P .ge. DABS(0.5d0*E*Q)) GO TO 70
65    E = D
      D = P/Q
      GO TO 80
70    D = XM
      E = D
80    A = B
      FA = FB
      if (DABS(D) .gt. TOL1) B = B + D
      if (DABS(D) .le. TOL1) B = B + DSIGN(TOL1, XM)
      FB = F(B)
      if ((FB*(FC/DABS(FC))) .gt. 0d0) GO TO 20
85    GO TO 30
90    ZEROIN = B
      return
end function

!***********************************************************************

double precision function F(R1)
      implicit real*8(A - H, O - Z)
      common A, B
      R2 = (1d0 + B)*2d0*A - R1
      F = (R2 - R1)/DLOG(R2/R1) - A
      return
end function

!**********************************************************************
!    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!    N - NUMBER OF POINTS                                             *
!    Z - DIVISION POINTS                                              *
!    W - WEIGHTS                                                      *
!**********************************************************************

subroutine GAUSS(N, IND1, IND2, Z, W)
      implicit real*8(A - H, P - Z)
      real*8 Z(N), W(N)
      A = 1d0
      B = 2d0
      C = 3d0
      IND = mod(N, 2)
      K = N/2 + IND
      F = DFLOAT(N)
      do 100 I = 1, K
          M = N + 1 - I
          if (I .eq. 1) X = A - B/((F + A)*F)
          if (I .eq. 2) X = (Z(N) - A)*4d0 + Z(N)
          if (I .eq. 3) X = (Z(N - 1) - Z(N))*1.6d0 + Z(N - 1)
          if (I .gt. 3) X = (Z(M + 1) - Z(M + 2))*C + Z(M + 3)
          if (I .eq. K .and. IND .eq. 1) X = 0d0
          NITER = 0
          CHECK = 1d-16
10        PB = 1d0
          NITER = NITER + 1
          if (NITER .le. 100) GO TO 15
          CHECK = CHECK*10d0
15        PC = X
          DJ = A
          do 20 J = 2, N
              DJ = DJ + A
              PA = PB
              PB = PC
              PC = X*PB + (X*PB - PA)*(DJ - A)/DJ
20        continue
          PA = A/((PB - X*PC)*F)
          PB = PA*PC*(A - X*X)
          X = X - PB
          if (DABS(PB) .gt. CHECK*DABS(X)) GO TO 10
          Z(M) = X
          W(M) = PA*PA*(A - X*X)
          if (IND1 .eq. 0) W(M) = B*W(M)
          if (I .eq. K .and. IND .eq. 1) GO TO 100
          Z(I) = -Z(M)
          W(I) = W(M)
100   end do
      if (IND2 .ne. 1) GO TO 110
      print 1100, N
1100  format(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',&
          & ' OF ', I4, '-TH ORDER')
      do 105 I = 1, K
          ZZ = -Z(I)
          print 1200, I, ZZ, I, W(I)
105   continue
1200  format(' ', 4x, 'X(', I4, ') = ', F17.14, 5x, 'W(', I4, ') = ', F17.14)
      GO TO 115
110   continue
!     PRINT 1300,N
1300  format(' GAUSSIAN QUADRATURE FORMULA OF ', I4, '-TH ORDER IS USED')
115   continue
      if (IND1 .eq. 0) GO TO 140
      do 120 I = 1, N
          Z(I) = (A + Z(I))/B
120   continue
140   continue
      return
end subroutine

!****************************************************************

subroutine DISTRB(NNK, YY, WY, NDISTR, AA, BB, GAM, R1, R2, REFF,&
      & VEFF, PI)
      implicit real*8(A - H, O - Z)
      real*8 YY(NNK), WY(NNK)
      if (NDISTR .eq. 2) GO TO 100
      if (NDISTR .eq. 3) GO TO 200
      if (NDISTR .eq. 4) GO TO 300
      if (NDISTR .eq. 5) GO TO 360
      print 1001, AA, BB, GAM
1001  format('MODIFIED GAMMA DISTRIBUTION, alpha=', F6.4, '  r_c=',&
          & F6.4, '  gamma=', F6.4)
      A2 = AA/GAM
      DB = 1d0/BB
      do 50 I = 1, NNK
          X = YY(I)
          Y = X**AA
          X = X*DB
          Y = Y*DEXP(-A2*(X**GAM))
          WY(I) = WY(I)*Y
50    end do
      GO TO 400
100   print 1002, AA, BB
1002  format('LOG-NORMAL DISTRIBUTION, r_g=', F8.4,&
          & '  [ln(sigma_g)]**2=', F6.4)
      DA = 1d0/AA
      do 150 I = 1, NNK
          X = YY(I)
          Y = DLOG(X*DA)
          Y = DEXP(-Y*Y*0.5d0/BB)/X
          WY(I) = WY(I)*Y
150   end do
      GO TO 400
200   print 1003
1003  format('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')
      do 250 I = 1, NNK
          X = YY(I)
          WY(I) = WY(I)/(X*X*X)
250   end do
      GO TO 400
300   continue
      ! print 1004,AA,BB
1004  format('GAMMA DISTRIBUTION,  a=', F6.3, '  b=', F6.4)
      B2 = (1d0 - 3d0*BB)/BB
      DAB = 1d0/(AA*BB)
      do 350 I = 1, NNK
          X = YY(I)
          X = (X**B2)*DEXP(-X*DAB)
          WY(I) = WY(I)*X
350   end do
      GO TO 400
360   print 1005, BB
1005  format('MODIFIED POWER LAW DISTRIBUTION,  alpha=', D10.4)
      do 370 I = 1, NNK
          X = YY(I)
          if (X .le. R1) WY(I) = WY(I)
          if (X .gt. R1) WY(I) = WY(I)*(X/R1)**BB
370   end do
400   continue
      SUM = 0d0
      do 450 I = 1, NNK
          SUM = SUM + WY(I)
450   end do
      SUM = 1d0/SUM
      do 500 I = 1, NNK
          WY(I) = WY(I)*SUM
500   end do
      G = 0d0
      do 550 I = 1, NNK
          X = YY(I)
          G = G + X*X*WY(I)
550   end do
      REFF = 0d0
      do 600 I = 1, NNK
          X = YY(I)
          REFF = REFF + X*X*X*WY(I)
600   end do
      REFF = REFF/G
      VEFF = 0d0
      do 650 I = 1, NNK
          X = YY(I)
          XI = X - REFF
          VEFF = VEFF + XI*XI*X*X*WY(I)
650   end do
      VEFF = VEFF/(G*REFF*REFF)
      return
end subroutine

!*************************************************************

subroutine HOVENR(L1, A1, A2, A3, A4, B1, B2)
      implicit real*8(A - H, O - Z)
      real*8 A1(L1), A2(L1), A3(L1), A4(L1), B1(L1), B2(L1)
      do 100 L = 1, L1
          KONTR = 1
          LL = L - 1
          DL = DFLOAT(LL)*2d0 + 1d0
          DDL = DL*0.48d0
          AA1 = A1(L)
          AA2 = A2(L)
          AA3 = A3(L)
          AA4 = A4(L)
          BB1 = B1(L)
          BB2 = B2(L)
          if (LL .ge. 1 .and. DABS(AA1) .ge. DL) KONTR = 2
          if (DABS(AA2) .ge. DL) KONTR = 2
          if (DABS(AA3) .ge. DL) KONTR = 2
          if (DABS(AA4) .ge. DL) KONTR = 2
          if (DABS(BB1) .ge. DDL) KONTR = 2
          if (DABS(BB2) .ge. DDL) KONTR = 2
          if (KONTR .eq. 2) print 3000, LL
          C = -0.1d0
          do 50 I = 1, 11
              C = C + 0.1d0
              CC = C*C
              C1 = CC*BB2*BB2
              C2 = C*AA4
              C3 = C*AA3
              if ((DL - C*AA1)*(DL - C*AA2) - CC*BB1*BB1 .le. -1d-4) KONTR = 2
              if ((DL - C2)*(DL - C3) + C1 .le. -1d-4) KONTR = 2
              if ((DL + C2)*(DL - C3) - C1 .le. -1d-4) KONTR = 2
              if ((DL - C2)*(DL + C3) - C1 .le. -1d-4) KONTR = 2
              if (KONTR .eq. 2) print 4000, LL, C
50        continue
100   end do
      ! if(KONTR.EQ.1) PRINT 2000
2000  format('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
3000  format('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=', I3)
4000  format('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=', I3,&
          & '   A=', D9.2)
      return
end subroutine

!****************************************************************
!
!    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
!    COEFFICIENTS
!
!    A1,...,B2 - EXPANSION COEFFICIENTS
!    LMAX   - NUMBER OF COEFFICIENTS MINUS 1
!    Nang   - NUMBER OF SCATTERING ANGLES
!    muang  - array of cosine of the angles
!    scatmat - six scattering matrix elements (F11,F12,F33,F34,F22,F44)
!               for each angle.

subroutine MATR(A1, A2, A3, A4, B1, B2, LMAX, Nang, muang, scatmat)

      include 'shdScat_tmx_par.f90'

      implicit real*8(A - H, O - Z)
      integer Nang
      real*8 muang(Nang), scatmat(6, Nang)
      real*8 A1(NPL), A2(NPL), A3(NPL), A4(NPL), B1(NPL), B2(NPL)

      ! N=NPNA
      ! DN=1D0/DFLOAT(N-1)
      ! DA=DACOS(-1D0)*DN
      ! DB=180D0*DN
      L1MAX = LMAX + 1
      ! PRINT 1000
! 1000  FORMAT(' ')
      ! PRINT 1001
! 1001  FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',&
          ! & 6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
      ! DO 10 L1=1,L1MAX
          ! L=L1-1
          ! PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
! 10   CONTINUE
! 1002  FORMAT(' ',I3,6F12.5)
      ! TB=-DB
      ! TAA=-DA
      ! PRINT 1000
      ! PRINT 1003
! 1003  FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',&
         ! & 8X,'F44',8X,'F12',8X,'F34')

      D6 = DSQRT(6d0)*0.25d0
      do 500 I1 = 1, Nang
          ! TAA=TAA+DA
          ! TB=TB+DB
          ! U=DCOS(TAA)
          U = muang(I1)
          F11 = 0d0
          F2 = 0d0
          F3 = 0d0
          F44 = 0d0
          F12 = 0d0
          F34 = 0d0
          P1 = 0d0
          P2 = 0d0
          P3 = 0d0
          P4 = 0d0
          PP1 = 1d0
          PP2 = 0.25d0*(1d0 + U)*(1d0 + U)
          PP3 = 0.25d0*(1d0 - U)*(1d0 - U)
          PP4 = D6*(U*U - 1d0)
          do 400 L1 = 1, L1MAX
              L = L1 - 1
              DL = DFLOAT(L)
              DL1 = DFLOAT(L1)
              F11 = F11 + A1(L1)*PP1
              F44 = F44 + A4(L1)*PP1
              if (L .eq. LMAX) GO TO 350
              PL1 = DFLOAT(2*L + 1)
              P = (PL1*U*PP1 - DL*P1)/DL1
              P1 = PP1
              PP1 = P
350           if (L .lt. 2) GO TO 400
              F2 = F2 + (A2(L1) + A3(L1))*PP2
              F3 = F3 + (A2(L1) - A3(L1))*PP3
              F12 = F12 + B1(L1)*PP4
              F34 = F34 + B2(L1)*PP4
              if (L .eq. LMAX) GO TO 400
              PL2 = DFLOAT(L*L1)*U
              PL3 = DFLOAT(L1*(L*L - 4))
              PL4 = 1d0/DFLOAT(L*(L1*L1 - 4))
              P = (PL1*(PL2 - 4d0)*PP2 - PL3*P2)*PL4
              P2 = PP2
              PP2 = P
              P = (PL1*(PL2 + 4d0)*PP3 - PL3*P3)*PL4
              P3 = PP3
              PP3 = P
              P = (PL1*U*PP4 - DSQRT(DFLOAT(L*L - 4))*P4)/DSQRT(DFLOAT(L1*L1 - 4))
              P4 = PP4
              PP4 = P
400       continue
          F22 = (F2 + F3)*0.5d0
          F33 = (F2 - F3)*0.5d0
          ! print 1004,TB,F11,F22,F33,F44,F12,F34
          scatmat(1, I1) = F11
          scatmat(2, I1) = F12
          scatmat(3, I1) = F33
          scatmat(4, I1) = F34
          scatmat(5, I1) = F22
          scatmat(6, I1) = F44
500   end do
      ! print 1000
! 1004  format(' ',F6.2,6F11.4)
      return
end subroutine

subroutine WIGNERFCT(X, NRANK, M, M1, DMM1)
      ! The routine computes the DMM1N vector coefficients for
      ! M >= 0, M1 >= 0,  N = 0,...,NRANK, and -1 < X= COS(BETA) < 1
      implicit none
      integer NRANK, M, M1
      double precision X
      double precision DMM1(0:NRANK)
      integer N0, N
      double precision FACT1, FACT2, DMM1_N0

      N0 = max(M, M1)
      do N = 0, NRANK
          DMM1(n) = 0.0d0
      end do
      if (N0 .eq. 0) then
          DMM1(0) = 1.d0
          DMM1(1) = X
          do N = 1, NRANK - 1
              FACT1 = dble(2*N + 1)*X/dble(N + 1)
              FACT2 = dble(N)/dble(N + 1)
              DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
          end do
      else
          DMM1(N0) = DMM1_N0(X, M, M1)
          do N = N0, NRANK - 1
              FACT1 = dble(N*(N + 1))*X - dble(M*M1)
              FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M**2))
              FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M1**2))
              FACT1 = FACT1*dble(2*N + 1)/dble(N)
              FACT2 = DSQRT(dble(N**2 - M**2))&
                  & *DSQRT(dble(N**2 - M1**2))
              FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
              FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M1**2))
              FACT2 = FACT2*dble(N + 1)/dble(N)
              DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
          end do
      end if
      return
end subroutine

double precision function DMM1_N0(X, M, M1)
      ! The routine computes the Wigner functions for
      ! M >= 0, M1 >= 0, AND N0 = MAX(M,M1)
      implicit none
      integer M, M1
      double precision X
      integer P, MAXM, MINM
      double precision PROD, FACT, CMM1, DMM1N0

      if (M .eq. M1) then
          FACT = ((1.d0 + X)/2.d0)**M
          DMM1N0 = FACT
      else
          if (M1 .gt. M) then
              CMM1 = 1.d0
          else
              CMM1 = (-1.d0)**(M - M1)
          end if
          MAXM = max(M, M1)
          MINM = min(M, M1)
          PROD = 1.d0
          do P = 1, MAXM - MINM
              FACT = DSQRT(dble(M + M1 + P)/dble(P))
              PROD = PROD*FACT
          end do
          FACT = DSQRT((1.d0 - X)/2.d0)
          DMM1N0 = CMM1*PROD*FACT**(MAXM - MINM)
          FACT = DSQRT((1.d0 + X)/2.d0)
          DMM1N0 = DMM1N0*FACT**(MAXM + MINM)
      end if
      DMM1_N0 = DMM1N0
      return
end function
