!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!   --
!
!   Copyright (C) 2025-
!   EaR³T Community: https://discord.gg/ntqsguwaWv
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
!   SHDOM was originally developed by K. Franklin Evans and is available
!   at https://coloradolinux.com/shdom.
!
!   This file was adapted from HPARX by
!   Hironobu Iwabuchi(hironobu.iwabuchi.c4@tohoku.ac.jp), Souichiro Hioki,
!   and Rintaro Okamura.
!
! History:
!   Created by: Hong Chen
!               hong.chen-1@colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

module globals_const

  implicit none
  public            ! all in this module are public!

  ! Programmer-dependent kinds of real/complex variables
  !// Note: The R_ is default and RD_ should be used only when needed.
  !    One's choice of precision for R_ and RD_ may depend on required accuracy and used memory size.
  !integer,   parameter :: R_  = selected_real_kind(6)  ! default precision
  integer,   parameter :: R_  = selected_real_kind(13) ! default precision
  integer,   parameter :: RD_ = selected_real_kind(13) ! higher  precision

  ! System-specific kinds of variables (no freedom for programmers)
  integer,   parameter :: R4_ = selected_real_kind(6)  ! 4-byte real
  integer,   parameter :: R8_ = selected_real_kind(13) ! 8-byte real

  integer,   parameter :: I_  = selected_int_kind(9)   ! 4-byte integer
  integer,   parameter :: I4_ = selected_int_kind(9)   ! 4-byte integer
  integer,   parameter :: I8_ = selected_int_kind(18)  ! 8-byte integer

  ! Numerical parameters for real variables
  real(R_),  parameter :: REPS_   = epsilon(1.0_R_)    ! 1.0e-7, typically
  real(R_),  parameter :: RTINY_  = tiny(1.0_R_)       ! 1.0e-37
  real(R_),  parameter :: RHUGE_  = huge(1.0_R_)       ! 1.0e+37
  real(R_),  parameter :: RSML_   = RTINY_  * 3.0_R_   ! 1.0e-37
  real(R_),  parameter :: RLRG_   = RHUGE_  / 3.0_R_   ! 1.0e+37
  real(RD_), parameter :: RDEPS_  = epsilon(1.0_RD_)   ! 1.0e-15, typically
  real(RD_), parameter :: RDTINY_ = tiny(1.0_RD_)      ! 1.0e-307
  real(RD_), parameter :: RDHUGE_ = huge(1.0_RD_)      ! 1.0e+307
  real(RD_), parameter :: RDSML_  = RDTINY_ * 3.0_RD_  ! 1.0e-307
  real(RD_), parameter :: RDLRG_  = RDHUGE_ / 3.0_RD_  ! 1.0e+307

  ! Numerical parameters for integer variables
  integer,   parameter :: IMAX_ =  2147483647          ! max integer of default kind
  integer,   parameter :: IMIN_ = -2147483647 - 1      ! min integer of default kind

  ! Mathematical constants
  real(R_),  parameter :: PI_     = 3.141592653589793238462643383279502884197_R_   ! pi
  real(R_),  parameter :: PIH_    = 1.570796326794896619231321691639751442098_R_   ! pi/2
  real(R_),  parameter :: PI2_    = 6.283185307179586476925286766559005768394_R_   ! pi*2
  real(RD_), parameter :: DPI_    = 3.141592653589793238462643383279502884197_RD_  ! pi
  real(RD_), parameter :: DPIH_   = 1.570796326794896619231321691639751442098_RD_  ! pi/2
  real(RD_), parameter :: DPI2_   = 6.283185307179586476925286766559005768394_RD_  ! pi*2
  real(R_) , parameter :: OPI_    = 1.0_R_  / PI_      ! 1/PI
  real(RD_), parameter :: DOPI_   = 1.0_RD_ / DPI_     ! 1/PI
  real(R_) , parameter :: O4PI_   = 1.0_R_  / (4.0_R_*PI_)      ! 1/(4*PI)
  real(RD_), parameter :: DO4PI_  = 1.0_RD_ / (4.0_RD_*DPI_)    ! 1/(4*PI)
  real(R_),  parameter :: DTOR_   = PI_ / 180.0_R_     ! degree-to-radian factor

  real(R_),  parameter :: SQRT2_  = 1.414213562373095048801688724209698078569_R_   ! sqrt(2)
  real(R_),  parameter :: SQRT3_  = 1.732050807568877293527446341505872366942_R_   ! sqrt(3)
  real(R_),  parameter :: SQRT2A_ = 0.7071067811865475244008443621048490392848_R_  ! 1/sqrt(2)
  real(R_),  parameter :: SQRT3A_ = 0.5773502691896257645091487805019574556476_R_  ! 1/sqrt(3)

  real(R_),  parameter :: FRAC12_  = 1.0_R_ / 2.0_R_    ! 1/2
  real(R_),  parameter :: FRAC13_  = 1.0_R_ / 3.0_R_    ! 1/3
  real(R_),  parameter :: FRAC23_  = 2.0_R_ / 3.0_R_    ! 2/3
  real(R_),  parameter :: FRAC43_  = 4.0_R_ / 3.0_R_    ! 4/3
  real(R_),  parameter :: FRAC53_  = 5.0_R_ / 3.0_R_    ! 5/3
  real(R_),  parameter :: FRAC16_  = 1.0_R_ / 6.0_R_    ! 1/6
  real(R_),  parameter :: FRAC112_ = 1.0_R_ / 12.0_R_   ! 1/12

  ! Physical constants (SI unit)
  real(R_),  parameter :: CSPEED_ = 2.99792458e+8_R_     ! speed of light (m/s)
  real(R_),  parameter :: PLANCK_ = 6.6260689633e-34_R_  ! Planck's constant (J*s)
  real(R_),  parameter :: BOLTZ_  = 1.380650424e-23_R_   ! Boltzmann's constant (J/K)
  real(R_),  parameter :: STEBOL_ = 5.67040040e-8_R_     ! Stefan-Boltzmann's constant (W/m^2/K^4)
  real(R_),  parameter :: AVOGAD_ = 6.0221417930e+23_R_  ! Avogadro number (/mol)
  real(R_),  parameter :: GASCON_ = 8.31447215_R_        ! universal gas constant (J/mol/K)
  real(R_),  parameter :: VOLMOL_ = 2.241399639e-2_R_    ! molar volume (m^3/mol) of ideal gas (0 deg.C, 1 atm)
  real(R_),  parameter :: EGRAVE_ = 9.7803267715_R_      ! earth equatorial gravity (m/s^2)
  real(R_),  parameter :: EGRAV0_ = 9.80665_R_           ! earth standard   gravity (m/s^2)
  real(R_),  parameter :: WGTVAP_ = 18.016e-3_R_         ! molar mass of water vapor    (H2O) (kg/mol)
  real(R_),  parameter :: WGTCO2_ = 44.01e-3_R_          ! molar mass of carbon dioxide (CO2) (kg/mol)
  real(R_),  parameter :: WGTOZN_ = 48.00e-3_R_          ! molar mass of ozone          (O3)  (kg/mol)
  real(R_),  parameter :: WGTDRY_ = 28.96454e-3_R_       ! molar mass of dry air (kg/mol)
  real(R_),  parameter :: RHOWAT0_ = 0.99984e+3_R_       ! density (kg/m^3) of liquid water (0 deg.C)
  real(R_),  parameter :: RHOICE0_ = 0.917e+3_R_         ! density (kg/m^3) of water ice (0 deg.C)

end module globals_const


module globals_cfg

    implicit none
    public

    integer, parameter :: charMaxLen_ = 256 ! file name max length (number of characters)

    real(kind=8), save :: accSolv_
    real(kind=4), save :: memTotSysMB_ = 16000.0
    real(kind=4), save :: memPerCpuMB_
    real(kind=4), save :: planetRadius_ = 6371.009

    integer(kind=1), save :: verbose_ = 5
    !Severity    Level  Description
    !Emergency     0      System unusable
    !Alert         1      Action must be taken immediately
    !Critical      2      Critical conditions
    !Error         3      Error conditions
    !Warning       4      Warning conditions
    !Notice        5      Normal but significant condition
    !Informational 6      Informational messages
    !All           7      All-level messages

    ! Namelist input
    logical, save :: doNmlInp_ = .true.

    ! High order radiance (use more SH terms)
    logical, save :: doHighOrder_

    ! Boundary conditions (BCFLAG)
    logical, save :: isOpenX_, isOpenYx_, isMpiX_, isMpiY_, isMpi_

    ! Indepedent pixel conditions (IPFLAG)
    logical, save :: is1DX_, is1DYx_, is1D_, is3Ddir1Ddif_

    ! Pure 1D calculation
    logical, save :: is1D1_

    ! Polarization
    logical, save :: isPol_, isPol4_

    ! Radiation sources (SRCTYPE)
    logical, save :: isSolar_, isThermal_

end module globals_cfg


module globals_type

    use globals_cfg, only: charMaxLen_

    implicit none
    public

    integer, parameter, private :: I__ = selected_int_kind(9)
    integer, parameter, private :: R__ = selected_real_kind(6)

    type :: shd_type_rtm

        character(len=charMaxLen_) :: runName = 'NONE'
        character(len=1) :: type = ' '

        integer(kind=I__) :: flag_bc = -999
        integer(kind=I__) :: flag_ip = -999

        integer(kind=I__) :: Nx = -999
        integer(kind=I__) :: Ny = -999
        integer(kind=I__) :: Nz = -999

        integer(kind=I__) :: maxNz   = -999
        integer(kind=I__) :: maxNpts = -999
        integer(kind=I__) :: maxNptsLeg = -999

        real(kind=R__) :: dx = -999.0
        real(kind=R__) :: dy = -999.0
        real(kind=R__) :: xStart  = -999.0
        real(kind=R__) :: yStart  = -999.0
        real(kind=R__) :: xDomain = -999.0
        real(kind=R__) :: yDomain = -999.0

        real(kind=R__),  allocatable :: x(:)
        real(kind=R__),  allocatable :: y(:)
        real(kind=R__),  allocatable :: z(:)

        logical :: isAllocated  = .false.

    end type shd_type_rtm

    type :: shd_type_grid

        integer(kind=I__) :: Nx = -999
        integer(kind=I__) :: Ny = -999
        integer(kind=I__) :: Nz = -999

        integer(kind=I__) :: Npts0  = -999
        integer(kind=I__) :: Ncell0 = -999

        integer(kind=I__) :: Npts  = -999
        integer(kind=I__) :: Ncell = -999

        integer(kind=I__) :: maxNpts  = -999
        integer(kind=I__) :: maxNcell = -999

        integer(kind=I__), pointer:: i6Cell(:, :) => null()
        integer(kind=I__), pointer:: i2Tree(:, :) => null()
        integer(kind=2),   pointer:: flag(:) => null()

        logical :: isAllocated  = .false.

    end type shd_type_grid

    ! SHDO
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_shdo

        integer(kind=I__) :: Nmu  = -999
        integer(kind=I__) :: Nphi = -999
        integer(kind=I__) :: Nang = -999
        integer(kind=I__) :: Nlm  = -999
        integer(kind=I__) :: Nsh  = -999
        integer(kind=I__) :: Nstokes = -999
        integer(kind=I__) :: NlegSet = -999

        integer(kind=I__) :: maxNl = -999
        integer(kind=I__) :: maxNm = -999
        integer(kind=I__) :: maxNphi0  = -999
        integer(kind=I__) :: maxNptsSH = -999
        integer(kind=I__) :: maxNptsDO = -999

        integer(kind=I__), allocatable :: Nphi0(:)
        real(kind=R__),    allocatable :: mu(:)
        real(kind=R__),    allocatable :: phi(:, :)
        real(kind=R__),    allocatable :: wgt(:, :)

        logical :: isAllocated  = .false.

    end type shd_type_shdo
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SRC (radiation source, e.g., solar, thermal, or both)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_src

        character(len=1) :: type = ' '
        character(len=1) :: unit = ' '

        real(kind=R__), dimension(2) :: waveNum = [-999.0, -999.0]
        real(kind=R__) :: avgWaveNum = -999.0

        real(kind=R__), dimension(2) :: waveLen = [-999.0, -999.0]
        real(kind=R__) :: avgWaveLen = -999.0

        real(kind=R__) :: avgMu = -999.0
        real(kind=R__) :: avgAz = -999.0

        real(kind=R__) :: inpSolar  = -999.0 ! solar read from input file
        real(kind=R__) :: avgSolar0 = -999.0 ! solar at every "g", scaled with avgMu and solarFac
        real(kind=R__) :: avgSolar  = -999.0 ! averaged solar from all "g"s of correlated-k

        real(kind=R__) :: skyRad = -999.0

    end type shd_type_src
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! PRP
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_prp

        character(len=charMaxLen_) :: fname = 'NONE'
        character(len=1) :: type = ' '

        integer(kind=I__) :: Nx = -999
        integer(kind=I__) :: Ny = -999
        integer(kind=I__) :: Nz = -999
        integer(kind=I__) :: Npha = -999
        integer(kind=I__) :: Nleg = -999
        integer(kind=I__) :: NlegSet = -999

        integer(kind=I__) :: maxNz   = -999
        integer(kind=I__) :: maxNpts = -999
        integer(kind=I__) :: maxNleg = -999
        integer(kind=I__) :: maxNptsLeg = -999

        real(kind=R__) :: maxAsym = -999.0

        real(kind=R__) :: dx = -999.0
        real(kind=R__) :: dy = -999.0
        real(kind=R__) :: xStart  = -999.0
        real(kind=R__) :: yStart  = -999.0
        real(kind=R__) :: xDomain = -999.0
        real(kind=R__) :: yDomain = -999.0

        real(kind=R__), allocatable :: z(:)

        logical :: isAllocated  = .false.

    end type shd_type_prp
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SFC
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_sfc

        character(len=charMaxLen_) :: fname = 'NONE'
        character(len=2) :: type = '  '

        integer(kind=I__) :: Nx = -999
        integer(kind=I__) :: Ny = -999
        integer(kind=I__) :: Npar = -999

        integer(kind=I__) :: maxNpts = -999
        integer(kind=I__) :: maxNpar = -999

        real(kind=R__) :: dx = -999.0
        real(kind=R__) :: dy = -999.0
        real(kind=R__) :: avgTemp = -999.0
        real(kind=R__) :: avgAlb  = -999.0

        logical :: isUniform    = .false.
        logical :: isLambertian = .false.
        logical :: isAllocated  = .false.

    end type shd_type_sfc
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! CKD
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_ckd

        character(len=charMaxLen_) :: fname = 'NONE'

        integer(kind=I__) :: Ng = -999
        integer(kind=I__) :: Nz = -999
        integer(kind=I__) :: ig = -999
        integer(kind=I__) :: maxNg = -999
        integer(kind=I__) :: maxNz = -999

        real(kind=R__), dimension(2) :: waveNum = [-999.0, -999.0]
        real(kind=R__) :: avgWaveNum = -999.0

        real(kind=R__), dimension(2) :: waveLen = [-999.0, -999.0]
        real(kind=R__) :: avgWaveLen = -999.0

        real(kind=R__), allocatable :: z(:)

        logical :: isAllocated  = .false.

    end type shd_type_ckd
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SEN
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type :: shd_type_sen

        character(len=charMaxLen_) :: fname = 'NONE'
        integer(kind=I__) :: Npts = -999
        integer(kind=I__) :: Npar = -999

        logical :: isAllocated  = .false.

    end type shd_type_sen
    !╰────────────────────────────────────────────────────────────────────────────╯!

end module globals_type


module globals_obj

    use globals_type, only: shd_type_rtm, &
                            shd_type_grid, &
                            shd_type_shdo, &
                            shd_type_src, &
                            shd_type_prp, &
                            shd_type_sfc, &
                            shd_type_sfc, &
                            shd_type_ckd, &
                            shd_type_sen

    implicit none
    public

    integer, parameter, private :: I__ = selected_int_kind(9)
    integer, parameter, private :: R__ = selected_real_kind(6)


    ! RTM
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_rtm),  save :: rtm_obj_

    real(kind=R__),  allocatable, save ::   temp_Rtm_(:), &    ! (grid_obj%maxNpts)
                                             ext_Rtm_(:), &    ! (grid_obj%maxNpts)
                                             ssa_Rtm_(:), &    ! (grid_obj%maxNpts)
                                             leg_Rtm_(:, :, :) ! (prp_obj%NlegSet, prp_obj%maxNleg, prp_obj%Npha)
    integer(kind=2), allocatable, save ::   ipha_Rtm_(:)       ! (grid_obj%maxNpts)

    real(kind=R__),  allocatable, save :: planck_Rtm_(:)       ! (grid_obj%maxNpts)
    real(kind=R__),  allocatable, save ::   brdf_Rtm_(:, :)    ! (sfc_obj%maxNpar, xxx_obj%xxx)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Grid (and cells)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_grid), save :: grid_obj_

    real(kind=R__),    allocatable, save ::   xyz_Grid_(:, :)     ! (3, grid_obj%maxNpts)
    integer(kind=I__), allocatable, save :: i8Pts_Grid_(:, :) ! (8, grid_obj%maxNpts)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SHDO
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_shdo), save :: shdo_obj_
    integer(kind=I__), allocatable, save :: iradI_Shdo_(:)    ! (grid_obj%maxNpts+2)
    real(kind=R__),    allocatable, save ::  radI_Shdo_(:, :) ! (shdo_obj_%Nstokes, shdo_obj_%maxNptsSH+grid_obj%maxNpts)
    integer(kind=I__), allocatable, save :: isrcJ_Shdo_(:)    ! (grid_obj%maxNpts+1)
    real(kind=R__),    allocatable, save ::  srcJ_Shdo_(:, :) ! (shdo_obj%Nstokes, shdo_obj%maxNptsSH)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SRC (radiation source, e.g., solar, thermal, or both)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_src),  save :: src_obj_
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! PRP
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_prp), save :: prp_obj_

    real(kind=R__),  allocatable, save :: temp_Prp_(:),  &   ! (prp_obj%maxNpts)
                                           ext_Prp_(:),  &   ! (prp_obj%maxNpts)
                                           ssa_Prp_(:),  &   ! (prp_obj%maxNpts)
                                           leg_Prp_(:, :, :) ! (prp_obj%NlegSet, prp_obj%maxNleg, prp_obj%Npha)
    integer(kind=2), allocatable, save :: ipha_Prp_(:)       ! (prp_obj%maxNpts)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SFC
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sfc), save :: sfc_obj_
    real(kind=R__), allocatable, save :: brdf_Sfc_(:, :)  ! (sfc_obj%maxNpar, sfc_obj%maxNpts)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! CKD
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_ckd), save :: ckd_obj_
    real(kind=R__), allocatable, save :: solar_Ckd_(:), & ! (ckd_obj%Ng)
                                           wgt_Ckd_(:), & ! (ckd_obj%Ng)
                                          kabs_Ckd_(:, :) ! (ckd_obj%Nz, ckd_obj%Ng)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! SEN
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sen), save :: sen_obj_
    real(kind=R__), allocatable, save :: view_Sen_(:, :)  ! (sen_obj%Npar, sen_obj%Npts)
    !╰────────────────────────────────────────────────────────────────────────────╯!

end module globals_obj
