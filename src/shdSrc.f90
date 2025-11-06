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
! History:
!   Created by:  K. Franklin Evans
!                evans@nit.colorado.edu
!
!   Modified by: Hong Chen
!                hong.chen-1@colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

module shd_src

    use globals_const, only: RD_, DPI_, O4PI_, DO4PI_, &
                             FRAC12_, FRAC13_, FRAC112_
    use globals_cfg,   only: accSolv_, isOpenX_, isOpenYx_, &
                             isPol_, isPol4_, &
                             isSolar_, isThermal_
    use globals_type,  only: shd_type_grid

    use globals_obj,   only: rtm_obj_, xyz_Rtm_, planck_Rtm_, ext_Rtm_, ssa_Rtm_, leg_Rtm_, ipha_Rtm_, &
                             i8Grid_Cell_, i6Cell_Cell_, flag_Cell_, &
                            shdo_obj_, &
                             src_obj_, &
                             prp_obj_, &
                             sfc_obj_

    use shd_scat, only: shd_scat__computeLegendre, &
                        shd_scat__computeWigner

    use shd_shdo, only: shd_shdo__computeYlm

    use shd_util, only: shd_util__findCell, &
                        shd_util__findNextCell, &
                        shd_util__interpTriLin

    use shd_init, only: shd_init__makeBotRadV

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_src__computeSrcOnePnt
    public :: shd_src__computeSrcOnePnt_
    public :: shd_src__computeSrcOneCell
    public :: shd_src__computeSrcOneCell_
    public :: shd_src__computeSrcOneDir
    public :: shd_src__computeSrc

    public :: shd_src__integrateSrcOneRay
    public :: shd_src__integrateSrc

    private :: shd_src__computeSolarStokesVect
    private :: shd_src__rotatePolarPlane
    private :: shd_src__findBndryRad

contains

subroutine shd_src__computeSrcOneCell(&
    &icell, GRIDPTR,&
    &NSTOKES, ML, MM, NLM,&
    &NPTS, DELTAM, SRCTYPE, SOLARMU,&
    &DIRFLUX,&
    &isrcJ_Shdo, srcJ_Shdo, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
    &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
    &EXTINCT8, SRCEXT8)

    !   Computes the source function times extinction for gridpoints
    ! belonging to cell icell in the direction (MU,PHI). The results
    ! are returned in SRCEXT8 and EXTINCT8.
    ! The spherical harmonic source function series is input in srcJ_Shdo.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer icell, NSTOKES, NPTS, ML, MM, NLM
    integer GRIDPTR(8, *), isrcJ_Shdo(NPTS + 1)
    integer DONETHIS(8), OLDIPTS(8)
    logical DELTAM
    real SOLARMU
    real DIRFLUX(NPTS), srcJ_Shdo(NSTOKES, *)
    real YLMDIR(prp_obj_%NlegSet, NLM), YLMSUN(prp_obj_%NlegSet, NLM)
    real SINGSCAT(NSTOKES, prp_obj_%Npha)
    real OEXTINCT8(8), OSRCEXT8(NSTOKES, 8)
    real EXTINCT8(8), SRCEXT8(NSTOKES, 8)
    real(RD_) :: SUNDIRLEG(0:prp_obj_%Nleg)
    character(len=1) :: SRCTYPE

    integer IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
    real SECMU0, F, DA, A, A1, B1

    SECMU0 = 1.0d0/abs(SOLARMU)

    ! Loop over the grid points, computing the source function
    !   at the viewing angle from the spherical harmonics source fun
    do N = 1, 8
        IP = GRIDPTR(N, icell)
        I = DONETHIS(N)
        if (I .gt. 0 .and. IP .eq. OLDIPTS(N)) then
            EXTINCT8(N) = OEXTINCT8(I)
            SRCEXT8(:, N) = OSRCEXT8(:, I)
        else if (I .lt. 0) then
            EXTINCT8(N) = EXTINCT8(abs(I))
            SRCEXT8(:, N) = SRCEXT8(:, abs(I))
        else

            OLDIPTS(N) = IP
            IS = isrcJ_Shdo(IP)
            NS = isrcJ_Shdo(IP + 1) - IS
            ! Sum over the real generalized spherical harmonic series
            !   of the source function
            SRCEXT8(:, N) = 0.0
            do J = 1, NS
                SRCEXT8(1, N) = SRCEXT8(1, N) + srcJ_Shdo(1, IS + J)*YLMDIR(1, J)
            end do
            if (isPol_) then
                do J = 1, NS
                    SRCEXT8(2, N) = SRCEXT8(2, N) + srcJ_Shdo(2, IS + J)*YLMDIR(2, J)&
                        &+ srcJ_Shdo(3, IS + J)*YLMDIR(5, J)
                    SRCEXT8(3, N) = SRCEXT8(3, N) + srcJ_Shdo(2, IS + J)*YLMDIR(6, J)&
                        &+ srcJ_Shdo(3, IS + J)*YLMDIR(3, J)
                end do
            end if
            if (isPol4_) then
                do J = 1, NS
                    SRCEXT8(4, N) = SRCEXT8(4, N) + srcJ_Shdo(4, IS + J)*YLMDIR(4, J)
                end do
            end if

            ! Special case for solar source and Delta-M
            if (isSolar_ .and. DELTAM) then
                if (prp_obj_%Npha .gt. 0) then
                    IPH = ipha_Rtm_(IP)
                else
                    IPH = IP
                end if
                ! First subtract off the truncated single scattering
                DA = ssa_Rtm_(IP)*DIRFLUX(IP)*SECMU0
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    MS = -ME
                    A1 = DA*leg_Rtm_(1, L, IPH)
                    B1 = DA*leg_Rtm_(5, L, IPH)
                    if (J .le. NS) then
                        JT = J
                        do M = MS, ME
                            SRCEXT8(1, N) = SRCEXT8(1, N) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                        if (isPol_) then
                            J = JT
                            do M = MS, ME
                                SRCEXT8(2, N) = SRCEXT8(2, N) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                                SRCEXT8(3, N) = SRCEXT8(3, N) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                                J = J + 1
                            end do
                        end if
                    end if
                end do
                ! Then add in the single scattering contribution for the
                !   original unscaled phase function.
                if (prp_obj_%Npha .gt. 0) then
                    SRCEXT8(:, N) = SRCEXT8(:, N) + DA*SINGSCAT(:, IPH)
                else if (.not. isPol_) then
                    F = leg_Rtm_(1, ML + 1, IPH)
                    do L = 0, prp_obj_%Nleg
                        if (L .le. ML) then
                            A = DA*(leg_Rtm_(1, L, IPH) + F/(1 - F))
                        else
                            A = DA*leg_Rtm_(1, L, IPH)/(1 - F)
                        end if
                        SRCEXT8(1, N) = SRCEXT8(1, N) + A*SUNDIRLEG(L)
                    end do
                end if
            end if

            SRCEXT8(:, N) = SRCEXT8(:, N)*ext_Rtm_(IP)
            EXTINCT8(N) = ext_Rtm_(IP)
        end if
    end do

    return
end subroutine shd_src__computeSrcOneCell

subroutine shd_src__computeSrcOneCell_(&
    &icell, GRIDPTR,&
    &ML, MM, NLM,&
    &NPTS, DELTAM, SRCTYPE, SOLARMU,&
    &DIRFLUX,&
    &isrcJ_Shdo, srcJ_Shdo, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
    &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
    &EXTINCT8, SRCEXT8)

    !   Computes the source function times extinction for gridpoints
    ! belonging to cell icell in the direction (MU,PHI).  The results
    ! are returned in SRCEXT8 and EXTINCT8.
    ! The spherical harmonic source function series is input in srcJ_Shdo.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer icell, NPTS, ML, MM, NLM
    integer GRIDPTR(8, *), isrcJ_Shdo(NPTS + 1)
    integer DONETHIS(8), OLDIPTS(8)
    logical DELTAM
    real SOLARMU
    real DIRFLUX(NPTS), srcJ_Shdo(*)
    real YLMDIR(NLM), YLMSUN(NLM)
    real SINGSCAT(prp_obj_%Npha)
    real OEXTINCT8(8), OSRCEXT8(8)
    real EXTINCT8(8), SRCEXT8(8)
    real(RD_) :: SUNDIRLEG(0:prp_obj_%Nleg)
    character(len=1) :: SRCTYPE

    integer IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
    real SECMU0, F, DA, A, A1, B1

    SECMU0 = 1.0d0/abs(SOLARMU)

    ! Loop over the grid points, computing the source function
    !   at the viewing angle from the spherical harmonics source fun
    do N = 1, 8
        IP = GRIDPTR(N, icell)
        I = DONETHIS(N)
        if (I .gt. 0 .and. IP .eq. OLDIPTS(N)) then
            EXTINCT8(N) = OEXTINCT8(I)
            SRCEXT8(N) = OSRCEXT8(I)
        else if (I .lt. 0) then
            EXTINCT8(N) = EXTINCT8(abs(I))
            SRCEXT8(N) = SRCEXT8(abs(I))
        else

            OLDIPTS(N) = IP
            IS = isrcJ_Shdo(IP)
            NS = isrcJ_Shdo(IP + 1) - IS
            ! Sum over the real generalized spherical harmonic series
            !   of the source function
            SRCEXT8(N) = 0.0
            do J = 1, NS
                SRCEXT8(N) = SRCEXT8(N) + srcJ_Shdo(IS + J)*YLMDIR(J)
            end do

            ! Special case for solar source and Delta-M
            if (isSolar_ .and. DELTAM) then
                if (prp_obj_%Npha .gt. 0) then
                    IPH = ipha_Rtm_(IP)
                else
                    IPH = IP
                end if
                ! First subtract off the truncated single scattering
                DA = ssa_Rtm_(IP)*DIRFLUX(IP)*SECMU0
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    MS = -ME
                    A1 = DA*leg_Rtm_(1, L, IPH)
                    if (J .le. NS) then
                        JT = J
                        do M = MS, ME
                            SRCEXT8(N) = SRCEXT8(N) - A1*YLMDIR(J)*YLMSUN(J)
                            J = J + 1
                        end do
                    end if
                end do
                ! Then add in the single scattering contribution for the
                !   original unscaled phase function.
                if (prp_obj_%Npha .gt. 0) then
                    SRCEXT8(N) = SRCEXT8(N) + DA*SINGSCAT(IPH)
                else
                    F = leg_Rtm_(1, ML + 1, IPH)
                    do L = 0, prp_obj_%Nleg
                        if (L .le. ML) then
                            A = DA*(leg_Rtm_(1, L, IPH) + F/(1 - F))
                        else
                            A = DA*leg_Rtm_(1, L, IPH)/(1 - F)
                        end if
                        SRCEXT8(N) = SRCEXT8(N) + A*SUNDIRLEG(L)
                    end do
                end if
            end if

            SRCEXT8(N) = SRCEXT8(N)*ext_Rtm_(IP)
            EXTINCT8(N) = ext_Rtm_(IP)
        end if
    end do

    return
end subroutine shd_src__computeSrcOneCell_

subroutine shd_src__computeSrcOneDir(&
    &DELTAM, MU, PHI,&
    &DIRFLUX,&
    &src1_Shdo,&
    &grid_obj, isrcJ_Shdo, srcJ_Shdo)

    ! Computes the source function (src1_Shdo) in the direction (MU,PHI)
    ! for the whole domain (rtm_obj_%Nx,rtm_obj_%Nz).  The spherical harmonic source function
    ! series is input in srcJ_Shdo.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    type(shd_type_grid), intent(in) :: grid_obj

    integer(I__), intent(inout) :: isrcJ_Shdo(*)
    real(R__),    intent(inout) ::  srcJ_Shdo(shdo_obj_%Nstokes, *)

    logical DELTAM
    real MU, PHI
    real DIRFLUX(grid_obj%Npts)
    real src1_Shdo(shdo_obj_%Nstokes, grid_obj%Npts)
    integer I, J, JT, K, L, M, N, IPH, IS, NS, ME, MS
    real DA, F, A, A1, B1, SECMU0, COSSCAT
    real, allocatable :: YLMDIR(:, :), YLMSUN(:, :)
    real(RD_), allocatable :: SUNDIRLEG(:)
    real, allocatable :: UNSCLEGEN(:, :), SINGSCAT(:, :)
    logical TMS

    allocate(YLMDIR(prp_obj_%NlegSet, shdo_obj_%Nlm), YLMSUN(prp_obj_%NlegSet, shdo_obj_%Nlm))
    allocate(SUNDIRLEG(0:prp_obj_%Nleg))
    allocate(UNSCLEGEN(2, 0:prp_obj_%Nleg), SINGSCAT(shdo_obj_%Nstokes, prp_obj_%Npha))
    if (abs(MU) .gt. 1.0) stop 'Error [_computeOneSrc]: Bad mu.'

    TMS = DELTAM
    ! TMS = .FALSE.

    ! Precompute Ylm's for output direction and solar direction
    call shd_shdo__computeYlm(.false., MU, PHI, shdo_obj_%maxNl, shdo_obj_%maxNm, prp_obj_%NlegSet, YLMDIR)
    if (isSolar_) then
        call shd_shdo__computeYlm(.true., src_obj_%avgMu, src_obj_%avgAz, shdo_obj_%maxNl, shdo_obj_%maxNm, prp_obj_%NlegSet,&
            &YLMSUN)
        SECMU0 = 1.0/abs(src_obj_%avgMu)
        COSSCAT = src_obj_%avgMu*MU + sqrt((1.0d0 - src_obj_%avgMu**2)*(1.0d0 - MU**2))&
            &*cos(src_obj_%avgAz - PHI)
        if (prp_obj_%Npha .eq. 0 .and. TMS) then
            ! If non-tabulate phase functions (unpolarized) then
            !   compute the Legendre polynomials for the scattering angle
            !   for the untruncated solar single scattering computation.
            call shd_scat__computeLegendre(dble(COSSCAT), prp_obj_%Nleg, SUNDIRLEG)
            do L = 0, prp_obj_%Nleg
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)*O4PI_
            end do
        end if
    end if

    ! If solar source and tabulated phase functions then precompute
    !   the single scattered solar Stokes vector to outgoing direction
    if ((src_obj_%type .ne. 'T') .and. TMS .and. (prp_obj_%Npha .gt. 0)) then
        do IPH = 1, prp_obj_%Npha
            ! Unscaled the needed Legendre coefficients, but divide by 1-f
            !   because extinction is still scaled (for TMS method)
            F = leg_Rtm_(1, shdo_obj_%maxNl + 1, IPH)
            do L = 0, prp_obj_%Nleg
                if (L .le. shdo_obj_%maxNl) then
                    UNSCLEGEN(1, L) = leg_Rtm_(1, L, IPH) + F/(1 - F)
                else
                    UNSCLEGEN(1, L) = leg_Rtm_(1, L, IPH)/(1 - F)
                end if
                if (prp_obj_%NlegSet .gt. 1) then
                    if (L .le. shdo_obj_%maxNl) then
                        UNSCLEGEN(2, L) = leg_Rtm_(5, L, IPH)
                    else
                        UNSCLEGEN(2, L) = leg_Rtm_(5, L, IPH)/(1 - F)
                    end if
                end if
            end do
            ! Sum the Wigner series (for A1 and B1), rotate the plane of
            !   polarization to the outgoing direction, and return the Stoke
            call shd_src__computeSolarStokesVect(UNSCLEGEN, COSSCAT,&
                &src_obj_%avgMu, MU, src_obj_%avgAz - PHI, shdo_obj_%Nstokes,&
                &SINGSCAT(1, IPH))
            if (SINGSCAT(1, IPH) .le. 0.0) then
                write (6, *) 'Error [_computeOneSrc]: negative source function',&
                    &' for tabulated phase function:', IPH, MU, PHI, SINGSCAT(1, IPH)
                stop
            end if
        end do
    end if

    ! Loop over all the grid points, computing the source function
    !   at the viewing angle from the spherical harmonic source func
    do I = 1, grid_obj%Npts
        IS = isrcJ_Shdo(I)
        NS = isrcJ_Shdo(I + 1) - IS
        ! Sum over the real generalized spherical harmonic series
        !   of the source function
        src1_Shdo(:, I) = 0.0
        do J = 1, NS
            src1_Shdo(1, I) = src1_Shdo(1, I) + srcJ_Shdo(1, IS + J)*YLMDIR(1, J)
        end do
        if (isPol_) then
            do J = 1, NS
                src1_Shdo(2, I) = src1_Shdo(2, I) + srcJ_Shdo(2, IS + J)*YLMDIR(2, J)&
                    &+ srcJ_Shdo(3, IS + J)*YLMDIR(5, J)
                src1_Shdo(3, I) = src1_Shdo(3, I) + srcJ_Shdo(2, IS + J)*YLMDIR(6, J)&
                    &+ srcJ_Shdo(3, IS + J)*YLMDIR(3, J)
            end do
        end if
        if (isPol4_) then
            do J = 1, NS
                src1_Shdo(4, I) = src1_Shdo(4, I) + srcJ_Shdo(4, IS + J)*YLMDIR(4, J)
            end do
        end if

        ! Special case for solar source and Delta-M
        if ((src_obj_%type .ne. 'T') .and. TMS) then
            if (prp_obj_%Npha .gt. 0) then
                IPH = ipha_Rtm_(I)
            else
                IPH = I
            end if
            ! First subtract off the truncated single scattering
            DA = ssa_Rtm_(I)*DIRFLUX(I)*SECMU0
            J = 1
            do L = 0, shdo_obj_%maxNl
                ME = min(L, shdo_obj_%maxNm)
                MS = -ME
                A1 = DA*leg_Rtm_(1, L, IPH)
                B1 = DA*leg_Rtm_(5, L, IPH)
                if (J .le. NS) then
                    JT = J
                    do M = MS, ME
                        src1_Shdo(1, I) = src1_Shdo(1, I) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                        J = J + 1
                    end do
                    if (isPol_) then
                        J = JT
                        do M = MS, ME
                            src1_Shdo(2, I) = src1_Shdo(2, I) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                            src1_Shdo(3, I) = src1_Shdo(3, I) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                    end if
                end if
            end do

            ! Then add in the single scattering contribution for the
            !   original unscaled phase matrix.
            if (prp_obj_%Npha .gt. 0) then
                src1_Shdo(:, I) = src1_Shdo(:, I) + DA*SINGSCAT(:, IPH)
            else if (.not. isPol_) then
                ! If we have a phase function for each grid point (which
                !   should be the unpolarized case) then unscale the Legendre
                !   coefficients for L<=shdo_obj_%maxNl and sum over the Legendre series
                !   in cosine of the scattering angle.
                F = leg_Rtm_(1, shdo_obj_%maxNl + 1, IPH)
                do L = 0, prp_obj_%Nleg
                    if (L .le. shdo_obj_%maxNl) then
                        A = DA*(leg_Rtm_(1, L, IPH) + F/(1 - F))
                    else
                        A = DA*leg_Rtm_(1, L, IPH)/(1 - F)
                    end if
                    src1_Shdo(1, I) = src1_Shdo(1, I) + A*SUNDIRLEG(L)
                end do
            end if
        end if
    end do

    deallocate(YLMDIR, YLMSUN, SUNDIRLEG, UNSCLEGEN, SINGSCAT)

    return

end subroutine shd_src__computeSrcOneDir

subroutine shd_src__computeSrc(&
    FIXSH, &
    YLMSUN, &
    DIRFLUX, SHACC, &
    isrcD_Shdo, srcD_Shdo, &
    FIRST, ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM, &
    grid_obj, &
    iradI_Shdo, radI_Shdo, &
    isrcJ_Shdo, srcJ_Shdo)

    ! Computes the source function (srcJ_Shdo) in spherical harmonic space
    !   for all the grid points.  The thermal source and/or solar
    !   pseudo-source (in planck_Rtm_ or DIRFLUX) is added to the scattering source
    !   (computed from leg_Rtm_ and the spherical harmonic expansion in radI_Shdo).
    !   Also computes the adaptive spherical harmonic truncation for
    !   the source function, unless FIXSH is true.  The SH truncation is set
    !   so that the higher terms of the source function are smaller in
    !   absolute value than SHACC. The truncation is to the nearest l,
    !   so that a triangular (for M=L) truncation is preserved.  If there
    !   is little or no scattering then the SH truncation can be set to
    !   zero terms.  The source function SH truncation is output in terms
    !   of the pointer isrcJ_Shdo (note that iradI_Shdo refers to the radI_Shdo array
    !   and isrcD_Shdo to the srcD_Shdo array).
    ! To save memory this routine also computes some things for the
    !   series acceleration: the dot products of the difference in source
    !   function between successive iterations and the new difference source
    !   function vector (srcD_Shdo).  To avoid using another array and
    !   because of the adaptive truncation pointer stuff, this requires
    !   two extra calculations of the source function.
    ! Computes the sum of squares of the new source function for calculation
    !   of the solution criterion (SOLCRIT, the RMS difference in succesive
    !   source function fields normalized by the RMS of the field).

    type(shd_type_grid), intent(in) :: grid_obj

    integer(I__), intent(inout) :: iradI_Shdo(*)
    real(R__),    intent(inout) ::  radI_Shdo(shdo_obj_%Nstokes, *)

    integer(I__), intent(inout) :: isrcJ_Shdo(*)
    real(R__),    intent(inout) ::  srcJ_Shdo(shdo_obj_%Nstokes, *)


    integer isrcD_Shdo(grid_obj%Npts + 1)
    logical FIXSH, FIRST, ACCELFLAG
    real SHACC, YLMSUN(prp_obj_%NlegSet, shdo_obj_%Nlm)
    real DELJDOT, DELJOLD, DELJNEW, JNORM
    real DIRFLUX(grid_obj%Npts)
    real srcD_Shdo(shdo_obj_%Nstokes, *)
    integer IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
    integer, allocatable :: LOFJ(:)
    real, allocatable :: srcT_Shdo(:, :)
    real SRCMIN, C, SECMU0, D

    allocate(LOFJ(shdo_obj_%Nlm), srcT_Shdo(shdo_obj_%Nstokes, shdo_obj_%Nlm))

    SRCMIN = SHACC
    ! Make source function from scattering integral and real sources:
    !   Solar source: DIRFLUX has F0*exp(-tau)
    !   Thermal source: planck_Rtm_ has (1-omega)*B(T)
    SECMU0 = 1.0/abs(src_obj_%avgMu)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, shdo_obj_%maxNl
        ME = min(L, shdo_obj_%maxNm)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    ! Compute the dot products of the successive source differences
    if (.not. FIRST) then
        DELJDOT = 0.0
        DELJOLD = 0.0
        DELJNEW = 0.0
        JNORM = 0.0
        do I = 1, grid_obj%Npts
            ! Compute the temporary source function for this point
            if (prp_obj_%Npha .gt. 0) then
                IPH = ipha_Rtm_(I)
            else
                IPH = I
            end if
            IR = iradI_Shdo(I)
            NR = iradI_Shdo(I + 1) - IR
            if (NR .gt. shdo_obj_%Nlm) then
                write (6, *) 'Error [_computeSrc]: NR>NLM 1', I, IR, iradI_Shdo(I + 1)
                stop
            end if
            if (isPol_) then
                call shd_src__computeSrcOnePnt(shdo_obj_%Nstokes, shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                    &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                    !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
            else
                call shd_src__computeSrcOnePnt_(shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                    &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                    !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
            end if

            ! Compute the dot product of the new and old source difference,
            !   and the solution criterion (don't use new points for Jnorm).
            if (ACCELFLAG) then
                IS = isrcJ_Shdo(I)
                ISO = isrcD_Shdo(I)
                NS = min(isrcJ_Shdo(I + 1) - IS, isrcD_Shdo(I + 1) - ISO)
                do K = 1, shdo_obj_%Nstokes
                    do J = 1, NS
                        DELJDOT = DELJDOT +&
                            &(srcT_Shdo(K, J) - srcJ_Shdo(K, IS + J))*srcD_Shdo(K, ISO + J)
                        DELJOLD = DELJOLD + srcD_Shdo(K, ISO + J)**2
                        DELJNEW = DELJNEW + (srcT_Shdo(K, J) - srcJ_Shdo(K, IS + J))**2
                        JNORM = JNORM + srcJ_Shdo(K, IS + J)**2
                    end do
                end do
            else
                IS = isrcJ_Shdo(I)
                NS = isrcJ_Shdo(I + 1) - IS
                do K = 1, shdo_obj_%Nstokes
                    do J = 1, NS
                        DELJNEW = DELJNEW + (srcT_Shdo(K, J) - srcJ_Shdo(K, IS + J))**2
                        JNORM = JNORM + srcJ_Shdo(K, IS + J)**2
                    end do
                end do
            end if
        end do
    end if

    ! Compute and store the new source function difference vector
    if (.not. FIRST .and. ACCELFLAG) then
        do I = 1, grid_obj%Npts
            ! Compute the temporary source function for this point
            if (prp_obj_%Npha .gt. 0) then
                IPH = ipha_Rtm_(I)
            else
                IPH = I
            end if
            IR = iradI_Shdo(I)
            NR = iradI_Shdo(I + 1) - IR
            if (isPol_) then
                call shd_src__computeSrcOnePnt(shdo_obj_%Nstokes, shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                    &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                    !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
            else
                call shd_src__computeSrcOnePnt_(shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                    &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                    !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
            end if

            ! Make the source difference vector
            IS = isrcJ_Shdo(I)
            NS = isrcJ_Shdo(I + 1) - IS
            isrcD_Shdo(I) = IS
            do J = 1, NS
                srcD_Shdo(:, IS + J) = srcT_Shdo(:, J) - srcJ_Shdo(:, IS + J)
            end do
        end do
        isrcD_Shdo(grid_obj%Npts + 1) = isrcJ_Shdo(grid_obj%Npts + 1)
    end if

    ! Compute the source function in the temporary array.
    IS = 0
    do I = 1, grid_obj%Npts
        if (prp_obj_%Npha .gt. 0) then
            IPH = ipha_Rtm_(I)
        else
            IPH = I
        end if
        IR = iradI_Shdo(I)
        NR = iradI_Shdo(I + 1) - IR
        if (NR .gt. shdo_obj_%Nlm) then
            write (6, *) 'Error [_computeSrc]: NR>NLM 3', I, IR, iradI_Shdo(I + 1)
            stop
        end if
        if (isPol_) then
            call shd_src__computeSrcOnePnt(shdo_obj_%Nstokes, shdo_obj_%Nlm, LOFJ,&
                &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
        else
            call shd_src__computeSrcOnePnt_(shdo_obj_%Nlm, LOFJ,&
                &src_obj_%type, DIRFLUX(I)*SECMU0, YLMSUN,&
                &planck_Rtm_(I), ssa_Rtm_(I), IPH, leg_Rtm_,&
                &NR, radI_Shdo(1, IR+1), srcT_Shdo)
                !&NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), srcT_Shdo)
        end if
        if (FIXSH) then
            NS = isrcJ_Shdo(I + 1) - IS
        else
            ! Then do the adaptive spherical harmonics:
            !   Find how many source function terms to use, and how many
            !   terms wanted for the radiance.
            if (src_obj_%type .eq. 'S') then
                JS = 0
            else
                JS = 1
            end if
            do J = 1, shdo_obj_%Nlm
                do K = 1, shdo_obj_%Nstokes
                    if (abs(srcT_Shdo(K, J)) .gt. SRCMIN) JS = J
                end do
            end do
            ! Use all the m's for the last l for the source function.
            if (JS .eq. 0) then
                NS = 0
            else
                LS = LOFJ(JS)
                if (LS .le. shdo_obj_%maxNm) then
                    NS = (LS*(LS + 1)) + LS + 1
                else
                    NS = (2*shdo_obj_%maxNm + 1)*LS - (shdo_obj_%maxNm*(1 + (shdo_obj_%maxNm - 1))) + shdo_obj_%maxNm + 1
                end if
            end if
            isrcJ_Shdo(I) = IS
        end if
        if (IS + NS .gt. shdo_obj_%maxNptsSH) then
            write (6, *) 'Error [_computeSrc]: shdo_obj_%maxNptsSH exceeded.', shdo_obj_%maxNptsSH
            write (6, *) '  Out of memory for more spherical harmonic terms.'
            stop
        end if
        ! Transfer the new source function
        do J = 1, NS
            srcJ_Shdo(:, IS + J) = srcT_Shdo(:, J)
        end do
        IS = IS + NS
    end do
    isrcJ_Shdo(grid_obj%Npts + 1) = IS

    deallocate(LOFJ, srcT_Shdo)

    return

end subroutine shd_src__computeSrc

subroutine shd_src__integrateSrcOneRay(&
    &NSTOKES,&
    &NSTPHASE, NSCATANGLE, PHASETAB,&
    &ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &DELTAM, SRCTYPE, WAVELEN, SOLARMU, SOLARAZ,&
    &DIRFLUX, isrcJ_Shdo, srcJ_Shdo, YLMSUN,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &MU2, PHI2, X0, Y0, Z0,&
    &XE, YE, ZE, SIDE, TRANSMIT, radI_Shdo, VALIDRAD,&
    &grid_obj)

    !   Integrates the source function through the extinction field
    ! (ext_Rtm_) backward from the outgoing direction (MU2,PHI2) to find the
    ! radiance (radI_Shdo) at the point X0,Y0,Z0.
    ! The transmission and radiance of the ray so far (TRANSMIT, radI_Shdo)
    ! are input and returned after the integration along with the exitting
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
    ! 5=-Z,6=+Z).

    type(shd_type_grid), intent(in) :: grid_obj

    integer NSTOKES, NSTPHASE, NSCATANGLE
    integer SIDE
    integer ML, MM, NLM
    integer MAXNBC, NTOPPTS, NBOTPTS
    integer NMU, NPHI0MAX, NPHI0(NMU)
    integer isrcJ_Shdo(grid_obj%Npts + 1)
    integer BCPTR(MAXNBC, 2)
    logical DELTAM, VALIDRAD
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ
    real BCRAD(*)
    real DIRFLUX(grid_obj%Npts), srcJ_Shdo(NSTOKES, *)
    real YLMSUN(prp_obj_%NlegSet, NLM)
    real PHASETAB(NSTPHASE, prp_obj_%Npha, NSCATANGLE)
    real(RD_) :: MU2, PHI2, X0, Y0, Z0, XE, YE, ZE
    real(RD_) :: TRANSMIT, radI_Shdo(NSTOKES)
    character(len=1) :: SRCTYPE

    integer BITX, BITY, BITZ, IOCT, icell, INEXTCELL, IFACE
    integer IOPP, NTAU, IT, I, IPT1, IPT2, J, L
    logical DONE, IPINX, IPINY, OPENBCFACE
    integer JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
    integer OPPFACE(6), OLDIPTS(8), DONETHIS(8)
    integer DONEFACE(8, 7), ONEX(8), ONEY(8)
    real OEXTINCT8(8), OSRCEXT8(NSTOKES, 8)
    real EXTINCT8(8), SRCEXT8(NSTOKES, 8)
    real EXT0, EXT1, EXTN
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RADBND(NSTOKES)
    real XM, YM, F
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XN, YN, ZN, XI, YI, ZI
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: TAUTOL, TAUGRID, S, DELS, TRANSCUT
    real(RD_) :: EXT, TAU, TRANSCELL, ABSCELL
    real(RD_) :: SRC(NSTOKES)
    real(RD_) :: U, V, W, DELX, DELY, DELZ, INVDELX, INVDELY, INVDELZ
    real(RD_) :: COSSCAT
    real, allocatable :: YLMDIR(:, :), SINGSCAT(:, :)
    real(RD_), allocatable :: SUNDIRLEG(:)

    data OPPFACE/2, 1, 4, 3, 6, 5/
    data ONEY/0, 0, -1, -2, 0, 0, -5, -6/, ONEX/0, -1, 0, -3, 0, -5, 0, -7/
    data DONEFACE/0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 5, 0, 7, 2, 0, 4, 0, 6, 0, 8, 0,&
        &0, 0, 1, 2, 0, 0, 5, 6, 3, 4, 0, 0, 5, 6, 0, 0,&
        &0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0/

    ! TRANSCUT is the transmission to stop the integration at
    !TRANSCUT = 5.0e-5
    TRANSCUT = accSolv_
    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.2

    EPS = 1.0e-5*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))
    MAXCELLSCROSS = 50*max(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz)

    ! Calculate the generalized spherical harmonics for this direction
    allocate(YLMDIR(prp_obj_%NlegSet, NLM))
    call shd_shdo__computeYlm(.false., sngl(MU2), sngl(PHI2), ML, MM, prp_obj_%NlegSet, YLMDIR)

    if ((SRCTYPE .ne. 'T') .and. (DELTAM)) then
        ! Get the solar single scattering Stokes vector for the outgoin
        !   direction my interpolating in scattering angle in the PHASETA
        !   table and then rotating the Q/U polarization to the outgoing
        COSSCAT = SOLARMU*MU2 + dsqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MU2**2))&
            &*dcos(SOLARAZ - PHI2)
        if (prp_obj_%Npha .gt. 0) then
            allocate(SINGSCAT(NSTOKES, prp_obj_%Npha))
            F = (NSCATANGLE - 1)*(dacos(COSSCAT)/DPI_) + 1
            J = min(NSCATANGLE - 1, int(F))
            F = F - J
            do I = 1, prp_obj_%Npha
                SINGSCAT(1:NSTPHASE, I)&
                    &= (1 - F)*PHASETAB(:, I, J) + F*PHASETAB(:, I, J + 1)
                if (isPol_) then
                    call shd_src__rotatePolarPlane(NSTOKES, COSSCAT, SOLARMU,&
                        &sngl(MU2), SOLARAZ - sngl(PHI2), SINGSCAT(:, I))
                end if
            end do
        else
            allocate(SUNDIRLEG(0:prp_obj_%Nleg))
            call shd_scat__computeLegendre(COSSCAT, prp_obj_%Nleg, SUNDIRLEG)
            do L = 0, prp_obj_%Nleg
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)/(4*DPI_)
            end do
        end if
    end if

    ! Make the ray direction (opposite to the outgoing direction)
    CX = dsqrt(1.0d0 - MU2**2)*dcos(PHI2 - DPI_)
    CY = dsqrt(1.0d0 - MU2**2)*dsin(PHI2 - DPI_)
    CZ = -MU2
    if (abs(CX) .gt. accSolv_) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. accSolv_) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    if (abs(CZ) .gt. accSolv_) then
        CZINV = 1.0d0/CZ
    else
        CZ = 0.0
        CZINV = 1.0e6
    end if

    ! Setup for the ray path direction
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. 0.0) then
        BITZ = 1
    else
        BITZ = 0
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5d0*(rtm_obj_%x(1) + rtm_obj_%x(rtm_obj_%Nx))
    YM = 0.5d0*(rtm_obj_%y(1) + rtm_obj_%y(rtm_obj_%Ny))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    ! ZE = MAX(MIN(Z0,DBLE(rtm_obj_%z(rtm_obj_%Nz))),DBLE(rtm_obj_%z(1)))
    call shd_util__findCell(XE, YE, ZE, icell, grid_obj)
    ! print '(A,5(1X,F8.5),1X,I6)',&
    !     &'_integrateSrcOneRay:',X0,Y0,Z0,MU2,PHI2,icell
    IFACE = 0
    NGRID = 0

    ! Loop until reach a Z boundary or transmission is very small
    VALIDRAD = .false.
    do while ((.not. VALIDRAD) .and. (icell .gt. 0))
        ! Make sure current cell is valid
        if (icell .le. 0) then
            write (6, *) 'Error [_integrateSrcOneRay]: icell=', icell,&
                &MU2, PHI2, XE, YE, ZE
            stop
        end if
        NGRID = NGRID + 1

        ! Decide which of the eight grid points we need the source fun
        do I = 1, 8
            DONETHIS(I) = DONEFACE(I, IFACE + 1)
            if (rtm_obj_%Nx .eq. 1 .and. ONEX(I) .lt. 0) DONETHIS(I) = ONEX(I)
            if (rtm_obj_%Ny .eq. 1 .and. ONEY(I) .lt. 0) DONETHIS(I) = ONEY(I)
            OEXTINCT8(I) = EXTINCT8(I)
            OSRCEXT8(:, I) = SRCEXT8(:, I)
        end do
        ! Compute the source function times extinction in direction (MU2
        if (isPol_) then
            call shd_src__computeSrcOneCell(icell, i8Grid_Cell_,&
                &NSTOKES, ML, MM, NLM,&
                &grid_obj%Npts, DELTAM, SRCTYPE, SOLARMU,&
                &DIRFLUX,&
                &isrcJ_Shdo, srcJ_Shdo, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        else
            call shd_src__computeSrcOneCell_(icell, i8Grid_Cell_,&
                &ML, MM, NLM,&
                &grid_obj%Npts, DELTAM, SRCTYPE, SOLARMU,&
                &DIRFLUX,&
                &isrcJ_Shdo, srcJ_Shdo, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        end if

        ! Interpolate the source and extinction to the current point
        IPT1 = i8Grid_Cell_(1, icell) ! (-x, -y, -z)
        IPT2 = i8Grid_Cell_(8, icell) ! (+x, +y, +z)
        DELX = xyz_Rtm_(1, IPT2) - xyz_Rtm_(1, IPT1)
        if (DELX .le. 0.0) then
            INVDELX = 1.0
        else
            INVDELX = 1.0d0/DELX
        end if
        DELY = xyz_Rtm_(2, IPT2) - xyz_Rtm_(2, IPT1)
        if (DELY .le. 0.0) then
            INVDELY = 1.0
        else
            INVDELY = 1.0d0/DELY
        end if
        DELZ = xyz_Rtm_(3, IPT2) - xyz_Rtm_(3, IPT1)
        INVDELZ = 1.0d0/DELZ
        U = (XE - xyz_Rtm_(1, IPT1))*INVDELX
        V = (YE - xyz_Rtm_(2, IPT1))*INVDELY
        W = (ZE - xyz_Rtm_(3, IPT1))*INVDELZ
        SRCEXT1(:) = (1 - W)*((1 - V)*((1 - U)*SRCEXT8(:, 1) + U*SRCEXT8(:, 2))&
            &+ V*((1 - U)*SRCEXT8(:, 3) + U*SRCEXT8(:, 4)))&
            &+ W*((1 - V)*((1 - U)*SRCEXT8(:, 5) + U*SRCEXT8(:, 6))&
            &+ V*((1 - U)*SRCEXT8(:, 7) + U*SRCEXT8(:, 8)))
        SRCEXT1(1) = max(0.0, SRCEXT1(1))
        EXT1 = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
            &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
            &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
            &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))

        ! This cell is independent pixel if IP mode or open boundary
        !   conditions and ray is leaving domain (i.e. not entering)
        IPINX = btest(int(flag_Cell_(icell)), 0) .and.&
            &.not. (isOpenX_ .and.&
            &((CX .gt. 0 .and. XE .lt. XM) .or. (CX .lt. 0 .and. XE .gt. XM)))
        IPINY = btest(int(flag_Cell_(icell)), 1) .and.&
            &.not. (isOpenYx_ .and.&
            &((CY .gt. 0 .and. YE .lt. YM) .or. (CY .lt. 0 .and. YE .gt. YM)))

        ! Find boundaries of the current cell
        ! Find the three possible intersection planes (X,Y,Z)
        !   from the coordinates of the opposite corner grid point
        IOPP = i8Grid_Cell_(9 - IOCT, icell)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX) then
            SOX = 1.0e20
        else
            SOX = (xyz_Rtm_(1, IOPP) - XE)*CXINV
        end if
        if (IPINY) then
            SOY = 1.0e20
        else
            SOY = (xyz_Rtm_(2, IOPP) - YE)*CYINV
        end if
        SOZ = (xyz_Rtm_(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [_integrateSrcOneRay]: SO<0 ',&
                &MU2, PHI2, XE, YE, ZE, SO, icell
            stop
        end if
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        U = (XN - xyz_Rtm_(1, IPT1))*INVDELX
        V = (YN - xyz_Rtm_(2, IPT1))*INVDELY
        W = (ZN - xyz_Rtm_(3, IPT1))*INVDELZ
        EXTN = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
            &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
            &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
            &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))
        TAUGRID = SO*0.5*(EXT1 + EXTN)
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU

        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            XI = XE + S*CX
            YI = YE + S*CY
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            U = (XI - xyz_Rtm_(1, IPT1))*INVDELX
            V = (YI - xyz_Rtm_(2, IPT1))*INVDELY
            W = (ZI - xyz_Rtm_(3, IPT1))*INVDELZ
            EXT0 = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
                &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
                &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
                &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))
            SRCEXT0(:) = (1 - W)*((1 - V)*((1 - U)*SRCEXT8(:, 1) + U*SRCEXT8(:, 2))&
                &+ V*((1 - U)*SRCEXT8(:, 3) + U*SRCEXT8(:, 4)))&
                &+ W*((1 - V)*((1 - U)*SRCEXT8(:, 5) + U*SRCEXT8(:, 6))&
                &+ V*((1 - U)*SRCEXT8(:, 7) + U*SRCEXT8(:, 8)))
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0 + EXT1)
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0d0 - 0.5d0*TAU*(1.0d0 - FRAC13_*TAU))
                TRANSCELL = 1.0d0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5d0*(SRCEXT0(:) + SRCEXT1(:))&
                    &+ FRAC112_*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*DELS&
                    &*(1.0d0 - 0.05d0*(EXT1 - EXT0)*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            radI_Shdo(:) = radI_Shdo(:) + TRANSMIT*SRC(:)*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            EXT1 = EXT0
            SRCEXT1(:) = SRCEXT0(:)
        ! End of sub grid cell loop
        end do

        ! Get the intersection face number (i.e. neighptr index)
        if ((SOX .le. SOZ) .and. (SOX .le. SOY)) then
            IFACE = 2 - BITX
            JFACE = 1
            OPENBCFACE = btest(int(flag_Cell_(icell)), 0) .and. isOpenX_
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
            OPENBCFACE = btest(int(flag_Cell_(icell)), 1) .and. isOpenYx_
        else
            IFACE = 6 - BITZ
            JFACE = 3
            OPENBCFACE = .false.
        end if
        ! Get the next cell to go to
        INEXTCELL = i6Cell_Cell_(IFACE, icell)
        if (INEXTCELL .lt. 0) then
            call shd_util__findNextCell(XN, YN, ZN, IFACE, JFACE, icell,&
                &INEXTCELL, grid_obj)
        end if
        ! If going to same or larger face then use previous face
        if ((i6Cell_Cell_(IFACE, icell) .ge. 0) .and. (.not. OPENBCFACE)) then
            KFACE = IFACE
            IC = icell
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
            IFACE = 0
        end if
        ! Get the location coordinate
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XN = xyz_Rtm_(1, i8Grid_Cell_(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YN = xyz_Rtm_(2, i8Grid_Cell_(IOCT, INEXTCELL))
            else
                ZN = xyz_Rtm_(3, i8Grid_Cell_(IOCT, INEXTCELL))
            end if
        end if

        ! If the transmission is greater than zero and not at a
        !   boundary then prepare for next cell
        if ((TRANSMIT .lt. TRANSCUT) .or. (NGRID .gt. MAXCELLSCROSS)) then
            VALIDRAD = .true.
        else if ((INEXTCELL .eq. 0) .and. (IFACE .ge. 5)) then
            VALIDRAD = .true.
            call shd_src__findBndryRad(NSTOKES, XN, YN,&
                &sngl(MU2), sngl(PHI2),&
                &IC, KFACE, i8Grid_Cell_, xyz_Rtm_,&
                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                &RADBND)
            radI_Shdo(:) = radI_Shdo(:) + TRANSMIT*RADBND(:)
        else
            icell = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do

    SIDE = IFACE

    if (isSolar_ .and. DELTAM) then
        if (prp_obj_%Npha .gt. 0) then
            deallocate(SINGSCAT)
        else
            deallocate(SUNDIRLEG)
        end if
    end if

    deallocate(YLMDIR)

    return
end subroutine shd_src__integrateSrcOneRay

subroutine shd_src__integrateSrc(NSTOKES,&
    &MU, PHI, GRIDRAD, srcJ_Shdo,&
    &X0, Y0, Z0, XE, YE, ZE, SIDE,&
    &TRANSMIT, radI_Shdo, VALIDRAD,&
    &grid_obj)

    ! Integrates the source function (srcJ_Shdo) through the extinction
    ! field (ext_Rtm_) backward in the direction opposite to (MU,PHI)
    ! to find the outgoing radiance (radI_Shdo) at the point X0,Y0,Z0.
    ! The radiances at the opposite boundary are input in GRIDRAD.
    ! The transmission and radiance of the ray so far (TRANSMIT, radI_Shdo)
    ! are input and returned after the integration along with the exitting
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
    ! 5=-Z,6=+Z).

    type(shd_type_grid), intent(in) :: grid_obj

    integer NSTOKES, SIDE
    logical VALIDRAD
    real MU, PHI
    real(RD_) :: X0, Y0, Z0, XE, YE, ZE, TRANSMIT, radI_Shdo(NSTOKES)
    real srcJ_Shdo(NSTOKES, grid_obj%Npts), GRIDRAD(NSTOKES, grid_obj%Npts)
    integer BITX, BITY, BITZ, IOCT, icell, INEXTCELL, IFACE
    integer I1, I2, I3, I4, IOPP, NTAU, IT, K
    logical IPINX, IPINY, OPENBCFACE
    integer GRIDFACE(4, 6), OPPFACE(6), JFACE, KFACE, IC
    real EXT0(1), EXT1(1), EXTN(1)
    real XM, YM
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RAD0(NSTOKES)
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XN, YN, ZN, XI, YI, ZI
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: TAUTOL, TAUGRID, S, DELS
    real(RD_) :: EXT, TAU, TRANSCELL, ABSCELL
    real(RD_) :: SRC(NSTOKES), RAD(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/

    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.1
    EPS = 1.0e-4*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    CX = sqrt(1.0 - MU**2)*cos(PHI + DPI_)
    CY = sqrt(1.0 - MU**2)*sin(PHI + DPI_)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup for the ray path direction
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [_integrateSrc]: Bad MU.'
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5*(rtm_obj_%x(1) + rtm_obj_%x(rtm_obj_%Nx))
    YM = 0.5*(rtm_obj_%y(1) + rtm_obj_%y(rtm_obj_%Ny))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    call shd_util__findCell(XE, YE, ZE, icell, grid_obj)
    ! print '(A,5(1X,F8.5),1X,I6)', &
    !     &'_integrateSrc:', X0, Y0, Z0, MU, PHI, icell
    call shd_util__interpTriLin(icell,&
        &XE, YE, ZE, 1, ext_Rtm_, 1, 1, EXT1, grid_obj)
    call shd_util__interpTriLin(icell,&
        &XE, YE, ZE, NSTOKES, srcJ_Shdo, 1, NSTOKES, SRCEXT1, grid_obj)
    SRCEXT1(1) = max(0.0, SRCEXT1(1))

    ! Loop until finding a face with known radiances
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. icell .gt. 0)
        ! This cell is independent pixel if IP mode or open boundary
        !   conditions and ray is leaving domain (i.e. not entering)
        IPINX = btest(int(flag_Cell_(icell)), 0) .and.&
            &.not. (isOpenX_ .and.&
            &((CX .gt. 0 .and. XE .lt. XM) .or. (CX .lt. 0 .and. XE .gt. XM)))
        IPINY = btest(int(flag_Cell_(icell)), 1) .and.&
            &.not. (isOpenYx_ .and.&
            &((CY .gt. 0 .and. YE .lt. YM) .or. (CY .lt. 0 .and. YE .gt. YM)))

        ! Find boundaries of the current cell
        ! Find the three possible intersection planes (X,Y,Z)
        ! from the coordinates of the opposite corner grid point
        IOPP = i8Grid_Cell_(9 - IOCT, icell)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX .or. CX .eq. 0.0) then
            SOX = 1.0e20
        else
            SOX = (xyz_Rtm_(1, IOPP) - XE)*CXINV
        end if
        if (IPINY .or. CY .eq. 0.0) then
            SOY = 1.0e20
        else
            SOY = (xyz_Rtm_(2, IOPP) - YE)*CYINV
        end if
        SOZ = (xyz_Rtm_(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [_integrateSrc]: SO<0 ',&
                &MU, PHI, X0, Y0, Z0, XE, YE, ZE, SO, SOX, SOY, SOZ, rtm_obj_%x(1), rtm_obj_%y(1),&
                &xyz_Rtm_(1, IOPP), xyz_Rtm_(2, IOPP), xyz_Rtm_(3, IOPP)
            stop
        end if
        if (IPINX) then
            XN = XE
        else
            XN = XE + SO*CX
        end if
        if (IPINY) then
            YN = YE
        else
            YN = YE + SO*CY
        end if
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        call shd_util__interpTriLin(icell,&
            &XN, YN, ZN, 1, ext_Rtm_, 1, 1, EXTN, grid_obj)
        TAUGRID = SO*0.5*(EXT1(1) + EXTN(1))
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU
        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            if (IPINX) then
                XI = XE
            else
                XI = XE + S*CX
            end if
            if (IPINY) then
                YI = YE
            else
                YI = YE + S*CY
            end if
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            if (IT .eq. NTAU) then
                EXT0(1) = EXTN(1)
            else
                call shd_util__interpTriLin(icell,&
                    &XI, YI, ZI, 1, ext_Rtm_, 1, 1, EXT0, grid_obj)
            end if
            call shd_util__interpTriLin(icell,&
                &XI, YI, ZI, NSTOKES, srcJ_Shdo, 1, NSTOKES, SRCEXT0, grid_obj)
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0(1) + EXT1(1))
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - FRAC13_*TAU))
                TRANSCELL = 1.0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_*&
                    &(EXT0(1)*SRCEXT1(:) - EXT1(1)*SRCEXT0(:))*DELS&
                    &*(1.0 - 0.05*(EXT1(1) - EXT0(1))*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            radI_Shdo(:) = radI_Shdo(:) + TRANSMIT*SRC(:)*ABSCELL
            SRCEXT1(:) = SRCEXT0(:)
            EXT1(1) = EXT0(1)
            TRANSMIT = TRANSMIT*TRANSCELL
        ! End of sub grid cell loop
        end do

        ! Get the intersection face number (i.e. neighptr index)
        if (SOX .le. SOZ .and. SOX .le. SOY) then
            IFACE = 2 - BITX
            JFACE = 1
            OPENBCFACE = btest(int(flag_Cell_(icell)), 0) .and. isOpenX_
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
            OPENBCFACE = btest(int(flag_Cell_(icell)), 1) .and. isOpenYx_
        else
            IFACE = 6 - BITZ
            JFACE = 3
            OPENBCFACE = .false.
        end if
        ! Get the next cell to go to
        INEXTCELL = i6Cell_Cell_(IFACE, icell)
        if (INEXTCELL .lt. 0) then
            call shd_util__findNextCell(XN, YN, ZN, IFACE, JFACE, icell,&
                &INEXTCELL, grid_obj)
        end if
        ! Get the four grid points on the intersection face for
        !   telling if at the boundary
        !   If going to same or larger face then use previous face
        if (i6Cell_Cell_(IFACE, icell) .ge. 0 .and. .not. OPENBCFACE) then
            KFACE = IFACE
            IC = icell
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = i8Grid_Cell_(GRIDFACE(1, KFACE), IC)
        I2 = i8Grid_Cell_(GRIDFACE(2, KFACE), IC)
        I3 = i8Grid_Cell_(GRIDFACE(3, KFACE), IC)
        I4 = i8Grid_Cell_(GRIDFACE(4, KFACE), IC)
        ! Get the location coordinate
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XN = xyz_Rtm_(1, i8Grid_Cell_(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YN = xyz_Rtm_(2, i8Grid_Cell_(IOCT, INEXTCELL))
            else
                ZN = xyz_Rtm_(3, i8Grid_Cell_(IOCT, INEXTCELL))
            end if
        end if

        ! See if there are valid radiances for the intersection face:
        !   If so, set flag, interpolate radiance, and add in
        !   If not, only prepare for next cell

        !if (TRANSMIT .lt. 1.0e-5) then
        if (TRANSMIT .lt. accSolv_) then
            VALIDRAD = .true.
        else if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0 .and.&
            &GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            call shd_util__interpTriLin(icell, XI, YI, ZI,&
                &NSTOKES, GRIDRAD, 1, NSTOKES, RAD0, grid_obj)
            radI_Shdo(:) = radI_Shdo(:) + TRANSMIT*RAD0(:)
        else
            icell = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do
    SIDE = IFACE
    return
end subroutine shd_src__integrateSrc

subroutine shd_src__computeSolarStokesVect(&
    &SINGLEGEN, COSSCAT,&
    &SOLARMU, MU, DELPHI, NSTOKES,&
    &SCATVECT)

    ! Calculates the Stokes vector for a single scattering of unpolarized
    ! sunlight from the solar direction (SOLARMU) to the outgoing
    ! direction (MU) and with cosine scattering angle (COSSCAT).
    ! The A1 and B1 elements of the phase matrix are calculated from
    ! the Wigner function expansion coefficients (SINGLEGEN).
    ! The plane of polarization is rotated from the scattering plane to
    ! the outgoing meridional plane.  Also divides the output vector by 4*pi.

    integer NSTOKES
    real SINGLEGEN(2, 0:prp_obj_%Nleg), COSSCAT, SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    integer I, J, K, L
    real(RD_) :: X, FCT, A1, B1
    real(RD_), allocatable :: DMM1(:)
    real(RD_) :: SIN_SCAT, SIN_THETA1, SIN_THETA2
    real(RD_) :: SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

    allocate(DMM1(0:prp_obj_%Nleg))

    ! Sum the first Wigner function series, which is actually a Legend
    X = dble(COSSCAT)
    call shd_scat__computeWigner(X, prp_obj_%Nleg, 0, 0, DMM1)
    A1 = 0.0d0
    do L = 0, prp_obj_%Nleg
        FCT = 2.0d0*L + 1.0d0
        A1 = A1 + FCT*dble(SINGLEGEN(1, L))*DMM1(L)
    end do
    SCATVECT(1) = A1*DO4PI_

    if (isPol_) then
        ! If doing polarization, sum the second Wigner function series
        call shd_scat__computeWigner(X, prp_obj_%Nleg, 2, 0, DMM1)
        B1 = 0.0d0
        do L = 0, prp_obj_%Nleg
            FCT = 2.0d0*L + 1.0d0
            B1 = B1 - FCT*dble(SINGLEGEN(2, L))*DMM1(L)
        end do

        ! Rotate the plane of polarization from the scattering plane
        !   to the outgoing meridional plane.
        SIN_SCAT = dsqrt(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = dsqrt(1.d0 - SOLARMU**2)
        SIN_THETA2 = dsqrt(1.d0 - MU**2)
        SINPHI = dsin(dble(DELPHI))
        COSPHI = dcos(dble(DELPHI))
        if (SIN_SCAT .eq. 0.0) then
            SIN2 = 0.0d0
            COS2 = -1.0d0
        else
            SIN2 = SIN_THETA1*SINPHI/SIN_SCAT
            COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        end if
        SIN22 = 2.0d0*SIN2*COS2
        COS22 = 1.0d0 - 2.0d0*SIN2**2
        SCATVECT(2) = B1*COS22*DO4PI_
        SCATVECT(3) = B1*SIN22*DO4PI_
    end if
    if (isPol4_) then
        SCATVECT(4) = 0.0
    end if
    deallocate(DMM1)

    return

end subroutine shd_src__computeSolarStokesVect

subroutine shd_src__findBndryRad(&
    &NSTOKES, XB, YB, MU2, PHI2,&
    &icell, KFACE, GRIDPTR, GRIDPOS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    &RADBND)

    !   Returns the interpolated Stokes radiance at the boundary (RADBND)
    ! Inputs are the boundary location (XB,YB), ray direction away from
    ! boundary (MU2,PHI2), cell number (icell) and face (KFACE) at
    ! the boundary point.

    integer NSTOKES, icell, KFACE, MAXNBC, NTOPPTS, NBOTPTS
    integer GRIDPTR(8, *), BCPTR(MAXNBC, 2)
    integer NMU, NPHI0MAX, NPHI0(NMU)
    real MU2, PHI2, RADBND(NSTOKES)
    real(RD_) :: XB, YB
    real GRIDPOS(3, *)
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE

    integer IL, IM, IU, IP, IBC, J
    real X(4), Y(4), RAD(NSTOKES, 4), U, V
    integer GRIDFACE(4, 6)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/

    ! Loop over the four gridpoints of the cell face on the boundary
    do J = 1, 4
        IP = GRIDPTR(GRIDFACE(J, KFACE), icell)
        X(J) = GRIDPOS(1, IP)
        Y(J) = GRIDPOS(2, IP)

        if (MU2 .lt. 0.0) then
            ! Do a binary search to locate the top boundary point
            IL = 1
            IU = NTOPPTS
            do while (IU - IL .gt. 1)
                IM = (IU + IL)/2
                if (IP .ge. BCPTR(IM, 1)) then
                    IL = IM
                else
                    IU = IM
                end if
            end do
            IBC = IL
            if (BCPTR(IBC, 1) .ne. IP) IBC = IU
            if (BCPTR(IBC, 1) .ne. IP)&
                &stop 'Error [_findBndryRad]: Not at boundary.'
            RAD(:, J) = BCRAD(:, IBC)

        else

            ! Do a binary search to locate the bottom boundary point
            IL = 1
            IU = NBOTPTS
            do while (IU - IL .gt. 1)
                IM = (IU + IL)/2
                if (IP .ge. BCPTR(IM, 2)) then
                    IL = IM
                else
                    IU = IM
                end if
            end do
            IBC = IL
            if (BCPTR(IBC, 2) .ne. IP) IBC = IU
            if (BCPTR(IBC, 2) .ne. IP)&
                &stop 'Error [_findBndryRad]: Not at boundary.'

            if (.not. sfc_obj_%isLambertian) then
                call shd_init__makeBotRadV(NBOTPTS, IBC, IBC, BCPTR(1, 2),&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,&
                    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                    &NSTOKES,&
                    &BCRAD(:, 1 + NTOPPTS))
            end if
            RAD(:, J) = BCRAD(:, NTOPPTS + IBC)
        end if
    end do
    if (X(2) - X(1) .gt. 0.0) then
        U = (XB - X(1))/(X(2) - X(1))
    else
        U = 0.0
    end if
    if (Y(3) - Y(1) .gt. 0.0) then
        V = (YB - Y(1))/(Y(3) - Y(1))
    else
        V = 0.0
    end if
    RADBND(:) = (1 - U)*(1 - V)*RAD(:, 1) + U*(1 - V)*RAD(:, 2)&
        &+ (1 - U)*V*RAD(:, 3) + U*V*RAD(:, 4)
    return
end subroutine shd_src__findBndryRad

subroutine shd_src__rotatePolarPlane(NSTOKES, COSSCAT,&
    &SOLARMU, MU, DELPHI, SCATVECT)

    ! Rotates the plane of polarization from the scattering plane
    ! to the outgoing meridional plane for the solar single
    ! scattered Stokes vector in SINGSCAT.

    integer NSTOKES
    real(RD_) :: COSSCAT
    real SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    real B1
    real(RD_) :: SIN_SCAT, SIN_THETA1, SIN_THETA2
    real(RD_) :: SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

    if (isPol_) then
        B1 = SCATVECT(2)
        SIN_SCAT = dsqrt(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = dsqrt(1.d0 - SOLARMU**2)
        SIN_THETA2 = dsqrt(1.d0 - MU**2)
        SINPHI = dsin(dble(DELPHI))
        COSPHI = dcos(dble(DELPHI))
        if (SIN_SCAT .eq. 0.0) then
            SIN2 = 0.0d0
            COS2 = -1.0d0
        else
            SIN2 = SIN_THETA1*SINPHI/SIN_SCAT
            COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        end if
        SIN22 = 2.0d0*SIN2*COS2
        COS22 = 1.0d0 - 2.0d0*SIN2**2
        SCATVECT(2) = B1*COS22
        SCATVECT(3) = B1*SIN22
    end if
    if (isPol4_) then
        SCATVECT(4) = 0.0
    end if

    return

end subroutine shd_src__rotatePolarPlane

subroutine shd_src__computeSrcOnePnt(&
    &NSTOKES, NLM, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ssa0, ipha0, leg, NR, radI_Shdo,&
    &srcT_Shdo)

    ! Calculates the polarized source function in real generalized
    ! spherical harmonics at one grid point from the radiance in
    ! real generalized spherical harmonics and thermal and/or solar
    ! source information.

    integer, intent(in) :: NSTOKES, NLM, LOFJ(NLM), ipha0, NR
    real, intent(in) :: FLUX0, YLMSUN(prp_obj_%NlegSet, NLM), PLANCK, ssa0
    real, intent(in) :: leg(prp_obj_%NlegSet, 0:prp_obj_%Nleg, *)
    real, intent(in) :: radI_Shdo(NSTOKES, NR)
    character(len=1), intent(in) :: SRCTYPE

    real, intent(out) :: srcT_Shdo(NSTOKES, NLM)

    real, parameter :: C = 3.544907703
    integer J

    srcT_Shdo(:, :) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.  Do different parts of the
    ! matrix multiply depending on NSTOKES.
    do J = 1, NR
        srcT_Shdo(1, J) = srcT_Shdo(1, J) + leg(1, LOFJ(J), ipha0)*radI_Shdo(1, J)
    end do

    if (isPol_) then

        do J = 1, NR
            srcT_Shdo(1, J) = srcT_Shdo(1, J) + leg(5, LOFJ(J), ipha0)*radI_Shdo(2, J)
        end do

        do J = 5, NR
            srcT_Shdo(2, J) = srcT_Shdo(2, J) + leg(5, LOFJ(J), ipha0)*radI_Shdo(1, J)&
                &+ leg(2, LOFJ(J), ipha0)*radI_Shdo(2, J)
            srcT_Shdo(3, J) = srcT_Shdo(3, J) + leg(3, LOFJ(J), ipha0)*radI_Shdo(3, J)
        end do

    end if

    if (isPol4_) then

        do J = 1, NR
            srcT_Shdo(3, J) = srcT_Shdo(3, J) + leg(6, LOFJ(J), ipha0)*radI_Shdo(4, J)
            srcT_Shdo(4, J) = srcT_Shdo(4, J) - leg(6, LOFJ(J), ipha0)*radI_Shdo(3, J)&
                &+ leg(4, LOFJ(J), ipha0)*radI_Shdo(4, J)
        end do

    end if

    ! Do the final multiplication by the single scattering albedo
    srcT_Shdo(:, :) = ssa0*srcT_Shdo(:, :)

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (isSolar_) then

        do J = 1, NLM
            srcT_Shdo(1, J) = srcT_Shdo(1, J)&
                &+ FLUX0*ssa0*leg(1, LOFJ(J), ipha0)*YLMSUN(1, J)
        end do

        if (isPol_) then
            do J = 5, NLM
                srcT_Shdo(2, J) = srcT_Shdo(2, J)&
                    &+ FLUX0*ssa0*leg(5, LOFJ(J), ipha0)*YLMSUN(1, J)
            end do
        end if

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (isThermal_) then
        srcT_Shdo(1, 1) = srcT_Shdo(1, 1) + C*PLANCK
    end if

    return

end subroutine shd_src__computeSrcOnePnt

subroutine shd_src__computeSrcOnePnt_(&
    &NLM, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ssa0, ipha0, leg, NR, radI_Shdo,&
    &srcT_Shdo)

    ! Calculates the scalar source function in spherical harmonics at
    ! one grid point from the radiance in spherical harmonics and
    ! thermal and/or solar source information.

    integer, intent(in) :: NLM, LOFJ(NLM), ipha0, NR
    character(len=1), intent(in) :: SRCTYPE
    real, intent(in) :: FLUX0, YLMSUN(NLM), PLANCK, ssa0
    real, intent(in) :: leg(0:prp_obj_%Nleg, *)
    real, intent(in) :: radI_Shdo(NR)

    real, intent(out) :: srcT_Shdo(NLM)

    real, parameter :: C = 3.544907703
    integer J

    srcT_Shdo(:) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (isSolar_) then

        do J = 1, NLM
            srcT_Shdo(J) = FLUX0*ssa0*leg(LOFJ(J), ipha0)*YLMSUN(J)
        end do

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (isThermal_) then
        srcT_Shdo(1) = srcT_Shdo(1) + C*PLANCK
    end if

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.
    do J = 1, NR
        srcT_Shdo(J) = srcT_Shdo(J) + ssa0*leg(LOFJ(J), ipha0)*radI_Shdo(J)
    end do

    return

end subroutine shd_src__computeSrcOnePnt_

end module shd_src
