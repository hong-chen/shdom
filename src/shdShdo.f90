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

module shd_shdo

    use globals_const, only: RD_, PI_, PI2_, DPI2_, OPI_, FRAC23_, FRAC43_
    use globals_cfg,   only: isPol_, isPol4_

    use globals_obj,   only: src_obj_

    use shd_util, only: shd_util__computeGausQuad, &
                        shd_util__computeGausQuadDble

    use shd_scat, only: shd_scat__computeWignerDM0, &
                        shd_scat__computeWignerDM0DM2PDM2M

    use shd_sfc,  only: shd_sfc__brdfSurface

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_shdo__computeSH

    public :: shd_shdo__makeSHDO
    public :: shd_shdo__makeSfcDO

    ! unpolarized
    public :: shd_shdo__transSHtoDO_
    public :: shd_shdo__transDOtoSH_

    ! polarized
    public :: shd_shdo__transSHtoDO
    public :: shd_shdo__transDOtoSH

    public :: shd_shdo__makeAngle

    public :: shd_shdo__computeYlm
    public :: shd_shdo__computePlm

    private :: shd_shdo__computeYlm_


contains

subroutine shd_shdo__computeSH(&
    &NSHOUT, NSTOKES, NPTS, SRCTYPE,&
    &SOLARMU, SOLARAZ, DIRFLUX,&
    &iradI_Shdo, ML, MM, radI_Shdo, SHTERMS)

    !   Computes the quantities for the spherical harmonic output.
    ! At each grid point the mean radiance and net flux (Fx, Fy, Fz)
    ! are computed.  The output includes the diffuse and direct solar
    ! components.  If NSHOUT=5 then the isHighOrder flag is set and the
    ! rms of the higher order terms in the radiance series is computed.

    integer NSTOKES, NPTS, iradI_Shdo(NPTS + 1), ML, MM, NSHOUT
    real SOLARMU, SOLARAZ, DIRFLUX(NPTS)
    real radI_Shdo(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    character(len=1) :: SRCTYPE
    integer I, IR, J, NLM
    real A, C, C0, IMEAN, FX, FY, FZ, HORMS
    real SECSOL, SINSOL, SOLX, SOLY, SOLM, F0
    real, allocatable :: YLMSUN(:)

    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate(YLMSUN(NLM))

    ! Spherical Harmonic output: individual SH terms
    A = sqrt(PI_*FRAC23_)
    C = sqrt(PI_*FRAC43_)
    C0 = 1.0/sqrt(4.0*PI_)
    SECSOL = 1.0/abs(SOLARMU)
    SINSOL = sqrt(1.0 - SOLARMU**2)
    SOLX = SINSOL*cos(SOLARAZ)*SECSOL
    SOLY = SINSOL*sin(SOLARAZ)*SECSOL
    SOLM = SECSOL/(4.0*PI_)
    F0 = 0.0
    if (src_obj_%doSolar) then
        call shd_shdo__computeYlm(.true., SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
    end if
    do I = 1, NPTS
        IR = iradI_Shdo(I)
        IMEAN = C0*radI_Shdo(1, IR + 1)
        FZ = C*radI_Shdo(1, IR + 3)
        if (MM .gt. 0) then
            ! FX = -C*radI_Shdo(1,IR+4)   ! old spherical harmonics basis
            ! real generaliz
            FX = -A*(radI_Shdo(1, IR + 4) + radI_Shdo(1, IR + 2))
        else
            FX = 0.0
        end if
        if (MM .gt. 0) then
            ! FY = -C*radI_Shdo(1,IR+2)
            ! real generaliz
            FY = -A*(radI_Shdo(1, IR + 2) - radI_Shdo(1, IR + 4))
        else
            FY = 0.0
        end if
        if (src_obj_%doSolar) then
            IMEAN = IMEAN + DIRFLUX(I)*SOLM
            FX = FX + DIRFLUX(I)*SOLX
            FY = FY + DIRFLUX(I)*SOLY
            FZ = FZ - DIRFLUX(I)
        end if
        SHTERMS(1, I) = IMEAN
        SHTERMS(2, I) = FX
        SHTERMS(3, I) = FY
        SHTERMS(4, I) = FZ
        if (NSHOUT .eq. 5) then
            HORMS = 0.0
            if (src_obj_%doSolar) then
                F0 = SECSOL*DIRFLUX(I)
            end if
            do J = 5, iradI_Shdo(I + 1) - IR
                HORMS = HORMS + (F0*YLMSUN(J) + radI_Shdo(1, IR + J))**2
            end do
            SHTERMS(5, I) = sqrt(HORMS)
        end if
    end do
    deallocate(YLMSUN)

    return

end subroutine shd_shdo__computeSH

subroutine shd_shdo__makeSHDO(NSTLEG, ML, MM, NLM, NMU, NPHI0,&
    &NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,&
    &CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    ! Makes the transformation coefficients for the real generalized
    ! spherical harmonic transform.  The main coefficients are output in
    ! four arrays: 1 is for the SH_TO_DO forward transform, 2 is for the
    ! DO_TO_SH back transform, which contains the discrete angle integration
    ! weights.  The coefficients are divided up into the mu dependent set
    ! CMUn (function of l, m, mu_j), and the phi dependent set CPHIn
    ! (function of m, phi_k) for each mu_j.
    ! The FFTPACK phase coefficients for the FFT in azimuth are also
    ! output in WPHISAVE.

    integer NSTLEG, ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
    logical FFTFLAG(NMU)
    real MU(NMU), PHI(NMU, *), WTMU(NMU), WTDO(NMU, *)
    real CMU1(NSTLEG, NLM, NMU), CMU2(NSTLEG, NMU, NLM)
    real CPHI1(-16:16, 32, NMU), CPHI2(32, -16:16, NMU)
    real WPHISAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, M, Q
    real X, W
    real, allocatable :: PRC(:, :)

    ! Make the mu part of the real generalized spherical harmonics
    allocate(PRC(6, NLM))
    do I = 1, NMU
        X = MU(I)
        PRC(:, :) = 0.0
        call shd_shdo__computePlm(.false., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU1(Q, J, I) = PRC(Q, J)
            end do
        end do
        PRC(:, :) = 0.0
        call shd_shdo__computePlm(.true., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU2(Q, I, J) = PRC(Q, J)*WTMU(I)
            end do
        end do
    end do
    deallocate(PRC)

    ! Make the to and from Fourier coefficients for each mu
    do I = 1, NMU
        if (.not. FFTFLAG(I)) then
            ! If not doing an FFT for this mu then make the DFT coeffici
            W = 1.0/WTMU(I)
            do K = 1, NPHI0(I)
                if (NPHI0(I) .gt. 32 .or. MM .gt. 16) stop &
                    &'Error [_makeSHDO]: Fourier coefficients array exceeded.'
                CPHI1(0, K, I) = 1.0
                CPHI2(K, 0, I) = WTDO(I, K)*W
                do M = 1, MM
                    CPHI1(M, K, I) = cos(M*PHI(I, K))
                    CPHI1(-M, K, I) = sin(M*PHI(I, K))
                    CPHI2(K, M, I) = CPHI1(M, K, I)*WTDO(I, K)*W
                    CPHI2(K, -M, I) = CPHI1(-M, K, I)*WTDO(I, K)*W
                end do
            end do
        else
            ! Precompute the phase factors for the FFTs
            call DFFTI(NPHI0(I), WPHISAVE(1, I))
        end if
    end do
    return
end subroutine shd_shdo__makeSHDO

subroutine shd_shdo__makeSfcDO(SFCTYPE, WAVELEN,&
    &NXSFC, NYSFC, NSFCPAR, SFCPARMS, NSTOKES,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &SFC_BRDF_DO)

    !   Makes the precomputed BRDF surface reflection matrix for
    ! discrete ordinate incident and outgoing directions.  For each
    ! pair of discrete ordinates shd_sfc__brdfSurface is called and the
    ! integration weights are included.

    integer NXSFC, NYSFC, NSFCPAR
    integer NSTOKES, NMU, NPHI0MAX, NPHI0(NMU)
    real SFCPARMS(NSFCPAR, NYSFC+1, NXSFC+1), WAVELEN
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    character(len=2) :: SFCTYPE
    integer IMU, IMU1, IPHI, JMU, JPHI, K1, K2
    real :: REFLECT(4, 4), W

    ! Loop over outgoing discrete ordinates (IMU,IPHI)
    do IMU = NMU/2 + 1, NMU
        IMU1 = IMU - NMU/2
        do IPHI = 1, NPHI0(IMU)
            ! Loop over incident discrete ordinates (JMU,JPHI)
            do JMU = 1, NMU/2
                do JPHI = 1, NPHI0(JMU)
                    ! Get the BRDF for this
                    call shd_sfc__brdfSurface(SFCTYPE(2:2), SFCPARMS(2, 1, 1),&
                        &WAVELEN, MU(IMU), PHI(IMU, IPHI),&
                        &MU(JMU), PHI(JMU, JPHI), NSTOKES, REFLECT)
                    W = OPI_*abs(MU(JMU))*WTDO(JMU, JPHI)
                    SFC_BRDF_DO(:, IMU1, IPHI, :, JMU, JPHI)&
                        &= W*REFLECT(1:NSTOKES, 1:NSTOKES)
                end do
            end do
        end do
    end do

    return

end subroutine shd_shdo__makeSfcDO

subroutine shd_shdo__transSHtoDO(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    &SHPTR, INDATA, OUTDATA)

    !   Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for one zenith angle (IMU).  Successive transforms are
    ! done in zenith and then azimuth angle.  The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The spherical harmonic truncation of the
    ! input array is adaptive (with pointer SHPTR); the last term for
    ! each grid point must be the maximum m for the particular l term.

    integer NPTS, NSTOKES, NSTLEG
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, SHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NSTOKES, *), OUTDATA(NSTOKES, NPHI0MAX, NPTS)
    real CMU1(NSTLEG, NLM, NMU), CPHI1(-16:16, 32, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, L, M, M2, N, IPHI, MS, ME, MR, IS, NSH
    real SUM1
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:, :), SUMCS(:, :), TMP(:)

    allocate(MOFJ(NLM), TMP(2*MM + 2))
    allocate(SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        IS = SHPTR(I)
        NSH = SHPTR(I + 1) - IS
        if (NSH .eq. 0) then
            do K = 1, NPHI0
                OUTDATA(:, K, I) = 0.0
            end do
        else
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
            MS = -ME
            ! First do mu part of real generalized spherical harmonics
            !   transform by summing over l for each m.
            SUMUV(:, :) = 0.0
            do J = 1, NSH
                M = MOFJ(J)
                SUMUV(1, M) = SUMUV(1, M) + CMU1(1, J, IMU)*INDATA(1, IS + J)
                if (isPol_) then
                    SUMUV(2, M) = SUMUV(2, M)&
                        &+ CMU1(2, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(5, J, IMU)*INDATA(3, IS + J)
                    SUMUV(3, M) = SUMUV(3, M)&
                        &+ CMU1(6, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(3, J, IMU)*INDATA(3, IS + J)
                end if
                if (isPol4_) then
                    SUMUV(4, M) = SUMUV(4, M) + CMU1(4, J, IMU)*INDATA(4, IS + J)
                end if
            end do
            SUMCS(:, 0) = SUMUV(:, 0)
            do M = 1, ME
                SUMCS(1, M) = SUMUV(1, M) + SUMUV(1, -M)
                SUMCS(1, -M) = SUMUV(1, -M) - SUMUV(1, M)
            end do
            if (isPol_) then
                do M = 1, ME
                    SUMCS(2, M) = SUMUV(2, M) + SUMUV(2, -M)
                    SUMCS(2, -M) = SUMUV(2, -M) - SUMUV(2, M)
                    SUMCS(3, -M) = SUMUV(3, M) - SUMUV(3, -M)
                    SUMCS(3, M) = SUMUV(3, M) + SUMUV(3, -M)
                end do
            end if
            if (isPol4_) then
                do M = 1, ME
                    SUMCS(4, -M) = SUMUV(4, M) - SUMUV(4, -M)
                    SUMCS(4, M) = SUMUV(4, M) + SUMUV(4, -M)
                end do
            end if

            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                do N = 1, NSTOKES
                    ! For sines and cosines use a real FFT, and rearrange terms
                    TMP(1) = SUMCS(N, 0)
                    MR = 2
                    do M = 1, ME
                        TMP(MR) = 0.5*SUMCS(N, M)
                        TMP(MR + 1) = -0.5*SUMCS(N, -M)
                        MR = MR + 2
                    end do
                    do M2 = MR, NPHI0
                        TMP(M2) = 0.0
                    end do
                    call DFFTB(NPHI0, TMP, WSAVE(1, IMU))
                    do K = 1, NPHI0
                        OUTDATA(N, K, I) = TMP(K)
                    end do
                end do

            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0
                    do N = 1, NSTOKES
                        SUM1 = 0.0
                        do M = MS, ME
                            SUM1 = SUM1 + CPHI1(M, IPHI, IMU)*SUMCS(N, M)
                        end do
                        OUTDATA(N, IPHI, I) = SUM1
                    end do
                end do
            end if
        end if
    end do
    deallocate(MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transSHtoDO

subroutine shd_shdo__transSHtoDO_(&
    &NPTS, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    &SHPTR, INDATA, OUTDATA)

    !   Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for one zenith angle (IMU) for the unpolarized case
    ! Successive transforms are done in zenith and then azimuth angle.
    ! The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The spherical harmonic truncation of the
    ! input array is adaptive (with pointer SHPTR); the last term for
    ! each grid point must be the maximum m for the particular l term.
    ! The number of floating point operations for the DFT is 2*Nmu*Nphi0*Nm.
    ! The number of floating point operations for the zenith angle transform
    ! is 2*Nmu*Nlm.

    integer NPTS
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, SHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(*), OUTDATA(NPHI0MAX, NPTS)
    real CMU1(NLM, NMU), CPHI1(-16:16, 32, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, L, M, M2, N, IPHI, MS, ME, MR, IS, NSH
    real SUM1
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:), SUMCS(:), TMP(:)

    allocate(MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(2*MM + 2))

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        IS = SHPTR(I)
        NSH = SHPTR(I + 1) - IS
        if (NSH .eq. 0) then
            do K = 1, NPHI0
                OUTDATA(K, I) = 0.0
            end do
        else
            ! Figure the max Fourier azimuthal mode we can do for this Nph0
            ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
            MS = -ME
            ! First do mu part of real generalized spherical harmonics
            !   transform by summing over l for each m.
            SUMUV(:) = 0.0
            do J = 1, NSH
                M = MOFJ(J)
                SUMUV(M) = SUMUV(M) + CMU1(J, IMU)*INDATA(IS + J)
            end do
            SUMCS(0) = SUMUV(0)
            do M = 1, ME
                SUMCS(M) = SUMUV(M) + SUMUV(-M)
                SUMCS(-M) = SUMUV(-M) - SUMUV(M)
            end do

            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                ! For sines and cosines use a real FFT, and rearrange terms
                TMP(1) = SUMCS(0)
                MR = 2
                do M = 1, ME
                    TMP(MR) = 0.5*SUMCS(M)
                    TMP(MR + 1) = -0.5*SUMCS(-M)
                    MR = MR + 2
                end do
                do M2 = MR, NPHI0
                    TMP(M2) = 0.0
                end do
                call DFFTB(NPHI0, TMP, WSAVE(1, IMU))
                do K = 1, NPHI0
                    OUTDATA(K, I) = TMP(K)
                end do
            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0
                    SUM1 = 0.0
                    do M = MS, ME
                        SUM1 = SUM1 + CPHI1(M, IPHI, IMU)*SUMCS(M)
                    end do
                    OUTDATA(IPHI, I) = SUM1
                end do
            end if
        end if
    end do
    deallocate(MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transSHtoDO_

subroutine shd_shdo__transDOtoSH(&
    &NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &iradI_Shdo, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_shdo__transSHtoDO for more).

    integer NPTS, NSTOKES, NSTLEG
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, iradI_Shdo(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NSTOKES, NPHI0MAX, NPTS), OUTDATA(NSTOKES, *)
    real CMU2(NSTLEG, NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:, :), SUMCS(:, :), TMP(:)

    allocate(MOFJ(NLM), TMP(NPHI0))
    allocate(SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

    DELPHI = PI2_/NPHI0

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    IS = 0
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        NSH = iradI_Shdo(I + 1) - iradI_Shdo(I)
        ! Figure the max Fourier azimuthal mode we can do for this Nphi0
        ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
        MS = -ME
        ! First do Fourier transform from phi to m for each mu
        if (FFTFLAG(IMU)) then
            do N = 1, NSTOKES
                do K = 1, NPHI0
                    TMP(K) = INDATA(N, K, I)
                end do
                ! For sines and cosines use a real FFT, and rearrange terms
                call DFFTF(NPHI0, TMP, WSAVE(1, IMU))
                SUMCS(N, 0) = TMP(1)*DELPHI
                MR = 2
                do M = 1, ME
                    SUMCS(N, M) = TMP(MR)*DELPHI
                    SUMCS(N, -M) = -TMP(MR + 1)*DELPHI
                    MR = MR + 2
                end do
            end do
        else
            ! Else do a slow DFT
            do M = MS, ME
                do N = 1, NSTOKES
                    SUM1 = 0.0
                    do IPHI = 1, NPHI0
                        SUM1 = SUM1 + CPHI2(IPHI, M, IMU)*INDATA(N, IPHI, I)
                    end do
                    SUMCS(N, M) = SUM1
                end do
            end do
        end if
        ! Then do mu part of real generalized spherical harmonics
        !   transform by summing over l for each m.
        SUMUV(:, ME + 1:MM) = 0.0
        SUMUV(:, -MM:MS - 1) = 0.0
        SUMUV(:, 0) = SUMCS(:, 0)
        do M = 1, ME
            SUMUV(1, M) = SUMCS(1, M) - SUMCS(1, -M)
            SUMUV(1, -M) = SUMCS(1, M) + SUMCS(1, -M)
        end do
        if (isPol_) then
            do M = 1, ME
                SUMUV(2, M) = SUMCS(2, M) - SUMCS(2, -M)
                SUMUV(2, -M) = SUMCS(2, M) + SUMCS(2, -M)
                SUMUV(3, M) = SUMCS(3, M) + SUMCS(3, -M)
                SUMUV(3, -M) = SUMCS(3, M) - SUMCS(3, -M)
            end do
        end if
        if (isPol4_) then
            do M = 1, ME
                SUMUV(4, M) = SUMCS(4, M) + SUMCS(4, -M)
                SUMUV(4, -M) = SUMCS(4, M) - SUMCS(4, -M)
            end do
        end if
        do J = 1, NSH
            M = MOFJ(J)
            OUTDATA(1, IS + J) = OUTDATA(1, IS + J) + CMU2(1, IMU, J)*SUMUV(1, M)
            if (isPol_) then
                OUTDATA(2, IS + J) = OUTDATA(2, IS + J) + CMU2(2, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(5, IMU, J)*SUMUV(3, M)
                OUTDATA(3, IS + J) = OUTDATA(3, IS + J) + CMU2(6, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(3, IMU, J)*SUMUV(3, M)
            end if
            if (isPol4_) then
                OUTDATA(4, IS + J) = OUTDATA(4, IS + J) + CMU2(4, IMU, J)*SUMUV(4, M)
            end if
        end do
        IS = IS + NSH
    end do
    deallocate(MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transDOtoSH

subroutine shd_shdo__transDOtoSH_(NPTS, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &iradI_Shdo, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_shdo__transSHtoDO for more).

    integer NPTS
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, iradI_Shdo(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NPHI0MAX, NPTS), OUTDATA(*)
    real CMU2(NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:), SUMCS(:), TMP(:)

    allocate(MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(NPHI0))

    DELPHI = PI2_/NPHI0
    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    IS = 0
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        NSH = iradI_Shdo(I + 1) - iradI_Shdo(I)
        SUMCS(:) = 0.0
        ! Figure the max Fourier azimuthal mode we can do for this Nphi0
        ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
        MS = -ME
        ! First do Fourier transform from phi to m for each mu
        if (FFTFLAG(IMU)) then
            do K = 1, NPHI0
                TMP(K) = INDATA(K, I)
            end do
            ! For sines and cosines use a real FFT, and rearrange terms
            call DFFTF(NPHI0, TMP, WSAVE(1, IMU))
            SUMCS(0) = TMP(1)*DELPHI
            MR = 2
            do M = 1, ME
                SUMCS(M) = TMP(MR)*DELPHI
                SUMCS(-M) = -TMP(MR + 1)*DELPHI
                MR = MR + 2
            end do
        else
            ! Else do a slow DFT
            do M = MS, ME
                SUM1 = 0.0
                do IPHI = 1, NPHI0
                    SUM1 = SUM1 + CPHI2(IPHI, M, IMU)*INDATA(IPHI, I)
                end do
                SUMCS(M) = SUM1
            end do
        end if
        ! Then do mu part of real generalized spherical harmonics
        !   transform by summing over l for each m.
        SUMUV(0) = SUMCS(0)
        do M = 1, MM
            SUMUV(M) = SUMCS(M) - SUMCS(-M)
            SUMUV(-M) = SUMCS(M) + SUMCS(-M)
        end do
        do J = 1, NSH
            M = MOFJ(J)
            OUTDATA(IS + J) = OUTDATA(IS + J) + CMU2(IMU, J)*SUMUV(M)
        end do
        IS = IS + NSH
    end do
    deallocate(MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transDOtoSH_

subroutine shd_shdo__makeAngle(&
    &NMU, NPHI, ITYPE, NPHI0MAX,&
    &NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the set of angles for the discrete space representation.
    ! The number of mu's (cosine zenith angles) and maximum number of
    ! phi's (azimuth angles) is input.  The actual number of azimuthal
    ! angles for each mu is output (NPHI0).  There are three types of
    ! discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
    ! gaussian grid, and 3 is a reduced double gaussian set.
    ! The output is the NMU mu values, the NPHI0 phi values for each mu,
    ! and the integration weight for each ordinate. The first NMU/2 mu
    ! angles are the downwelling (mu<0) angles.  Also output are the
    ! maximum number of azimuthal angles (NPHI0MAX), the flags for doing
    ! an azimuthal FFT for each mu and the total number of angles.

    integer(kind=I__), intent(in)    :: NMU, NPHI, ITYPE
    integer(kind=I__), intent(inout) :: NPHI0MAX

    integer(kind=I__), intent(out) :: NPHI0(NMU), NANG
    real(kind=R__)   , intent(out) :: MU(NMU), PHI(NMU, *), WTMU(NMU), WTDO(NMU, *)
    logical          , intent(out) :: FFTFLAG(NMU)

    integer(kind=I__) :: J, K, MM
    real(kind=R__)    :: DELPHI
    integer(kind=I__), parameter :: MAXNPHI = 256

    integer(kind=I__) :: GOODNFFT(MAXNPHI)
    data GOODNFFT/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 16,&
        &18, 18, 20, 20, 24, 24, 24, 24, 25, 27, 27, 30, 30, 30, 32, 32,&
        &36, 36, 36, 36, 40, 40, 40, 40, 45, 45, 45, 45, 45, 48, 48, 48,&
        &50, 50, 54, 54, 54, 54, 60, 60, 60, 60, 60, 60, 64, 64, 64, 64,&
        &72, 72, 72, 72, 72, 72, 72, 72, 80, 80, 80, 80, 80, 80, 80, 80,&
        &81, 90, 90, 90, 90, 90, 90, 90, 90, 90, 96, 96, 96, 96, 96, 96,&
        &100, 100, 100, 100, 108, 108, 108, 108,&
        &108, 108, 108, 108, 120, 120, 120, 120,&
        &120, 120, 120, 120, 120, 120, 120, 120,&
        &128, 128, 128, 128, 128, 128, 128, 128,&
        &144, 144, 144, 144, 144, 144, 144, 144,&
        &144, 144, 144, 144, 144, 144, 144, 144,&
        &160, 160, 160, 160, 160, 160, 160, 160,&
        &160, 160, 160, 160, 160, 160, 160, 160,&
        &180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
        &180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
        &192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192,&
        &200, 200, 200, 200, 200, 200, 200, 200,&
        &216, 216, 216, 216, 216, 216, 216, 216,&
        &216, 216, 216, 216, 216, 216, 216, 216,&
        &240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240,&
        &240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240,&
        &256, 256, 256, 256, 256, 256, 256, 256,&
        &256, 256, 256, 256, 256, 256, 256, 256/

    MM = max(0, int(NPHI/2) - 1)
    NANG = 0
    if (ITYPE .ge. 1 .or. ITYPE .le. 3) then
        if (ITYPE .le. 2) then
            call shd_util__computeGausQuad(NMU, MU, WTMU)
            !call shd_util__computeGausQuad_new(NMU, MU, WTMU)
        else if (ITYPE .eq. 3) then
            call shd_util__computeGausQuadDble(NMU, MU, WTMU)
        end if
        do J = 1, NMU
            NPHI0MAX = NPHI
            if (ITYPE .eq. 1) then
                NPHI0(J) = NPHI0MAX
            else
                if (NPHI0MAX .gt. MAXNPHI)&
                    &stop 'Error [_makeAngle]: Exceeded GOODNFFT.'
                NPHI0(J) = int(0.9 + NPHI0MAX*sqrt(1 - MU(J)**2))
                NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J)))
            end if
            DELPHI = PI2_/NPHI0(J)
            do K = 1, NPHI0(J)
                PHI(J, K) = (K - 1)*DELPHI
                WTDO(J, K) = DELPHI*WTMU(J)
            end do
            NANG = NANG + NPHI0(J)
            FFTFLAG(J) = (NPHI0(J) .gt. 14) .or. (MM .gt. 15)
        end do

    else
        stop 'Error [_makeAngle]: Invalid discrete ordinate type.'
    end if

    return

end subroutine shd_shdo__makeAngle

subroutine shd_shdo__computeYlm(TRANSPOSED, MU, PHI, ML, MM, NSTLEG, YR)

    !   This subroutine computes a set of normalized real generalized
    ! spherical harmonic functions (YR) for a particular direction mu,phi.
    ! (Doicu et al. 2013; http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mode.
    ! NSTLEG is the number of matrix elements returned, either 1 for the
    ! 1,1 term or 6 for the 6 unique elements of the 4x4 matrix.
    ! J = 2*(L*(L+1))/2 + M+1  for L<=MM
    ! J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM

    integer ML, MM, NSTLEG
    real MU, PHI
    logical TRANSPOSED
    real YR(NSTLEG, *)
    integer J, M, MABS, L
    real(RD_) :: X, FCT, P1, P2, P3, COSM, SINM, SIGN
    real(RD_), allocatable :: DM0(:), DM2P(:), DM2M(:)

    if (NSTLEG .eq. 1) then
        call shd_shdo__computeYlm_(MU, PHI, ML, MM, YR)
        return
    end if

    allocate(DM0(0:ML), DM2P(0:ML), DM2M(0:ML))
    X = dble(MU)

    FCT = 1.d0/dsqrt(DPI2_)

    ! sign accounting for transpose matrices
    if (.not. TRANSPOSED) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if

    ! ..................................................................
    !                               M = 0
    M = 0
    call shd_scat__computeWignerDM0DM2PDM2M(X, ML, M, DM0, DM2P, DM2M)
    do L = 0, ML
        if (L .le. MM) then
            J = L*(L + 1) + M + 1
        else
            J = (2*MM + 1)*L - MM**2 + M + 1
        end if
        P1 = FCT*DM0(L)
        P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
        P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

    ! --- Unm = Pnm*cos(m*phi) ---
        YR(1, J) = P1
        if (NSTLEG .eq. 6) then
            YR(2, J) = P2
            YR(3, J) = P2
            YR(4, J) = P1
            YR(5, J) = P3
            YR(6, J) = P3
        end if
    end do

    ! ..................................................................
    !                               M =/ 0
    do MABS = 1, MM
        call shd_scat__computeWignerDM0DM2PDM2M(X, ML, MABS, DM0, DM2P, DM2M)
        COSM = cos(MABS*PHI)
        SINM = sin(MABS*PHI)
        do L = MABS, ML
    ! ..............................................................
    !                               M > 0
            M = MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if
            P1 = FCT*DM0(L)
            P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
            P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

            ! --- Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi) ---
            YR(1, J) = P1*COSM - P1*SINM
            if (NSTLEG .eq. 6) then
                YR(2, J) = P2*COSM - P2*SINM
                YR(3, J) = P2*COSM + P2*SINM
                YR(4, J) = P1*COSM + P1*SINM
                YR(5, J) = P3*COSM - SIGN*P3*SINM
                YR(6, J) = P3*COSM + SIGN*P3*SINM
            end if
    ! ..............................................................
    !                               M < 0
            M = -MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! --- Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi) ---
            YR(1, J) = P1*SINM + P1*COSM
            if (NSTLEG .eq. 6) then
                YR(2, J) = P2*SINM + P2*COSM
                YR(3, J) = P2*SINM - P2*COSM
                YR(4, J) = P1*SINM - P1*COSM
                YR(5, J) = P3*SINM + SIGN*P3*COSM
                YR(6, J) = P3*SINM - SIGN*P3*COSM
            end if
        end do
    end do

    deallocate(DM0, DM2P, DM2M)

    return

end subroutine shd_shdo__computeYlm

subroutine shd_shdo__computeYlm_(MU, PHI, ML, MM, YR)

    !   This subroutine computes a set of the unpolarized (I-I) element
    ! normalized real generalized spherical harmonic functions (YR) for
    ! a particular direction mu,phi.  (Doicu et al. 2013;
    ! http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mod
    !   J = 2*(L*(L+1))/2 + M+1  for L<=MM
    !   J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM

    integer ML, MM
    real MU, PHI
    real YR(*)
    integer J, M, L
    real(RD_) :: X, FCT, COSM, SINM
    real(RD_), allocatable :: DM0(:)

    allocate(DM0(0:ML))
    X = dble(MU)

    FCT = 1.d0/dsqrt(DPI2_)

    do M = 0, MM
        call shd_scat__computeWignerDM0(X, ML, M, DM0)
        DM0(:) = FCT*DM0(:)
        if (M .gt. 0) then
            COSM = cos(M*PHI)
            SINM = sin(M*PHI)
        else
            COSM = 1.0d0
            SINM = 0.0d0
        end if
        do L = M, ML
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! M>=0:  Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi)
            YR(J) = (COSM - SINM)*DM0(L)

            J = J - 2*M
            ! M<0:  Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi)
            YR(J) = (COSM + SINM)*DM0(L)
        end do
    end do

    deallocate(DM0)

    return

end subroutine shd_shdo__computeYlm_

subroutine shd_shdo__computePlm(TRANSPOSED, MU, ML, MM, PRC)

    integer ML, MM
    real MU
    logical TRANSPOSED
    real PRC(6, *)
    ! --- local variables ---
    integer J, M, MABS, L
    real(RD_) :: X, FCT, P1, P2, P3, SIGN
    real(RD_), allocatable :: DM0(:), DM2P(:),&
        &DM2M(:)

    ! --- allocate ---
    allocate(DM0(0:ML), DM2P(0:ML), DM2M(0:ML))

    X = dble(MU)

    ! --- normalization constant ---
    FCT = 1.d0/dsqrt(DPI2_)

    ! --- sign accounting for transpose matrices ---
    if (.not. TRANSPOSED) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if
    ! ..................................................................
    !                            M = 0
    ! ..................................................................
    M = 0
    call shd_scat__computeWignerDM0DM2PDM2M&
        &(X, ML, M,&
        &DM0, DM2P, DM2M)
    do L = 0, ML
        if (L .le. MM) then
            J = L*(L + 1) + M + 1
        else
            J = (2*MM + 1)*L - MM**2 + M + 1
        end if
        P1 = FCT*DM0(L)
        P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
        P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

        ! --- PRCnm = Pnm ---
        PRC(1, J) = P1
        PRC(2, J) = P2
        PRC(3, J) = P2
        PRC(4, J) = P1
        PRC(5, J) = P3
        PRC(6, J) = P3
    end do
    ! ..................................................................
    !                               M =/ 0
    ! ..................................................................
    do MABS = 1, MM
        call shd_scat__computeWignerDM0DM2PDM2M&
            &(X, ML, MABS,&
            &DM0, DM2P, DM2M)
        do L = MABS, ML
    ! ..................................................................
    !                               M > 0
    ! ..................................................................
            M = MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if
            P1 = FCT*DM0(L)
            P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
            P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

            ! --- PRCnm = Pnm ---
            PRC(1, J) = P1
            PRC(2, J) = P2
            PRC(3, J) = P2
            PRC(4, J) = P1
            PRC(5, J) = P3
            PRC(6, J) = P3
    ! ..................................................................
    !                               M < 0
    ! ..................................................................
            M = -MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! --- PRCnm = D*Pnm ---
            PRC(1, J) = P1
            PRC(2, J) = P2
            PRC(3, J) = -P2
            PRC(4, J) = -P1
            PRC(5, J) = SIGN*P3
            PRC(6, J) = -SIGN*P3
        end do
    end do
    deallocate(DM0, DM2P, DM2M)

    return

end subroutine shd_shdo__computePlm
 
end module shd_shdo
