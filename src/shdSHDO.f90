module shd_shdo

    use globals
    use shd_pha
    use shd_sfc

    implicit none

    public :: shd_shdo__computeSH

    public :: shd_shdo__makeSHDO
    public :: shd_shdo__makeSfcDO

    public :: shd_shdo__transSHtoDO
    public :: shd_shdo__transSHtoDO_

    public :: shd_shdo__transDOtoSH
    public :: shd_shdo__transDOtoSH_

contains

subroutine shd_shdo__computeSH(&
    &NSHOUT, NSTOKES, NPTS, SRCTYPE,&
    &SOLARMU, SOLARAZ, DIRFLUX,&
    &RSHPTR, ML, MM, RADIANCE, SHTERMS)

    !   Computes the quantities for the spherical harmonic output.
    ! At each grid point the mean radiance and net flux (Fx, Fy, Fz)
    ! are computed.  The output includes the diffuse and direct solar
    ! components.  If NSHOUT=5 then the HIGHORDERRAD flag is set and the
    ! rms of the higher order terms in the radiance series is computed.

    integer NSTOKES, NPTS, RSHPTR(NPTS + 1), ML, MM, NSHOUT
    real SOLARMU, SOLARAZ, DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    character(len=1) :: SRCTYPE
    integer I, IR, J, NLM
    real PI, A, C, C0, IMEAN, FX, FY, FZ, HORMS
    real SECSOL, SINSOL, SOLX, SOLY, SOLM, F0
    real, allocatable :: YLMSUN(:)

    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate (YLMSUN(NLM))

    ! Spherical Harmonic output: individual SH terms
    PI = acos(-1.0)
    A = sqrt(2.0*PI/3.0)
    C = sqrt(4.0*PI/3.0)
    C0 = 1.0/sqrt(4.0*PI)
    SECSOL = 1.0/abs(SOLARMU)
    SINSOL = sqrt(1.0 - SOLARMU**2)
    SOLX = SINSOL*cos(SOLARAZ)*SECSOL
    SOLY = SINSOL*sin(SOLARAZ)*SECSOL
    SOLM = SECSOL/(4.0*PI)
    F0 = 0.0
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
    end if
    do I = 1, NPTS
        IR = RSHPTR(I)
        IMEAN = C0*RADIANCE(1, IR + 1)
        FZ = C*RADIANCE(1, IR + 3)
        if (MM .gt. 0) then
            ! FX = -C*RADIANCE(1,IR+4)   ! old spherical harmonics basis
            ! real generaliz
            FX = -A*(RADIANCE(1, IR + 4) + RADIANCE(1, IR + 2))
        else
            FX = 0.0
        end if
        if (MM .gt. 0) then
            ! FY = -C*RADIANCE(1,IR+2)
            ! real generaliz
            FY = -A*(RADIANCE(1, IR + 2) - RADIANCE(1, IR + 4))
        else
            FY = 0.0
        end if
        if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
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
            if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
                F0 = SECSOL*DIRFLUX(I)
            end if
            do J = 5, RSHPTR(I + 1) - IR
                HORMS = HORMS + (F0*YLMSUN(J) + RADIANCE(1, IR + J))**2
            end do
            SHTERMS(5, I) = sqrt(HORMS)
        end if
    end do
    deallocate (YLMSUN)

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
    allocate (PRC(6, NLM))
    do I = 1, NMU
        X = MU(I)
        PRC(:, :) = 0.0
        call shd_pha__computePLM(.false., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU1(Q, J, I) = PRC(Q, J)
            end do
        end do
        PRC(:, :) = 0.0
        call shd_pha__computePLM(.true., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU2(Q, I, J) = PRC(Q, J)*WTMU(I)
            end do
        end do
    end do
    deallocate (PRC)

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
    real SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1), WAVELEN
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    character(len=2) :: SFCTYPE
    integer IMU, IMU1, IPHI, JMU, JPHI, K1, K2
    real OPI, REFLECT(4, 4), W

    OPI = 1.0/acos(-1.0)
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
                    W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
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

    allocate (MOFJ(NLM), TMP(2*MM + 2))
    allocate (SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

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
                if (NSTOKES .gt. 1) then
                    SUMUV(2, M) = SUMUV(2, M)&
                        &+ CMU1(2, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(5, J, IMU)*INDATA(3, IS + J)
                    SUMUV(3, M) = SUMUV(3, M)&
                        &+ CMU1(6, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(3, J, IMU)*INDATA(3, IS + J)
                end if
                if (NSTOKES .eq. 4) then
                    SUMUV(4, M) = SUMUV(4, M) + CMU1(4, J, IMU)*INDATA(4, IS + J)
                end if
            end do
            SUMCS(:, 0) = SUMUV(:, 0)
            do M = 1, ME
                SUMCS(1, M) = SUMUV(1, M) + SUMUV(1, -M)
                SUMCS(1, -M) = SUMUV(1, -M) - SUMUV(1, M)
            end do
            if (NSTOKES .gt. 1) then
                do M = 1, ME
                    SUMCS(2, M) = SUMUV(2, M) + SUMUV(2, -M)
                    SUMCS(2, -M) = SUMUV(2, -M) - SUMUV(2, M)
                    SUMCS(3, -M) = SUMUV(3, M) - SUMUV(3, -M)
                    SUMCS(3, M) = SUMUV(3, M) + SUMUV(3, -M)
                end do
            end if
            if (NSTOKES .eq. 4) then
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
    deallocate (MOFJ, SUMUV, SUMCS, TMP)

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

    allocate (MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(2*MM + 2))

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
    deallocate (MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transSHtoDO_

subroutine shd_shdo__transDOtoSH(&
    &NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &RSHPTR, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_shdo__transSHtoDO for more).

    integer NPTS, NSTOKES, NSTLEG
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, RSHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NSTOKES, NPHI0MAX, NPTS), OUTDATA(NSTOKES, *)
    real CMU2(NSTLEG, NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:, :), SUMCS(:, :), TMP(:)

    allocate (MOFJ(NLM), TMP(NPHI0))
    allocate (SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

    DELPHI = 2.0*acos(-1.0)/NPHI0

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
        NSH = RSHPTR(I + 1) - RSHPTR(I)
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
        if (NSTOKES .gt. 1) then
            do M = 1, ME
                SUMUV(2, M) = SUMCS(2, M) - SUMCS(2, -M)
                SUMUV(2, -M) = SUMCS(2, M) + SUMCS(2, -M)
                SUMUV(3, M) = SUMCS(3, M) + SUMCS(3, -M)
                SUMUV(3, -M) = SUMCS(3, M) - SUMCS(3, -M)
            end do
        end if
        if (NSTOKES .eq. 4) then
            do M = 1, ME
                SUMUV(4, M) = SUMCS(4, M) + SUMCS(4, -M)
                SUMUV(4, -M) = SUMCS(4, M) - SUMCS(4, -M)
            end do
        end if
        do J = 1, NSH
            M = MOFJ(J)
            OUTDATA(1, IS + J) = OUTDATA(1, IS + J) + CMU2(1, IMU, J)*SUMUV(1, M)
            if (NSTOKES .gt. 1) then
                OUTDATA(2, IS + J) = OUTDATA(2, IS + J) + CMU2(2, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(5, IMU, J)*SUMUV(3, M)
                OUTDATA(3, IS + J) = OUTDATA(3, IS + J) + CMU2(6, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(3, IMU, J)*SUMUV(3, M)
            end if
            if (NSTOKES .eq. 4) then
                OUTDATA(4, IS + J) = OUTDATA(4, IS + J) + CMU2(4, IMU, J)*SUMUV(4, M)
            end if
        end do
        IS = IS + NSH
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)

    return

end subroutine shd_shdo__transDOtoSH

subroutine shd_shdo__transDOtoSH_(NPTS, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &RSHPTR, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_shdo__transSHtoDO for more).

    integer NPTS
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, RSHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NPHI0MAX, NPTS), OUTDATA(*)
    real CMU2(NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:), SUMCS(:), TMP(:)

    allocate (MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(NPHI0))

    DELPHI = 2.0*acos(-1.0)/NPHI0
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
        NSH = RSHPTR(I + 1) - RSHPTR(I)
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
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine shd_shdo__transDOtoSH_

end module shd_shdo
