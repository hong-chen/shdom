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

module shd_scat

    use globals_const, only: RD_, DPI_, DOPI_, DO4PI_
    use globals_cfg,   only: isPol_

    implicit none
    private

    public :: shd_scat__computePhase
    public :: shd_scat__computeLegendre

    public :: shd_scat__computeWigner
    public :: shd_scat__computeWignerDM0
    public :: shd_scat__computeWignerDM0DM2PDM2M

    private :: wignerDMM1N0
    private :: wignerDMM0N0

contains

subroutine shd_scat__computePhase(&
    &NSCATANGLE, NUMPHASE, NSTPHASE,&
    &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)

    ! Description:
    !   Compute the phase matrix elements I-I and I-Q as a function
    !   of scattering angle for solar direct scattering for all the
    !   tabulated phase functions. Output is in PHASETAB.
    !   NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.

    integer NSCATANGLE, NUMPHASE, NSTPHASE
    integer NSTOKES, NSTLEG, ML, NLEG
    real LEGEN(NSTLEG, 0:NLEG, NUMPHASE)
    real PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE)
    integer IPH, I, J, L
    real(RD_) :: COSSCAT, FCT, F, X, A1, B1
    real(RD_), allocatable :: UNSCLEGEN(:, :)
    real(RD_), allocatable :: DMM1(:), DMM2(:)

    allocate(UNSCLEGEN(2, 0:NLEG), DMM1(0:NLEG), DMM2(0:NLEG))

    do J = 1, NSCATANGLE
        COSSCAT = cos(DPI_*DFLOAT(J - 1)/(NSCATANGLE - 1))
        X = dble(COSSCAT)
        call shd_scat__computeWigner(X, NLEG, 0, 0, DMM1)
        if (isPol_) then
            call shd_scat__computeWigner(X, NLEG, 2, 0, DMM2)
        end if

        ! Unscaled the needed Legendre coefficients, but divide by 1-f
        !   because extinction is still scaled (for TMS method)
        do IPH = 1, NUMPHASE
            F = LEGEN(1, ML + 1, IPH)
            do L = 0, NLEG
                if (L .le. ML) then
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH) + F/(1 - F)
                else
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH)/(1 - F)
                end if
                if (NSTLEG .gt. 1) then
                    if (L .le. ML) then
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)
                    else
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)/(1 - F)
                    end if
                end if
            end do

            ! Sum the first Wigner function series, which is actually a Le
            A1 = 0.0d0
            do L = 0, NLEG
                FCT = 2.0d0*L + 1.0d0
                A1 = A1 + FCT*dble(UNSCLEGEN(1, L))*DMM1(L)
            end do
            if (A1 .le. 0.0) then
                write(6, *) 'Error [_computePhase]: negative phase function',&
                    &' for tabulated phase function: ', IPH, J, A1
                stop
            end if
            PHASETAB(1, IPH, J) = sngl(A1*DO4PI_)

            if (isPol_) then
                ! If doing polarization, sum the second Wigner function series
                B1 = 0.0d0
                do L = 0, NLEG
                    FCT = 2.0d0*L + 1.0d0
                    B1 = B1 - FCT*dble(UNSCLEGEN(2, L))*DMM2(L)
                end do
                PHASETAB(2, IPH, J) = sngl(B1*DO4PI_)
            end if
        end do
    end do

    deallocate(UNSCLEGEN, DMM1, DMM2)

    return

end subroutine shd_scat__computePhase

subroutine shd_scat__computeLegendre(COSSCAT, NLEG, P)

    !   This subroutine computes a set of Legendre polynomials for
    ! a particular scattering angle COSSCAT.  NLEG is the maximum term.
    ! The Legendre functions evaluated at COSSCAT are returned in
    ! P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).

    integer, intent(in) :: NLEG
    real(RD_) :: COSSCAT, P(0:NLEG)
    integer L
    real(RD_) :: X, PL, PL1, PL2

    X = dble(COSSCAT)
    if (X*X .gt. 1.) stop 'Error [_computeLegendre]: <COSSCAT> is greater than 1.'

    ! Use the stable upward recursion on l, starting from P_0
    PL2 = 1.0d0
    P(0) = PL2
    if (NLEG .gt. 1) then
        PL1 = X
        P(1) = X
    end if

    do L = 2, NLEG
        PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        P(L) = PL
        PL2 = PL1
        PL1 = PL
    end do

    return

end subroutine shd_scat__computeLegendre

subroutine shd_scat__computeWigner(X, NRANK, M, M1, DMM1)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0, M1 >= 0,  N = 0,...,NRANK, and
    ! |               -1 < X= COS(BETA) < 1
    ! ----------------------------------------------------------------

    integer NRANK, M, M1
    real(RD_) :: X, DMM1(0:NRANK)

    integer N0, N
    real(RD_) :: FACT1, FACT2

    N0 = max(M, M1)
    DMM1 = 0.d0
    if (N0 .eq. 0) then
        DMM1(0) = 1.d0
        DMM1(1) = X
        do N = 1, NRANK - 1
            FACT1 = dble(2*N + 1)*X/dble(N + 1)
            FACT2 = dble(N)/dble(N + 1)
            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    else
        DMM1(N0) = wignerDMM1N0(X, M, M1)
        do N = N0, NRANK - 1
            FACT1 = dble(N*(N + 1))*X - dble(M*M1)
            FACT1 = FACT1/dsqrt(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/dsqrt(dble((N + 1)**2 - M1**2))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)

            FACT2 = dsqrt(dble(N**2 - M**2))&
                  * dsqrt(dble(N**2 - M1**2))
            FACT2 = FACT2/dsqrt(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/dsqrt(dble((N + 1)**2 - M1**2))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    end if
    return
end subroutine shd_scat__computeWigner

subroutine shd_scat__computeWignerDM0DM2PDM2M(&
    &MIU, NRANK, M, DM0, DM2P, DM2M)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0; M1 = 0, 2, -2; N = 0,...,NRANK, and
    ! |              -1 <= MIU <= 1
    ! |  For M1 = -2, we use the symmetry relation
    ! |               DM-M1(MIU) = (-1)**(N+M) * DMM1(-MIU)
    ! ----------------------------------------------------------------

    integer NRANK, M
    real(RD_) :: MIU
    real(RD_) :: DM0(0:NRANK), DM2P(0:NRANK), DM2M(0:NRANK)
    integer N0, N
    real(RD_) :: XP, XM, FACT1, FACT2, FACTP, FACTM

    !for M1 = 0 and M1 = 2
    XP = MIU
    !for M1 = - 2
    XM = -MIU
    ! ................................................................
    ! DM0 for M1 = 0
    ! ................................................................
    DM0 = 0.d0
    if (M .le. NRANK) then
        if (M .eq. 0) then
            DM0(0) = 1.d0
            DM0(1) = XP
            do N = 1, NRANK - 1
                FACT1 = dble(2*N + 1)*XP/dble(N + 1)
                FACT2 = dble(N)/dble(N + 1)
                DM0(N + 1) = FACT1*DM0(N) - FACT2*DM0(N - 1)
            end do
        else
            DM0(M) = wignerDMM1N0(XP, M, 0)
            do N = M, NRANK - 1
                FACT1 = dble(N*(N + 1))*XP
                FACT1 = FACT1/dsqrt(dble((N + 1)**2 - M**2))
                FACT1 = FACT1/dble(N + 1)
                FACT1 = FACT1*dble(2*N + 1)/dble(N)

                FACT2 = dsqrt(dble(N**2 - M**2))*dble(N)
                FACT2 = FACT2/dsqrt(dble((N + 1)**2 - M**2))
                FACT2 = FACT2/dble(N + 1)
                FACT2 = FACT2*dble(N + 1)/dble(N)

                DM0(N + 1) = FACT1*DM0(N) - FACT2*DM0(N - 1)
            end do
        end if
        do N = 0, NRANK
            DM0(N) = dsqrt(N + 0.5d0)*DM0(N)
        end do
    end if
    ! ................................................................
    ! DM2P and DM2M for M1 = 2 and M1 = -2
    ! ................................................................
    N0 = max(M, 2)
    DM2P = 0.d0
    DM2M = 0.d0
    if (M .le. NRANK .and. NRANK .ge. 2) then
        DM2P(N0) = wignerDMM1N0(XP, M, 2)
        DM2M(N0) = wignerDMM1N0(XM, M, 2)
        do N = N0, NRANK - 1
            FACTP = dble(N*(N + 1))*XP - dble(2*M)
            FACTM = dble(N*(N + 1))*XM - dble(2*M)

            FACT1 = 1.d0/dsqrt(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/dsqrt(dble((N + 1)**2 - 4))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)

            FACT2 = dsqrt(dble(N**2 - M**2))&
                  * dsqrt(dble(N**2 - 4))
            FACT2 = FACT2/dsqrt(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/dsqrt(dble((N + 1)**2 - 4))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DM2P(N + 1) = FACTP*FACT1*DM2P(N) - FACT2*DM2P(N - 1)
            DM2M(N + 1) = FACTM*FACT1*DM2M(N) - FACT2*DM2M(N - 1)
        end do
        do N = 0, NRANK
            DM2M(N) = (-1.d0)**(N + M)*DM2M(N)

            DM2P(N) = dsqrt(N + 0.5d0)*DM2P(N)
            DM2M(N) = dsqrt(N + 0.5d0)*DM2M(N)
        end do
    end if

    return

end subroutine shd_scat__computeWignerDM0DM2PDM2M

subroutine shd_scat__computeWignerDM0(X, NRANK, M, DM0)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for      |
    ! |               M >= 0; M1 = 0; N = 0,...,NRANK, and           |
    ! |              -1 <= X <= 1                                    |
    ! ----------------------------------------------------------------

    integer NRANK, M
    real(RD_) :: X
    real(RD_) :: DM0(0:NRANK)
    integer N

    ! DM0 for M1 = 0
    DM0(:) = 0.d0
    if (M .le. NRANK) then
        if (M .eq. 0) then
            DM0(0) = 1.0d0
            DM0(1) = X
            do N = 1, NRANK - 1
                DM0(N + 1) = ((2*N + 1)*X*DM0(N) - N*DM0(N - 1))/(N + 1)
            end do
        else
            DM0(M) = wignerDMM0N0(X, M)
            do N = M, NRANK - 1
                DM0(N + 1) = ((2*N + 1)*X*DM0(N)&
                    &- dsqrt(dble(N**2 - M**2))*DM0(N - 1))&
                    &/dsqrt(dble((N + 1)**2 - M**2))
            end do
        end if
        do N = 0, NRANK
            DM0(N) = dsqrt(N + 0.5d0)*DM0(N)
        end do
    end if
    return
end subroutine shd_scat__computeWignerDM0

function wignerDMM1N0(X, M, M1) result(DMM1N0)

    ! ---------------------------------------------------------------
    ! | THE ROUTINE COMPUTES THE WIGNER FUNCTIONS FOR               |
    ! |     M >= 0, M1 >= 0, AND N0 = MAX(M,M1)                     |
    ! ---------------------------------------------------------------

    integer M, M1
    real(RD_) :: X, DMM1N0

    integer P, MAXM, MINM
    real(RD_) :: PROD, FACT, CMM1

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
            FACT = dsqrt(dble(M + M1 + P)/dble(P))
            PROD = PROD*FACT
        end do
        FACT = dsqrt((1.d0 - X)/2.d0)
        DMM1N0 = CMM1*PROD*FACT**(MAXM - MINM)
        FACT = dsqrt((1.d0 + X)/2.d0)
        DMM1N0 = DMM1N0*FACT**(MAXM + MINM)
    end if

    return

end function wignerDMM1N0

function wignerDMM0N0(X, M) result(DMM0N0)

    ! The routine computes the Wigner functions for M >= 0, M1 = 0, a
    ! N0 = M

    integer M
    real(RD_) :: X, DMM0N0
    integer P
    real(RD_) :: CM, PROD

    if (M .eq. 0) then
        DMM0N0 = 1.d0
    else
        CM = (-1.d0)**M
        PROD = 1.d0
        do P = 1, M
            PROD = PROD*dsqrt((M + P)/dble(P))
        end do
        DMM0N0 = CM*PROD*(0.5d0*dsqrt((1.d0 - X)*(1.d0 + X)))**M
    end if

    return

end function wignerDMM0N0

end module shd_scat
