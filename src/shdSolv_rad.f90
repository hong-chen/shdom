subroutine shd_comp__interpRad(NSTOKES, OLDNPTS, NPTS, RSHPTR,&
    &RADIANCE, NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)

    ! Interpolates the spherical harmonic radiance array to the new
    !   grid points (>OLDNPTS) from the old grid points.

    implicit none
    integer NSTOKES, OLDNPTS, NPTS, NBCELLS, NCELLS
    integer RSHPTR(NPTS + 1), GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    real RADIANCE(NSTOKES, *), GRIDPOS(3, *)
    integer IP, IC, ICELL, IPARENT, IDIR, I, J, K, DIR1, DIR2, IP1, IP2
    integer IR, IR1, IR2, NR, NR1, NR2
    integer GRIDCORNER(2, 4, 3)
    real RAD1(NSTOKES), RAD2(NSTOKES)
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/

    if (OLDNPTS .le. 0) return

    ! Go through all the non-base grid cells
    do ICELL = NBCELLS + 1, NCELLS
        do IC = 1, 8
            IP = GRIDPTR(IC, ICELL)
            if (IP .gt. OLDNPTS) then
                IPARENT = TREEPTR(1, ICELL)
                ! Go through all 12 edges of parent cell to find which one h
                do IDIR = 1, 3
                    DIR1 = mod(IDIR - 1 + 1, 3) + 1
                    DIR2 = mod(IDIR - 1 + 2, 3) + 1
                    do I = 1, 4
                        IP1 = GRIDPTR(GRIDCORNER(1, I, IDIR), IPARENT)
                        IP2 = GRIDPTR(GRIDCORNER(2, I, IDIR), IPARENT)
                        if (GRIDPOS(DIR1, IP1) .eq. GRIDPOS(DIR1, IP) .and.&
                            &GRIDPOS(DIR2, IP1) .eq. GRIDPOS(DIR2, IP)) goto 120
                    end do
                end do
                stop 'Error [_interpRad]: Point not on an edge.'
120             continue
                ! Interpolate the data from the two parent points at the
                !   ends of the edge; assume the cell was split in half.
                IR1 = RSHPTR(IP1)
                IR2 = RSHPTR(IP2)
                NR1 = RSHPTR(IP1 + 1) - RSHPTR(IP1)
                NR2 = RSHPTR(IP2 + 1) - RSHPTR(IP2)
                NR = max(NR1, NR2)
                IR = RSHPTR(IP)
                RSHPTR(IP + 1) = IR + NR
                do J = 1, NR
                    if (J .le. NR1) then
                        RAD1(:) = RADIANCE(:, IR1 + J)
                    else
                        RAD1(:) = 0.0
                    end if
                    if (J .le. NR2) then
                        RAD2(:) = RADIANCE(:, IR2 + J)
                    else
                        RAD2(:) = 0.0
                    end if
                    RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
                end do
            end if
        end do
    end do

    return
end subroutine shd_comp__interpRad

subroutine shd_solv__truncateRad(HIGHORDERRAD,&
    &NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE,&
    &NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
    &MAXIR, FIXSH, SHACC, RSHPTR)

    ! Computes the next radiance spherical harmonic series truncation
    ! to use based on the current source function truncation and the
    ! scattering properties.  The radiance truncation is the smaller
    ! of that based on a minimum SH term radiance and that needed for
    ! the computation of the source function.  The truncation is always
    ! made for a complete set in the L index.  The radiance SH truncation
    ! is set to a few higher L index than the source function, so that
    ! the source function truncation may grow (since initialized with
    ! a low order truncation).  The minimum radiance truncation is
    ! L=1 so that the mean radiance and net flux terms may be computed,
    ! unless HIGHORDERRAD is on, in which case all terms are kept.

    implicit none
    integer NPTS, NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE, MAXIR
    integer SHPTR(NPTS + 1), RSHPTR(NPTS + 2)
    integer*2 IPHASE(NPTS)
    logical HIGHORDERRAD, FIXSH
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *), RADIANCE(NSTOKES, *)
    real SHACC
    integer IRO, IR, NRO, NR, I, J, K, L, LR, LS, M, MS, ME, NS, NLM
    integer, allocatable :: LOFJ(:)
    logical NOTEND
    real RADMIN

    if (FIXSH) goto 190

    RADMIN = SHACC
    NOTEND = .true.
    ! Make the l index as a function of SH term (J)
    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate (LOFJ(NLM))
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    IRO = 0
    IR = 0
    RSHPTR(1) = 0
    do I = 1, NPTS
        ! Find how many terms wanted for the radiance based on the
        !   radiance threshold (RADMIN) (use the mean radiance as an
        !   upper limit to the SH term radiance).
        !   We don't have radiances for any new points, so don't use them
        NRO = RSHPTR(I + 1) - IRO
        if (NRO .eq. 0) NOTEND = .false.
        if (NOTEND) then
            if (NUMPHASE .gt. 0) then
                K = IPHASE(I)
            else
                K = I
            end if
            LR = 1
            do L = 1, ML
                if (ALBEDO(I)*LEGEN(1, L, K)*RADIANCE(1, IRO + 1) .gt. RADMIN) LR = L
            end do
        else
            LR = ML
        end if
        IRO = IRO + NRO

        ! The radiance L truncation is the smaller of the radiance cutoff
        !   from above and a few more than the source truncation, but
        !   it is at least 1 (so low order terms are always available).
        !   For higher order radiance option, use all terms.
        !   But if out of memory, use the source function truncation.
        NS = max(1, SHPTR(I + 1) - SHPTR(I))
        LS = LOFJ(NS)
        LR = min(LR, LS + ML/8 + 2)
        if (HIGHORDERRAD) LR = ML
        if (LR .le. MM) then
            NR = (LR*(LR + 1)) + LR + 1
        else
            NR = (2*MM + 1)*LR - (MM*(1 + (MM - 1))) + MM + 1
        end if
        IR = IR + NR
        RSHPTR(I + 1) = IR
        if (IR .gt. MAXIR) then
            write (6, *)&
                &'Message [_truncateRad]: Out of memory for more radiance terms.'
            goto 190
        end if
    end do
    RSHPTR(NPTS + 2) = IR
    deallocate (LOFJ)
    return

    ! If out of memory or the SH truncation is fixed then
    ! just set the radiance truncation to that of the source.
190 continue

    RSHPTR(1) = 0
    IR = 0
    do I = 1, NPTS
        NR = max(4, SHPTR(I + 1) - SHPTR(I))
        if (HIGHORDERRAD) then
            if (ML .le. MM) then
                NR = (ML*(ML + 1)) + ML + 1
            else
                NR = (2*MM + 1)*ML - (MM*(1 + (MM - 1))) + MM + 1
            end if
        end if
        IR = IR + NR
        if (IR .gt. MAXIR) then
            write (6, *) 'Error [_truncateRad]: Really out of memory ',&
                &'for more radiance terms.'
            write (6, *) '  Increase MAXIV.'
            stop
        end if
        RSHPTR(I + 1) = IR
    end do
    RSHPTR(NPTS + 2) = IR
    return
end subroutine shd_solv__truncateRad
