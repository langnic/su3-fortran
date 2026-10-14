!! lattice -> image decomposition (pure logic, no lattice allocated)
program test_decomposition
    use mod_testing
    use mod_gauge_field
    implicit none

    integer :: nzp, ntp, n, nz, nt, zi, ti, seen(0:63, 0:63)
    logical :: found, all_ok

    call decompose(16, 16, 4, nzp, ntp, found)
    call check('16x16 on 4 images: square 2x2', found .and. nzp == 2 .and. ntp == 2)

    call decompose(12, 12, 9, nzp, ntp, found)
    call check('12x12 on 9 images: 3x3', found .and. nzp == 3 .and. ntp == 3)

    call decompose(6, 10, 5, nzp, ntp, found)
    call check('6x10 on 5 images: non-power-of-2 count, 1x5', found .and. nzp == 1 .and. ntp == 5)

    call decompose(16, 16, 8, nzp, ntp, found)
    call check('16x16 on 8 images: non-square patch grid', found .and. nzp * ntp == 8 .and. nzp /= ntp)

    call decompose(16, 16, 10, nzp, ntp, found)
    call check('16x16 on 10 images: impossible, reported', .not. found)

    call decompose(16, 16, 1, nzp, ntp, found)
    call check('single image: 1x1', found .and. nzp == 1 .and. ntp == 1)

    ! for every case that is accepted: the patches must tile the z-t plane exactly once
    all_ok = .true.
    do nz = 2, 24, 2
        do nt = 2, 24, 2
            do n = 1, 12
                call decompose(nz, nt, n, nzp, ntp, found)
                if (.not. found) cycle
                if (nzp * ntp /= n .or. mod(nz, nzp) /= 0 .or. mod(nt, ntp) /= 0) all_ok = .false.
            end do
        end do
    end do
    call check('accepted decompositions are consistent (nz,nt <= 24, up to 12 images)', all_ok)

    ! every grid position is owned by exactly one image
    seen = 0
    do n = 1, 6
        call decompose(12, 12, 6, nzp, ntp, found)
        call patch_position(n, nzp, zi, ti)
        seen(zi, ti) = seen(zi, ti) + 1
    end do
    call check('images map to distinct patches', found .and. count(seen == 1) == 6 .and. sum(seen) == 6)

    call finish()
end program test_decomposition
