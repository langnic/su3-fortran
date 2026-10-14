!! the distributed lattice must be the same lattice for any number of images:
!!  - the link index array is a bijection onto the (image, local index) storage
!!  - the plaquette of a deterministic, non-trivial configuration matches a serial reference
!! The lattice size scales with the image count so that every count can be distributed.
program test_distribution
    use mod_testing
    use mod_gauge_field
    use mod_measure
    implicit none

    type(gauge_field) :: u
    integer, parameter :: nx = 4, ny = 4
    integer :: nz, nt, n, c, l, k, j, i, m, owner, idx
    logical, allocatable :: used(:,:)
    logical :: bijective
    real :: pref

    n = num_images()
    nz = 2 * n
    nt = 4 * n
    call u%initialize(nx, ny, nz, nt, .false.)

    ! bijection: each image owns exactly size(data,3) links, each local slot used once
    allocate(used(size(data, 3), n))
    used = .false.
    bijective = .true.
    do m = 1, 4
        do i = 1, nt
            do j = 1, nz
                do k = 1, ny
                    do l = 1, nx
                        idx = u%indices(0,l,k,j,i,m)
                        owner = u%indices(1,l,k,j,i,m)
                        if (owner < 1 .or. owner > n .or. idx < 1 .or. idx > size(data, 3)) then
                            bijective = .false.
                        else if (used(idx, owner)) then
                            bijective = .false.
                        else
                            used(idx, owner) = .true.
                        end if
                    end do
                end do
            end do
        end do
    end do
    call check('index array is a bijection onto the storage', bijective .and. all(used))

    ! set U_m(x) = diag(e^{i phi}, e^{-i phi}, 1) on the links owned by this image
    do m = 1, 4
        do i = 1, nt
            do j = 1, nz
                do k = 1, ny
                    do l = 1, nx
                        if (u%indices(1,l,k,j,i,m) == this_image()) then
                            data(:,:,u%indices(0,l,k,j,i,m)) = diag_link(phi(l, k, j, i, m))
                        end if
                    end do
                end do
            end do
        end do
    end do
    sync all

    pref = reference_plaquette()
    call check_range('configuration is not trivial', pref, 0.0, 0.95)
    call check_range('plaquette matches serial reference', plaquette(u), pref - 2e-5, pref + 2e-5)
    call check_range('links are SU(3)', max_unitarity_violation(), 0.0, 1e-6)

    call finish()

contains

    !! smooth-ish pseudo-random phase of the link in direction m at site (x,y,z,t)
    pure function phi(x, y, z, t, m)
        integer, intent(in) :: x, y, z, t, m
        real :: phi
        phi = 0.9 * sin(1.3 * x + 2.1 * y + 0.7 * z + 1.9 * t + 0.5 * m)
    end function phi

    pure function diag_link(a) result(mtx)
        real, intent(in) :: a
        complex :: mtx(3,3)
        mtx = 0
        mtx(1,1) = cmplx(cos(a), sin(a))
        mtx(2,2) = cmplx(cos(a), -sin(a))
        mtx(3,3) = 1
    end function diag_link

    !! average plaquette computed directly from the phases (periodic boundaries):
    !! for diagonal links tr(plaquette) = 1 + 2 cos(Phi), Phi = phi_a(x) + phi_b(x+a) - phi_a(x+b) - phi_b(x)
    function reference_plaquette() result(res)
        real :: res
        real(kind=8) :: total, ph
        integer :: x(4), x1, x2, x3, x4, a, b, e_a(4), e_b(4), xa(4), xb(4)

        total = 0
        do x4 = 1, nt
            do x3 = 1, nz
                do x2 = 1, ny
                    do x1 = 1, nx
                        x = [x1, x2, x3, x4]
                        do a = 1, 3
                            do b = a + 1, 4
                                e_a = 0
                                e_a(a) = 1
                                e_b = 0
                                e_b(b) = 1
                                xa = wrap(x + e_a)
                                xb = wrap(x + e_b)
                                ph = phi(x(1), x(2), x(3), x(4), a) + phi(xa(1), xa(2), xa(3), xa(4), b) &
                                    - phi(xb(1), xb(2), xb(3), xb(4), a) - phi(x(1), x(2), x(3), x(4), b)
                                total = total + (1 + 2 * cos(ph)) / 3
                            end do
                        end do
                    end do
                end do
            end do
        end do
        res = real(total / (6.0d0 * nx * ny * nz * nt))
    end function reference_plaquette

    pure function wrap(x) result(y)
        integer, intent(in) :: x(4)
        integer :: y(4)
        y(1) = modulo(x(1) - 1, nx) + 1
        y(2) = modulo(x(2) - 1, ny) + 1
        y(3) = modulo(x(3) - 1, nz) + 1
        y(4) = modulo(x(4) - 1, nt) + 1
    end function wrap

end program test_distribution
