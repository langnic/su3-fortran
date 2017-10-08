program main
    use metropolis
    use utils
    use su3_matrix_module
    use gauge_field_module
    implicit none


    INTEGER, parameter :: Nx = 4, Ny = 4, Nz = 4, Nt = 4
    INTEGER, parameter :: nMeasure = 100, nDiscard = 10, nThermalize = 100
    INTEGER :: i, j
    COMPLEX :: U(4, Nx, Ny, Nz, Nt, 3, 3)
    REAL, parameter :: beta = 0.1
    REAL :: acc = 0
    REAL :: re(size(U)), im(size(U))

    type(gauge_field) :: uf, vf, wf
    !type(gauge_field), dimension(5) :: af
    type(su3_matrix) :: su, sv
    type(su3_matrix), allocatable :: test(:)
    complex, allocatable, target :: data2(:,:,:)[:], data3(:,:,:)[:]
    complex, allocatable :: tmp(:,:,:,:)[:]
    complex, dimension(:,:,:), pointer :: ptr
    complex, dimension(:,:,:), pointer :: ptr2
    complex, allocatable :: datatest(:,:,:)
    complex, allocatable :: datatest2(:,:,:)[:]
    !type(gauge_field), pointer :: ptr_gauge ! geht nich

    ! observables
    REAL, dimension(nMeasure) :: sp  ! plaquette
    REAL, dimension(nMeasure, Nt) :: spt, sps    ! spatial and temporal plaquettes (time sliced)


    !init uf
    !uf = gauge_field(4,4,4,4)
    call uf%initialize(4,4,4,4)
    call vf%initialize(4,4,4,4)
    call wf%initialize(4,4,4,4)

    if (allocated(uf%data)) then
        write(*,*) 'allocated'
    else
        write(*,*) 'unallocated'
    end if
    su%m = reshape((/ 1, 1, 1,  1, 1, 1,  1, 1, 1 /), shape(su%m))
    sv = su
    sv = 3

    su = uf + wf

    allocate(data3(1,3,3)[*])
    data3(1,2,2) = 2.0
    allocate(data2(1,3,3)[*], source=data3 )
    call move_alloc(data2, uf%data)
    write(*,*) 'uf data(1,2,2): ', uf%data(1,2,2)

    allocate(tmp(3,1,3,3)[*])
    tmp(1,:,:,:) = data3

    ptr => data3(1:1,:,:)
    ptr2 => ptr
    write(*,*) ptr(1,2,2)
    write(*,*) ptr2(1,2,2)

    allocate(datatest(1:1,3,3), source=ptr)
    write(*,*) datatest(1,2,2)


    !uf%data = data2 + data3
    !uf%data(0) = su
    !write(*,*) uf%data(0)%m

    !allocate(test(3))
    !test(0) = sv

    !write(*,*) sv%m
    !uf%data(0) = su
    !uf%data(0) = type(
    !su = uf%get(0)
    !uf%data(0)%elements(0,0) = 3.0
    !write(*,*) 'bla: ', uf%data(0)%elements(0,0)

    ! initialize (random) - FIXME: not SU, just random complex!
    call random_seed()
    call random_number(re)
    call random_number(im)
    U = cmplx(reshape(re, SHAPE(U)), reshape(im, SHAPE(U)))

    ! thermalize
    !do i = 0, nThermalize
    !    call metro(U, acc, beta)
    !end do
    acc = 0

    ! measure
    do i = 1, nMeasure

        !call wilson(U, spt(i,:), sps(i,:), sp(i))
        !do j = 1, nDiscard    ! discard nine configurations
        !    call metro(U, acc, beta)
        !end do

    end do

    ! output mean plaquette
    write(*,*) 'Mean plaquette: ', SUM(sp) / nMeasure
    write(*,*) 'Acceptance: ', acc



end program main
