program main
    use metropolis
    use utils
    implicit none


    INTEGER, parameter :: Nx = 4, Ny = 4, Nz = 4, Nt = 4
    INTEGER, parameter :: nMeasure = 100, nDiscard = 10, nThermalize = 100
    INTEGER :: i, j
    COMPLEX :: U(4, Nx, Ny, Nz, Nt, 3, 3)
    REAL, parameter :: beta = 0.1
    REAL :: acc = 0
    REAL :: re(size(U)), im(size(U))

    ! observables
    REAL, dimension(nMeasure) :: sp  ! plaquette
    REAL, dimension(nMeasure, Nt) :: spt, sps    ! spatial and temporal plaquettes (time sliced)


    ! initialize (random)
    call random_seed()
    call random_number(re)
    call random_number(im)
    U = cmplx(reshape(re, SHAPE(U)), reshape(im, SHAPE(U)))

    ! thermalize
    do i = 0, nThermalize
        call metro(U, acc, beta)
    end do
    acc = 0

    ! measure
    do i = 1, nMeasure

        call wilson(U, spt(i,:), sps(i,:), sp(i))
        do j = 1, nDiscard    ! discard nine configurations
            call metro(U, acc, beta)
        end do

    end do

    ! output mean plaquette
    write(*,*) 'Mean plaquette: ', SUM(sp) / nMeasure
    write(*,*) 'Acceptance: ', acc



end program main
