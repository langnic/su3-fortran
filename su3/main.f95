program main
    use metropolis
    use utils
    implicit none


    !real :: sdot
    INTEGER, parameter :: Nx = 4
    INTEGER, parameter :: Ny = 4
    INTEGER, parameter :: Nz = 4
    INTEGER, parameter :: Nt = 4
    COMPLEX :: U(4, Nx, Ny, Nz, Nt, 3, 3)
    REAL :: acc
    COMPLEX :: su3(3,3), a(3,3), b(3,3), c(3,3)
    !REAL :: test(2,3) = reshape( (/7, 7, 8, 8, 2, 6/), (/2, 3/) )

    REAL :: re(size(U)), im(size(U)), test(5)


    ! initialize
    call random_seed()
    call random_number(re)
    call random_number(im)
    U = cmplx(reshape(re, SHAPE(U)), reshape(im, SHAPE(U)))

    !su3 = get_small_su3(0.001)

    !write(*,*) 'before'
    !write(*,*) ' ', su3(1,:)
    !write(*,*) ' ', su3(2,:)
    !write(*,*) ' ', su3(3,:)

    !call unitarize(su3)

    !write(*,*) 'after'
    !write(*,*) ' ', su3(1,:)
    !write(*,*) ' ', su3(2,:)
    !write(*,*) ' ', su3(3,:)

    !a = reshape( (/ (4.,1.), (0.,0.), (0.,0.),  (0.,0.), (4.,2.), (0.,0.),  (0.,0.), (0.,0.), (4.,3.) /), (/ 3, 3 /) )
    !b = reshape( (/ (0.,0.), (0.,0.), (2.,0.),  (0.,0.), (2.,0.), (0.,0.),  (2.,0.), (0.,0.), (0.,0.) /), (/ 3, 3 /) )


    !write(*,*) 'a ='
    !write(*,*) ' ', a(1,:)
    !write(*,*) ' ', a(2,:)
    !write(*,*) ' ', a(3,:)

    !write(*,*) 'b ='
    !write(*,*) ' ', b(1,:)
    !write(*,*) ' ', b(2,:)
    !write(*,*) ' ', b(3,:)

    !c = multSu3(a, .False., b, .FALSE.)

    !write(*,*) 'c ='
    !write(*,*) ' ', c(1,:)
    !rite(*,*) ' ', c(2,:)
    !write(*,*) ' ', c(3,:)

    !c = staple(U, 1, 1, 1, 1, 1)

    !write(*,*) 'c ='
    !write(*,*) ' ', c(1,:)
    !write(*,*) ' ', c(2,:)
    !write(*,*) ' ', c(3,:)


    acc = 0
    call metro(U, acc, 0.1)
    write(*,*) 'Done. Acceptance: ', acc



end program main
