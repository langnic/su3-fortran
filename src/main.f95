program main

    use mod_su3_field_var
    use mod_gauge_field
    use mod_real_field_var
    use mod_interfaces
    use mod_metropolis
    use mod_measure
    implicit none


    type(gauge_field) :: u
    real :: acc = 0
    real :: beta, beta_stop, beta_step, eps, mean, sigma
    real, allocatable :: sp(:)
    integer :: i, st = 0, n_measure, n_thermal
    logical :: output
    real, dimension(2) :: tarray
    real :: result


    call random_seed()
    call u%initialize(16,16,16,16,.false.)


    ! file output
    open(20, file='plaquette_test.txt', iostat=st)
    if (st /= 0) error stop 'Failed to open file'

    ! initial value for beta and epsilon
    beta = 0.2
    beta_stop = 0.6
    beta_step = 0.2
    eps = 0.3
    n_measure = 100
    n_thermal = 100
    output = .true.

    ! execution time
    call ETIME(tarray, result)
    print *, 'Elapsed time: ', result
    print *, 'User time: ', tarray(1)
    print *, 'System time: ', tarray(2)

    allocate(sp(n_measure))
    print *, 'Start run'

    !! Perform beta scan
    do
        print *, 'Beta = ', beta

        ! Thermalize
        print *, 'Thermalize...'
        do i = 1, n_thermal
            call metro(u, eps, beta, acc)
            !if (this_image() == 1) print *, 'Iteration: ', i, 'Acc: ', acc
        end do

        ! Measure
        print *, 'Measure...'
        do i = 1, n_measure
            call metro(u, eps, beta, acc)
            call metro(u, eps, beta, acc)
            sp(i) = plaquette(u)
            !if (this_image() == 1) print *, 'Iteration: ', i, 'Acc: ', acc, ', s_wilson: ', sp(i)
        end do
        print *, 'Acceptance: ', acc

        ! Calc. mean and stat. error, write to file
        sync all
        if (this_image() == 1) then
            print *, 'Calculate mean and sigma...'
            mean = sum(sp)/size(sp)
            sigma = sqrt(sum( (sp - mean)**2 )/size(sp))
            write(20, *) beta, mean, sigma

            if (output) then
                print *, 'Mean: ', mean, ' sigma: ', sigma
                print *, 'Orthogonal: ', dot_product(data(1,:,1), data(2,:,1)), &
                    dot_product(data(1,:,1), data(3,:,1)), &
                    dot_product(data(2,:,1), data(3,:,1))
                print *, 'Normal: ', norm2(abs(data(1,:,1))), norm2(abs(data(2,:,1))), &
                    norm2(abs(data(3,:,1)))
            end if

        end if

        call unitarize_field()

        ! increase beta
        beta = beta + beta_step

        ! break loop when beta_stop is reached
        if (beta > beta_stop) exit
    end do

    ! execution time
    call ETIME(tarray, result)
    print *, 'Elapsed time: ', result
    print *, 'User time: ', tarray(1)
    print *, 'System time: ', tarray(2)

    print *, 'Done'

    close(20)

end program main
