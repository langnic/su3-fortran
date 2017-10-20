program main

    use mod_su3_field_var
    use mod_gauge_field
    use mod_real_field_var
    use mod_interfaces
    use mod_metropolis
    use mod_measure
    implicit none


    type(gauge_field) :: u
    type(su3_field_var) :: v
    type(real_field_var) :: r, er
    !write(*,*) num_images()
    !write(*,*) this_image()
    complex :: res(3,3), res2(3,3)
    real :: acc = 0
    real :: sw
    integer :: i


    call random_seed()
    call u%initialize(16,16,16,16,.false.)


    !! Some tests

    if (this_image() == 1) then
        !data(:,:,u%indices(0,1,1,1,1,1))[u%indices(1,1,1,1,1,1)] = reshape((/ 2, 0, 1, 0, 2, 0, 1, 0, 2 /), (/ 3, 3 /))
        !data(:,:,u%indices(0,2,2,3,3,3))[u%indices(1,2,2,3,3,3)] = reshape((/ 4, 0, 3, 0, 4, 0, 3, 0, 4 /), (/ 3, 3 /))

        !data(:,:,u%indices(0,3,1,1,1,1))[u%indices(1,3,1,1,1,1)] = reshape((/ 2, 0, 3, 0, 2, 0, 1, 0, 2 /), (/ 3, 3 /))
        !data(:,:,u%indices(0,4,2,3,3,3))[u%indices(1,4,2,3,3,3)] = reshape((/ 4, 0, 5, 0, 4, 0, 6, 0, 4 /), (/ 3, 3 /))

        !print *, 'idx1: ', u%indices(0,1,1,1,1,1), ', co-idx1: ', u%indices(1,1,1,1,1,1)
        !print *, 'idx2: ', u%indices(0,2,2,3,3,3), ', co-idx2: ', u%indices(1,2,2,3,3,3)
    end if


    !v = ((u.sel.[0,1]) - (u.sel.[0,2])) * ((u.shift.[1,1].shift.[2,1].shift.[4,2].shift.[3,2].sel.[0,-3]) - (u.sel.[0,2]))
    !r = v%real_trace()
    !er = exp(-r)

    !v = (u.sel.[0,3]) * (u.sel.[1,-4].shift.[1,1]) * (u.sel.[0,2])
    if (this_image() == 1) then
    !    res = v%data(:,:,1)![1]
    !    res2 = v%data(:,:,2)![1]
    !    print *, res(1,:)
    !    print *, res(2,:)
    !    print *, res(3,:)
    !    print *
    !    print *, res2(1,:)
    !    print *, res2(2,:)
    !    print *, res2(3,:)
    !    print *
    !    print *, r%data(1)
    !    print *, r%data(2)
    !    print *
    !    print *, er%data(1)
    !    print *, er%data(2)
    end if

    ! 10 metropolis sweeps

    do i = 1, 10
        call metro(u, 0.01, 6.0, acc)
        print *, 'Iteration: ', i, 'Acc: ', acc
    end do

    sw = (s_wilson(u)-1)*(-1)

    print *, 'Acc: ', acc, ', s_wilson: ', sw

end program main
