module mod_interfaces

    implicit none

    interface operator(*)
        module procedure multiply_loc_fv_gf
        module procedure mult_scalar_rf
        module procedure mult_scalar_loc_rf
    end interface

    interface operator(+)
        module procedure add_loc_fv_gf
    end interface

    interface operator(-)
        module procedure subtract_loc_fv_gf
    end interface

    interface exp
        module procedure exp_rf
        module procedure exp_loc_rf
    end interface

contains

    !! multiply loc_su3_field_var and gauge_field
    !! return loc_su3_field_var
    function multiply_loc_fv_gf(b, a) result (res)
        use mod_su3_field_var
        use mod_gauge_field
        implicit none
        type(loc_su3_field_var), intent(in) :: b
        type(gauge_field), intent(in) :: a
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, la, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim_fv(a, size(b%data)/9)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l)
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                            else                ! all directions
                                da = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                            else                        ! undefined parity
                                la = l
                            end if

                            res%data(:,:,idx) = multSu3( &
                            b%data(:,:,idx), &
                            b%conj, &
                            data( :,:,a%indices(0,la,k,j,i,da)) &
                            [a%indices(1,la,k,j,i,da)], &
                            (a%dir < 0))
                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! add gauge_field and loc_su3_field_var
    !! return loc_su3_field_var
    function add_loc_fv_gf(b, a) result (res)
        use mod_gauge_field
        type(gauge_field), intent(in) :: a
        type(loc_su3_field_var), intent(in) :: b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, la, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim_fv(a, size(b%data)/9)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l)
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                            else                ! all directions
                                da = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                            else                        ! undefined parity
                                la = l
                            end if

                            if (a%dir >= 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                + b%data(:,:,idx)
                            else if (a%dir >= 0 .and. b%conj) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                + conjg(transpose(b%data(:,:,idx)))
                            else if (a%dir < 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                + b%data(:,:,idx)
                            else
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                + conjg(transpose(b%data(:,:,idx)))
                            end if

                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! add gauge_field and su3_field_var
    !! return loc_su3_field_var
    function subtract_loc_fv_gf(b, a) result (res)
        use mod_gauge_field
        type(gauge_field), intent(in) :: a
        type(loc_su3_field_var), intent(in) :: b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, la, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim_fv(a, size(b%data)/9)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l)
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                            else                ! all directions
                                da = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                            else                        ! undefined parity
                                la = l
                            end if

                            if (a%dir >= 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = b%data(:,:,idx) &
                                - data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)]
                            else if (a%dir >= 0 .and. b%conj) then
                                res%data(:,:,idx) = conjg(transpose(b%data(:,:,idx))) &
                                - data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)]
                            else if (a%dir < 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = b%data(:,:,idx) &
                                - conjg(transpose(data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)]))
                            else
                                res%data(:,:,idx) = conjg(transpose(b%data(:,:,idx))) &
                                - conjg(transpose(data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)]))
                            end if

                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! multiply loc_real_field_var with scalar
    !! returns loc_real_field_var
    function mult_scalar_loc_rf(b,a) result (res)
        use mod_real_field_var
        type(loc_real_field_var), intent(in) :: a
        real, intent(in) :: b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*b)
    end function

    !! multiply real_field_var with scalar
    !! returns loc_real_field_var
    function mult_scalar_rf(b,a) result (res)
        use mod_real_field_var
        type(real_field_var), intent(in) :: a
        real, intent(in) :: b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*b)
    end function

    !! take exp() of loc_real_field_var
    !! returns loc_real_field_var
    function exp_loc_rf(a) result (res)
        use mod_real_field_var
        type(loc_real_field_var), intent(in) :: a
        type(loc_real_field_var) :: res
        integer :: idx

        allocate(res%data(size(a%data)), source=exp(a%data))
        !do idx = 1, size(a%data)
        !    res%data(idx) = exp(a%data(idx))
        !end do
    end function

    !! take exp() of real_field_var
    !! returns loc_real_field_var
    function exp_rf(a) result (res)
        use mod_real_field_var
        type(real_field_var), intent(in) :: a
        type(loc_real_field_var) :: res
        integer :: idx

        allocate(res%data(size(a%data)), source=exp(a%data))
        !do idx = 1, size(a%data)
        !    res%data(idx) = exp(a%data(idx))
        !end do
    end function

end module mod_interfaces
