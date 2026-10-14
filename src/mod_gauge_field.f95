module mod_gauge_field
    use mod_su3_matrix
    use mod_su3_field_var
    use mod_logical_field_var
    use mod_utils
    implicit none

    !private
    ! global data variable (for now: only one field instance allowed)
    ! (has to be global because co-array pointers are not supported)
    complex, allocatable :: data(:,:,:)[:]

    !! gauge field type - contains index arrays, dimensions and settings
    !! does not contain any data
    type, public :: gauge_field
        !private
        integer, allocatable :: indices(:,:,:,:,:,:)        ! index/co-index array
        !integer, pointer :: coordinates(:,:,:)      ! TODO: maybe this is not needed at all...
        logical :: tmp  ! indicates whether the instance is a temporary variable
        integer :: dir  ! 0 means all directions
        integer :: par  ! -1 indicates undefined parity
        integer :: x_low, x_up, y_low, y_up, z_low, z_up, t_low, t_up   ! boundaries for this image
        integer :: nx, ny, nz, nt   ! dimensions on this image
    contains
        procedure :: get => get_element
        procedure :: copy_shift => copy_shift
        procedure :: initialize
        procedure :: update
        procedure :: add
        procedure :: add_fv
        procedure :: add_loc_fv
        procedure :: subtract
        procedure :: subtract_fv
        procedure :: subtract_loc_fv
        procedure :: multiply
        procedure :: multiply_fv
        procedure :: multiply_loc_fv
        procedure :: select
        procedure :: shift
        generic :: operator(+) => add, add_fv, add_loc_fv
        generic :: operator(-) => subtract, subtract_fv, subtract_loc_fv
        generic :: operator(*) => multiply, multiply_fv, multiply_loc_fv
        generic :: operator(.shift.) => shift
        generic :: operator(.sel.) => select
    end type gauge_field

contains

    function get_element(self, i, j, k, l, d)
        class(gauge_field), intent(in) :: self
        integer :: i, j, k, l, d
        complex, dimension(3,3) :: get_element
        get_element = data(:,:,self%indices(0,l,k,j,i,d))[self%indices(1,l,k,j,i,d)]
    end function


    !! return range of directions and step size for iteration
    !! allocate memory of result variable accordingly
    subroutine comp_dimensions(a, res, d_low, d_up, step_size)
        class(gauge_field), intent(in) :: a
        type(loc_su3_field_var), intent(inout) :: res
        integer, intent(inout) :: d_low, d_up, step_size
        integer :: nx, ny, nz, nt

        ! TODO: compare dimensions - error if different
        ! for now: assume all fields have same dimensions
        ! TODO: if x_low, x_up etc. differ calculate deltaX, ... between the two fields
        nx = a%x_up - a%x_low + 1
        ny = a%y_up - a%y_low + 1
        nz = a%z_up - a%z_low + 1
        nt = a%t_up - a%t_low + 1

        if (a%par >= 0) then        ! defined parity
            if (a%dir /= 0) then    ! defined direction
                if (.not. allocated(res%data)) allocate(res%data(3,3,nx*ny*nz*nt/2))
                d_low = 0
                d_up = 0
            else                    ! undefined direction -> iterate over all directions
                if (.not. allocated(res%data)) allocate(res%data(3,3,nx*ny*nz*nt*2))
                d_low = 1
                d_up = 4
            end if
            step_size = 2
        else                        ! undefined parity -> iterate over even AND odd sites
            if (a%dir /= 0) then    ! defined direction
                if (.not. allocated(res%data)) allocate(res%data(3,3,nx*ny*nz*nt))
                d_low = 0
                d_up = 0
            else                    ! undefined direction -> iterate over all directions
                if (.not. allocated(res%data)) allocate(res%data(3,3,nx*ny*nz*nt*4))
                d_low = 1
                d_up = 4
            end if
            step_size = 1
        end if

    end subroutine


    !! compare parity and selected directions for two gauge_fields
    !! cartesian dimensions of gauge fields assumed to be equal
    subroutine check_dim(a, b)
        class(gauge_field), intent(in) :: a, b
        ! parity has to be set or unset on both fields
        if (sign(1, a%par) /= sign(1, b%par)) then
            error stop "Multiply: field dimensions don't match"
        end if

        ! direction has to be set or unset on both fields
        if (a%dir == 0 .neqv. b%dir == 0) then
            error stop "Multiply: field dimensions don't match"
        end if
    end subroutine

    !! compare dimensions and selected directions for a gauge_fields and a local field var
    subroutine check_dim_fv(a, sb)
        class(gauge_field), intent(in) :: a
        integer :: sb

        !print *, 'size data: ', size(data)/72, ', size fv: ', sb

        if (a%par >= 0) then        ! defined parity
            if (a%dir /= 0) then    ! defined direction
                if (size(data)/72 /= sb) error stop "Multiply: dimensions don't match"
            else                    ! undefined direction
                if (size(data)/18 /= sb) error stop "Multiply: dimensions don't match"
            end if
        else                        ! undefined parity
            if (a%dir /= 0) then    ! defined direction
                if (size(data)/36 /= sb) error stop "Multiply: dimensions don't match"
            else                    ! undefined direction
                if (size(data)/9 /= sb) error stop "Multiply: dimensions don't match"
            end if
        end if
    end subroutine

    !! add two gauge fields site by site
    !! return loc_su3_field_var
    function add(a, b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, db, la, lb, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim(a, b)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l) - seems to slow the program down rather than speed it up
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                                db = abs(b%dir)
                            else                ! all directions
                                da = d
                                db = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                                lb = l-abs(mod(i+j+k,2)-b%par)
                            else                        ! undefined parity
                                la = l
                                lb = l
                            end if

                            if (a%dir >= 0 .and. b%dir >= 0) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                + data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]
                            else if (a%dir >= 0 .and. b%dir < 0) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                + conjg(transpose(data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]))
                            else if (a%dir < 0 .and. b%dir >= 0) then
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                + data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]
                            else
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                + conjg(transpose(data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]))
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
    function add_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
        type(su3_field_var), intent(in) :: b
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

    !! add gauge_field and loc_su3_field_var
    !! return loc_su3_field_var
    function add_loc_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
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


    !! subtract two gauge fields site by site
    !! return loc_su3_field_var
    function subtract(a, b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, db, la, lb, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim(a, b)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l) - seems to slow the program down rather than speed it up
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                                db = abs(b%dir)
                            else                ! all directions
                                da = d
                                db = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                                lb = l-abs(mod(i+j+k,2)-b%par)
                            else                        ! undefined parity
                                la = l
                                lb = l
                            end if

                            if (a%dir >= 0 .and. b%dir >= 0) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                - data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]
                            else if (a%dir >= 0 .and. b%dir < 0) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                - conjg(transpose(data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]))
                            else if (a%dir < 0 .and. b%dir >= 0) then
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]
                            else
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - conjg(transpose(data( :,:,b%indices(0,lb,k,j,i,db))[b%indices(1,lb,k,j,i,db)]))
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
    function subtract_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
        type(su3_field_var), intent(in) :: b
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
                                - b%data(:,:,idx)
                            else if (a%dir >= 0 .and. b%conj) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                - conjg(transpose(b%data(:,:,idx)))
                            else if (a%dir < 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - b%data(:,:,idx)
                            else
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - conjg(transpose(b%data(:,:,idx)))
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
    function subtract_loc_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
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
                                - b%data(:,:,idx)
                            else if (a%dir >= 0 .and. b%conj) then
                                res%data(:,:,idx) = data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                - conjg(transpose(b%data(:,:,idx)))
                            else if (a%dir < 0 .and. .not. b%conj) then
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - b%data(:,:,idx)
                            else
                                res%data(:,:,idx) = conjg(transpose( &
                                data( :,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)])) &
                                - conjg(transpose(b%data(:,:,idx)))
                            end if

                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! multiply two gauge fields site by site
    !! return loc_su3_field_var
    function multiply(a, b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, db, la, lb, step_size

        call comp_dimensions(a, res, d_low, d_up, step_size)
        call check_dim(a, b)

        idx = 1
        do d = d_low, d_up              ! directions
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        !!$omp parallel do private(l) - seems to slow the program down rather than speed it up
                        do l = a%x_low-1+step_size, a%x_up, step_size        ! x
                            if (d == 0) then    ! specific direction
                                da = abs(a%dir)
                                db = abs(b%dir)
                            else                ! all directions
                                da = d
                                db = d
                            end if
                            if (step_size == 2) then    ! defined parity
                                la = l-abs(mod(i+j+k,2)-a%par)
                                lb = l-abs(mod(i+j+k,2)-b%par)
                            else                        ! undefined parity
                                la = l
                                lb = l
                            end if

                            res%data(:,:,idx) = multSu3( &
                            data( :,:,a%indices(0,la,k,j,i,da)) &
                            [a%indices(1,la,k,j,i,da)], &
                            (a%dir < 0), &
                            data( :,:,b%indices(0,lb,k,j,i,db)) &
                            [b%indices(1,lb,k,j,i,db)], &
                            (b%dir < 0) )
                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! multiply gauge_field and su3_field_var
    !! return loc_su3_field_var
    function multiply_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
        type(su3_field_var), intent(in) :: b
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
                            data( :,:,a%indices(0,la,k,j,i,da)) &
                            [a%indices(1,la,k,j,i,da)], &
                            (a%dir < 0), &
                            b%data(:,:,idx), &
                            b%conj )
                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! multiply gauge_field and loc_su3_field_var
    !! return loc_su3_field_var
    function multiply_loc_fv(a, b) result (res)
        class(gauge_field), intent(in) :: a
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

                            res%data(:,:,idx) = multSu3( &
                            data( :,:,a%indices(0,la,k,j,i,da)) &
                            [a%indices(1,la,k,j,i,da)], &
                            (a%dir < 0), &
                            b%data(:,:,idx), &
                            b%conj )
                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

    end function


    !! copies the gauge_field
    !! deep-copies index arrays depending on whether the source object is a temporary variable
    !! optionally applies a shift while copying
    subroutine copy_shift(self, c, shift_dim, shift_length)
        class(gauge_field), intent(in) :: self
        type(gauge_field), intent(inout) :: c
        integer, dimension(6) :: s
        integer, dimension(3) :: t
        integer, intent(in) :: shift_dim, shift_length
        s = shape(self%indices)
        !t = shape(self%coordinates)
        c = self
        if (shift_dim /= 0) then
            c%indices = cshift(c%indices, shift_length, shift_dim+1)
        end if


        ! only deep-copy the index arrays if self is not a temporary
        !if (.not. self%tmp) then
        !    if (shift_dim /= 0) then
        !        allocate(c%indices(0:1, s(2), s(3), s(4), s(5), s(6)), &
        !            source=cshift(self%indices, shift_length, shift_dim+1))
        !    else
        !        allocate(c%indices(0:1, s(2), s(3), s(4), s(5), s(6)), source=self%indices)
        !        !allocate(c%coordinates(t(1), t(2), t(3)), source=self%coordinates)
        !    end if
        !else if (shift_dim /= 0) then
        !    c%indices = cshift(c%indices, shift_length, shift_dim+1)
        !    !if (associated(c%coordinates)) deallocate(c%coordinates)
        !end if
    end subroutine


    !! Select a parity and direction (used by multiply, etc.)
    function select (a, b) result (tmp)
        class(gauge_field), intent(in) :: a
        integer, dimension(2), intent(in) :: b
        type(gauge_field) :: tmp

        call a%copy_shift(tmp, 0, 0)
        tmp%par = b(1)
        tmp%dir = b(2)
        tmp%tmp = .true.
    end function


    !! applies a shift in direction b(1) by b(2) nodes
    function shift (a, b) result (tmp)
        class(gauge_field), intent(in) :: a
        integer, dimension(2), intent(in) :: b  !(1): direction, (2): distance
        type(gauge_field) :: tmp

        call a%copy_shift(tmp, b(1), b(2))
        tmp%tmp = .true.
        !print *, 'par: ', tmp%par, ', dir: ', tmp%dir
    end function


    !! updates the lattice sites with the value of fv if lfv is true
    subroutine update(a, pin, din, fv, lfv, acc)
        class(gauge_field), intent(inout) :: a
        integer, intent(in) :: pin, din
        type(loc_su3_field_var), intent(inout) :: fv
        type(loc_logical_field_var), intent(in) :: lfv
        integer, intent(inout) :: acc
        integer :: i,j,k,l,d, idx
        integer :: d_low, d_up, da, la, step_size

        a%par = pin
        a%dir = din

        call check_dim_fv(a, size(fv%data)/9)
        call check_dim_fv(a, size(lfv%data))
        call comp_dimensions(a, fv, d_low, d_up, step_size)

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

                            if (lfv%data(idx)) then
                                if (fv%conj) then
                                    data(:,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                    = conjg(transpose(fv%data(:,:,idx)))
                                else
                                    data(:,:,a%indices(0,la,k,j,i,da))[a%indices(1,la,k,j,i,da)] &
                                    = fv%data(:,:,idx)
                                end if
                                acc = acc + 1
                            end if
                            idx = idx + 1

                        end do
                        !!$omp end parallel do
                    end do
               end do
            end do
        end do

        ! reset parity and direction again
        a%par = -1
        a%dir = 0

    end subroutine

    !! project matrices back to unitarity
    subroutine unitarize_field()
        integer :: i

        do i = 1, size(data)/9
            call unitarize(data(:,:,i))
        end do
    end subroutine


    !! split the z-t plane into nz_patches x nt_patches patches, one per image.
    !! Any factorization n_images = nz_patches * nt_patches with nz_patches | nz and
    !! nt_patches | nt is allowed; the one with the shortest patch boundary
    !! (smallest nz/nz_patches + nt/nt_patches) is chosen to minimize communication.
    subroutine decompose(nz, nt, n_images, nz_patches, nt_patches, found)
        integer, intent(in) :: nz, nt, n_images
        integer, intent(out) :: nz_patches, nt_patches
        logical, intent(out) :: found
        integer :: d, best_cost

        found = .false.
        nz_patches = 0
        nt_patches = 0
        best_cost = huge(best_cost)
        do d = 1, n_images
            if (mod(n_images, d) /= 0) cycle
            if (mod(nz, d) /= 0 .or. mod(nt, n_images/d) /= 0) cycle
            if (nz/d + nt/(n_images/d) < best_cost) then
                best_cost = nz/d + nt/(n_images/d)
                nz_patches = d
                nt_patches = n_images/d
                found = .true.
            end if
        end do
    end subroutine decompose


    !! position (0-based) of the patch owned by an image in the nz_patches x nt_patches grid
    !! images are numbered with the z patch running fastest
    pure subroutine patch_position(image, nz_patches, zi, ti)
        integer, intent(in) :: image, nz_patches
        integer, intent(out) :: zi, ti
        zi = mod(image - 1, nz_patches)
        ti = (image - 1) / nz_patches
    end subroutine patch_position


    !! initialize the gauge field - set dimensions
    !! the number of images has to factor into patches along z and t (see decompose)
    !! random indicates hot or cold start
    subroutine initialize(self, nx, ny, nz, nt, random)
        class(gauge_field), intent(inout) :: self
        integer, intent(in) :: nx, ny, nz, nt
        logical, intent(in) :: random
        integer :: alloc, i, j, k, l, lp, m, p, idx, co_idx
        integer :: nz_patches, nt_patches, nzp, ntp, zi, ti
        logical :: found

        ! enforce even lattice dimensions
        if ( mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0 .or. mod(nz, 2) /= 0 .or. mod(nt, 2) /= 0) then
            error stop "Latice dimensions should be even"
        end if

        ! for simplicity we only divide along z and t
        call decompose(nz, nt, num_images(), nz_patches, nt_patches, found)
        if (.not. found) then
            if (this_image() == 1) then
                write (*,*) 'Cannot distribute a lattice with nz =', nz, ', nt =', nt, ' over', num_images(), ' images:'
                write (*,*) 'the image count must factor into nz_patches * nt_patches with nz_patches | nz and nt_patches | nt'
            end if
            error stop "Lattice cannot be distributed over this number of images"
        end if
        allocate(self%indices(0:1, nx, ny, nz, nt, 4))

        ! allocate data co-array
        nzp = nz / nz_patches   ! nz per patch
        ntp = nt / nt_patches   ! nt per patch
        if (allocated(data)) then
            error stop "For now only one lattice instance is allowed"
        end if
        allocate(data(3,3,nx*ny*nzp*ntp*4)[*], stat=alloc)
        if (alloc /= 0) then
            error stop "Co-array allocation failed"
        end if
        !allocate(self%coordinates(nx*ny*nzp*ntp*4,num_images(),5))

        ! set boundaries on this image
        self%x_low = 1
        self%x_up = nx
        self%y_low = 1
        self%y_up = ny
        call patch_position(this_image(), nz_patches, zi, ti)
        self%z_low = zi*nzp+1
        self%z_up = zi*nzp+nzp
        self%t_low = ti*ntp+1
        self%t_up = ti*ntp+ntp

        ! (local) dimensions - for convenience
        self%nx = self%x_up - self%x_low + 1
        self%ny = self%y_up - self%y_low + 1
        self%nz = self%z_up - self%z_low + 1
        self%nt = self%t_up - self%t_low + 1

        ! set partiy and selected direction to undefined
        self%dir = 0
        self%par = -1
        self%tmp = .false.

        if (this_image() == 1) then
            write (*,*) 'Configuration:'
            write (*,*) num_images(), ' nodes'
            write (*,*) nz_patches, ' x ', nt_patches, ' distributed patches'
            write (*,*) nx, 'x', ny, 'x', nzp, 'x', ntp, ' lattice sites per patch'
        end if

        ! fill index array and initialize matrices
        do co_idx = 1, num_images()     ! patches
            call patch_position(co_idx, nz_patches, zi, ti)
            idx = 0
            do m = 1, 4                 ! link direction
                do p = 0, 1             ! parity
                    do i = ti*ntp+1, ti*ntp+ntp          ! t
                        do j = zi*nzp+1, zi*nzp+nzp          ! z
                            do k = 1, ny                                                             ! y
                                do l = 2-abs(mod(i+j+k, 2)-p), nx-abs(mod(i+j+k, 2)-p), 2            ! x
                                    idx = idx + 1
                                    !write (*,*) l, ', ', k, ', ', j, ', ', i
                                    !write (*,*) 'Idx: ', idx, 'co-index: ', co_idx
                                    self%indices(:,l,k,j,i,m) = (/ idx, co_idx /)
                                    !self%coordinates(idx,co_idx,:) = (/ l, k, j, i, m /)
                                    if (co_idx == this_image()) then
                                        if (random) then
                                            data(:,:,idx) = get_random_su3(1.0)
                                        else
                                            data(:,:,idx) = id_su3()
                                        end if
                                    end if
                                end do
                            end do
                       end do
                    end do
                end do
            end do
        end do
        write (*,*) 'Image ', this_image(), ': Field initialized'
        sync all()

    end subroutine



end module mod_gauge_field
