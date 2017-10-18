module mod_real_field_var
    use mod_logical_field_var
    implicit none
    private

    !! complex valued field variable
    !! contains a distributed data co-array
    type, public :: real_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        real, allocatable :: data(:)[:]
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary

    contains
        procedure :: add
        procedure :: assign
        procedure :: neg
        procedure :: init_rnd
        generic :: operator (+) => add
        generic :: operator (-) => neg
        generic :: assignment (=) => assign
    end type real_field_var


    !! local real valued field variable
    !! same as real_field_var but with a local data store
    type, public :: loc_real_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        real, allocatable :: data(:)
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary
    contains
        procedure :: add_loc
        procedure :: subtract_loc
        procedure :: mult_loc
        procedure :: mult_loc_scalar
        procedure :: lt_loc
        procedure :: le_loc
        procedure :: eq_loc
        procedure :: ge_loc
        procedure :: gt_loc
        procedure :: neg_loc
        procedure :: init_rnd => init_rnd_loc
        generic :: operator (+) => add_loc
        generic :: operator (-) => subtract_loc, neg_loc
        generic :: operator (*) => mult_loc, mult_loc_scalar
        generic :: operator (<) => lt_loc
        generic :: operator (<=) => le_loc
        generic :: operator (==) => eq_loc
        generic :: operator (>=) => ge_loc
        generic :: operator (>) => gt_loc
    end type loc_real_field_var

contains
    
    !! add routine for real field var
    !! returns loc_real_field_var
    function add(a,b) result (res)
        class(real_field_var), intent(in) :: a, b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data+b%data)
    end function add

    !! add routine for local real field var
    !! returns loc_real_field_var
    function add_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data+b%data)  ! TODO: parallel.
    end function add_loc

    !! subtract routine for local real field var
    !! returns loc_real_field_var
    function subtract_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data-b%data)  ! TODO: parallel.
    end function

    !! multiply routine for local real field var
    !! returns loc_real_field_var
    function mult_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*b%data)  ! TODO: parallel.
    end function

    !! multiply with scalar
    !! returns loc_real_field_var
    function mult_loc_scalar(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a
        real, intent(in) :: b
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*b)  ! TODO: parallel.
    end function

    !! assigns a loc_real_field_var to a real_field_var
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(loc_real_field_var), intent(in) :: x
        class(real_field_var), intent(inout) :: y

        if (allocated(y%data)) then
            deallocate(y%data)
        end if
        allocate(y%data(size(x%data))[*])
        y%data = x%data
        return
    end subroutine

    !! initialize real_field_var with random numbers between 0 and 1
    subroutine init_rnd(self, s)
        class(real_field_var), intent(inout) :: self
        integer, intent(in) :: s

        if (mod(s,num_images()) /= 0) &
            error stop "Init field var: size has to be multiple of number of nodes"
        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(int(s/num_images()))[*])
        !call random_seed()          ! FIXME: random_seed should be called just once in the program
        call random_number(self%data)
    end subroutine

    !! initialize loc_real_field_var with random numbers between 0 and 1
    subroutine init_rnd_loc(self, s)
        integer, intent(in) :: s
        class(loc_real_field_var), intent(inout) :: self

        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(s))
        !call random_seed()          ! FIXME: random_seed should be called just once in the program
        call random_number(self%data)
    end subroutine

    !! compare two scalar fields
    !! returns loc_logical_field_var
    function lt_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_logical_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data<b%data)
    end function

    !! compare two scalar fields
    !! returns loc_logical_field_var
    function le_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_logical_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data<=b%data)
    end function

    !! compare two scalar fields
    !! returns loc_logical_field_var
    function eq_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_logical_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data==b%data)
    end function

    !! compare two scalar fields
    !! returns loc_logical_field_var
    function ge_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_logical_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data>=b%data)
    end function

    !! compare two scalar fields
    !! returns loc_logical_field_var
    function gt_loc(a,b) result (res)
        class(loc_real_field_var), intent(in) :: a, b
        type(loc_logical_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data>b%data)
    end function

    !! negate loc_real_field_var
    !! returns loc_real_field_var
    function neg_loc(a) result (res)
        class(loc_real_field_var), intent(in) :: a
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*(-1.0))
    end function

    !! negate loc_real_field_var
    !! returns loc_real_field_var
    function neg(a) result (res)
        class(real_field_var), intent(in) :: a
        type(loc_real_field_var) :: res
        allocate(res%data(size(a%data)), source=a%data*(-1.0))
    end function

end module mod_real_field_var
