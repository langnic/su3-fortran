module mod_real_field_var
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
        generic :: operator (+) => add
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
        generic :: operator (+) => add_loc
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

    !! assigns a loc_real_field_var to a real_field_var
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(loc_real_field_var), intent(in) :: x
        class(real_field_var), intent(inout) :: y

        y%data = x%data
        return
    end subroutine

end module mod_real_field_var
