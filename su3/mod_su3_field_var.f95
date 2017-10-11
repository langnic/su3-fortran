module mod_su3_field_var
    implicit none
    private

    !! su3-matrix valued field variable
    !! contains a distributed data co-array
    type, public :: su3_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:,:,:)[:]
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary

    contains
        procedure :: add
        procedure :: assign
        generic :: operator (+) => add
        generic :: assignment (=) => assign
    end type su3_field_var


    !! local su3-matrix valued field variable
    !! same as su3_field_var but with a local data store
    type, public :: loc_su3_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:,:,:)
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary
    contains
        procedure :: add_loc
        generic :: operator (+) => add_loc
    end type loc_su3_field_var

contains
    
    !! add routine for su3 field var
    !! returns loc_su3_field_var
    function add(a,b) result (res)
        class(su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)
    end function add

    !! add routine for local su3 field var
    !! returns loc_su3_field_var
    function add_loc(a,b) result (res)
        class(loc_su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)  ! TODO: parallel.
    end function add_loc

    !! assigns a loc_su3_field_var to a su3_field_var
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(loc_su3_field_var), intent(in) :: x
        class(su3_field_var), intent(inout) :: y

        y%data = x%data
        return
    end subroutine

end module mod_su3_field_var
