module mod_logical_field_var
    implicit none
    private

    !! complex valued field variable
    !! contains a distributed data co-array
    type, public :: logical_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        logical, allocatable :: data(:)[:]
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary

    contains
        procedure :: assign
        generic :: assignment (=) => assign
    end type logical_field_var


    !! local logical valued field variable
    !! same as logical_field_var but with a local data store
    type, public :: loc_logical_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        logical, allocatable :: data(:)
        integer, allocatable :: indices(:,:,:,:)    ! might be unnecessary
    contains

    end type loc_logical_field_var

contains

    !! assigns a loc_logical_field_var to a logical_field_var
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(loc_logical_field_var), intent(in) :: x
        class(logical_field_var), intent(inout) :: y

        y%data = x%data
        return
    end subroutine



end module mod_logical_field_var
