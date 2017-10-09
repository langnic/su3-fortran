module mod_gauge_field_component
    implicit none
    private

    !! gauge field component reflecting links of definite parity and orientation
    !! contains a distributed data co-array
    type, public :: gauge_field_component
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:)[:]      ! change back to data(:,:,:)[:]
        integer, allocatable :: indices(:,:,:,:)

    contains
        procedure :: add
        procedure :: assign
        generic :: operator (+) => add
        generic :: assignment (=) => assign
    end type gauge_field_component


    !! auxiliary gauge field component - same as gauge_field_component but with a local data store
    type, public :: aux_gauge_field_component
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:)     ! change back to data(:,:,:)
        integer, allocatable :: indices(:,:,:,:)
    contains
        procedure :: add_aux
        generic :: operator (+) => add_aux
    end type aux_gauge_field_component

contains
    

    !! add routine for gauge field components
    !! returns aux_gauge_field_component
    function add(a,b) result (res)
        class(gauge_field_component), intent(in) :: a, b
        type(aux_gauge_field_component) :: res
        allocate(res%data(size(a%data)), source=a%data+b%data)
    end function add

    !! add routine for auxiliary gauge field components
    !! returns aux_gauge_field_component
    function add_aux(a,b) result (res)
        class(aux_gauge_field_component), intent(in) :: a, b
        type(aux_gauge_field_component) :: res
        allocate(res%data(size(a%data)), source=a%data+b%data)
    end function add_aux

    !! assigns an aux_gauge_field_component to a gauge_field_component
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(aux_gauge_field_component), intent(in) :: x
        class(gauge_field_component), intent(inout) :: y

        y%data = x%data
        return
    end subroutine

end module mod_gauge_field_component
