module mod_su3_matrix
    implicit none
    private

    !public :: assignment(=)

    !interface assignment(=)
        !module procedure su3_assign
    !    module procedure bogus_assign
    !end interface

    type, public :: su3_matrix
        complex :: m(3,3)
    contains
        final :: destructor
        procedure :: su3_assign2
        procedure :: bogus_assign
        generic :: assignment(=) => su3_assign2, bogus_assign
    end type su3_matrix

contains
    
    subroutine destructor(self)
        type(su3_matrix), intent(in) :: self
    end subroutine

    subroutine su3_assign(y, x)
        type(su3_matrix), intent(in) :: x
        type(su3_matrix), intent(inout) :: y

        y%m = x%m
        return
    end subroutine

    subroutine bogus_assign(y, x)
        integer, intent(in) :: x
        class(su3_matrix), intent(inout) :: y
    end subroutine

    subroutine su3_assign2(y, x)
        type(su3_matrix), intent(in) :: x
        class(su3_matrix), intent(inout) :: y

        y%m = x%m
        return
    end subroutine

end module mod_su3_matrix
