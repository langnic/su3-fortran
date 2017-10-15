module mod_su3_matrix
    implicit none
    private


    type, public :: su3_matrix
        complex :: m(3,3)
    contains
        final :: destructor
        procedure :: su3_assign
        generic :: assignment(=) => su3_assign
    end type su3_matrix

contains
    
    subroutine destructor(self)
        type(su3_matrix), intent(in) :: self
    end subroutine

    subroutine su3_assign(y, x)
        type(su3_matrix), intent(in) :: x
        class(su3_matrix), intent(inout) :: y

        y%m = x%m
        return
    end subroutine

    !TODO: add, subtract, multiply, conjugate

end module mod_su3_matrix
