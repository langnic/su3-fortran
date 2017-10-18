module mod_su3_matrix
    implicit none
    private


    type, public :: su3_matrix
        complex :: m(3,3)
    contains

    end type su3_matrix

contains

    !TODO: add, subtract, multiply, conjugate

end module mod_su3_matrix
