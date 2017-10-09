module mod_gauge_field
    use mod_su3_matrix
    implicit none
    private
    complex, allocatable :: tmp(:,:,:)[:]


    !interface operator(+)
    !    module procedure add
    !end interface

    type, public :: gauge_field
        !private
        complex, allocatable :: data(:,:,:)[:]
        integer, allocatable :: indices(:,:,:,:)
    contains
        procedure :: get => get_site
        procedure :: initialize
        procedure :: add
        final :: destructor
        generic :: operator(+) => add
    end type gauge_field


contains

    function add(a,b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(su3_matrix) :: res
        integer :: shp(3)
        shp = shape(a%data)
        allocate(tmp(shp(1), shp(2), shp(3))[*], source=a%data+b%data)
    end function add

    function get_site(self, i)
        class(gauge_field), intent(in) :: self
        integer :: i
        complex, dimension(:,:), allocatable :: get_site
    end function

    subroutine destructor(self)
        type(gauge_field), intent(in) :: self
    end subroutine

    subroutine initialize(self, nx, ny, nz, nt)
        class(gauge_field), intent(inout) :: self
        integer, intent(in) :: nx, ny, nz, nt
        integer :: alloc
        allocate(self%indices(nx, ny, nz, nt))
        allocate(self%data(nx*ny*nz*nt,3,3)[*], stat=alloc)
        write(*,*) 'Alloc: ', alloc
    end subroutine



end module mod_gauge_field
