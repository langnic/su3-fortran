module mod_su3_field_var
    use mod_su3_matrix
    use mod_real_field_var
    use mod_utils
    implicit none
    private

    !! su3-matrix valued field variable
    !! contains a distributed data co-array
    type, public :: su3_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:,:,:)[:]
        logical :: conj = .false. ! compl. conj.

    contains
        procedure :: add
        procedure :: add_glob_loc
        procedure :: assign
        procedure :: increment
        procedure :: real_trace
        procedure :: init_rnd
        procedure :: init_zero
        generic :: operator (+) => add, add_glob_loc
        generic :: assignment (=) => assign
        generic :: operator (.rtrc.) => real_trace
        !generic :: assignment (.inc.) => increment
    end type su3_field_var


    !! local su3-matrix valued field variable
    !! same as su3_field_var but with a local data store
    type, public :: loc_su3_field_var
        !private    !TODO: make private again eventually - should be hidden from the user
        complex, allocatable :: data(:,:,:)
        logical :: conj = .false. ! compl. conj.
    contains
        procedure :: add_loc
        procedure :: add_loc_glob
        procedure :: subtract_loc
        procedure :: multiply_loc
        procedure :: real_trace => real_trace_loc
        procedure :: init_rnd => init_rnd_loc
        procedure :: init_zero => init_zero_loc
        procedure :: increment => increment_loc
        generic :: operator (+) => add_loc, add_loc_glob
        generic :: operator (-) => subtract_loc
        generic :: operator (*) => multiply_loc
        generic :: operator (.rtrc.) => real_trace
        !generic :: assignment (.inc.) => increment
    end type loc_su3_field_var

contains
    
    !! add routine for su3 field var
    !! returns loc_su3_field_var
    function add(a,b) result (res)
        class(su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)    ! FIXME: conj.transp. not considered!!! (same below)
    end function add

    !! routine for adding loc_su3_field_var to su3_field_var
    !! returns loc_su3_field_var
    function add_glob_loc(a,b) result (res)
        class(su3_field_var), intent(in) :: a
        type(loc_su3_field_var), intent(in) :: b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)
    end function add_glob_loc

    !! routine for adding loc_su3_field_var to su3_field_var
    !! returns loc_su3_field_var
    function add_loc_glob(a,b) result (res)
        class(loc_su3_field_var), intent(in) :: a
        type(su3_field_var), intent(in) :: b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)
    end function add_loc_glob

    !! add routine for local su3 field var
    !! returns loc_su3_field_var
    function add_loc(a,b) result (res)
        class(loc_su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data+b%data)  ! TODO: parallel.
    end function add_loc

    !! subtract routine for local su3 field var
    !! returns loc_su3_field_var
    function subtract_loc(a,b) result (res)
        class(loc_su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: s(3)
        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)), source=a%data-b%data)  ! TODO: parallel.
    end function

    !! assigns a loc_su3_field_var to a su3_field_var
    !! the local data is copied into the local part of the co-array on every image
    subroutine assign(y, x)
        type(loc_su3_field_var), intent(in) :: x
        class(su3_field_var), intent(inout) :: y
        integer :: s(3)
        s = shape(x%data)

        if (allocated(y%data)) then
            deallocate(y%data)
        end if
        allocate(y%data(3,3,s(3))[*])
        y%data = x%data
        return
    end subroutine


    !! multiply two loc_su3_field_vars
    function multiply_loc(a,b) result (res)
        class(loc_su3_field_var), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: idx
        integer :: s(3)

        s = shape(a%data)
        allocate(res%data(s(1),s(2),s(3)))

        if (size(a%data) /= size(b%data)) error stop "Multiply: dimensions don't match"

        do idx = 1, s(3)
            res%data(:,:,idx) = multSu3( &
            a%data(:,:,idx), &
            a%conj, &
            b%data(:,:,idx), &
            b%conj )
        end do

    end function


    !! increment assignment
    subroutine increment(y, x)
        type(loc_su3_field_var), intent(in) :: x
        class(su3_field_var), intent(inout) :: y
        integer :: s(3), t(3)
        s = shape(x%data)
        t = shape(y%data)

        if (.not. allocated(y%data) .or. .not. allocated(x%data)) then
            error stop "Increment: variable has to be allocated"
        end if

        if (s(3) /= t(3)) then
            error stop "Increment: dimensions don't match"
        end if

        y%data = y%data + x%data
        return
    end subroutine

    !! increment assignment (local)
    subroutine increment_loc(y, x)
        type(loc_su3_field_var), intent(in) :: x
        class(loc_su3_field_var), intent(inout) :: y
        integer :: s(3), t(3)
        s = shape(x%data)
        t = shape(y%data)

        if (.not. allocated(y%data) .or. .not. allocated(x%data)) then
            error stop "Increment: variable has to be allocated"
        end if

        if (s(3) /= t(3)) then
            error stop "Increment: dimensions don't match"
        end if

        y%data = y%data + x%data
        return
    end subroutine


    !! get real trace from loc_su3_field_var
    function real_trace_loc(self) result (res)
        class(loc_su3_field_var), intent(in) :: self
        type(loc_real_field_var) :: res
        integer :: idx
        integer :: s(3)

        s = shape(self%data)
        allocate(res%data(s(3)))

        do idx = 1, s(3)
            res%data(idx) = REAL(self%data(1,1,idx) + self%data(2,2,idx) + self%data(3,3,idx))
        end do

    end function

    !! get real trace from su3_field_var
    function real_trace(self) result (res)
        class(su3_field_var), intent(in) :: self
        type(loc_real_field_var) :: res
        integer :: idx
        integer :: s(3)

        s = shape(self%data)
        allocate(res%data(s(3)))

        do idx = 1, s(3)
            res%data(idx) = REAL(self%data(1,1,idx) + self%data(2,2,idx) + self%data(3,3,idx))
        end do

    end function

    !! initialize su3_field_var with random matrices
    !! epsilon controls proximity to unity
    subroutine init_rnd(self, s, eps)
        class(su3_field_var), intent(inout) :: self
        real, intent(in) :: eps
        integer, intent(in) :: s
        integer :: idx

        if (mod(s,num_images()) /= 0) &
            error stop "Init field var: size has to be multiple of number of nodes"
        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(3,3,int(s/num_images()))[*])
        do idx = 1, s
            self%data(:,:,idx) = get_random_su3(eps)
        end do

    end subroutine

    !! initialize loc_su3_field_var with random matrices
    !! epsilon controls proximity to unity
    subroutine init_rnd_loc(self, s, eps)
        class(loc_su3_field_var), intent(inout) :: self
        real, intent(in) :: eps
        integer, intent(in) :: s
        integer :: idx

        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(3,3,s))
        do idx = 1, s
            self%data(:,:,idx) = get_random_su3(eps)
        end do

    end subroutine

    !! initialize su3_field_var with zeros
    subroutine init_zero(self, s)
        class(su3_field_var), intent(inout) :: self
        integer, intent(in) :: s
        integer :: idx

        if (mod(s,num_images()) /= 0) &
            error stop "Init field var: size has to be multiple of number of nodes"
        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(3,3,int(s/num_images()))[*])
        self%data = 0

    end subroutine

    !! initialize loc_su3_field_var with zeros
    subroutine init_zero_loc(self, s)
        class(loc_su3_field_var), intent(inout) :: self
        integer, intent(in) :: s
        integer :: idx

        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(3,3,s))
        self%data = 0

    end subroutine

end module mod_su3_field_var
