module mod_gauge_field
    use mod_su3_matrix
    use mod_su3_field_var
    use mod_utils
    implicit none

    private
    ! global data variable (for now: only one field instance allowed)
    complex, allocatable :: data(:,:,:)[:]

    !! gauge field type - contains index arrays, dimensions and settings
    !! does not contain any data
    type, public :: gauge_field
        !private
    integer, pointer :: indices(:,:,:,:,:,:)        ! index/co-index array
        integer, pointer :: coordinates(:,:,:)      ! TODO: maybe this is not needed at all...
        logical :: tmp  ! indicates whether the instance is a temporary variable
        integer :: dir  ! 0 means all directions
        integer :: par  ! -1 means undefined parity
        integer :: x_low, x_up, y_low, y_up, z_low, z_up, t_low, t_up   ! boundaries for this image
    contains
        procedure :: get => get_site
        procedure :: copy_shift => copy_shift
        procedure :: initialize
        procedure :: add
        procedure :: multiply
        procedure :: select => select
        procedure :: shift => shift
        generic :: operator(+) => add
        generic :: operator(*) => multiply
    end type gauge_field


contains

    function add(a,b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(su3_matrix) :: res
        !integer :: shp(3)
        !shp = shape(a%data)
        !allocate(tmp(shp(1), shp(2), shp(3))[*], source=a%data+b%data)
    end function add

    function get_site(self, i)
        class(gauge_field), intent(in) :: self
        integer :: i
        complex, dimension(:,:), allocatable :: get_site
    end function


    function multiply(a, b) result (res)
        class(gauge_field), intent(in) :: a, b
        type(loc_su3_field_var) :: res
        integer :: i,j,k,l,d, idx
        integer :: nx, ny, nz, nt, d_low, d_up, a_dir, b_dir

        ! TODO: compare dimensions - error if different
        ! for now: assume all fields have same dimensions
        ! TODO: if x_low, x_up etc. differ calculate deltaX, ... between the two fields
        nx = a%x_up - a%x_low + 1
        ny = a%y_up - a%y_low + 1
        nz = a%z_up - a%z_low + 1
        nt = a%t_up - a%t_low + 1

        ! parity has to be set or unset on both fields
        if (sign(1, a%par) /= sign(1, b%par)) then
            error stop "Multiply: field dimensions don't match"
        end if

        ! direction has to be set or unset on both fields
        if (a%dir == 0 .neqv. b%dir == 0) then
            error stop "Multiply: field dimensions don't match"
        end if

        if (a%par >= 0) then
            if (a%dir /= 0) then
                allocate(res%data(nx*ny*nz*nt/2,3,3))
                d_low = 0
                d_up = 0
            else
                allocate(res%data(nx*ny*nz*nt*2,3,3))
                d_low = 1
                d_up = 4
            end if
        else
            if (a%dir /= 0) then
                allocate(res%data(nx*ny*nz*nt,3,3))
                d_low = 0
                d_up = 0
            else
                allocate(res%data(nx*ny*nz*nt*4,3,3))
                d_low = 1
                d_up = 4
            end if
        end if

        !TODO: put the above in subroutine

        idx = 1
        do d = d_low, d_up
            do i = a%t_low, a%t_up              ! t
                do j = a%z_low, a%z_up          ! z
                    do k = a%y_low, a%y_up                 ! y
                        do l = a%x_low+1, a%x_up, 2        ! x  !FIXME: assumed parity!
                            if (d == 0) then
                                a_dir = abs(a%dir)
                                b_dir = abs(b%dir)
                            else
                                a_dir = d
                                b_dir = d
                            end if
                            res%data(idx,:,:) = multSu3( &
                            data( a%indices(l-abs(mod(i+j+k,2)-a%par),k,j,i,a_dir,0), :,:) &
                            [a%indices(l-abs(mod(i+j+k,2)-a%par),k,j,i,a_dir,1)], &
                            (a%dir < 0), &
                            data( b%indices(l-abs(mod(i+j+k,2)-b%par),k,j,i,b_dir,0), :,:) &
                            [b%indices(l-abs(mod(i+j+k,2)-b%par),k,j,i,b_dir,1)], &
                            (b%dir < 0) )
                            idx = idx + 1

                        end do
                    end do
               end do
            end do
        end do

    end function


    !! copies the gauge_field
    !! deep-copies index arrays depending on whether the source object is a temporary variable
    !! optionally applies a shift while copying
    subroutine copy_shift(self, c, shift_dim, shift_length)
        class(gauge_field), intent(in) :: self
        type(gauge_field), intent(inout) :: c
        integer, dimension(6) :: s
        integer, dimension(3) :: t
        integer, intent(in) :: shift_dim, shift_length
        s = shape(self%indices)
        t = shape(self%coordinates)
        c = self
        ! only deep-copy the index arrays if self is not a temporary
        if (.not. self%tmp) then
            if (shift_dim /= 0) then
                allocate(c%indices(s(1), s(2), s(3), s(4), s(5), s(6)), &
                    source=cshift(self%indices, shift_length, shift_dim))
            else
                allocate(c%indices(s(1), s(2), s(3), s(4), s(5), s(6)), source=self%indices)
                allocate(c%coordinates(t(1), t(2), t(3)), source=self%coordinates)
            end if
        else if (shift_dim /= 0) then
            c%indices = cshift(c%indices, shift_length, shift_dim)
            deallocate(c%coordinates)
        end if
    end subroutine


    !! Select a parity and direction (used by multiply, etc.)
    function select (self, p, d) result (tmp)
        class(gauge_field), intent(in) :: self
        type(gauge_field) :: tmp
        integer, intent(in) :: p, d

        call self%copy_shift(tmp, 0, 0)
        tmp%par = p
        tmp%dir = d
        tmp%tmp = .true.
    end function


    !! applies a shift in direction d by x nodes
    function shift (self, d, x) result (tmp)
        class(gauge_field), intent(in) :: self
        type(gauge_field) :: tmp
        integer, intent(in) :: d, x

        call self%copy_shift(tmp, d, x)
        tmp%tmp = .true.
    end function


    !! initialize the gauge field - set dimensions
    !! nz*nt has to be a multiple of the number of nodes
    !! random indicates hot or cold start
    subroutine initialize(self, nx, ny, nz, nt, random)
        class(gauge_field), intent(inout) :: self
        integer, intent(in) :: nx, ny, nz, nt
        logical, intent(in) :: random
        integer :: alloc, i, j, k, l, lp, m, p, idx, co_idx
        integer :: nz_patches, nt_patches, nzp, ntp

        ! enforce even lattice dimensions
        if ( mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0 .or. mod(nz, 2) /= 0 .or. mod(nt, 2) /= 0) then
            error stop "Latice dimensions should be even"
        end if

        ! for simplicity we only divide along z and t
        ! nz * nt has to be a multiple of the number of nodes
        allocate(self%indices(nx, ny, nz, nt, 4, 0:1))
        if ( mod(real(nz), sqrt(real(num_images()))) == 0.0 .and. mod(real(nt), sqrt(real(num_images()))) == 0.0) then
            nz_patches = sqrt(real(num_images()))
            nt_patches = sqrt(real(num_images()))
        else if (mod(nz, num_images()) == 0) then
            nz_patches = num_images()
            nt_patches = 1
        else if (mod(nt, num_images()) == 0) then
            nt_patches = num_images()
            nz_patches = 1
        else
            error stop "Number of lattice sites Nz*Nt has to be a multiple of the number of nodes"
        end if

        ! allocate data co-array and coordinates array
        nzp = nz / nz_patches   ! nz per patch
        ntp = nt / nt_patches   ! nt per patch
        if (allocated(data)) then
            error stop "For now only one lattice instance is allowed"
        end if
        allocate(data(nx*ny*nzp*ntp*4,3,3)[*], stat=alloc)
        if (alloc /= 0) then
            error stop "Co-array allocation failed"
        end if
        allocate(self%coordinates(nx*ny*nzp*ntp*4,num_images(),5))

        ! set boundaries on this image
        self%x_low = 1
        self%x_up = nx
        self%y_low = 1
        self%y_up = ny
        self%z_low = mod(this_image()-1, nt_patches)*nzp+1
        self%z_up = mod(this_image()-1, nt_patches)*nzp+nzp
        self%t_low = int((this_image()-1)/nz_patches)*ntp+1
        self%t_up = int((this_image()-1)/nz_patches)*ntp+ntp

        ! set partiy and selected direction to undefined
        self%dir = 0
        self%par = -1
        self%tmp = .false.

        if (this_image() == 1) then
            write (*,*) 'Configuration:'
            write (*,*) num_images(), ' nodes'
            write (*,*) nz_patches, ' x ', nt_patches, ' distributed patches'
            write (*,*) nx, 'x', ny, 'x', nzp, 'x', ntp, ' lattice sites per patch'
        end if

        ! fill index array and initialize matrices
        do co_idx = 1, num_images()     ! patches
            idx = 0
            do m = 1, 4                 ! link direction
                do p = 0, 1             ! parity
                    do i = int((co_idx-1)/nz_patches)*ntp+1, int((co_idx-1)/nz_patches)*ntp+ntp      ! t
                        do j = mod(co_idx-1, nt_patches)*nzp+1, mod(co_idx-1, nt_patches)*nzp+nzp    ! z
                            do k = 1, ny                                                             ! y
                                do l = 2-abs(mod(i+j+k, 2)-p), nx-abs(mod(i+j+k, 2)-p), 2            ! x
                                    idx = idx + 1
                                    !write (*,*) l, ', ', k, ', ', j, ', ', i
                                    !write (*,*) 'Idx: ', idx, 'co-index: ', co_idx
                                    self%indices(l,k,j,i,m,:) = (/ idx, co_idx /)
                                    self%coordinates(idx,co_idx,:) = (/ l, k, j, i, m /)
                                    if (co_idx == this_image()) then
                                        if (random) then
                                            data(idx,:,:) = get_random_su3(1.0)
                                        else
                                            data(idx,:,:) = id_su3()
                                        end if
                                    end if
                                end do
                            end do
                       end do
                    end do
                end do
            end do
        end do
        write (*,*) 'Image ', this_image(), ': Field initialized'
        sync all()

    end subroutine



end module mod_gauge_field
