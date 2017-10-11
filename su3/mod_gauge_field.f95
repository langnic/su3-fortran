module mod_gauge_field
    use mod_su3_matrix
    use mod_utils
    implicit none
    private
    complex, allocatable :: data(:,:,:)[:]

    type, public :: gauge_field
        !private
        integer, allocatable :: indices(:,:,:,:,:,:)
        integer, allocatable :: coordinates(:,:,:)
        integer :: dir
        integer :: par
        integer :: nx, ny, nz, nt
    contains
        procedure :: get => get_site
        procedure :: initialize
        procedure :: add
        generic :: operator(+) => add
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


    !! initialize the gauge field - set dimensions
    !! nx*ny has to be a multiple of the number of nodes
    !! random indicates hot or cold start
    subroutine initialize(self, nx, ny, nz, nt, random)
        class(gauge_field), intent(inout) :: self
        integer, intent(in) :: nx, ny, nz, nt
        logical, intent(in) :: random
        integer :: alloc, i, j, k, l, lp, m, p, idx, co_idx
        integer :: nz_patches, nt_patches, nzp, ntp

        if ( mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0 .or. mod(nz, 2) /= 0 .or. mod(nt, 2) /= 0) then
            error stop "Latice dimensions should be even"
        end if

        ! fill index array and initialize matrices
        ! for simplicity we only divide along x and y
        allocate(self%indices(nx, ny, nz, nt, 4, 0:1))
        if ( mod(nz, int(sqrt(real(num_images())))) == 0 .and. mod(nt, int(sqrt(real(num_images())))) == 0) then
            nz_patches = sqrt(real(num_images()))
            nt_patches = sqrt(real(num_images()))
        else if (mod(nz, num_images()) == 0) then
            nz_patches = num_images()
            nt_patches = 1
        else
            error stop "Number of lattice sites Nz*Nt has to be a multiple of the number of nodes"
        end if

        nzp = nz / nz_patches
        ntp = nt / nt_patches
        if (allocated(data)) then
            error stop "For now only one lattice instance is allowed"
        end if
        allocate(data(nx*ny*nzp*ntp,3,3)[*], stat=alloc)
        if (alloc /= 0) then
            error stop "Co-array allocation failed"
        end if
        allocate(self%coordinates(nx*ny*nzp*ntp,num_images(),5))

        if (this_image() == 1) then
            write (*,*) 'Configuration:'
            write (*,*) nz_patches, ' x ', nt_patches, ' distributed patches'
            write (*,*) nx, 'x', ny, 'x', nzp, 'x', ntp, ' lattice sites per patch'
        end if

        do co_idx = 0, num_images()-1   ! patches
            idx = 0
            do m = 1, 4                 ! link direction
                do p = 0, 1             ! parity
                    do i = int(co_idx/nz_patches)*ntp+1, int(co_idx/nz_patches)*ntp+ntp          ! t
                        do j = mod(co_idx, nt_patches)*nzp+1, mod(co_idx, nt_patches)*nzp+nzp    ! z
                            do k = 1, ny                                                            ! y
                                do l = 2-abs(mod(i+j+k, 2)-p), nx-abs(mod(i+j+k, 2)-p), 2           ! x
                                    idx = idx + 1
                                    !write (*,*) l, ', ', k, ', ', j, ', ', i
                                    !write (*,*) 'Idx: ', idx, 'co-index: ', co_idx+1
                                    self%indices(l,k,j,i,m,:) = (/ idx, co_idx+1 /)
                                    self%coordinates(idx,co_idx+1,:) = (/ l, k, j, i, m /)
                                    if (co_idx+1 == this_image()) then
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

    end subroutine



end module mod_gauge_field
