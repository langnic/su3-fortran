module mod_utils
    implicit none

    contains

        !! generate random su(2) matrix close to unity
        function get_small_su2(eps)
            complex, dimension(2,2) :: get_small_su2
            real, intent(in) :: eps
            real :: r(4), x(4)
            integer :: seed
            ! Pauli matrices
            complex, parameter, dimension(2,2) :: sig0 = reshape( (/ (1,0), (0,0), (0,0), (1,0) /), (/ 2, 2 /) ), &
            sig1 = reshape( (/ (0,0), (1,0), (1,0), (0,0) /), (/ 2, 2 /) ), &
            sig2 = reshape( (/ (0,0), (0,1), (0,-1), (0,0) /), (/ 2, 2 /) ), &
            sig3 = reshape( (/ (1,0), (0,0), (0,0), (-1,0) /), (/ 2, 2 /) )

            ! get random numbers
            !call random_seed(seed)
            call random_number(r)

            r = r - 0.5
            !eps = 0.01
            x(2:4) = eps * r(2:4) / NORM2(r(2:4))
            !x(1) = sign(1.0, r(1)) * sqrt(1 - eps**2)
            x(1) = 1!sqrt(1 - eps**2)

            get_small_su2 = x(1) * sig0 + (0,1) * ( x(2) * sig1 + x(3) * sig2 + x(4) * sig3 )

            !write(*,*) 'su2: ', get_small_su2

        end function get_small_su2

        !! generate random su(3) matrix close to unity
        function get_random_su3(eps)
            complex, dimension(3,3) :: get_random_su3
            complex, dimension(2,2) :: rp, sp, tp
            complex, dimension(3,3) :: r, s, t
            real, intent(in) :: eps

            rp = get_small_su2(eps)
            sp = get_small_su2(eps)
            tp = get_small_su2(eps)

            r(1,1:2) = rp(1,1:2)
            r(2,1:2) = rp(2,1:2)
            r(:,3) = (/ 0, 0, 1 /)
            r(3,:) = (/ 0, 0, 1 /)

            s(1,1:3:2) = sp(1,1:2)
            s(3,1:3:2) = sp(2,1:2)
            s(:,2) = (/ 0, 1, 0 /)
            s(2,:) = (/ 0, 1, 0 /)

            t(2,2:3) = tp(1,1:2)
            t(3,2:3) = tp(2,1:2)
            t(:,1) = (/ 1, 0, 0 /)
            t(1,:) = (/ 1, 0, 0 /)

            !get_random_su3 = matmul(r, matmul(s, t)) ! TODO: replace with blas function
            get_random_su3 = multSu3(r, .false., &
                multSu3(s, .false., t, .false.), .false.)

            call unitarize(get_random_su3)

            !write(*,*) ' ', get_small_su3(1,:)
            !write(*,*) ' ', get_small_su3(2,:)
            !write(*,*) ' ', get_small_su3(3,:)

        end function get_random_su3

        !! vector product
        pure function cross(a, b)
          complex, dimension(3) :: cross
          complex, dimension(3), INTENT(IN) :: a, b

          cross(1) = a(2) * b(3) - a(3) * b(2)
          cross(2) = a(3) * b(1) - a(1) * b(3)
          cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross


        !! project complex square matrix to unitarity
        subroutine unitarize(m)
            real :: scnrm2
            complex, dimension(3,3), intent(inout) :: m
            complex, dimension(3) :: v
            complex :: cdotc

            ! normalize first row
            m(1,1:3) = m(1,1:3) / scnrm2(3,m(1,1:3),1)

            ! normalize / orthogonalize second row
            v = m(2,1:3) - m(1,1:3) * conjg(cdotc(3, m(2,1:3), 1, m(1,1:3), 1)) ! attention: cdotc conj. the FIRST vect.
            m(2,1:3) = v / scnrm2(3,v,1)

            ! construct third row from first two
            m(3,1:3) = cross(conjg(m(1,1:3)), conjg(m(2,1:3)))

            ! check scalar product 0
            !write(*,*) 'check dot prod. 1st * 2nd: ', cdotc(3, m(2,1:3), 1, m(1,1:3), 1)
            !write(*,*) 'check dot prod. 1st * 3rd: ', cdotc(3, m(3,1:3), 1, m(1,1:3), 1)
            !write(*,*) 'check dot prod. 2nd * 3rd: ', cdotc(3, m(3,1:3), 1, m(2,1:3), 1)
            ! check norm
            !write(*,*) 'check nrm 1st row: ', scnrm2(3,m(1,1:3),1)
            !write(*,*) 'check nrm 2nd row: ', scnrm2(3,m(2,1:3),1)
            !write(*,*) 'check nrm 3rd row: ', scnrm2(3,m(3,1:3),1)
            ! final matrix
            !write(*,*) ' ', m(1,:)
            !write(*,*) ' ', m(2,:)
            !write(*,*) ' ', m(3,:)

        end subroutine unitarize

        !! multiply complex 3x3 matrix using blas
        function multSu3(a, conjA, b, conjB)

            complex, intent(in), dimension(3,3) :: a, b
            complex, dimension(3,3) :: multSu3
            logical, intent(in) :: conjA, conjB

            if (conjA) then
                if (conjB) then
                    call cgemm('c', 'c', 3, 3, 3, (1.0, 0.0), a, 3, b, 3, (0, 0), multSu3, 3)
                else
                    call cgemm('c', 'n', 3, 3, 3, (1.0, 0.0), a, 3, b, 3, (0, 0), multSu3, 3)
                endif
            else
                if (conjB) then
                    call cgemm('n', 'c', 3, 3, 3, (1.0, 0.0), a, 3, b, 3, (0, 0), multSu3, 3)
                else
                    call cgemm('n', 'n', 3, 3, 3, (1.0, 0.0), a, 3, b, 3, (0, 0), multSu3, 3)
                endif
            endif


        end function multSu3

        !! calculate trace of square complex matrix
        function trace(U)
            complex :: U(:,:)
            INTEGER, dimension(2) :: dm(2)
            INTEGER :: i
            complex :: trace

            trace = 0

            dm = shape(U)
            if (dm(1) == dm(2)) then
                do i = 1, dm(1)
                    trace = trace + U(i,i)
                end do
            endif
        end function trace

        !! get the identity matrix
        function id_su3()
            complex, dimension(3,3) :: id_su3
            id_su3 = reshape( (/ 1, 0, 0,  0, 1, 0,  0, 0, 1 /), (/ 3, 3 /) )
        end function id_su3

end module mod_utils
