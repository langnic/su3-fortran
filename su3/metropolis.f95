module metropolis
    use utils
    implicit none

    contains

        function staple(U, d, x, y, z, t)
            COMPLEX, intent(in) :: U(:,:,:,:,:,:,:)
            INTEGER, dimension(7) :: dm(7)
            INTEGER, intent(in) :: d, x, y, z, t
            INTEGER :: Nx, Ny, Nz, Nt, xp, xm, yp, ym, zp, zm, tp, tm;
            complex, dimension(3,3) :: staple

            dm = shape(U)

            Nx = dm(2)
            Ny = dm(3)
            Nz = dm(4)
            Nt = dm(5)

            xp = modulo(x, Nx)+1
            xm = modulo(x-2, Nx)+1
            yp = modulo(y, Ny)+1
            ym = modulo(y-2, Ny)+1
            zp = modulo(z, Nz)+1
            zm = modulo(z-2, Nz)+1
            tp = modulo(t, Nt)+1
            tm = modulo(t-2, Nt)+1

            select case (d)
                case (1)
                    staple = multSu3( U(2, x, y, z, t, :, :), .false., &
                                multSu3( U(1, x, yp, z, t, :, :), .false., U(2, xp, y, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(2, x, ym, z, t, :, :), .true., &
                                multSu3( U(1, x, ym, z, t, :, :), .false., U(2, xp, ym, z, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(3, x, y, z, t, :, :), .false., &
                                multSu3( U(1, x, y, zp, t, :, :), .false., U(3, xp, y, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(3, x, y, zm, t, :, :), .true., &
                                multSu3( U(1, x, y, zm, t, :, :), .false., U(3, xp, y, zm, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, t, :, :), .false., &
                                multSu3( U(1, x, y, z, tp, :, :), .false., U(4, xp, y, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, tm, :, :), .true., &
                                multSu3( U(1, x, y, z, tm, :, :), .false., U(3, xp, y, z, tm, :, :), .false. ), &
                                .false.)
                case (2)
                    staple = multSu3( U(1, x, y, z, t, :, :), .false., &
                                multSu3( U(2, xp, y, z, t, :, :), .false., U(1, x, yp, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(1, xm, y, z, t, :, :), .true., &
                                multSu3( U(2, xm, y, z, t, :, :), .false., U(1, xm, yp, z, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(3, x, y, z, t, :, :), .false., &
                                multSu3( U(2, x, y, zp, t, :, :), .false., U(3, x, yp, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(3, x, y, zm, t, :, :), .true., &
                                multSu3( U(2, x, y, zm, t, :, :), .false., U(3, x, yp, zm, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, t, :, :), .false., &
                                multSu3( U(2, x, y, z, tp, :, :), .false., U(4, x, yp, z, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, tm, :, :), .true., &
                                multSu3( U(2, x, y, z, tm, :, :), .false., U(3, x, yp, z, tm, :, :), .false. ), &
                                .false.)
                case (3)
                    staple = multSu3( U(1, x, y, z, t, :, :), .false., &
                                multSu3( U(3, xp, y, z, t, :, :), .false., U(1, x, y, zp, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(1, xm, y, z, t, :, :), .true., &
                                multSu3( U(3, xm, y, z, t, :, :), .false., U(1, xm, y, zp, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(2, x, y, z, t, :, :), .false., &
                                multSu3( U(3, x, yp, z, t, :, :), .false., U(2, x, y, zp, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(2, x, ym, z, t, :, :), .true., &
                                multSu3( U(3, x, ym, z, t, :, :), .false., U(2, x, ym, zp, t, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, t, :, :), .false., &
                                multSu3( U(3, x, y, z, tp, :, :), .false., U(4, x, y, zp, t, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(4, x, y, z, tm, :, :), .true., &
                                multSu3( U(3, x, y, z, tm, :, :), .false., U(4, x, y, zp, tm, :, :), .false. ), &
                                .false.)
                case (4)
                    staple = multSu3( U(1, x, y, z, t, :, :), .false., &
                                multSu3( U(4, xp, y, z, t, :, :), .false., U(1, x, y, z, tp, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(1, xm, y, z, t, :, :), .true., &
                                multSu3( U(4, xm, y, z, t, :, :), .false., U(1, xm, y, z, tp, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(2, x, y, z, t, :, :), .false., &
                                multSu3( U(4, x, yp, z, t, :, :), .false., U(2, x, y, z, tp, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(2, x, ym, z, t, :, :), .true., &
                                multSu3( U(4, x, ym, z, t, :, :), .false., U(2, x, ym, z, tp, :, :), .false. ), &
                                .false.) &
                           + multSu3( U(3, x, y, z, t, :, :), .false., &
                                multSu3( U(4, x, y, zp, t, :, :), .false., U(3, x, y, z, tp, :, :), .true. ), &
                                .false.) &
                           + multSu3( U(3, x, y, zm, t, :, :), .true., &
                                multSu3( U(4, x, y, zm, t, :, :), .false., U(3, x, y, zm, tp, :, :), .false. ), &
                                .false.)

                case default
                    write(*,*) 'Dimension out of range'
            end select

        end function staple


        subroutine metro(U, acc, beta)

            COMPLEX, intent(inout) :: U(:,:,:,:,:,:,:)
            INTEGER, dimension(5) :: dm(5)
            INTEGER :: N, i, d, x, y, z, t
            REAL, intent(inout) :: acc
            REAL, intent(in) :: beta
            REAL :: rnd, deltaS
            complex, dimension(3,3) :: stp, oldLink, newLink, wloop
            REAL, dimension(size(U(:,:,:,:,:,1,1))) :: rndTab
            INTEGER, dimension(size(U(:,:,:,:,:,1,1))) :: dTab, xTab, yTab, zTab, tTab
            COMPLEX, dimension(1000, 3, 3) :: updateTab

            dm = shape(U(:,:,:,:,:,1,1))
            N = size(U(:,:,:,:,:,1,1))   ! Nx*Ny*Nz*Nt*Nd
            acc = 0

            ! generate random integer tables for link access
            call random_seed()
            call random_number(rndTab)
            dTab = floor( rndTab * dm(1) ) + 1
            call random_number(rndTab)
            xTab = floor( rndTab * dm(2) ) + 1
            call random_number(rndTab)
            yTab = floor( rndTab * dm(3) ) + 1
            call random_number(rndTab)
            zTab = floor( rndTab * dm(4) ) + 1
            call random_number(rndTab)
            tTab = floor( rndTab * dm(5) ) + 1

            ! generate link updates
            do i = 1, size(updateTab(:,1,1))
                updateTab(i,:,:) = get_small_su3(0.001)
            end do

            ! check int tabs
            !write(*,*) 'xTab(1:5): ', xTab(1:5)
            !write(*,*) 'yTab(1:5): ', yTab(1:5)

            do i = 1, N
                d = dTab(i)
                x = xTab(i)
                y = yTab(i)
                z = zTab(i)
                t = tTab(i)

                stp = staple(U, d, x, y, z, t)
                oldLink = U(d,x,y,z,t,:,:)

                ! generate new link - multiply with compl. conj. update mtx. with equal probability
                call random_number(rnd)
                if (rnd >= 0.5) then
                    newLink = multSu3(U(d,x,y,z,t,:,:), .false., updateTab(modulo(i, 1000)+1,:,:), .false.)
                else
                    newLink = multSu3(U(d,x,y,z,t,:,:), .false., updateTab(modulo(i, 1000)+1,:,:), .true.)
                end if

                ! calculate action change
                wloop = multSu3((newLink - oldLink), .false., stp, .false.)
                deltaS = beta / 3 * REAL(wloop(1,1) + wloop(2,2) + wloop(3,3))

                ! accept with some prob. dep. on action
                call random_number(rnd)
                if (rnd <= exp(-deltaS)) then
                    U(d,x,y,z,t,:,:) = newLink
                    acc = acc + 1
                end if

            end do

            acc = acc / N

        end subroutine metro

end module metropolis
