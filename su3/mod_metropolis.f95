module mod_metropolis
    use mod_utils
    use mod_gauge_field
    use mod_real_field_var
    use mod_interfaces
    implicit none

    type(loc_su3_field_var), allocatable, private :: update_tab(:)

    contains

        !! calculate staple for parity p and link direction d
        function staple(u, p, d)
            type(gauge_field), intent(in) :: u
            integer, intent(in) :: d, p
            integer :: i, s
            type(loc_su3_field_var) :: staple

            s = u%nx * u%ny * u%nz * u%nt / 2

            call staple%init_zero(s)
            do i = 1, 4
                if (i /= d) then
                    call staple%increment( &
                        (u.sel.[p,i]) * (u.shift.[i,1].sel.[p,d]) * (u.shift.[d,1].sel.[p,-i]) &
                        + (u.shift.[i,-1].sel.[p,-i]) * (u.shift.[i,-1].sel.[p,d]) * (u.shift.[i,-1].shift.[d,1].sel.[p,i]) &
                        )
                end if
            end do

        end function staple


        !! allocate and initialize the update tab
        !! st: tab size, su: field_var size, eps: epsilon
        subroutine init_update_tab(st, su, eps)
            integer, intent(in) :: st, su
            real, intent(in) :: eps
            integer :: i
            if (.not. allocated(update_tab)) then
                allocate(update_tab(st))
                do i = 1, st
                    call update_tab(i)%init_rnd(su, eps)
                end do
            end if
        end subroutine


        !! perform metropolis sweep
        subroutine metro(u, eps, beta, acc)
            type(gauge_field), intent(inout) :: u
            real, intent(inout) :: acc
            real, intent(in) :: eps, beta
            type(loc_su3_field_var) :: stp, new_link
            type(loc_real_field_var) :: delta_s, rnd_field
            integer :: p, d, s, acc_count
            real :: rnd1, rnd2
            real :: scnrm2

            s = u%nx * u%ny * u%nz * u%nt / 2
            acc_count = 0

            ! generate link updates
            call init_update_tab(100, s, eps)

            ! iterate over both parities and all directions
            do p = 1,2
                do d = 1,4

                ! TODO: multi-hit

                ! get staple
                stp = staple(u,p,d)

                ! generate new link - multiply with compl. conj. update mtx. with equal probability
                call random_number(rnd1)
                call random_number(rnd2)
                if (rnd1 >= 0.5) then
                    update_tab(ceiling(rnd2*100))%conj = .true.     ! TODO: some kind of operator would be nicer for this
                else
                    update_tab(ceiling(rnd2*100))%conj = .false.
                end if

                new_link = (u.sel.[p,d]) * update_tab(ceiling(rnd2*100))

                !! Tests
                if (this_image() == 1) then
                    print *, 'Rnd2: ', rnd2
                    print *

                    print *, 'Update:'
                    print *, update_tab(ceiling(rnd2*100))%data(1,:,1)
                    print *, update_tab(ceiling(rnd2*100))%data(2,:,1)
                    print *, update_tab(ceiling(rnd2*100))%data(3,:,1)
                    print *

                    print *, 'New link:'
                    print *, new_link%data(1,:,1)
                    print *, new_link%data(2,:,1)
                    print *, new_link%data(3,:,1)

                    print *
                    print *, 'Orthogonal: ', dot_product(new_link%data(1,:,1), new_link%data(2,:,1)), &
                        dot_product(new_link%data(1,:,1), new_link%data(3,:,1)), &
                        dot_product(new_link%data(2,:,1), new_link%data(3,:,1))
                    print *, 'Normal: ', scnrm2(3,new_link%data(1,:,1),1), scnrm2(3,new_link%data(2,:,1),1), &
                        scnrm2(3,new_link%data(3,:,1),1)
                end if

                ! calculate action change
                delta_s = beta / 3.0 * .rtrc.((new_link - (u.sel.[p,d])) * stp)

                if (this_image() == 1) print *, 'Delta S: ', delta_s%data(1)

                ! accept with some prob. dep. on action
                call rnd_field%init_rnd(s)
                call u%update(p, d, new_link, (rnd_field <= exp(-delta_s)), acc_count)
                sync all

                end do
            end do

            acc = real(acc_count) / real((s*8))

        end subroutine metro

end module mod_metropolis
