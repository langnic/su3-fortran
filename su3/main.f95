program main

    use mod_su3_field_var
    use mod_gauge_field
    implicit none


    type(gauge_field) :: u, ut
    type(su3_field_var) :: v
    !write(*,*) num_images()
    !write(*,*) this_image()

    call u%initialize(4,4,4,4,.false.)

    !write(*,*) 'Elem: ', u%get(4,4,4,4,4)
    !write(*,*) 'Elem: ', minval(u%indices)

    ut = u%shift(1,1)
    !v = u%shift(1,1) * u
    if (this_image() == 1) then
        write(*,*) 'Elem: ', v%data(1,:,:)[4]
    end if

end program main
