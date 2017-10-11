program main
    !use mod_metropolis
    !use mod_utils
    !use mod_su3_matrix
    use mod_su3_field_var
    use mod_gauge_field
    implicit none


    type(gauge_field) :: u, v
    !write(*,*) num_images()
    !write(*,*) this_image()

    call u%initialize(4,4,4,4,.true.)




end program main
