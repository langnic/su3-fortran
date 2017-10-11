program main
    !use mod_metropolis
    !use mod_utils
    !use mod_su3_matrix
    use mod_gauge_field
    use mod_gauge_field_component
    implicit none


    type(gauge_field) :: u, v
    !write(*,*) num_images()
    !write(*,*) this_image()

    call u%initialize(32,32,32,32,.true.)




end program main
