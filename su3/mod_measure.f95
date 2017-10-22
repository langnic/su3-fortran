module mod_measure
    use mod_gauge_field
    use mod_su3_field_var
    use mod_real_field_var
    use mod_interfaces
    implicit none


    contains

    !! calculates the average plaquette of the current lattice
    function plaquette(u)
        type(gauge_field), intent(in) :: u
        real :: plaquette
        type(real_field_var) :: sp
        integer :: i,j

        call sp%init_zero(u%nx * u%ny * u%nz * u%nt * num_images())

        do i = 1, 3
            do j = i+1, 4
                call sp%increment( .rtrc.( &
                    (u.sel.[-1,i]) * (u.shift.[i,1].sel.[-1,j]) &
                    * (u.shift.[j,1].sel.[-1,-i]) * (u.sel.[-1,-j]) ))
            end do
        end do
        sync all

        plaquette = (sp%field_avg()/3.0) / 6.0

    end function

    !! calculates the wilson action
    function s_wilson(u, beta)
        type(gauge_field), intent(in) :: u
        real, intent(in) :: beta
        real :: s_wilson
        s_wilson = beta * (1 - plaquette(u))
    end function

end module mod_measure
