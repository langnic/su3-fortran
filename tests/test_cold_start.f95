!! cold start: all links = 1, so every plaquette is 1
program test_cold_start
    use mod_testing
    use mod_gauge_field
    use mod_measure
    implicit none

    type(gauge_field) :: u

    call u%initialize(4, 4, 4, 4, .false.)

    call check_range('cold start: plaquette = 1', plaquette(u), 1.0 - 1e-6, 1.0 + 1e-6)
    call check_range('cold start: links are SU(3)', max_unitarity_violation(), 0.0, 1e-6)
    call check_range('cold start: Wilson action = 0', s_wilson(u, 5.0), -1e-5, 1e-5)

    call finish()
end program test_cold_start
