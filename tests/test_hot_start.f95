!! hot start: random links must be proper SU(3) matrices with a small average plaquette
program test_hot_start
    use mod_testing
    use mod_gauge_field
    use mod_measure
    use mod_utils
    implicit none

    type(gauge_field) :: u

    call seed_rng(1234, this_image())
    call u%initialize(4, 4, 4, 4, .true.)

    call check_range('hot start: links are SU(3)', max_unitarity_violation(), 0.0, 1e-5)
    call check_range('hot start: |plaquette| small', abs(plaquette(u)), 0.0, 0.2)

    call finish()
end program test_hot_start
