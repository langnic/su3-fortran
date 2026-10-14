!! Metropolis on a 4^4 lattice, started cold:
!! links stay SU(3), plaquette follows the known strong/weak coupling behaviour
program test_metropolis
    use mod_testing
    use mod_gauge_field
    use mod_metropolis
    use mod_measure
    use mod_utils
    implicit none

    type(gauge_field) :: u
    real :: acc, p_strong, p_weak
    real, parameter :: beta_strong = 0.6, beta_weak = 6.0

    call seed_rng(2024, this_image())
    call u%initialize(4, 4, 4, 4, .false.)

    ! strong coupling: <P> = beta/18 + O(beta^3)
    call run_beta(beta_strong, p_strong)
    call check_range('beta=0.6: plaquette ~ beta/18', p_strong, 0.8 * beta_strong / 18, 1.25 * beta_strong / 18)
    call check_range('beta=0.6: acceptance in (0,1)', acc, 0.05, 0.999)
    call check_range('beta=0.6: links are SU(3)', max_unitarity_violation(), 0.0, 1e-4)

    ! weak coupling, beta=6: <P> ~ 0.59 on large lattices
    call run_beta(beta_weak, p_weak)
    call check_range('beta=6: plaquette in expected range', p_weak, 0.5, 0.7)
    call check_range('beta=6: acceptance in (0,1)', acc, 0.05, 0.999)
    call check_range('beta=6: links are SU(3)', max_unitarity_violation(), 0.0, 1e-4)

    call check('plaquette increases with beta', p_weak > p_strong)

    call finish()

contains

    !! thermalize, then average the plaquette over 100 measurements
    subroutine run_beta(beta, p_avg)
        real, intent(in) :: beta
        real, intent(out) :: p_avg
        integer :: i

        do i = 1, 100
            call metro(u, 0.3, beta, acc)
            if (mod(i, 10) == 0) call unitarize_field()
        end do
        p_avg = 0
        do i = 1, 100
            call metro(u, 0.3, beta, acc)
            call metro(u, 0.3, beta, acc)
            if (mod(i, 5) == 0) call unitarize_field()
            p_avg = p_avg + plaquette(u)
        end do
        p_avg = p_avg / 100
    end subroutine run_beta

end program test_metropolis
