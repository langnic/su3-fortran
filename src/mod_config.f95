!! Typed run configuration, read from a TOML file (see config/default.toml)
module mod_config
    use, intrinsic :: iso_fortran_env, only: error_unit
    use mod_toml
    implicit none
    private

    public :: sim_params, load_params, n_beta, beta_at

    type :: sim_params
        ! [lattice]
        integer :: nx, ny, nz, nt
        logical :: hot_start
        ! [algorithm]
        real :: eps
        integer :: reunitarize_interval      ! in sweeps, 0 = never
        ! [scan]
        real :: beta_min, beta_max, beta_step
        ! [sampling]
        integer :: n_thermal                 ! thermalization sweeps per beta
        integer :: n_measure                 ! measurements per beta
        integer :: sweeps_per_measurement
        ! [random]
        integer :: seed                      ! 0 = pick one from the clock
        ! [output]
        character(len=256) :: output_dir
        character(len=256) :: output_name
    end type sim_params

    character(len=*), parameter :: known_keys(*) = [character(len=32) :: &
        'lattice.nx', 'lattice.ny', 'lattice.nz', 'lattice.nt', 'lattice.start', &
        'algorithm.epsilon', 'algorithm.reunitarize_interval', &
        'scan.beta_min', 'scan.beta_max', 'scan.beta_step', &
        'sampling.thermalization_sweeps', 'sampling.measurements', 'sampling.sweeps_per_measurement', &
        'random.seed', &
        'output.dir', 'output.name' ]

contains

    subroutine config_error(msg)
        character(len=*), intent(in) :: msg
        write (error_unit, '(a)') 'Config error: ' // msg
        error stop 1
    end subroutine config_error


    !! read, validate and return the run parameters
    function load_params(path) result(p)
        character(len=*), intent(in) :: path
        type(sim_params) :: p
        type(toml_table) :: t
        character(len=256) :: start
        integer :: i

        t = toml_load(path)

        ! typos in key names should not silently fall back to defaults
        do i = 1, t%n
            if (.not. any(trim(t%keys(i)) == known_keys)) call config_error('unknown key "' // trim(t%keys(i)) // '"')
        end do

        p%nx = t%get_int('lattice.nx', 16)
        p%ny = t%get_int('lattice.ny', 16)
        p%nz = t%get_int('lattice.nz', 16)
        p%nt = t%get_int('lattice.nt', 16)
        start = t%get_string('lattice.start', 'cold')

        p%eps = t%get_real('algorithm.epsilon', 0.3)
        p%reunitarize_interval = t%get_int('algorithm.reunitarize_interval', 10)

        p%beta_min = t%get_real('scan.beta_min')
        p%beta_max = t%get_real('scan.beta_max', p%beta_min)
        p%beta_step = t%get_real('scan.beta_step', 1.0)

        p%n_thermal = t%get_int('sampling.thermalization_sweeps', 100)
        p%n_measure = t%get_int('sampling.measurements', 100)
        p%sweeps_per_measurement = t%get_int('sampling.sweeps_per_measurement', 2)

        p%seed = t%get_int('random.seed', 0)

        p%output_dir = t%get_string('output.dir', 'output')
        p%output_name = t%get_string('output.name', 'plaquette_scan')

        select case (trim(start))
        case ('cold')
            p%hot_start = .false.
        case ('hot')
            p%hot_start = .true.
        case default
            call config_error('lattice.start must be "cold" or "hot", got "' // trim(start) // '"')
        end select

        call validate(p)
    end function load_params


    subroutine validate(p)
        type(sim_params), intent(in) :: p

        if (any([p%nx, p%ny, p%nz, p%nt] < 2)) call config_error('lattice dimensions must be >= 2')
        if (any(mod([p%nx, p%ny, p%nz, p%nt], 2) /= 0)) call config_error('lattice dimensions must be even')
        if (p%eps <= 0) call config_error('algorithm.epsilon must be positive')
        if (p%reunitarize_interval < 0) call config_error('algorithm.reunitarize_interval must be >= 0')
        if (p%beta_step <= 0) call config_error('scan.beta_step must be positive')
        if (p%beta_max < p%beta_min) call config_error('scan.beta_max must be >= scan.beta_min')
        if (p%n_thermal < 0) call config_error('sampling.thermalization_sweeps must be >= 0')
        if (p%n_measure < 1) call config_error('sampling.measurements must be >= 1')
        if (p%sweeps_per_measurement < 1) call config_error('sampling.sweeps_per_measurement must be >= 1')
        if (p%seed < 0) call config_error('random.seed must be >= 0')
        if (len_trim(p%output_dir) == 0 .or. len_trim(p%output_name) == 0) &
            call config_error('output.dir and output.name must not be empty')
    end subroutine validate


    !! number of beta values in the scan (beta_max is included if it lies on the grid)
    function n_beta(p)
        type(sim_params), intent(in) :: p
        integer :: n_beta
        n_beta = floor((p%beta_max - p%beta_min) / p%beta_step + 1.0e-3) + 1
    end function n_beta


    function beta_at(p, i)
        type(sim_params), intent(in) :: p
        integer, intent(in) :: i
        real :: beta_at
        beta_at = p%beta_min + (i - 1) * p%beta_step
    end function beta_at

end module mod_config
