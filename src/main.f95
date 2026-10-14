!! SU(3) lattice gauge theory: Metropolis beta scan
!! usage: su3 [config.toml]      (default: config/default.toml)
!! Writes the raw plaquette measurements (no averaging) to <output.dir>/<output.name>.csv
program main
    use mod_gauge_field
    use mod_utils
    use mod_config
    use mod_metropolis
    use mod_measure
    implicit none

    type(gauge_field) :: u
    type(sim_params) :: p
    character(len=256) :: config_path
    integer :: ib, i, j, seed, sweeps_done = 0, out = -1
    integer(kind=8) :: t_start, t_end, rate
    real :: beta, acc, acc_avg, plaq, dev


    call get_command_argument(1, config_path)
    if (len_trim(config_path) == 0) config_path = 'config/default.toml'
    p = load_params(trim(config_path))

    ! seed: the same on all images (each derives its own stream from it)
    seed = p%seed
    if (seed == 0) then
        if (this_image() == 1) then
            call system_clock(t_start)
            seed = int(mod(t_start, 2147483647_8)) + 1
        end if
        call co_broadcast(seed, 1)
    end if
    call seed_rng(seed, this_image())

    call u%initialize(p%nx, p%ny, p%nz, p%nt, p%hot_start)

    if (this_image() == 1) then
        call execute_command_line('mkdir -p ' // trim(p%output_dir))
        call open_output()
    end if
    sync all

    call system_clock(t_start, rate)
    if (this_image() == 1) print '(a,i0,a)', 'Start run: ', n_beta(p), ' beta values'

    ! beta scan
    do ib = 1, n_beta(p)
        beta = beta_at(p, ib)
        if (this_image() == 1) print '(a,f8.4)', 'Beta = ', beta

        ! thermalize
        do i = 1, p%n_thermal
            call sweep()
        end do

        ! measure
        do i = 1, p%n_measure
            do j = 1, p%sweeps_per_measurement
                call sweep()
            end do
            plaq = plaquette(u)
            acc_avg = acc
            call co_sum(acc_avg)
            acc_avg = acc_avg / num_images()
            if (this_image() == 1) then
                write (out, '(f8.4,",",i6,",",es16.8,",",f8.5)') beta, i, plaq, acc_avg
            end if
        end do

        if (this_image() == 1) then
            flush(out)
            print '(a,f8.5,a,f8.5)', '  last plaquette: ', plaq, '   acceptance: ', acc_avg
        end if
    end do

    call system_clock(t_end)
    dev = max_unitarity_violation()     ! collective: all images must call it
    if (this_image() == 1) then
        close(out)
        print '(a,f8.2,a)', 'Done in ', real(t_end - t_start) / real(rate), ' s'
        print '(a,es9.2)', 'Max link deviation from SU(3): ', dev
    end if

contains

    !! one Metropolis sweep over the whole lattice; re-projects links to SU(3) periodically
    subroutine sweep()
        call metro(u, p%eps, beta, acc)
        sweeps_done = sweeps_done + 1
        if (p%reunitarize_interval > 0) then
            if (mod(sweeps_done, p%reunitarize_interval) == 0) call unitarize_field()
        end if
    end subroutine sweep

    !! open the csv file and write the header with all run parameters
    subroutine open_output()
        character(len=512) :: path
        path = trim(p%output_dir) // '/' // trim(p%output_name) // '.csv'
        open(newunit=out, file=trim(path), status='replace', action='write')
        write (out, '(a)') '# su3-fortran raw plaquette measurements'
        write (out, '(a,a)') '# config = ', trim(config_path)
        write (out, '(a,i0)') '# seed = ', seed
        write (out, '(a,i0)') '# images = ', num_images()
        write (out, '(a,4(i0,1x))') '# lattice (nx ny nz nt) = ', p%nx, p%ny, p%nz, p%nt
        write (out, '(a,l1)') '# hot_start = ', p%hot_start
        write (out, '(a,f0.5)') '# epsilon = ', p%eps
        write (out, '(a,i0)') '# reunitarize_interval = ', p%reunitarize_interval
        write (out, '(a,i0)') '# thermalization_sweeps = ', p%n_thermal
        write (out, '(a,i0)') '# measurements = ', p%n_measure
        write (out, '(a,i0)') '# sweeps_per_measurement = ', p%sweeps_per_measurement
        write (out, '(a)') '# acceptance: Metropolis acceptance of the last sweep before the measurement, averaged over images'
        write (out, '(a)') 'beta,measurement,plaquette,acceptance'
    end subroutine open_output

end program main
