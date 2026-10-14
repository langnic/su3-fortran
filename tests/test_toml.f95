!! config parsing and seeding (no lattice needed)
program test_toml
    use mod_testing
    use mod_toml
    use mod_config
    use mod_utils
    implicit none

    type(toml_table) :: t
    type(sim_params) :: p
    integer :: unit
    character(len=64) :: path
    real :: a(3), b(3)

    ! write a small config file (one per image: all images run this program)
    write (path, '(a,i0,a)') 'build/test_config_', this_image(), '.toml'
    open(newunit=unit, file=trim(path), status='replace')
    write (unit, '(a)') '# comment line'
    write (unit, '(a)') 'top = 7'
    write (unit, '(a)') '[sec]'
    write (unit, '(a)') 'n = 12            # trailing comment'
    write (unit, '(a)') 'x = 1.5e-1'
    write (unit, '(a)') 'flag = true'
    write (unit, '(a)') 's = "a # not a comment"'
    close(unit)

    t = toml_load(trim(path))
    call check('top-level key', t%get_int('top') == 7)
    call check('int in section', t%get_int('sec.n') == 12)
    call check('real with exponent', abs(t%get_real('sec.x') - 0.15) < 1e-6)
    call check('logical', t%get_logical('sec.flag'))
    call check('string keeps #', trim(t%get_string('sec.s')) == 'a # not a comment')
    call check('default for missing key', t%get_int('sec.missing', 5) == 5)
    call check('has()', t%has('sec.n') .and. .not. t%has('sec.zzz'))

    ! the shipped configs must be valid
    p = load_params('config/default.toml')
    call check('default.toml: 16^4 lattice', p%nx == 16 .and. p%nt == 16)
    call check('default.toml: beta grid', n_beta(p) == 3 .and. abs(beta_at(p, 3) - 0.6) < 1e-6)
    p = load_params('config/small.toml')
    call check('small.toml: 4^4 lattice', p%nx == 4 .and. p%nt == 4)
    call check('small.toml: beta grid includes beta_max', n_beta(p) == 12 .and. abs(beta_at(p, n_beta(p)) - 6.0) < 1e-5)
    call check('small.toml: fixed seed', p%seed == 42)

    ! seeding: same (seed, image) -> same stream; different image -> different stream
    call seed_rng(42, 1)
    call random_number(a)
    call seed_rng(42, 1)
    call random_number(b)
    call check('same seed reproduces random numbers', all(a == b))
    call seed_rng(42, 2)
    call random_number(b)
    call check('different image gives different stream', any(a /= b))
    call seed_rng(43, 1)
    call random_number(b)
    call check('different seed gives different stream', any(a /= b))

    call finish()
end program test_toml
