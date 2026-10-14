!! Minimal TOML reader - supports the subset needed for the run configuration:
!!   [section] headers, key = value pairs, # comments,
!!   integers, floats, booleans and "double quoted" strings.
!! Not supported: arrays, inline tables, nested tables, multi-line strings.
!! Keys are stored flat as "section.key" ("key" for entries before any header).
module mod_toml
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private

    public :: toml_table, toml_load

    integer, parameter :: max_entries = 128
    integer, parameter :: key_len = 64
    integer, parameter :: val_len = 256

    type :: toml_table
        integer :: n = 0
        character(len=key_len) :: keys(max_entries) = ''
        character(len=val_len) :: vals(max_entries) = ''
    contains
        procedure :: has
        procedure :: get_int
        procedure :: get_real
        procedure :: get_logical
        procedure :: get_string
    end type toml_table

contains

    !! print message and abort all images
    subroutine fail(msg)
        character(len=*), intent(in) :: msg
        write (error_unit, '(a)') 'TOML error: ' // msg
        error stop 1
    end subroutine fail


    !! read and parse a TOML file
    function toml_load(path) result(tbl)
        character(len=*), intent(in) :: path
        type(toml_table) :: tbl
        character(len=1024) :: raw, line
        character(len=key_len) :: section
        character(len=key_len) :: key
        integer :: unit, st, lineno, eq, k, last

        open(newunit=unit, file=path, status='old', action='read', iostat=st)
        if (st /= 0) call fail('cannot open config file "' // trim(path) // '"')

        section = ''
        lineno = 0
        do
            read(unit, '(a)', iostat=st) raw
            if (st /= 0) exit
            lineno = lineno + 1

            line = adjustl(strip_comment(raw))
            if (len_trim(line) == 0) cycle

            if (line(1:1) == '[') then
                last = index(line, ']')
                if (last == 0) call fail(location(path, lineno) // 'unterminated section header')
                section = adjustl(line(2:last-1))
                if (len_trim(section) == 0) call fail(location(path, lineno) // 'empty section name')
                cycle
            end if

            eq = index(line, '=')
            if (eq < 2) call fail(location(path, lineno) // 'expected "key = value"')
            key = adjustl(line(:eq-1))
            if (len_trim(section) > 0) key = trim(section) // '.' // trim(key)
            if (tbl%has(key)) call fail(location(path, lineno) // 'duplicate key "' // trim(key) // '"')
            if (tbl%n == max_entries) call fail('too many entries in config file')

            k = tbl%n + 1
            tbl%keys(k) = key
            tbl%vals(k) = adjustl(line(eq+1:))
            if (len_trim(tbl%vals(k)) == 0) call fail(location(path, lineno) // 'missing value')
            tbl%n = k
        end do

        close(unit)
    end function toml_load


    function location(path, lineno) result(s)
        character(len=*), intent(in) :: path
        integer, intent(in) :: lineno
        character(len=:), allocatable :: s
        character(len=16) :: num
        write (num, '(i0)') lineno
        s = trim(path) // ':' // trim(num) // ': '
    end function location


    !! remove a trailing # comment (ignoring # inside quoted strings)
    function strip_comment(line) result(res)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: res
        logical :: in_quote
        integer :: i

        res = line
        in_quote = .false.
        do i = 1, len(line)
            if (line(i:i) == '"') in_quote = .not. in_quote
            if (line(i:i) == '#' .and. .not. in_quote) then
                res(i:) = ''
                exit
            end if
        end do
    end function strip_comment


    function has(self, key)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: has
        has = find(self, key) > 0
    end function has


    function find(self, key) result(idx)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: idx, i
        idx = 0
        do i = 1, self%n
            if (trim(self%keys(i)) == trim(key)) then
                idx = i
                return
            end if
        end do
    end function find


    !! fetch the raw value; aborts if missing and no default is available
    function raw_value(self, key, found) result(v)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        logical, intent(out) :: found
        character(len=val_len) :: v
        integer :: idx
        idx = find(self, key)
        found = idx > 0
        v = ''
        if (found) v = self%vals(idx)
    end function raw_value


    function get_int(self, key, default) result(res)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        integer, intent(in), optional :: default
        integer :: res, st
        character(len=val_len) :: v
        logical :: found

        v = raw_value(self, key, found)
        if (.not. found) then
            if (present(default)) then
                res = default
                return
            end if
            call fail('missing required key "' // trim(key) // '"')
        end if
        if (verify(trim(v), '+-0123456789') /= 0) call fail('"' // trim(key) // '" must be an integer, got ' // trim(v))
        read(v, *, iostat=st) res
        if (st /= 0) call fail('"' // trim(key) // '" must be an integer, got ' // trim(v))
    end function get_int


    function get_real(self, key, default) result(res)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        real, intent(in), optional :: default
        real :: res
        integer :: st
        character(len=val_len) :: v
        logical :: found

        v = raw_value(self, key, found)
        if (.not. found) then
            if (present(default)) then
                res = default
                return
            end if
            call fail('missing required key "' // trim(key) // '"')
        end if
        if (verify(trim(v), '+-0123456789.eE') /= 0) call fail('"' // trim(key) // '" must be a number, got ' // trim(v))
        read(v, *, iostat=st) res
        if (st /= 0) call fail('"' // trim(key) // '" must be a number, got ' // trim(v))
    end function get_real


    function get_logical(self, key, default) result(res)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        logical, intent(in), optional :: default
        logical :: res
        character(len=val_len) :: v
        logical :: found

        v = raw_value(self, key, found)
        if (.not. found) then
            if (present(default)) then
                res = default
                return
            end if
            call fail('missing required key "' // trim(key) // '"')
        end if
        select case (trim(v))
        case ('true')
            res = .true.
        case ('false')
            res = .false.
        case default
            call fail('"' // trim(key) // '" must be true or false, got ' // trim(v))
        end select
    end function get_logical


    !! string values must be written in double quotes
    function get_string(self, key, default) result(res)
        class(toml_table), intent(in) :: self
        character(len=*), intent(in) :: key
        character(len=*), intent(in), optional :: default
        character(len=val_len) :: res
        character(len=val_len) :: v
        integer :: n
        logical :: found

        v = raw_value(self, key, found)
        if (.not. found) then
            if (present(default)) then
                res = default
                return
            end if
            call fail('missing required key "' // trim(key) // '"')
        end if
        n = len_trim(v)
        if (n < 2 .or. v(1:1) /= '"' .or. v(n:n) /= '"') call fail('"' // trim(key) // '" must be a quoted string, got ' // trim(v))
        res = v(2:n-1)
    end function get_string

end module mod_toml
