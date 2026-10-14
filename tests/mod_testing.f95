!! Tiny test helper: check() records a result, finish() aborts if anything failed
module mod_testing
    implicit none
    integer, private :: n_checks = 0, n_failed = 0

contains

    subroutine check(name, ok)
        character(len=*), intent(in) :: name
        logical, intent(in) :: ok

        n_checks = n_checks + 1
        if (ok) then
            if (this_image() == 1) print '(2a)', '  ok    ', name
        else
            n_failed = n_failed + 1
            print '(a,i0,2a)', '  FAIL (image ', this_image(), ') ', name
        end if
    end subroutine check

    !! check a real value against an interval; prints the value on failure
    subroutine check_range(name, x, lo, hi)
        character(len=*), intent(in) :: name
        real, intent(in) :: x, lo, hi

        call check(name, x >= lo .and. x <= hi)
        if (x < lo .or. x > hi) print '(a,es12.4,a,es10.3,a,es10.3,a)', '        value ', x, ' not in [', lo, ',', hi, ']'
    end subroutine check_range

    subroutine finish()
        if (n_failed > 0) then
            print '(a,i0,a,i0,a)', 'FAILED: ', n_failed, ' of ', n_checks, ' checks'
            error stop 1
        end if
    end subroutine finish

end module mod_testing
