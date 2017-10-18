module mod_measure
    use mod_gauge_field
    use mod_su3_field_var
    use mod_real_field_var
    implicit none


    contains

    subroutine wilson(u, spt, sps, sp)
            type(gauge_field), intent(in) :: u
            INTEGER :: x, y, z, t
            REAL, dimension(6) :: plqs
            REAL, dimension(u%nt), intent(out) :: spt, sps
            REAL, intent(out) :: sp


            ! temporal and spatial plaquettes per time slice
            spt = 0
            sps = 0

            ! average plaquette (t + s)
            sp = 0

            plqs(1) = REAL(trace( &
                multSu3( U(1, x, y, z, t, :, :), .false., &
                multSu3( U(2, xp, y, z, t, :, :), .false., &
                multSu3( U(1, x, yp, z, t, :, :), .true., &
                U(2, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            plqs(2) = REAL(trace( &
                multSu3( U(1, x, y, z, t, :, :), .false., &
                multSu3( U(3, xp, y, z, t, :, :), .false., &
                multSu3( U(1, x, y, zp, t, :, :), .true., &
                U(3, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            plqs(3) = REAL(trace( &
                multSu3( U(1, x, y, z, t, :, :), .false., &
                multSu3( U(4, xp, y, z, t, :, :), .false., &
                multSu3( U(1, x, y, z, tp, :, :), .true., &
                U(4, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            plqs(4) = REAL(trace( &
                multSu3( U(2, x, y, z, t, :, :), .false., &
                multSu3( U(3, x, yp, z, t, :, :), .false., &
                multSu3( U(2, x, y, zp, t, :, :), .true., &
                U(3, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            plqs(5) = REAL(trace( &
                multSu3( U(2, x, y, z, t, :, :), .false., &
                multSu3( U(4, x, yp, z, t, :, :), .false., &
                multSu3( U(2, x, y, z, tp, :, :), .true., &
                U(4, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            plqs(6) = REAL(trace( &
                multSu3( U(3, x, y, z, t, :, :), .false., &
                multSu3( U(4, x, y, zp, t, :, :), .false., &
                multSu3( U(3, x, y, z, tp, :, :), .true., &
                U(4, x, y, z, t, :, :), .true.), &
                .false. ), .false.) &
                ))

            sps(t) = sps(t) + 1.0/3.0 * (plqs(1) + plqs(2) + plqs(4))
            spt(t) = spt(t) + 1.0/3.0 * (plqs(3) + plqs(5) + plqs(6))
            sp = sp + 1.0/3.0 * (plqs(1) + plqs(2) + plqs(3) + plqs(4) + plqs(5) + plqs(6))


            sps = sps / (3 * Nx * Ny * Nz)
            spt = spt / (3 * Nx * Ny * Nz)
            sp = sp / (6 * Nx * Ny * Nz * Nt)

        end subroutine wilson

end module mod_measure
