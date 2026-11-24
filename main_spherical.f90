program MainSpherical

    use utilities
    use grid
    use analytic_solution
    use cn_solver

    implicit none

    real(dp), allocatable :: r(:), v(:)
    real(dp) :: dr, dt, D, sigma0, r0
    real(dp) :: r_in, r_out, t_max, t
    integer :: N, n_steps, istep
    integer :: unit_diag, unit_snap
    integer :: snap_every, j
    character(len=64) :: snap_filename
    real(dp) :: v_exact_val, u_num, u_exact_val
    real(dp) :: max_relerr, relerr

    ! Generic validation setup: shell [5,15], Gaussian (in v=r*u space)
    ! centered mid-shell. Same D, sigma0 scale as the flat 1D test.
    ! Real crust scale: R_star = 10 km, crust thickness = 1 km
    r_in   = 9.0_dp
    r_out  = 10.0_dp
    r0     = 9.5_dp
    sigma0 = 0.05_dp
    D      = 0.1_dp    ! still placeholder -- not yet real magnetic diffusivity
    dr     = 0.001_dp
    dt     = 0.0003_dp
    t_max  = 0.03_dp

    print '(A,F10.4)', " r (stability diagnostic, CN unconditionally stable) = ", D*dt/dr**2

    call build_grid_fixed_dx(r_in, r_out, dr, r, N)
    print '(A,I0,A)', " Grid: N+1 = ", N+1, " points"
    print '(A,F6.2,A,F6.2,A)', " Shell: r in [", r_in, ", ", r_out, "]"

    allocate(v(0:N))
    ! v(r,0) = Gaussian centered at r0 (shifted analytic solution)
    v = u_exact(r - r0, 0.0_dp, sigma0, D)
    v(0) = 0.0_dp
    v(N) = 0.0_dp

    n_steps = nint(t_max/dt)
    snap_every = n_steps/10

    open(newunit=unit_diag, file="output/diagnostics_spherical.dat", status="replace", action="write")
    write(unit_diag, '(A)') "# t   max_relerr_u"

    t = 0.0_dp
    do istep = 1, n_steps

        call cn_step(v, N, D*dt/dr**2)
        t = real(istep, dp) * dt

        max_relerr = 0.0_dp
        block
            real(dp) :: sigma_t
            sigma_t = sqrt(sigma0**2 + 2.0_dp*D*t)
            do j = 1, N-1
                ! Only evaluate relative error within 3 sigma(t) of the
                ! peak -- outside that, u_exact is inherently small and
                ! relative error is dominated by that small denominator,
                ! not by real solver inaccuracy.
                if (abs(r(j) - r0) > 3.0_dp*sigma_t) cycle
                v_exact_val = u_exact(r(j) - r0, t, sigma0, D)
                u_num = v(j) / r(j)
                u_exact_val = v_exact_val / r(j)
                if (abs(u_exact_val) > 1.0e-12_dp) then
                    relerr = abs(u_num - u_exact_val) / abs(u_exact_val)
                    max_relerr = max(max_relerr, relerr)
                end if
            end do
        end block
        write(unit_diag, '(F10.4,3X,ES14.6)') t, max_relerr

        if (mod(istep, snap_every) == 0) then
            write(snap_filename, '(A,I0,A)') "output/snapshot_sph_t", istep/snap_every, ".dat"
            open(newunit=unit_snap, file=trim(snap_filename), status="replace", action="write")
            write(unit_snap, '(A,F8.3)') "# r   u_numeric   u_exact   at t = ", t
            do j = 0, N
                u_num = v(j) / r(j)
                u_exact_val = u_exact(r(j) - r0, t, sigma0, D) / r(j)
                write(unit_snap, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') r(j), u_num, u_exact_val
            end do
            close(unit_snap)
        end if

    end do

    close(unit_diag)

    print *
    print '(A,F10.4,A,ES14.6)', " Final t=", t, "  max relative error in u(r,t) = ", max_relerr
    print *
    print *, "Diagnostics written to output/diagnostics_spherical.dat"
    print *, "Snapshots written to output/snapshot_sph_t*.dat"
    print *

    deallocate(r, v)

end program MainSpherical
