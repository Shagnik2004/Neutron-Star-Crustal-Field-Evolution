program MainCN

    use utilities
    use grid
    use analytic_solution
    use cn_solver

    implicit none

    real(dp), allocatable :: x(:), u(:)
    real(dp) :: dx, dt, D, sigma0
    real(dp) :: x_min, x_max, t_max, t
    real(dp) :: r, m2, sigma2_exact, mass
    integer :: N, n_steps, istep
    integer :: unit_diag, unit_snap
    integer :: snap_every
    character(len=64) :: snap_filename
    real(dp) :: t_start, t_end

    ! Same physical parameters as the explicit run, for direct
    ! comparison -- but note Crank-Nicolson has NO stability limit
    ! on dt, so this dt is chosen for accuracy, not stability.
    sigma0 = 1.0_dp
    D      = 0.1_dp
    x_min  = -120.0_dp
    x_max  =  120.0_dp
    dx     = 0.02_dp
    dt     = 0.08_dp
    t_max  = 500.0_dp

    r = D*dt/dx**2
    print '(A,F10.6)', " r = D*dt/dx^2 = ", r
    print *, "(Crank-Nicolson is unconditionally stable -- this r is just"
    print *, " a diagnostic number now, not a hard limit.)"

    call cpu_time(t_start)

    call build_grid_fixed_dx(x_min, x_max, dx, x, N)
    print '(A,I0,A)', " Grid: N+1 = ", N+1, " points"

    allocate(u(0:N))
    u = u_exact(x, 0.0_dp, sigma0, D)

    n_steps = nint(t_max/dt)
    snap_every = n_steps/10

    open(newunit=unit_diag, file="output/diagnostics_cn.dat", status="replace", action="write")
    write(unit_diag, '(A)') "# t   mass_numeric   variance_numeric   variance_exact"

    t = 0.0_dp
    mass = sum(u(1:N-1))*dx + 0.5_dp*(u(0)+u(N))*dx
    m2 = sum(x(1:N-1)**2*u(1:N-1))*dx + 0.5_dp*(x(0)**2*u(0)+x(N)**2*u(N))*dx
    write(unit_diag, '(F10.4,3X,F14.10,3X,F14.10,3X,F14.10)') t, mass, m2, sigma0**2

    do istep = 1, n_steps

        call cn_step(u, N, r)
        t = real(istep, dp) * dt

        if (mod(istep, max(1, n_steps/50)) == 0) then
            print '(A,I0,A,I0,A,F6.1,A)', "  step ", istep, " / ", n_steps, &
                  "  (", 100.0_dp*real(istep,dp)/real(n_steps,dp), "%)"
        end if

        mass = sum(u(1:N-1))*dx + 0.5_dp*(u(0)+u(N))*dx
        m2 = sum(x(1:N-1)**2*u(1:N-1))*dx + 0.5_dp*(x(0)**2*u(0)+x(N)**2*u(N))*dx
        sigma2_exact = sigma0**2 + 2.0_dp*D*t
        write(unit_diag, '(F10.4,3X,F14.10,3X,F14.10,3X,F14.10)') t, mass, m2, sigma2_exact

        if (mod(istep, snap_every) == 0) then
            write(snap_filename, '(A,I0,A)') "output/snapshot_cn_t", istep/snap_every, ".dat"
            open(newunit=unit_snap, file=trim(snap_filename), status="replace", action="write")
            write(unit_snap, '(A,F8.3)') "# x   u_numeric   u_exact   at t = ", t
            block
                integer :: j
                do j = 0, N
                    write(unit_snap, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') &
                          x(j), u(j), u_exact(x(j), t, sigma0, D)
                end do
            end block
            close(unit_snap)
        end if

    end do

    close(unit_diag)
    call cpu_time(t_end)

    print *
    print '(A,F10.4,A,F14.10)', " Final t=", t, "  mass=", mass
    print '(A,F14.10,A,F14.10)', "   numeric variance=", m2, "   exact variance=", sigma0**2+2.0_dp*D*t
    print '(A,F10.3,A)', " Wall time: ", t_end-t_start, " seconds"
    print *
    print *, "Diagnostics written to output/diagnostics_cn.dat"
    print *, "Snapshots written to output/snapshot_cn_t*.dat"
    print *

    deallocate(x, u)

end program MainCN
