program MainMovie

    use utilities
    use grid
    use analytic_solution
    use cn_solver

    implicit none

    real(dp), allocatable :: x(:), u(:)
    real(dp) :: dx, dt, D, sigma0
    real(dp) :: x_min, x_max, t_max, t, r
    integer :: N, n_steps, istep
    integer :: n_frames, frame_every, frame_num
    integer :: stride
    integer :: unit_frame, j
    character(len=80) :: frame_filename
    real(dp) :: t_start, t_end

    ! Same validated physical setup as main_cn.f90
    sigma0 = 1.0_dp
    D      = 0.1_dp
    x_min  = -120.0_dp
    x_max  =  120.0_dp
    dx     = 0.02_dp
    dt     = 0.08_dp
    t_max  = 500.0_dp

    n_frames = 200
    stride   = 4   ! write every 4th grid point per frame (visual only)

    r = D*dt/dx**2
    print '(A,F10.4)', " r = ", r
    print *, "(Crank-Nicolson: unconditionally stable, this is just diagnostic.)"

    call cpu_time(t_start)

    call build_grid_fixed_dx(x_min, x_max, dx, x, N)
    print '(A,I0,A)', " Grid: N+1 = ", N+1, " points"

    allocate(u(0:N))
    u = u_exact(x, 0.0_dp, sigma0, D)

    n_steps = nint(t_max/dt)
    frame_every = max(1, n_steps/n_frames)

    print '(A,I0,A,I0,A)', " Will write ~", n_frames, " frames every ", frame_every, " steps"

    ! Frame 0: initial condition
    frame_num = 0
    write(frame_filename, '(A,I4.4,A)') "output/movie_frames/frame_", frame_num, ".dat"
    open(newunit=unit_frame, file=trim(frame_filename), status="replace", action="write")
    write(unit_frame, '(A,F10.4)') "# x   u_numeric   u_exact   t=", 0.0_dp
    do j = 0, N, stride
        write(unit_frame, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') x(j), u(j), u_exact(x(j), 0.0_dp, sigma0, D)
    end do
    close(unit_frame)

    t = 0.0_dp
    do istep = 1, n_steps

        call cn_step(u, N, r)
        t = real(istep, dp) * dt

        if (mod(istep, frame_every) == 0 .or. istep == n_steps) then
            frame_num = frame_num + 1
            write(frame_filename, '(A,I4.4,A)') "output/movie_frames/frame_", frame_num, ".dat"
            open(newunit=unit_frame, file=trim(frame_filename), status="replace", action="write")
            write(unit_frame, '(A,F10.4)') "# x   u_numeric   u_exact   t=", t
            do j = 0, N, stride
                write(unit_frame, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') x(j), u(j), u_exact(x(j), t, sigma0, D)
            end do
            close(unit_frame)
        end if

        if (mod(istep, max(1,n_steps/20)) == 0) then
            print '(A,F5.1,A)', "  progress: ", 100.0_dp*real(istep,dp)/real(n_steps,dp), "%"
        end if

    end do

    call cpu_time(t_end)
    print *
    print '(A,I0,A)', " Wrote ", frame_num+1, " frames to output/movie_frames/"
    print '(A,F10.3,A)', " Wall time: ", t_end-t_start, " seconds"
    print *

    deallocate(x, u)

end program MainMovie
