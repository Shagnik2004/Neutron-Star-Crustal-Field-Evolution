program MainMovieSpherical

    use utilities
    use grid
    use analytic_solution
    use cn_solver

    implicit none

    real(dp), allocatable :: r(:), v(:)
    real(dp) :: dr, dt, D, sigma0, r0
    real(dp) :: r_in, r_out, t_max, t
    integer :: N, n_steps, istep
    integer :: n_frames, frame_every, frame_num, j
    integer :: unit_frame
    character(len=80) :: frame_filename

    r_in   = 9.0_dp
    r_out  = 10.0_dp
    r0     = 9.5_dp
    sigma0 = 0.05_dp
    D      = 0.1_dp
    dr     = 0.001_dp
    dt     = 0.0003_dp
    t_max  = 0.03_dp

    n_frames = 150

    call build_grid_fixed_dx(r_in, r_out, dr, r, N)

    allocate(v(0:N))
    v = u_exact(r - r0, 0.0_dp, sigma0, D)
    v(0) = 0.0_dp
    v(N) = 0.0_dp

    n_steps = nint(t_max/dt)
    frame_every = max(1, n_steps/n_frames)

    print '(A,I0,A,I0,A)', " Writing ~", n_frames, " frames every ", frame_every, " steps"

    frame_num = 0
    write(frame_filename, '(A,I4.4,A)') "output/movie_frames_sph/frame_", frame_num, ".dat"
    open(newunit=unit_frame, file=trim(frame_filename), status="replace", action="write")
    write(unit_frame, '(A,F10.6)') "# r   u_numeric   u_exact   t=", 0.0_dp
    do j = 0, N
        write(unit_frame, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') &
              r(j), v(j)/r(j), u_exact(r(j)-r0, 0.0_dp, sigma0, D)/r(j)
    end do
    close(unit_frame)

    t = 0.0_dp
    do istep = 1, n_steps

        call cn_step(v, N, D*dt/dr**2)
        t = real(istep, dp) * dt

        if (mod(istep, frame_every) == 0 .or. istep == n_steps) then
            frame_num = frame_num + 1
            write(frame_filename, '(A,I4.4,A)') "output/movie_frames_sph/frame_", frame_num, ".dat"
            open(newunit=unit_frame, file=trim(frame_filename), status="replace", action="write")
            write(unit_frame, '(A,F10.6)') "# r   u_numeric   u_exact   t=", t
            do j = 0, N
                write(unit_frame, '(F10.4,3X,ES22.10E3,3X,ES22.10E3)') &
                      r(j), v(j)/r(j), u_exact(r(j)-r0, t, sigma0, D)/r(j)
            end do
            close(unit_frame)
        end if

    end do

    print '(A,I0,A)', " Wrote ", frame_num+1, " frames to output/movie_frames_sph/"

    deallocate(r, v)

end program MainMovieSpherical
