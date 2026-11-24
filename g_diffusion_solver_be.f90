module g_diffusion_solver_be

    use utilities
    use tridiagonal_solver
    implicit none

contains

    !--------------------------------------------------------------
    ! Fully implicit (backward Euler) step for
    !   dg/dt = S(x)*(g_xx - 2g/x^2)
    ! Unlike Crank-Nicolson, this is L-stable: the worst-case grid
    ! mode's amplification factor -> 0 as the step size grows,
    ! rather than CN's -> -1 (bounded but undamped). This actually
    ! damps the spurious oscillation that appears once the physical
    ! signal has decayed close to that mode's residual amplitude --
    ! CN cannot do this by construction, no matter how the step
    ! size or reaction-term treatment is tuned.
    ! Same boundary conditions as gcn_step (Dirichlet inner, Robin
    ! outer via ghost-point elimination).
    !--------------------------------------------------------------
    subroutine be_step(g, x, N, dx, x_s, S, dt)

        implicit none

        integer, intent(in)     :: N
        real(dp), intent(inout) :: g(N)
        real(dp), intent(in)    :: x(N), dx, x_s, dt
        real(dp), intent(in)    :: S(N)

        real(dp) :: a(N), b(N), c(N), d(N), sol(N)
        real(dp) :: r_j, Q_j
        integer :: j

        a(1) = 0.0_dp
        b(1) = 1.0_dp
        c(1) = 0.0_dp
        d(1) = 0.0_dp

        do j = 2, N-1
            r_j = dt*S(j)/(dx*dx)
            Q_j = dt*S(j)/x(j)**2

            a(j) = -r_j
            b(j) = 1.0_dp + 2.0_dp*r_j + Q_j
            c(j) = -r_j
            d(j) = g(j)
        end do

        r_j = dt*S(N)/(dx*dx)
        Q_j = dt*S(N)/x_s**2

        a(N) = -2.0_dp*r_j
        b(N) = 1.0_dp + 2.0_dp*r_j*(1.0_dp + dx/x_s) + Q_j
        c(N) = 0.0_dp
        d(N) = g(N)

        call thomas_solve(a, b, c, d, sol, N)
        g = sol

    end subroutine be_step

end module g_diffusion_solver_be
