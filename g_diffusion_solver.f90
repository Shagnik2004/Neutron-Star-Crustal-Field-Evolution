module g_diffusion_solver
    use utilities
    use tridiagonal_solver
    implicit none
contains
    !--------------------------------------------------------------
    ! Generalized theta-method step for
    !   dg/dt = S(x)*(g_xx - 2g/x^2)
    ! theta=0.5 reproduces the original Crank-Nicolson scheme
    ! exactly (default, backward-compatible with every prior call).
    ! theta slightly above 0.5 (e.g. 0.55-0.6) introduces mild
    ! numerical damping -- CN's known non-L-stability leaves its
    ! highest-frequency error mode with amplification factor -> -1
    ! (never fully decaying, just flipping sign each step) for
    ! stiff problems. A small theta shift stays in the same
    ! implicit-averaging family as CN (not a scheme swap to BE)
    ! while killing the persistent oscillation.
    !
    ! Both the diffusion term AND the reaction term (-2g/x^2) are
    ! theta-averaged (old/new split) -- a deliberate departure
    ! from the thesis's eq 5.25, which treats the reaction term
    ! fully explicitly.
    !
    ! Inner boundary (j=1): Dirichlet g=0 (eq 5.18, fixed -- no
    ! accretion). Outer boundary (j=N): Robin g_x + g/x_s = 0
    ! (eq 5.17), via ghost-point elimination.
    !--------------------------------------------------------------
    subroutine gcn_step(g, x, N, dx, x_s, S, dt, theta_in)
        implicit none
        integer, intent(in)     :: N
        real(dp), intent(inout) :: g(N)
        real(dp), intent(in)    :: x(N), dx, x_s, dt
        real(dp), intent(in)    :: S(N)
        real(dp), intent(in), optional :: theta_in
        real(dp) :: a(N), b(N), c(N), d(N), sol(N)
        real(dp) :: r_full, Q_full, theta_r, theta_Q, anti_r, anti_Q
        real(dp) :: theta
        integer :: j

        theta = 0.5_dp
        if (present(theta_in)) theta = theta_in

        ! Inner boundary: g=0
        a(1) = 0.0_dp
        b(1) = 1.0_dp
        c(1) = 0.0_dp
        d(1) = 0.0_dp

        ! Interior points
        do j = 2, N-1
            r_full = dt*S(j)/dx**2
            Q_full = 2.0_dp*dt*S(j)/x(j)**2
            theta_r = theta*r_full
            theta_Q = theta*Q_full
            anti_r = (1.0_dp-theta)*r_full
            anti_Q = (1.0_dp-theta)*Q_full

            a(j) = -theta_r
            b(j) = 1.0_dp + 2.0_dp*theta_r + theta_Q
            c(j) = -theta_r
            d(j) = anti_r*g(j-1) + (1.0_dp - 2.0_dp*anti_r - anti_Q)*g(j) + anti_r*g(j+1)
        end do

        ! Outer boundary: Robin, via ghost-point elimination,
        ! reaction term also theta-averaged
        r_full = dt*S(N)/dx**2
        Q_full = 2.0_dp*dt*S(N)/x_s**2
        theta_r = theta*r_full
        theta_Q = theta*Q_full
        anti_r = (1.0_dp-theta)*r_full
        anti_Q = (1.0_dp-theta)*Q_full

        a(N) = -2.0_dp*theta_r
        b(N) = 1.0_dp + 2.0_dp*theta_r*(1.0_dp + dx/x_s) + theta_Q
        c(N) = 0.0_dp
        d(N) = g(N) + anti_r*(2.0_dp*g(N-1) - 2.0_dp*g(N)*(1.0_dp+dx/x_s)) &
               - anti_Q*g(N)

        call thomas_solve(a, b, c, d, sol, N)
        g = sol
    end subroutine gcn_step
end module g_diffusion_solver
