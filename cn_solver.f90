module cn_solver

    use utilities
    use tridiagonal_solver
    implicit none

contains

    !--------------------------------------------------------------
    ! One Crank-Nicolson step, unconditionally stable for any r.
    ! Interior points j=1..N-1 solved as a tridiagonal system;
    ! u(0)=u(N)=0 held fixed (Dirichlet), consistent with the
    ! explicit FTCS solver, so the two are directly comparable.
    !--------------------------------------------------------------
    subroutine cn_step(u, N, r)

        implicit none

        integer, intent(in)     :: N
        real(dp), intent(inout) :: u(0:N)
        real(dp), intent(in)    :: r

        real(dp) :: a(N-1), b(N-1), c(N-1), d(N-1), x_sol(N-1)
        integer :: j, k

        ! Build tridiagonal system for interior points j=1..N-1
        ! (local index k=1..N-1 maps to global j=k)
        do k = 1, N-1
            j = k
            a(k) = -r/2.0_dp
            b(k) = 1.0_dp + r
            c(k) = -r/2.0_dp

            d(k) = (r/2.0_dp)*u(j-1) + (1.0_dp-r)*u(j) + (r/2.0_dp)*u(j+1)
        end do

        ! Boundary contributions: u(0)=u(N)=0, so the r/2*u(j-1) term
        ! at j=1 and r/2*u(j+1) term at j=N-1 are already correctly
        ! zero via the fixed boundary values -- no extra correction
        ! needed for this simple Dirichlet-zero case.

        call thomas_solve(a, b, c, d, x_sol, N-1)

        do k = 1, N-1
            u(k) = x_sol(k)
        end do

        u(0) = 0.0_dp
        u(N) = 0.0_dp

    end subroutine cn_step

end module cn_solver
