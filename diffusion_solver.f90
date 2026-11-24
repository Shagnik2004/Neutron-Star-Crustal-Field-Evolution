module diffusion_solver

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! One FTCS explicit time step:
    !   u_j^{n+1} = u_j^n + r*(u_{j+1}^n - 2*u_j^n + u_{j-1}^n)
    ! r = D*dt/dx^2 must satisfy r <= 0.5 for stability (checked
    ! by the caller via check_stability, not silently assumed here).
    ! Dirichlet boundaries: u(0) = u(N) = 0 held fixed every step.
    !--------------------------------------------------------------
    subroutine ftcs_step(u, N, r)

        implicit none

        integer, intent(in)     :: N
        real(dp), intent(inout) :: u(0:N)
        real(dp), intent(in)    :: r

        real(dp) :: u_new(0:N)
        integer :: j

        do j = 1, N-1
            u_new(j) = u(j) + r*(u(j+1) - 2.0_dp*u(j) + u(j-1))
        end do

        u_new(0) = 0.0_dp
        u_new(N) = 0.0_dp

        u = u_new

    end subroutine ftcs_step


    !--------------------------------------------------------------
    ! Stability check: r = D*dt/dx^2 must be <= 0.5 (von Neumann
    ! analysis). Stops the program rather than silently producing
    ! a blown-up, meaningless run.
    !--------------------------------------------------------------
    subroutine check_stability(D, dt, dx)

        implicit none

        real(dp), intent(in) :: D, dt, dx
        real(dp) :: r

        r = D*dt/dx**2

        print '(A,F10.6)', " Stability parameter r = D*dt/dx^2 = ", r

        if (r > 0.5_dp) then
            print *, "ERROR: r > 0.5 -- FTCS scheme is unstable with these parameters."
            print *, "Reduce dt or increase dx."
            stop 1
        else
            print '(A,F6.1,A)', " Stable (", 100.0_dp*r/0.5_dp, "% of the stability limit)."
        end if

    end subroutine check_stability


    !--------------------------------------------------------------
    ! Trapezoidal-rule total mass: integral of u(x) dx over the grid.
    ! Should stay ~1.0 for all time if boundaries are far enough out.
    !--------------------------------------------------------------
    function total_mass(u, N, dx) result(mass)
        implicit none
        integer, intent(in)  :: N
        real(dp), intent(in) :: u(0:N)
        real(dp), intent(in) :: dx
        real(dp) :: mass
        integer :: j

        mass = 0.5_dp*(u(0) + u(N))
        do j = 1, N-1
            mass = mass + u(j)
        end do
        mass = mass * dx

    end function total_mass


    !--------------------------------------------------------------
    ! Trapezoidal-rule second moment <x^2> = integral x^2 u(x) dx.
    ! Should grow as sigma0^2 + 2*D*t (slope 2D) if the scheme is
    ! correctly reproducing diffusive spreading.
    !--------------------------------------------------------------
    function second_moment(u, x, N, dx) result(m2)
        implicit none
        integer, intent(in)  :: N
        real(dp), intent(in) :: u(0:N), x(0:N)
        real(dp), intent(in) :: dx
        real(dp) :: m2
        integer :: j

        m2 = 0.5_dp*(x(0)**2*u(0) + x(N)**2*u(N))
        do j = 1, N-1
            m2 = m2 + x(j)**2*u(j)
        end do
        m2 = m2 * dx

    end function second_moment

end module diffusion_solver
