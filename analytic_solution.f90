module analytic_solution

    use utilities
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

    !--------------------------------------------------------------
    ! Exact solution: u(x,t) = 1/sqrt(2*pi*sigma(t)^2) * exp(-x^2/(2*sigma(t)^2))
    ! sigma(t)^2 = sigma0^2 + 2*D*t
    !--------------------------------------------------------------
    elemental function u_exact(x, t, sigma0, D) result(u)
        implicit none
        real(dp), intent(in) :: x, t, sigma0, D
        real(dp) :: u
        real(dp) :: sigma2

        sigma2 = sigma0**2 + 2.0_dp*D*t
        u = exp(-x**2/(2.0_dp*sigma2)) / sqrt(2.0_dp*pi*sigma2)

    end function u_exact

    elemental function sigma_of_t(t, sigma0, D) result(s)
        implicit none
        real(dp), intent(in) :: t, sigma0, D
        real(dp) :: s
        s = sqrt(sigma0**2 + 2.0_dp*D*t)
    end function sigma_of_t

end module analytic_solution
