module cubic_spline

    use utilities
    use binary_search
    implicit none

contains

    subroutine spline_coeffs(x, y, n, M)

        implicit none

        integer, intent(in)   :: n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(in)  :: y(n)
        real(dp), intent(out) :: M(n)

        real(dp) :: h(n-1)
        real(dp) :: alpha(n)
        real(dp) :: l(n)
        real(dp) :: mu(n)
        real(dp) :: z(n)

        integer :: i

        do i = 1, n-1
            h(i) = x(i+1) - x(i)
        end do

        alpha(1) = 0.0_dp
        do i = 2, n-1
            alpha(i) = (3.0_dp/h(i))*(y(i+1)-y(i)) - (3.0_dp/h(i-1))*(y(i)-y(i-1))
        end do

        l(1)  = 1.0_dp
        mu(1) = 0.0_dp
        z(1)  = 0.0_dp

        do i = 2, n-1
            l(i)  = 2.0_dp*(x(i+1)-x(i-1)) - h(i-1)*mu(i-1)
            mu(i) = h(i)/l(i)
            z(i)  = (alpha(i) - h(i-1)*z(i-1)) / l(i)
        end do

        l(n) = 1.0_dp
        z(n) = 0.0_dp
        M(n) = 0.0_dp

        do i = n-1, 2, -1
            M(i) = z(i) - mu(i)*M(i+1)
        end do
        M(1) = 0.0_dp

      
        ! IMPORTANT: the tridiagonal solve above computes
        ! c_j = S''(x_j)/2 (Burden-Faires convention), not the
        ! second derivative itself. Scale by 2 so M truly holds
        ! S''(x_j), matching what spline_eval expects.
     
        M = 2.0_dp * M

    end subroutine spline_coeffs


    subroutine spline_eval(x, y, M, n, x_query, y_query)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y(n)
        real(dp), intent(in) :: M(n)
        real(dp), intent(in) :: x_query
        real(dp), intent(out):: y_query

        integer  :: idx
        real(dp) :: h, a, b, c, d, dx

        call find_interval(x, n, x_query, idx)

        h = x(idx+1) - x(idx)

        a = y(idx)
        b = (y(idx+1)-y(idx))/h - h*(2.0_dp*M(idx)+M(idx+1))/6.0_dp
        c = M(idx)/2.0_dp
        d = (M(idx+1)-M(idx))/(6.0_dp*h)

        dx = x_query - x(idx)

        y_query = a + b*dx + c*dx**2 + d*dx**3

    end subroutine spline_eval

end module cubic_spline

