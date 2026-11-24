module pchip

    use utilities
    use binary_search
    implicit none

contains

  
    ! Computes derivative estimates d(1..n) using the Fritsch-Carlson
    ! method: weighted harmonic mean at interior points, with a zero
    ! derivative forced wherever the data changes monotonic direction
    ! (delta changes sign) or is locally flat.
   
    subroutine pchip_coeffs(x, y, n, d)

        implicit none

        integer, intent(in)   :: n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(in)  :: y(n)
        real(dp), intent(out) :: d(n)

        real(dp) :: h(n-1)
        real(dp) :: delta(n-1)
        real(dp) :: w1, w2

        integer :: i

        if (n < 3) then
            print *, "ERROR: pchip_coeffs requires n >= 3"
            stop 1
        end if

        do i = 1, n-1
            h(i) = x(i+1) - x(i)
            delta(i) = (y(i+1) - y(i)) / h(i)
        end do

        ! Interior points
        do i = 2, n-1
            if (delta(i-1)*delta(i) <= 0.0_dp) then
                ! Sign change (or a flat segment touching a non-flat one)
                ! -> zero derivative, this is what prevents overshoot
                d(i) = 0.0_dp
            else
                w1 = 2.0_dp*h(i) + h(i-1)
                w2 = h(i) + 2.0_dp*h(i-1)
                d(i) = (w1 + w2) / (w1/delta(i-1) + w2/delta(i))
            end if
        end do

        ! Left endpoint: non-centered three-point estimate, then clamp
        d(1) = ((2.0_dp*h(1) + h(2))*delta(1) - h(1)*delta(2)) / (h(1) + h(2))
        if (d(1)*delta(1) <= 0.0_dp) then
            d(1) = 0.0_dp
        else if ( (delta(1)*delta(2) < 0.0_dp) .and. (abs(d(1)) > abs(3.0_dp*delta(1))) ) then
            d(1) = 3.0_dp*delta(1)
        end if

        ! Right endpoint: mirror of the left-endpoint formula
        d(n) = ((2.0_dp*h(n-1) + h(n-2))*delta(n-1) - h(n-1)*delta(n-2)) / (h(n-1) + h(n-2))
        if (d(n)*delta(n-1) <= 0.0_dp) then
            d(n) = 0.0_dp
        else if ( (delta(n-1)*delta(n-2) < 0.0_dp) .and. (abs(d(n)) > abs(3.0_dp*delta(n-1))) ) then
            d(n) = 3.0_dp*delta(n-1)
        end if

    end subroutine pchip_coeffs



    ! Evaluates the PCHIP interpolant at x_query using cubic Hermite
    ! basis functions and the derivative estimates from pchip_coeffs.
   
    subroutine pchip_eval(x, y, d, n, x_query, y_query)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y(n)
        real(dp), intent(in) :: d(n)
        real(dp), intent(in) :: x_query
        real(dp), intent(out):: y_query

        integer  :: idx
        real(dp) :: h, t, t2, t3
        real(dp) :: h00, h10, h01, h11

        call find_interval(x, n, x_query, idx)

        h = x(idx+1) - x(idx)
        t = (x_query - x(idx)) / h

        t2 = t*t
        t3 = t2*t

        h00 =  2.0_dp*t3 - 3.0_dp*t2 + 1.0_dp
        h10 =  t3 - 2.0_dp*t2 + t
        h01 = -2.0_dp*t3 + 3.0_dp*t2
        h11 =  t3 - t2

        y_query = h00*y(idx) + h10*h*d(idx) + h01*y(idx+1) + h11*h*d(idx+1)

    end subroutine pchip_eval

end module pchip
