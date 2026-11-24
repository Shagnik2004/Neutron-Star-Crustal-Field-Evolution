module tridiagonal_solver

    use utilities
    implicit none

contains

   
    ! Solves a tridiagonal system A*x = d, where A has sub-diagonal
    ! a(2:n), main diagonal b(1:n), super-diagonal c(1:n-1).
    ! Thomas algorithm: forward elimination + back substitution,
    ! O(n) instead of O(n^3) for a general solve.
    ! a(1) and c(n) are unused (kept for uniform array indexing).
  
    subroutine thomas_solve(a, b, c, d, x, n)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: a(n), b(n), c(n), d(n)
        real(dp), intent(out):: x(n)

        real(dp) :: c_prime(n), d_prime(n)
        real(dp) :: m
        integer :: i

        c_prime(1) = c(1) / b(1)
        d_prime(1) = d(1) / b(1)

        do i = 2, n
            m = b(i) - a(i)*c_prime(i-1)
            c_prime(i) = c(i) / m
            d_prime(i) = (d(i) - a(i)*d_prime(i-1)) / m
        end do

        x(n) = d_prime(n)
        do i = n-1, 1, -1
            x(i) = d_prime(i) - c_prime(i)*x(i+1)
        end do

    end subroutine thomas_solve

end module tridiagonal_solver
