module g_profile_expelled

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Bhattacharya & Datta (1996), eq. (6):
    !   g(x,0) = 1/(sqrt(2*pi)*xw) * integral_{-inf}^{x}
    !            exp(-(x'-x0)^2 / (2*xw^2)) dx'
    !
    ! Note the prefactor here is 1/(sqrt(2*pi)*xw), NOT the
    ! 1/sqrt(2*pi) shown in the printed equation alone -- the
    ! paper's own text states g varies "from nearly unity to nearly
    ! zero," which is only true with the xw in the denominator
    ! (confirmed: integral of the exponential over all x' equals
    ! xw*sqrt(2*pi), so this prefactor is exactly what's needed for
    ! g to approach 1, matching the text).
    !
    ! Distinct from Konar's g_profile module: that formula uses
    ! exp(-(x'-xm)^2/gwidth^2) with NO factor of 2 in the exponent's
    ! denominator, a different (narrower, for the same gwidth value)
    ! Gaussian than eq. (6) specifies.
    !--------------------------------------------------------------
    subroutine g_profile_datta(x, g, npts, dx, x0, xw, x_cutoff)

        implicit none

        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: dx, x0, xw, x_cutoff

        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: anorm, tempg, gx
        integer :: i

        anorm = xw*sqrt(2.0_dp*pi)
        tempg = 0.0_dp

        do i = 1, npts
            if (x(i) > x_cutoff) then
                gx = x(i) - x0
                tempg = tempg + dx*exp(-gx*gx/(2.0_dp*xw*xw))/anorm
                g(i) = tempg
            else
                g(i) = 0.0_dp
            end if
        end do

    end subroutine g_profile_datta

end module g_profile_expelled
