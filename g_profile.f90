module g_profile_module

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Urpin-Muslimov initial field/current distribution.
    ! Direct F90 conversion of g_profile.f.
    !
    ! g(x) is built as a running (cumulative) integral of a Gaussian
    ! current distribution centered at xm, width gwidth -- NOT the
    ! Gaussian itself. This is why g(x) is an S-shaped, saturating
    ! curve (see Figure 5.2), and why J_phi ~ -(g'' - 2g/x^2) comes
    ! out two-lobed (g'' is the Gaussian's own derivative, odd and
    ! two-lobed).
    !
    ! x(1:npts): grid, must be UNIFORM (equal spacing dx) since the
    !            integral is a simple Riemann sum.
    ! xm       : center of the underlying Gaussian current
    ! gwidth   : width of the underlying Gaussian current
    ! xl       : below this x, g is forced to exactly 0
    !--------------------------------------------------------------
    subroutine g_profile(x, g, npts, dx, xm, gwidth, xl)

        implicit none

        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: dx, xm, gwidth, xl

        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: gwsqr, anorm, tempg, gx
        integer :: i

        gwsqr = gwidth*gwidth
        anorm = gwidth*sqrt(pi)
        tempg = 0.0_dp

        do i = 1, npts
            if (x(i) > xl) then
                gx = abs(xm - x(i))
                tempg = tempg + dx*exp(-gx*gx/gwsqr)/anorm
                g(i) = tempg
            else
                g(i) = 0.0_dp
            end if
        end do

    end subroutine g_profile

end module g_profile_module
