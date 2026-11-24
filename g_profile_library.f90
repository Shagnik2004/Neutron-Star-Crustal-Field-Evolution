module g_profile_library

    use utilities
    implicit none

contains

  
    ! Konar Ch.5 near-surface accreted-current profile (matches
    ! Fig 5.2): cumulative Gaussian integral, exp(-(x'-x0)^2/gwidth^2)
   
    subroutine g_profile_crustal(x, g, npts, dx, x0, gwidth, x_cutoff)
        implicit none
        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: dx, x0, gwidth, x_cutoff
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: gwsqr, anorm, tempg, gx
        integer :: i
        gwsqr = gwidth*gwidth
        anorm = gwidth*sqrt(pi)
        tempg = 0.0_dp
        do i = 1, npts
            if (x(i) > x_cutoff) then
                gx = abs(x0 - x(i))
                tempg = tempg + dx*exp(-gx*gx/gwsqr)/anorm
                g(i) = tempg
            else
                g(i) = 0.0_dp
            end if
        end do
    end subroutine g_profile_crustal


    ! Bhattacharya & Datta (1996) eq.(6) expelled-flux profile:
    ! cumulative Gaussian integral, exp(-(x'-x0)^2/(2*xw^2)), with
    ! the correct 1/(xw*sqrt(2*pi)) normalization.
    
    subroutine g_profile_expelled_datta(x, g, npts, dx, x0, xw, x_cutoff)
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
    end subroutine g_profile_expelled_datta


   
    ! Konar Ch.7 expelled-flux profile (matches Fig 7.1 and its
    ! predicted J_phi, Fig 7.2): baseline+Gaussian on the surface
    ! side, curvature-matched decaying Gaussian on the core side
    ! (continuous value, slope, AND curvature at x0 -- avoids a
    ! spurious spike in the predicted current).
    
    subroutine g_profile_expelled_ch7(x, g, npts, x0, sigma_right, amplitude)
        implicit none
        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: x0, sigma_right, amplitude
        real(dp) :: sigma_left, gx
        integer :: i
        sigma_left = sigma_right * sqrt((1.0_dp+amplitude)/amplitude)
        do i = 1, npts
            gx = x(i)-x0
            if (x(i) >= x0) then
                g(i) = 1.0_dp + amplitude*exp(-gx*gx/(2.0_dp*sigma_right*sigma_right))
            else
                g(i) = (1.0_dp+amplitude)*exp(-gx*gx/(2.0_dp*sigma_left*sigma_left))
            end if
        end do
    end subroutine g_profile_expelled_ch7

end module g_profile_library
