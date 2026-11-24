module g_profile_ch7_v2

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Confirmed Chapter 7 initial g-profile (matches Fig 7.1 shape
    ! and its predicted J_phi, Fig 7.2): baseline+Gaussian on the
    ! surface side, curvature-matched decaying Gaussian on the core
    ! side (so g'' is continuous across x0 -- avoids the spurious
    ! spike in J_phi that an equal-width join produces).
    !--------------------------------------------------------------
    subroutine g_profile_expelled_v2(x, g, npts, x0, sigma_right, amplitude)

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

    end subroutine g_profile_expelled_v2

end module g_profile_ch7_v2
