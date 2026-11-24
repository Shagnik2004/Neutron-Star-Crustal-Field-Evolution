module g_profile_ch7

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Chapter 7 initial g-profile: a flat-near-surface background
    ! (sigmoid, ~1 near x=1, declining toward 0 deeper in) plus a
    ! localized Gaussian bump at x0 (the freshly-deposited current).
    ! Confirmed against Fig 7.1's shape: flat baseline, sharp rise,
    ! peak, dip below baseline past the peak (NOT a plain Gaussian
    ! tail -- the background's own decline continues past the bump).
    !--------------------------------------------------------------
    subroutine g_profile_expelled_ch7(x, g, npts, x0, fwhm, amplitude, x_bg, w_bg)

        implicit none

        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: x0, fwhm, amplitude, x_bg, w_bg

        real(dp) :: xw, gx, background
        integer :: i

        xw = fwhm / 2.3548_dp

        do i = 1, npts
            background = 1.0_dp / (1.0_dp + exp((x_bg - x(i))/w_bg))
            gx = x(i)-x0
            g(i) = background + amplitude*exp(-gx*gx/(2.0_dp*xw*xw))
        end do

    end subroutine g_profile_expelled_ch7

end module g_profile_ch7
