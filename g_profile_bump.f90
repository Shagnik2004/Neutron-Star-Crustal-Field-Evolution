module g_profile_bump

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Direct (non-integrated) Gaussian bump, distinct from both
    ! Konar's g_profile.f (cumulative integral, S-curve) and eq.(6)
    ! of Bhattacharya & Datta (also a cumulative integral). Chapter
    ! 7's Figure 7.1 is visibly a peaked bump on a log scale (a
    ! downward parabola = log of a plain Gaussian), not an S-curve --
    ! so this uses the plain, un-integrated form directly:
    !   g(x,0) = A * exp(-(x-x0)^2 / (2*xw^2))
    !--------------------------------------------------------------
    subroutine g_profile_gaussian(x, g, npts, x0, xw, amplitude, x_cutoff)

        implicit none

        integer, intent(in)   :: npts
        real(dp), intent(in)  :: x(npts)
        real(dp), intent(out) :: g(npts)
        real(dp), intent(in)  :: x0, xw, amplitude, x_cutoff

        real(dp) :: gx
        integer :: i

        do i = 1, npts
            if (x(i) > x_cutoff) then
                gx = x(i) - x0
                g(i) = amplitude * exp(-gx*gx/(2.0_dp*xw*xw))
            else
                g(i) = 0.0_dp
            end if
        end do

    end subroutine g_profile_gaussian

end module g_profile_bump
