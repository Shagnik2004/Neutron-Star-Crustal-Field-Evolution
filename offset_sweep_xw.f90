program OffsetSweepXW

    use utilities
    use crust_model
    use g_profile_expelled

    implicit none

    integer, parameter :: npts = 1000
    real(dp) :: x(npts), g(npts)
    real(dp) :: x_c, x_s, dx, x0
    real(dp) :: xw_vals(3)
    real(dp) :: k_vals(5)
    integer :: i, j, k

    call init_crust_model("../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    x_c = get_x_crust_core()
    x_s = 1.0_dp
    print '(A,F10.6,A,F10.6)', " x_c = ", x_c, "   crust width = ", x_s-x_c
    print *

    dx = (x_s - x_c) / real(npts-1, dp)
    do i = 1, npts
        x(i) = x_c + real(i-1,dp)*dx
    end do

    xw_vals = (/ 0.010_dp, 0.005_dp, 0.003_dp /)
    k_vals  = (/ 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp /)

    do j = 1, 3
        print '(A,F6.4)', " --- xw = ", xw_vals(j)
        print '(A)', "   k    x0          g(x_s)     dist_to_xc/xw  dist_to_xs/xw"
        do k = 1, 5
            x0 = x_c + k_vals(k)*xw_vals(j)
            call g_profile_datta(x, g, npts, dx, x0, xw_vals(j), x_c)
            print '(F5.1,3X,F9.6,3X,F10.6,3X,F10.4,3X,F10.4)', &
                  k_vals(k), x0, g(npts), (x0-x_c)/xw_vals(j), (x_s-x0)/xw_vals(j)
        end do
        print *
    end do

end program OffsetSweepXW
