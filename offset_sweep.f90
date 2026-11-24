program OffsetSweep

    use utilities
    use crust_model
    use g_profile_module

    implicit none

    integer, parameter :: npts = 1000
    real(dp) :: x(npts), g(npts)
    real(dp) :: x_c, x_s, dx, x0
    real(dp), parameter :: xw = 0.01_dp
    real(dp) :: offsets(5)
    integer :: i, k

    call init_crust_model("../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    x_c = get_x_crust_core()
    x_s = 1.0_dp
    print '(A,F10.6,A,F10.6)', " x_c = ", x_c, "   crust width = ", x_s-x_c
    print '(A,F10.4)', " crust width in units of xw: ", (x_s-x_c)/xw
    print *

    dx = (x_s - x_c) / real(npts-1, dp)
    do i = 1, npts
        x(i) = x_c + real(i-1,dp)*dx
    end do

    offsets = (/ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp /)

    print '(A)', "  k    x0        g(x_c)        g(x_s)      dist_to_xs/xw"
    do k = 1, 5
        x0 = x_c + offsets(k)*xw
        call g_profile(x, g, npts, dx, x0, xw, 0.0_dp)
        print '(F4.1,3X,F8.6,3X,ES12.4,3X,F10.6,3X,F10.4)', &
              offsets(k), x0, g(1), g(npts), (x_s-x0)/xw
    end do

end program OffsetSweep
