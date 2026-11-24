program OffsetSweepDatta

    use utilities
    use crust_model
    use g_profile_expelled

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
    print *

    dx = (x_s - x_c) / real(npts-1, dp)
    do i = 1, npts
        x(i) = x_c + real(i-1,dp)*dx
    end do

    offsets = (/ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp /)

    print '(A)', "  k    x0        g(x_s)"
    do k = 1, 5
        x0 = x_c + offsets(k)*xw
        call g_profile_datta(x, g, npts, dx, x0, xw, x_c)
        print '(F4.1,3X,F8.6,3X,F10.6)', offsets(k), x0, g(npts)
    end do

end program OffsetSweepDatta
