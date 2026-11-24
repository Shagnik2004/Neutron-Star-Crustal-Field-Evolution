program ExpelledGProfile

    use utilities
    use crust_model
    use g_profile_expelled

    implicit none

    integer, parameter :: npts = 1000
    real(dp) :: x(npts), g(npts)
    real(dp) :: x_c, x_s, dx
    real(dp), parameter :: xw = 0.01_dp
    real(dp), parameter :: k_offset = 2.0_dp

    real(dp) :: x0, rho_at_x0, crust_thickness_km
    integer :: i, unit_out, unit_meta

    call init_crust_model("../../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    x_c = get_x_crust_core()
    x_s = 1.0_dp
    x0  = x_c + k_offset*xw
    crust_thickness_km = crust_R_star_km - get_crust_core_radius()

    print '(A,ES14.6,A)', " Crust-core boundary density used: ", rho_crust_core, " g/cm^3"
    print '(A,F10.6,A,F10.6)', " x_c (crust bottom) = ", x_c, "   x0 = ", x0
    print '(A,F10.4,A)', " crust thickness = ", crust_thickness_km, " km"

    rho_at_x0 = get_rho_at_r(x0*crust_R_star_km)
    print '(A,ES14.6)', " rho at x0 = ", rho_at_x0

    open(newunit=unit_meta, file="output/expelled_meta_gprofile.dat", status="replace", action="write")
    write(unit_meta, '(ES16.8)') rho_at_x0
    write(unit_meta, '(F12.6)') crust_thickness_km
    close(unit_meta)

    dx = (x_s - x_c) / real(npts-1, dp)
    do i = 1, npts
        x(i) = x_c + real(i-1,dp)*dx
    end do

    call g_profile_datta(x, g, npts, dx, x0, xw, x_c)

    print '(A,F12.8)', " g at surface (x=1): ", g(npts)
    print '(A,ES14.6)', " g at crust bottom (x=x_c): ", g(1)

    open(newunit=unit_out, file="output/expelled_gprofile.dat", status="replace", action="write")
    write(unit_out, '(A)') "# x   g(x)   log10(g(x)/g0)"
    do i = 1, npts
        write(unit_out, '(F12.8,3X,ES16.8,3X,F14.8)') x(i), g(i), &
              merge(log10(g(i)/g(npts)), -99.0_dp, g(i) > 1.0e-30_dp)
    end do
    close(unit_out)

    print *, "Written to output/expelled_gprofile.dat"

end program ExpelledGProfile
