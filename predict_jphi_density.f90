program PredictJphiDensity
    use utilities
    use crust_model
    use conductivity, only: uc

    implicit none

    integer, parameter :: npts = 2000
    real(dp) :: x(npts), g(npts), Jphi(npts)
    real(dp), parameter :: x0 = 0.925_dp, xw_x = 0.006_dp, amplitude = 10.5_dp
    real(dp) :: dx, rho0, rho_hi, u0, u_hi, dudx_at_x0, uw
    real(dp) :: rho_here, u_here, gu
    integer :: i, unit_out

    call init_crust_model("../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    ! Convert the stated x-space width into an equivalent width in
    ! u = (2*pi/9)*(log10(rho)-3), using the local slope du/dx at x0
    ! (from our real rho(x) profile) -- tests whether the profile is
    ! actually symmetric in a density-related coordinate rather than
    ! in x directly, which could explain the width/asymmetry issues.
    rho0   = get_rho_at_r(x0*crust_R_star_km)
    rho_hi = get_rho_at_r((x0+1.0e-4_dp)*crust_R_star_km)
    u0     = uc*(log10(rho0)-3.0_dp)
    u_hi   = uc*(log10(rho_hi)-3.0_dp)
    dudx_at_x0 = (u_hi-u0)/1.0e-4_dp
    uw = abs(dudx_at_x0) * xw_x

    print '(A,ES14.6)', " rho at x0: ", rho0
    print '(A,F10.4)', " u0: ", u0
    print '(A,F14.6)', " du/dx at x0: ", dudx_at_x0
    print '(A,F10.6,A,F10.6,A)', " x-space width ", xw_x, " -> u-space width ", uw, " (via local slope)"

    dx = (1.0_dp - 0.85_dp)/real(npts-1,dp)
    do i = 1, npts
        x(i) = 0.85_dp + real(i-1,dp)*dx
        rho_here = get_rho_at_r(x(i)*crust_R_star_km)
        u_here = uc*(log10(rho_here)-3.0_dp)
        gu = (u_here-u0)/uw
        g(i) = amplitude*exp(-0.5_dp*gu*gu)
    end do

    do i = 2, npts-1
        Jphi(i) = -(( (g(i+1)-2.0_dp*g(i)+g(i-1))/dx**2 ) - 2.0_dp*g(i)/x(i)**2)
    end do
    Jphi(1) = 0.0_dp
    Jphi(npts) = 0.0_dp

    open(newunit=unit_out, file="output/predicted_jphi_density.dat", status="replace", action="write")
    write(unit_out,'(A)') "# x   g(x)   Jphi(x)"
    do i = 1, npts
        write(unit_out, '(F12.8,3X,ES20.8E3,3X,ES20.8E3)') x(i), g(i), Jphi(i)
    end do
    close(unit_out)
    print *, "Written to output/predicted_jphi_density.dat"

end program PredictJphiDensity
