program KonarCh7AnalyticForm

    use utilities
    use crust_model
    use cubic_spline
    use conductivity, only: uc

    implicit none

    integer, parameter :: npts = 3000
    real(dp) :: x(npts), g(npts), Jphi(npts), logg(npts)
    real(dp), parameter :: x0 = 0.925_dp, xw_x = 0.006_dp, amplitude = 10.5_dp
    real(dp), parameter :: x_lo = 0.85_dp, x_hi = 1.0_dp

    ! Smooth (C2) representation of u(x) = (2*pi/9)*(log10(rho(x))-3),
    ! built from a plain cubic spline over our real density profile --
    ! NOT PCHIP, which only guarantees C1 continuity and produces
    ! spurious high-frequency noise once differentiated twice.
    integer, parameter :: n_prof_pts = 300
    real(dp) :: x_prof(n_prof_pts), u_prof(n_prof_pts), u2_prof(n_prof_pts)
    real(dp) :: rho0, u0, dx, gu

    integer :: i, unit_out

    call init_crust_model("../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    ! Build the smooth u(x) table
    do i = 1, n_prof_pts
        x_prof(i) = x_lo + (x_hi-x_lo)*real(i-1,dp)/real(n_prof_pts-1,dp)
        u_prof(i) = uc*(log10(get_rho_at_r(x_prof(i)*crust_R_star_km)) - 3.0_dp)
    end do
    call spline_coeffs(x_prof, u_prof, n_prof_pts, u2_prof)

    rho0 = get_rho_at_r(x0*crust_R_star_km)
    u0 = uc*(log10(rho0)-3.0_dp)

    print '(A,F10.6,A,F10.6)', " g-profile centred at x0 = ", x0, ", delta_x = ", xw_x
    print '(A,ES14.6,A,F10.4)', " rho at x0 = ", rho0, "   u0 = ", u0

    ! g(x) = A * exp(-0.5*((u(x)-u0)/uw)^2), where uw is set so that
    ! the FWHM in x matches xw_x, using the local slope du/dx at x0
    block
        real(dp) :: u_at_x0_plus, dudx0, uw
        call spline_eval(x_prof, u_prof, u2_prof, n_prof_pts, x0+1.0e-4_dp, u_at_x0_plus)
        dudx0 = (u_at_x0_plus - u0)/1.0e-4_dp
        uw = abs(dudx0) * xw_x
        print '(A,F14.6,A,F10.6)', " du/dx at x0 = ", dudx0, "   u-space width = ", uw

        dx = (x_hi - x_lo) / real(npts-1, dp)
        do i = 1, npts
            x(i) = x_lo + real(i-1,dp)*dx
            block
                real(dp) :: u_here
                call spline_eval(x_prof, u_prof, u2_prof, n_prof_pts, x(i), u_here)
                gu = (u_here-u0)/uw
                g(i) = amplitude*exp(-0.5_dp*gu*gu)
            end block
        end do
    end block

    do i = 2, npts-1
        Jphi(i) = -(( (g(i+1)-2.0_dp*g(i)+g(i-1))/dx**2 ) - 2.0_dp*g(i)/x(i)**2)
    end do
    Jphi(1) = 0.0_dp
    Jphi(npts) = 0.0_dp

    do i = 1, npts
        logg(i) = merge(log10(g(i)), -99.0_dp, g(i) > 1.0e-30_dp)
    end do

    open(newunit=unit_out, file="output/ch7_analytic_form.dat", status="replace", action="write")
    write(unit_out,'(A)') "# x   g(x)   log10(g(x))   Jphi(x)"
    do i = 1, npts
        write(unit_out, '(F12.8,3X,ES20.8E3,3X,F14.8,3X,ES20.8E3)') x(i), g(i), logg(i), Jphi(i)
    end do
    close(unit_out)
    print *, "Written to output/ch7_analytic_form.dat"

end program KonarCh7AnalyticForm
