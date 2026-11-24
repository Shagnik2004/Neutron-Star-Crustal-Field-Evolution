module magnetized_conductivity

    use utilities
    use cubic_spline
    use bps_conductivity_data
    use conductivity, only: uc, y2Z, y2A, sigma_ns

    implicit none

    ! Physical constants (Gaussian/cgs)
    real(dp), parameter :: e_charge = 4.803204e-10_dp   ! esu
    real(dp), parameter :: m_electron = 9.109384e-28_dp ! g
    real(dp), parameter :: m_u = 1.660539e-24_dp        ! atomic mass unit, g
    real(dp), parameter :: c_light_local = 2.99792458e10_dp

contains

    !--------------------------------------------------------------
    ! Perpendicular conductivity in a magnetized crust:
    !   sigma_perp = sigma_0 / (1 + (omega_B*tau)^2)
    ! where omega_B = eB/(m_e c) (electron cyclotron frequency),
    ! tau = sigma_0*m_e/(n_e*e^2) (collision time, from the
    ! FIELD-FREE conductivity sigma_0 already validated elsewhere),
    ! n_e = Z*rho/(A*m_u) from the existing BPS composition table.
    !--------------------------------------------------------------
    function sigma_perpendicular(rho, T8, Q, B_gauss) result(sigma_perp)
        implicit none
        real(dp), intent(in) :: rho, T8, Q, B_gauss
        real(dp) :: sigma_perp

        real(dp) :: rholog, u, Z_val, A_val, n_e, sigma_0, tau, omega_B, hall_param

        rholog = log10(rho)
        u = uc*(rholog - 3.0_dp)
        call spline_eval(u0, at_no, y2Z, n_tab, u, Z_val)
        call spline_eval(u0, at_wt, y2A, n_tab, u, A_val)

        n_e = Z_val*rho / (A_val*m_u)

        call sigma_ns(rholog, T8, Q, sigma_0)

        tau = sigma_0*m_electron / (n_e*e_charge**2)
        omega_B = e_charge*abs(B_gauss) / (m_electron*c_light_local)
        hall_param = omega_B*tau

        sigma_perp = sigma_0 / (1.0_dp + hall_param**2)

    end function sigma_perpendicular

end module magnetized_conductivity
