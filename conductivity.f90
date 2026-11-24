module conductivity

    use utilities
    use cubic_spline
    use bps_conductivity_data
    implicit none

    real(dp), parameter :: pi     = 3.14159265358979323846_dp
    real(dp), parameter :: third  = 1.0_dp/3.0_dp
    real(dp), parameter :: sixth  = 1.0_dp/6.0_dp
    real(dp), parameter :: uc     = 2.0_dp*pi/9.0_dp
    real(dp), parameter :: coulog_imp = 2.0_dp

    logical, save :: initialized = .false.
    real(dp), save :: y2Z(n_tab), y2A(n_tab), y2I171(n_tab), y2I5000(n_tab)

contains

    subroutine init_conductivity()
        implicit none
        if (.not. initialized) then
            call spline_coeffs(u0, at_no,    n_tab, y2Z)
            call spline_coeffs(u0, at_wt,    n_tab, y2A)
            call spline_coeffs(u0, I171_tab, n_tab, y2I171)
            call spline_coeffs(u0, I5000_tab,n_tab, y2I5000)
            initialized = .true.
        end if
    end subroutine init_conductivity


   
    ! rholog: log10(density, g/cm^3)
    ! T8    : temperature in units of 10^8 K
    ! Q     : impurity 
    ! sigma : electrical conductivity (s^-1)
    
    subroutine sigma_ns(rholog, T8, Q, sigma)

        implicit none

        real(dp), intent(in)  :: rholog, T8, Q
        real(dp), intent(out) :: sigma

        real(dp) :: rho6, u, Z, A, Ilow, Ihigh
        real(dp) :: rho6_by_A, x_YU, gamma0
        real(dp) :: xsqr, coulog_ei
        real(dp) :: v, Isig, Aprime, gamma, G0, Fsig
        real(dp) :: sigma_pure, sigma_imp
        real(dp) :: rhodrip

        call init_conductivity()

        rhodrip = log10(4.3e11_dp)

        rho6 = 10.0_dp**(rholog - 6.0_dp)
        u    = uc*(rholog - 3.0_dp)

        call spline_eval(u0, at_no,     y2Z,     n_tab, u, Z)
        call spline_eval(u0, at_wt,     y2A,     n_tab, u, A)
        call spline_eval(u0, I171_tab,  y2I171,  n_tab, u, Ilow)
        call spline_eval(u0, I5000_tab, y2I5000, n_tab, u, Ihigh)

        rho6_by_A = rho6/A
        x_YU      = (Z*rho6_by_A)**third
        gamma0    = 0.2275_dp*Z*Z*(rho6_by_A**third)/T8

        if (gamma0 < 171.0_dp) then
            ! Liquid layer (Yakovlev & Urpin 1980)
            xsqr      = x_YU*x_YU
            coulog_ei = log(((4.0_dp*Z)**third)*sqrt(1.0_dp+(2.0_dp/gamma0))) &
                        - 0.5_dp*xsqr/(1.0_dp+xsqr)
            sigma = 8.53e21_dp*xsqr*x_YU/(Z*coulog_ei*(1.0_dp+xsqr))
            return
        end if

        ! Solid: phonon conductivity (Itoh et al. 1984)
        call get_v(gamma0, v)
        Isig = (1.0_dp-v)*Ilow + v*Ihigh  

        if (rholog <= rhodrip) then
            Aprime = A
        else
            Aprime = 4.0_dp*Z
        end if

        gamma = 0.3443_dp*(rho6_by_A**sixth)*gamma0/(Z*sqrt(Aprime))
        G0    = 13.0_dp/sqrt(1.0_dp+0.0174_dp*gamma*gamma)
        Fsig  = Isig*G0

        sigma_pure = 1.0_dp + 1.018_dp*x_YU*x_YU
        sigma_pure = 1.61049e21_dp*(x_YU**4.0_dp)/sigma_pure
        sigma_pure = sigma_pure/(Fsig*T8)

        sigma = sigma_pure

        if (Q > 0.0_dp) then
            sigma_imp = 8.53e21_dp*x_YU*x_YU*x_YU*Z &
                        /(coulog_imp*Q*(1.0_dp+x_YU*x_YU))
            sigma = 1.0_dp/((1.0_dp/sigma_pure)+(1.0_dp/sigma_imp))
        end if

    end subroutine sigma_ns

   
    ! Electron number density of neutron star crustal matter.
    ! rholog: log10(density, g/cm^3)
    ! n_e   : electron number density (cm^-3)
    ! Reuses the SAME Z(rho), A(rho) spline lookups as sigma_ns,
    ! via n_e = Z*rho / (A*m_u), the standard fully-ionized-matter
    ! relation. Needed for the Hall coefficient H(x) = c/(4*pi*e*n_e),
    ! since our conductivity model never previously exposed n_e
    ! directly (it stays internal to sigma_ns otherwise).

    subroutine get_electron_density(rholog, n_e)
        implicit none
        real(dp), intent(in)  :: rholog
        real(dp), intent(out) :: n_e
        real(dp), parameter :: m_u = 1.66053906660e-24_dp   ! atomic mass unit, g
        real(dp) :: rho, u, Z, A
        call init_conductivity()
        rho = 10.0_dp**rholog
        u   = uc*(rholog - 3.0_dp)
        call spline_eval(u0, at_no, y2Z, n_tab, u, Z)
        call spline_eval(u0, at_wt, y2A, n_tab, u, A)
        n_e = Z*rho / (A*m_u)
    end subroutine get_electron_density


  
    ! tabulated I_sigma values (Itoh et al. 1984 fitting form).

    subroutine get_v(gamma0, v)
        implicit none
        real(dp), intent(in)  :: gamma0
        real(dp), intent(out) :: v
        real(dp) :: x

        if (gamma0 <= 171.0_dp) then
            v = 0.0_dp
            return
        end if
        if (gamma0 >= 5000.0_dp) then
            v = 1.0_dp
            return
        end if

        x = 1.0_dp/(gamma0**(1.0_dp/3.0_dp))
        v = 0.6898_dp + 449.1520_dp/gamma0
        v = v + 14.4310_dp*x - 182.2730_dp*x*x

    end subroutine get_v

end module conductivity
