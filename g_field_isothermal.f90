module g_field_isothermal

    use utilities
    use tov_constants, only: c_light
    use crust_model
    use conductivity, only: pi

    implicit none

contains

    !--------------------------------------------------------------
    ! S(x) for an ISOTHERMAL crust: single temperature T_i applied
    ! at every grid point (Bhattacharya & Datta 1996, eq. 9 -- same
    ! formula as Konar's eq. 2.11, but with NO position-dependent
    ! falloff -- unlike Konar's temperature_of_rho, which cools the
    ! outer crust faster than the interior).
    !--------------------------------------------------------------
    subroutine compute_S_isothermal(x, N, R_star_cm, T_i, Q_imp, S)
        implicit none
        integer, intent(in)   :: N
        real(dp), intent(in)  :: x(N), R_star_cm, T_i, Q_imp
        real(dp), intent(out) :: S(N)

        real(dp) :: rho_here, sigma_here, r_here
        integer :: j

        do j = 1, N
            r_here = x(j)*crust_R_star_km
            rho_here = get_rho_at_r(r_here)
            sigma_here = get_sigma_given_T(rho_here, T_i, Q_imp)
            S(j) = c_light**2 / (4.0_dp*pi*sigma_here*R_star_cm**2)
        end do

    end subroutine compute_S_isothermal

end module g_field_isothermal
