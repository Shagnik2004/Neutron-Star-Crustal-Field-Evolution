module stellar_structure

    use utilities
    use tov_constants
    use pchip
    use eos_reader
    use tov_solver

    implicit none

contains

    !--------------------------------------------------------------
    ! Solves for a star of mass M_target_Msun using the given EOS,
    ! writes its full radial density profile to output_profile_file.
    ! Wraps the already-validated bisection + integrate_star_profile
    ! machinery -- no new physics, just a clean callable interface.
    !--------------------------------------------------------------
    subroutine solve_star(eos_file, M_target_Msun, output_profile_file, &
                           R_star_km, M_achieved_Msun, rho_c_used)

        implicit none

        character(len=*), intent(in)  :: eos_file, output_profile_file
        real(dp), intent(in)          :: M_target_Msun
        real(dp), intent(out)         :: R_star_km, M_achieved_Msun, rho_c_used

        real(dp), allocatable :: rho(:), P(:), d(:)
        integer :: n

        real(dp) :: rho_lo, rho_hi, rho_mid, R_star, M_star, M_target
        integer :: n_steps, iter

        call read_eos_table(eos_file, rho, P, n)

        allocate(d(n))
        call pchip_coeffs(rho, P, n, d)

        M_target = M_target_Msun * Msun

        ! Bisection bracket: full EOS table range (minus the very
        ! endpoints, to stay safely inside PCHIP's valid domain).
        rho_lo = rho(2)
        rho_hi = rho(n-1)

        do iter = 1, 60
            rho_mid = 0.5_dp*(rho_lo+rho_hi)
            call integrate_star(rho, P, d, n, rho_mid, R_star, M_star, n_steps)

            if (M_star < M_target) then
                rho_lo = rho_mid
            else
                rho_hi = rho_mid
            end if

            if (abs(M_star-M_target)/M_target < 1.0e-6_dp) exit
        end do

        rho_c_used = rho_mid

        call integrate_star_profile(rho, P, d, n, rho_mid, R_star, M_star, &
                                     n_steps, output_profile_file)

        R_star_km = R_star/1.0e5_dp
        M_achieved_Msun = M_star/Msun

        deallocate(rho, P, d)

    end subroutine solve_star

end module stellar_structure
