program ProfileMain

    use utilities
    use tov_constants
    use pchip
    use eos_reader
    use tov_solver

    implicit none

    real(dp), allocatable :: rho(:), P(:), d(:)
    integer :: n

    real(dp) :: rho_lo, rho_hi, rho_mid
    real(dp) :: R_star, M_star, M_target
    integer  :: n_steps
    integer  :: iter

    call read_eos_table("../Interpolation_Project/data/eos_akmalpr.dat", rho, P, n)

    allocate(d(n))
    call pchip_coeffs(rho, P, n, d)

    M_target = 1.4_dp * Msun

    rho_lo = 5.0e14_dp
    rho_hi = 1.2e15_dp

    do iter = 1, 50
        rho_mid = 0.5_dp*(rho_lo+rho_hi)
        call integrate_star(rho, P, d, n, rho_mid, R_star, M_star, n_steps)

        if (M_star < M_target) then
            rho_lo = rho_mid
        else
            rho_hi = rho_mid
        end if

        if (abs(M_star - M_target)/M_target < 1.0e-6_dp) exit
    end do

    print *
    print '(A,ES14.6,A)', " Found rho_c = ", rho_mid, " g/cm^3  for M = 1.4 Msun"
    print '(A,F10.6,A)',  "   M = ", M_star/Msun, " Msun"
    print '(A,F10.4,A)',  "   R = ", R_star/1.0e5_dp, " km"

    call integrate_star_profile(rho, P, d, n, rho_mid, R_star, M_star, n_steps, &
                                 "output/density_profile_1p4Msun.dat")

    print *
    print *, "Profile written to output/density_profile_1p4Msun.dat"
    print *

    deallocate(rho, P, d)

end program ProfileMain
