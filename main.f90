program Main

    use utilities
    use tov_constants
    use pchip
    use eos_reader
    use tov_solver

    implicit none

    real(dp), allocatable :: rho(:), P(:), d(:)
    integer :: n

    !---------------- Mass-radius sweep ----------------
    integer, parameter :: n_sweep = 10000

    real(dp) :: rho_c(n_sweep)
    real(dp) :: R_star(n_sweep)
    real(dp) :: M_star(n_sweep)
    integer  :: steps_used(n_sweep)
    logical  :: extrapolated(n_sweep)

    real(dp) :: rho_c_min, rho_c_max, log_min, log_max
    real(dp) :: table_max
    integer  :: i, unit_out, idx_max

    !---------------- 1.4 Msun profile (bisection) ----------------
    real(dp) :: rho_lo, rho_hi, rho_mid
    real(dp) :: R_prof, M_prof, M_target
    integer  :: n_steps_prof, iter

    call read_eos_table("../Interpolation_Project/data/eos_akmalpr.dat", rho, P, n)

    allocate(d(n))
    call pchip_coeffs(rho, P, n, d)

    table_max = rho(n)

    !======================================================
    ! PART 1: Mass-radius sweep
    !======================================================
    rho_c_min = 1.0e14_dp
    rho_c_max = 1.0e16_dp

    log_min = log10(rho_c_min)
    log_max = log10(rho_c_max)

    print *
    print '(A,I0,A)', " Sweeping ", n_sweep, " central densities (stable branch range)"

    open(newunit=unit_out, file="output/mass_radius.dat", status="replace", action="write")
    write(unit_out, '(A)') "# rho_c[g/cm^3]   R[km]   M[Msun]   steps   extrapolated(1=yes)"

    do i = 1, n_sweep

        rho_c(i) = 10.0_dp**( log_min + (log_max-log_min)*real(i-1,dp)/real(n_sweep-1,dp) )
        extrapolated(i) = (rho_c(i) > table_max)

        call integrate_star(rho, P, d, n, rho_c(i), R_star(i), M_star(i), steps_used(i))

        write(unit_out, '(ES14.6,3X,F12.4,3X,F12.4,3X,I8,3X,I2)') &
              rho_c(i), R_star(i)/1.0e5_dp, M_star(i)/Msun, steps_used(i), merge(1,0,extrapolated(i))

    end do

    close(unit_out)

    idx_max = 1
    do i = 1, n_sweep
        if (.not. extrapolated(i)) then
            idx_max = i
            exit
        end if
    end do
    do i = idx_max+1, n_sweep
        if (.not. extrapolated(i) .and. M_star(i) > M_star(idx_max)) idx_max = i
    end do

    print '(A,ES14.6,A)', " rho_c(M_max) = ", rho_c(idx_max), " g/cm^3"
    print '(A,F10.4,A)',  " M_max        = ", M_star(idx_max)/Msun, " solar masses"
    print '(A,F10.4,A)',  " R(M_max)     = ", R_star(idx_max)/1.0e5_dp, " km"
    print *, "Sweep written to output/mass_radius.dat"

    !======================================================
    ! PART 2: Find rho_c for M = 1.4 Msun (bisection), write profile
    !======================================================
    M_target = 1.4_dp * Msun
    rho_lo = 5.0e14_dp
    rho_hi = 1.2e15_dp

    do iter = 1, 50
        rho_mid = 0.5_dp*(rho_lo+rho_hi)
        call integrate_star(rho, P, d, n, rho_mid, R_prof, M_prof, n_steps_prof)

        if (M_prof < M_target) then
            rho_lo = rho_mid
        else
            rho_hi = rho_mid
        end if

        if (abs(M_prof - M_target)/M_target < 1.0e-6_dp) exit
    end do

    print *
    print '(A,ES14.6,A)', " Found rho_c = ", rho_mid, " g/cm^3  for M = 1.4 Msun"
    print '(A,F10.6,A)',  "   M = ", M_prof/Msun, " Msun"
    print '(A,F10.4,A)',  "   R = ", R_prof/1.0e5_dp, " km"

    call integrate_star_profile(rho, P, d, n, rho_mid, R_prof, M_prof, n_steps_prof, &
                                 "output/density_profile_1p4Msun.dat")

    print *, "Profile written to output/density_profile_1p4Msun.dat"
    print *, "Test completed successfully."
    print *

    deallocate(rho, P, d)

end program Main
