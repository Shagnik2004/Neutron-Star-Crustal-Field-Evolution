program DemoStellarStructure

    use utilities
    use stellar_structure

    implicit none

    real(dp) :: R_star_km, M_achieved, rho_c

    call solve_star("../Interpolation_Project/data/eos_akmalpr.dat", 1.4_dp, &
                     "output/density_profile_1p4Msun.dat", &
                     R_star_km, M_achieved, rho_c)

    print '(A,F10.4,A)', " R_star = ", R_star_km, " km"
    print '(A,F10.6,A)', " M_achieved = ", M_achieved, " Msun"
    print '(A,ES14.6)',  " rho_c = ", rho_c

end program DemoStellarStructure
