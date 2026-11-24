module t_of_rho

    use utilities
    implicit none

contains

    function surface_gravity(M, R, G_grav, c_light) result(gs)
        implicit none
        real(dp), intent(in) :: M, R, G_grav, c_light
        real(dp) :: gs
        gs = G_grav*M / (R*R*sqrt(1.0_dp - 2.0_dp*G_grav*M/(R*c_light**2)))
    end function surface_gravity


    
    ! Isothermal interior temperature from surface temperature and
    ! surface gravity (eq. 2.11, Gudmundsson, Pethick & Epstein 1983).
    ! Ts in K, gs in cm/s^2. Returns Tb in K.
   
    function T_interior(Ts, gs) result(Tb)
        implicit none
        real(dp), intent(in) :: Ts, gs
        real(dp) :: Tb
        Tb = 1.288e8_dp * ( ((Ts/1.0e6_dp)**4) / (gs/1.0e14_dp) )**0.455_dp
    end function T_interior


  
    function temperature_of_rho(rho, Tb, Ts, gs) result(T)
        implicit none
        real(dp), intent(in) :: rho, Tb, Ts, gs
        real(dp) :: T
        real(dp) :: T6, rho6, rho_boundary6, T8, T8_in

        T6   = Ts/1.0e6_dp
        rho6 = rho/1.0e6_dp

        rho_boundary6 = 0.7_dp*(T6**3.28_dp) / ((gs/1.0e14_dp)**1.152_dp)

        T8_in = Tb/1.0e8_dp

        if (rho6 < rho_boundary6) then
            T8 = T8_in*((rho6/rho_boundary6)**0.25_dp)
        else
            T8 = T8_in
        end if

        T = T8*1.0e8_dp

    end function temperature_of_rho

end module t_of_rho
