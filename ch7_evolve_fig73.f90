program Ch7EvolveFig73

    use utilities
    use tov_constants, only: c_light, pi
    use crust_model
    use g_profile_ch7
    use g_diffusion_solver_be
    use t_of_rho, only: temperature_of_rho, T_interior
    use cooling_curve

    implicit none

    integer, parameter :: N = 2000
    integer, parameter :: n_snap = 9
    real(dp), parameter :: Q_imp = 0.1_dp

    real(dp) :: x(N), g(N)
    real(dp) :: x_c, x_s, dx
    real(dp), parameter :: x0 = 0.925_dp, fwhm = 0.006_dp, amplitude = 9.5_dp
    real(dp), parameter :: x_bg = 0.916_dp, w_bg = 0.003_dp

    real(dp) :: R_star_cm
    real(dp) :: t_yr, Ts_now, Tb_now
    real(dp) :: snap_times(n_snap)
    integer :: i, isnap, unit_out

    call init_crust_model("../../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)
    call init_cooling_curve("../crust_conductivity/FP_SF_B0.sparse")

    ! NOTE: x_c set to 0.90 pragmatically here, matching the visual
    ! extent of the published Fig 7.1/7.3 (not re-derived from first
    ! principles) -- the EOS-transition-sourced boundary (~0.924)
    ! would cut off the post-peak dip entirely.
    x_c = 0.90_dp
    x_s = 1.0_dp
    dx = (x_s - x_c) / real(N-1, dp)
    do i = 1, N
        x(i) = x_c + real(i-1,dp)*dx
    end do

    R_star_cm = crust_R_star_km * 1.0e5_dp

    call g_profile_expelled_ch7(x, g, N, x0, fwhm, amplitude, x_bg, w_bg)

    print '(A,F10.6,A,F10.6)', " x_c = ", x_c, "   x0 = ", x0
    print '(A,ES14.6)', " g(surface, t=0) = ", g(N)

    snap_times = (/ 10.0_dp, 1.0e2_dp, 1.0e3_dp, 1.0e4_dp, 1.0e5_dp, &
                    1.0e6_dp, 1.0e7_dp, 1.0e8_dp, 1.0e9_dp /)

    t_yr = 0.0_dp

    do isnap = 1, n_snap

        do while (t_yr < snap_times(isnap))
            block
                real(dp) :: t_sub_start, t_sub_end, dt_sub_yr, dt_sec, t_mid
                real(dp), parameter :: max_frac_step = 0.02_dp
                real(dp) :: S_local(N), g_trial(N), rho_here, T_here, sigma_here, r_here
                integer :: j

                t_sub_start = t_yr
                t_sub_end = min(snap_times(isnap), max(t_sub_start*(1.0_dp+max_frac_step), t_sub_start+1.0e-6_dp))
                dt_sub_yr = t_sub_end - t_sub_start
                if (dt_sub_yr <= 0.0_dp) exit
                dt_sec = dt_sub_yr * 365.25_dp * 86400.0_dp

                t_mid = sqrt(max(t_sub_start,1.0e-3_dp) * t_sub_end)
                Ts_now = Ts_at_age(t_mid)
                Tb_now = T_interior(Ts_now, crust_gs)

                ! Position-dependent "standard cooling" (matches
                ! Konar's own Ch.5 approach), NOT the flat isothermal
                ! temperature used for yesterday's Bhattacharya-Datta work.
                do j = 1, N
                    r_here = x(j)*crust_R_star_km
                    rho_here = get_rho_at_r(r_here)
                    T_here = temperature_of_rho(rho_here, Tb_now, Ts_now, crust_gs)
                    sigma_here = get_sigma_given_T(rho_here, T_here, Q_imp)
                    S_local(j) = c_light**2 / (4.0_dp*pi*sigma_here*R_star_cm**2)
                end do

                g_trial = g
                call be_step(g_trial, x, N, dx, x_s, S_local, dt_sec)
                g = g_trial
                t_yr = t_sub_end
                if (t_yr > 1.0e8_dp) then
                    print '(A,ES14.6,A,ES14.6)', "   t_yr=", t_yr, "  g(surface)=", g(N)
                end if
            end block
        end do

        block
            character(len=64) :: fname
            write(fname, '(A,I0,A)') "output/ch7_fig73_snap_", isnap, ".dat"
            open(newunit=unit_out, file=trim(fname), status="replace", action="write")
            write(unit_out, '(A,ES14.6)') "# x   g(x)   t_yr=", t_yr
            block
                integer :: j
                do j = 1, N
                    write(unit_out, '(F12.8,3X,ES16.8)') x(j), g(j)
                end do
            end block
            close(unit_out)
        end block

        print '(A,I0,A,ES12.4,A,ES14.6)', " Snapshot ", isnap, ": t=", t_yr, &
              " yr, g(surface)=", g(N)

    end do

end program Ch7EvolveFig73
