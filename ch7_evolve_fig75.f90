program Ch7EvolveFig75

    use utilities
    use tov_constants, only: c_light, pi
    use crust_model
    use g_profile_ch7_v2
    use g_diffusion_solver_be
    use t_of_rho, only: temperature_of_rho, T_interior
    use cooling_curve

    implicit none

    integer, parameter :: N = 2000
    integer, parameter :: n_Q = 5
    integer, parameter :: n_output = 400
    real(dp) :: Q_values(n_Q)

    real(dp) :: x(N), g(N)
    real(dp) :: x_c, x_s, dx
    real(dp), parameter :: x0 = 0.925_dp, sigma_right = 0.006_dp/2.3548_dp, amplitude = 9.5_dp

    real(dp) :: R_star_cm, g_surface_0
    real(dp) :: t_yr, log_t_min, log_t_max, Ts_now, Tb_now
    integer :: i, iQ, n_out, unit_out

    call init_crust_model("../../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)
    call init_cooling_curve("../crust_conductivity/FP_SF_B0.sparse")

    x_c = 0.90_dp
    x_s = 1.0_dp
    dx = (x_s - x_c) / real(N-1, dp)
    do i = 1, N
        x(i) = x_c + real(i-1,dp)*dx
    end do

    R_star_cm = crust_R_star_km * 1.0e5_dp

    Q_values = (/ 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp /)

    log_t_min = -1.0_dp
    log_t_max = 9.0_dp

    do iQ = 1, n_Q

        call g_profile_expelled_v2(x, g, N, x0, sigma_right, amplitude)
        g_surface_0 = g(N)

        block
            character(len=64) :: fname
            write(fname, '(A,I0,A)') "output/ch7_fig75_Q", iQ, ".dat"
            open(newunit=unit_out, file=trim(fname), status="replace", action="write")
        end block
        write(unit_out, '(A)') "# log10(t/yr)   g_surface/g_surface_0"

        t_yr = 10.0_dp**log_t_min
        write(unit_out, '(F10.4,3X,ES16.8)') log_t_min, 1.0_dp

        do n_out = 1, n_output
            block
                real(dp) :: t_target, t_sub_start, t_sub_end, dt_sub_yr, dt_sec, t_mid
                real(dp) :: max_frac_step, frac
                real(dp) :: S_local(N), g_trial(N), rho_here, T_here, sigma_here, r_here
                integer :: j

                t_target = 10.0_dp**(log_t_min + (log_t_max-log_t_min)*real(n_out,dp)/real(n_output,dp))
                if (t_target <= t_yr) cycle

                frac = g(N)/g_surface_0
                max_frac_step = merge(0.02_dp, 0.005_dp, frac > 0.3_dp .or. frac < 3.0_dp)

                do while (t_yr < t_target)
                    t_sub_start = t_yr
                    t_sub_end = min(t_target, t_sub_start*(1.0_dp+max_frac_step))
                    dt_sub_yr = t_sub_end - t_sub_start
                    if (dt_sub_yr <= 0.0_dp) exit
                    dt_sec = dt_sub_yr * 365.25_dp * 86400.0_dp

                    t_mid = sqrt(max(t_sub_start,1.0e-3_dp) * t_sub_end)
                    Ts_now = Ts_at_age(t_mid)
                    Tb_now = T_interior(Ts_now, crust_gs)

                    do j = 1, N
                        r_here = x(j)*crust_R_star_km
                        rho_here = get_rho_at_r(r_here)
                        T_here = temperature_of_rho(rho_here, Tb_now, Ts_now, crust_gs)
                        sigma_here = get_sigma_given_T(rho_here, T_here, Q_values(iQ))
                        S_local(j) = c_light**2 / (4.0_dp*pi*sigma_here*R_star_cm**2)
                    end do

                    g_trial = g
                    call be_step(g_trial, x, N, dx, x_s, S_local, dt_sec)
                    g = g_trial
                    t_yr = t_sub_end
                end do

                write(unit_out, '(F10.4,3X,ES16.8)') log10(t_yr), g(N)/g_surface_0
            end block
        end do

        close(unit_out)
        print '(A,F4.2,A,ES14.6)', " Q=", Q_values(iQ), "  final g_surface/g0 = ", g(N)/g_surface_0

    end do

end program Ch7EvolveFig75
