program ExpelledEvolveQ

    use utilities
    use crust_model
    use g_profile_expelled
    use g_diffusion_solver_be
    use g_field_isothermal
    use t_of_rho, only: T_interior
    use cooling_curve

    implicit none

    integer, parameter :: N = 2000
    integer, parameter :: n_Q = 4
    integer, parameter :: n_output = 500

    real(dp) :: x(N), g(N)
    real(dp) :: x_c, x_s, dx
    real(dp), parameter :: xw = 0.01_dp
    real(dp), parameter :: k_offset = 2.0_dp
    real(dp) :: x0

    real(dp) :: Q_values(n_Q)
    real(dp) :: R_star_cm, g_surface_0
    real(dp) :: t_yr, log_t_min, log_t_max
    real(dp) :: Ts_now, Ti_now
    real(dp) :: rho_at_x0, crust_thickness_km
    integer :: i, iQ, n_out, unit_out, unit_meta

    call init_crust_model("../../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)
    call init_cooling_curve("../crust_conductivity/FP_SF_B0.sparse")

    x_c = get_x_crust_core()
    x_s = 1.0_dp
    x0  = x_c + k_offset*xw
    dx  = (x_s - x_c) / real(N-1, dp)
    do i = 1, N
        x(i) = x_c + real(i-1,dp)*dx
    end do

    R_star_cm = crust_R_star_km * 1.0e5_dp
    crust_thickness_km = crust_R_star_km - get_crust_core_radius()

    print '(A,ES14.6,A)', " Crust-core boundary density used: ", rho_crust_core, " g/cm^3"
    print '(A,F10.6,A,F10.6)', " x_c (crust bottom) = ", x_c, "   x0 = ", x0
    print '(A,F10.4,A)', " crust thickness = ", crust_thickness_km, " km"

    rho_at_x0 = get_rho_at_r(x0*crust_R_star_km)
    print '(A,ES14.6)', " rho at x0 = ", rho_at_x0

    open(newunit=unit_meta, file="output/expelled_meta_Q.dat", status="replace", action="write")
    write(unit_meta, '(ES16.8)') rho_at_x0
    write(unit_meta, '(F12.6)') crust_thickness_km
    close(unit_meta)

    Q_values = (/ 0.001_dp, 0.003_dp, 0.01_dp, 0.03_dp /)

    log_t_min = -1.0_dp
    log_t_max = 10.0_dp

    do iQ = 1, n_Q

        call g_profile_datta(x, g, N, dx, x0, xw, x_c)
        g_surface_0 = g(N)

        block
            character(len=64) :: fname
            write(fname, '(A,I0,A)') "output/expelled_Q", iQ, ".dat"
            open(newunit=unit_out, file=trim(fname), status="replace", action="write")
        end block
        write(unit_out, '(A)') "# log10(t/yr)   Bs/B0"

        t_yr = 10.0_dp**log_t_min
        write(unit_out, '(F10.4,3X,ES16.8)') log_t_min, 1.0_dp

        do n_out = 1, n_output

            if (g(N)/g_surface_0 < 1.0e-9_dp) exit

            block
                real(dp) :: t_target, t_sub_start, t_sub_end, dt_sub_yr, dt_sec, t_mid
                real(dp) :: max_frac_step, frac
                real(dp) :: g_trial(N), S_local(N)

                t_target = 10.0_dp**(log_t_min + (log_t_max-log_t_min)*real(n_out,dp)/real(n_output,dp))
                if (t_target <= t_yr) cycle

                frac = g(N)/g_surface_0
                if (frac > 1.0e-3_dp) then
                    max_frac_step = 0.02_dp
                else
                    max_frac_step = 0.002_dp
                end if

                do while (t_yr < t_target)
                    t_sub_start = t_yr
                    t_sub_end = min(t_target, t_sub_start*(1.0_dp+max_frac_step))
                    dt_sub_yr = t_sub_end - t_sub_start
                    if (dt_sub_yr <= 0.0_dp) exit

                    t_mid = sqrt(max(t_sub_start,1.0e-3_dp) * t_sub_end)
                    Ts_now = Ts_at_age(t_mid)
                    Ti_now = T_interior(Ts_now, crust_gs)
                    call compute_S_isothermal(x, N, R_star_cm, Ti_now, Q_values(iQ), S_local)

                    dt_sec = dt_sub_yr * 365.25_dp * 86400.0_dp
                    g_trial = g
                    call be_step(g_trial, x, N, dx, x_s, S_local, dt_sec)
                    g = g_trial
                    t_yr = t_sub_end

                    if (g(N)/g_surface_0 < 1.0e-9_dp) exit
                end do

                write(unit_out, '(F10.4,3X,ES16.8)') log10(t_yr), g(N)/g_surface_0
            end block
        end do

        close(unit_out)
        print '(A,I0,A,F8.4,A,ES16.8)', " Q=", iQ, " (", Q_values(iQ), ")  final Bs/B0 = ", g(N)/g_surface_0

    end do

end program ExpelledEvolveQ
