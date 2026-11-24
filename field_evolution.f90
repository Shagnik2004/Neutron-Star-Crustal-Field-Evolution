module field_evolution

    use utilities
    use tov_constants, only: c_light, pi
    use crust_model
    use g_diffusion_solver
    use t_of_rho, only: temperature_of_rho, T_interior
    use cooling_curve

    implicit none

    integer, parameter :: TEMP_ISOTHERMAL = 1     ! flat T_i(t) everywhere (Bhattacharya & Datta style)
    integer, parameter :: TEMP_POSITION_DEPENDENT = 2   ! Konar "standard cooling" style

contains

    !--------------------------------------------------------------
    ! Master field-evolution subroutine. Given ANY initial g(x)
    ! array on a caller-supplied grid, evolves the surface value
    ! g(x_s,t)/g(x_s,0) via the validated backward-Euler solver.
    !
    ! Assumes init_crust_model and init_cooling_curve have ALREADY
    ! been called by the caller (so the same star/cooling curve can
    ! be reused across multiple evolutions without re-reading files).
    !
    ! Inputs:
    !   x(N), g0(N)     : grid and initial condition (caller-built,
    !                     any functional form -- this routine makes
    !                     no assumption about its shape)
    !   x_s             : outer boundary (Robin BC applied here)
    !   Q               : impurity parameter
    !   temp_mode       : TEMP_ISOTHERMAL or TEMP_POSITION_DEPENDENT
    !   log_t_min/max   : evolution range, log10(yr)
    !   n_output        : number of output points (log-spaced)
    !   decay_floor     : stop early if g_s/g_s0 drops below this
    !                     (avoids grinding through numerically
    !                     meaningless deep-decay regime)
    !
    ! Outputs:
    !   log_t_out(n_output), ratio_out(n_output), n_written
    !--------------------------------------------------------------
    subroutine evolve_surface_field(x, g0, N, x_s, Q, temp_mode, &
                                     log_t_min, log_t_max, n_output, &
                                     log_t_out, ratio_out, n_written, &
                                     decay_floor_in)

        implicit none

        integer, intent(in)  :: N, temp_mode, n_output
        real(dp), intent(in) :: x(N), g0(N), x_s, Q
        real(dp), intent(in) :: log_t_min, log_t_max
        real(dp), intent(out) :: log_t_out(n_output), ratio_out(n_output)
        integer, intent(out)  :: n_written
        real(dp), intent(in), optional :: decay_floor_in

        real(dp) :: g(N), dx, R_star_cm, g_surface_0, decay_floor
        real(dp) :: t_yr, Ts_now, Tb_now
        integer :: i, n_out

        decay_floor = 1.0e-9_dp
        if (present(decay_floor_in)) decay_floor = decay_floor_in

        g = g0
        g_surface_0 = g0(N)
        dx = x(2) - x(1)   ! assumes uniform grid, as built by every driver so far
        R_star_cm = crust_R_star_km * 1.0e5_dp

        t_yr = 10.0_dp**log_t_min
        n_written = 0

        do n_out = 1, n_output

            if (g(N)/g_surface_0 < decay_floor) exit

            block
                real(dp) :: t_target, t_sub_start, t_sub_end, dt_sub_yr, dt_sec, t_mid
                real(dp) :: max_frac_step, frac
                real(dp) :: S_local(N), g_trial(N), rho_here, T_here, sigma_here, r_here
                integer :: j

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

                    ! Second, physics-aware cap: shrink the step further
                    ! if Ts(t) itself is changing too fast within it --
                    ! the age-growth cap alone has no awareness of this,
                    ! and freezing conductivity across a step where Ts
                    ! collapses by orders of magnitude produces genuine
                    ! numerical error (found via instability_scan.f90,
                    ! confirmed independent of CN vs BE).
                    block
                        real(dp) :: Ts_step_start, Ts_step_end
                        real(dp), parameter :: max_frac_T_change = 0.05_dp
                        integer :: shrink_iter
                        Ts_step_start = Ts_at_age(max(t_sub_start,1.0e-3_dp))
                        do shrink_iter = 1, 40
                            Ts_step_end = Ts_at_age(t_sub_end)
                            if (abs(Ts_step_end-Ts_step_start) <= &
                                max_frac_T_change*Ts_step_start) exit
                            t_sub_end = t_sub_start + 0.5_dp*(t_sub_end-t_sub_start)
                            if (t_sub_end-t_sub_start < 1.0e-10_dp*max(t_sub_start,1.0e-3_dp)) exit
                        end do
                    end block

                    dt_sub_yr = t_sub_end - t_sub_start
                    if (dt_sub_yr <= 0.0_dp) exit
                    dt_sec = dt_sub_yr * 365.25_dp * 86400.0_dp

                    t_mid = sqrt(max(t_sub_start,1.0e-3_dp) * t_sub_end)
                    Ts_now = Ts_at_age(t_mid)
                    Tb_now = T_interior(Ts_now, crust_gs)

                    do j = 1, N
                        r_here = x(j)*crust_R_star_km
                        rho_here = get_rho_at_r(r_here)
                        if (temp_mode == TEMP_POSITION_DEPENDENT) then
                            T_here = temperature_of_rho(rho_here, Tb_now, Ts_now, crust_gs)
                        else
                            T_here = Tb_now
                        end if
                        sigma_here = get_sigma_given_T(rho_here, T_here, Q)
                        S_local(j) = c_light**2 / (4.0_dp*pi*sigma_here*R_star_cm**2)
                    end do

                    g_trial = g
                    call gcn_step(g_trial, x, N, dx, x_s, S_local, dt_sec, theta_in=0.55_dp)
                    g = g_trial
                    t_yr = t_sub_end

                    if (g(N)/g_surface_0 < decay_floor) exit
                end do

                n_written = n_written + 1
                log_t_out(n_written) = log10(t_yr)
                ratio_out(n_written) = g(N)/g_surface_0
            end block
        end do

    end subroutine evolve_surface_field

end module field_evolution
