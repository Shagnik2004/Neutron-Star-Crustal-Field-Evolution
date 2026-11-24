module tov_solver

    use utilities
    use tov_constants
    use pchip
    use log_linear_interp

    implicit none

contains

    function eos_P_of_rho(rho_table, P_table, d_table, n, rho_query) result(P_out)
        implicit none
        integer, intent(in)  :: n
        real(dp), intent(in) :: rho_table(n), P_table(n), d_table(n)
        real(dp), intent(in) :: rho_query
        real(dp) :: P_out
        call pchip_eval(rho_table, P_table, d_table, n, rho_query, P_out)
    end function eos_P_of_rho

    function eos_rho_of_P(rho_table, P_table, n, P_query) result(rho_out)
        implicit none
        integer, intent(in)  :: n
        real(dp), intent(in) :: rho_table(n), P_table(n)
        real(dp), intent(in) :: P_query
        real(dp) :: rho_out
        call loglinear_interpolate(P_table, rho_table, n, P_query, rho_out)
    end function eos_rho_of_P


    subroutine tov_rhs(r, P, m, rho_table, P_table, n, dPdr, dmdr)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: r, P, m
        real(dp), intent(in) :: rho_table(n), P_table(n)
        real(dp), intent(out):: dPdr, dmdr

        real(dp) :: rho

        rho = eos_rho_of_P(rho_table, P_table, n, P)

        dmdr = 4.0_dp * pi * r**2 * rho

        dPdr = -G_grav * (rho + P/c_light**2) * (m + 4.0_dp*pi*r**3*P/c_light**2) &
               / ( r**2 * (1.0_dp - 2.0_dp*G_grav*m/(r*c_light**2)) )

    end subroutine tov_rhs


    !--------------------------------------------------------------
    ! One Dormand-Prince RK45 (DOPRI5) trial step: 7 stages,
    ! producing both a 5th-order solution (P5,m5) and an embedded
    ! 4th-order solution (P4,m4). The difference between them is
    ! the local error estimate used for step-size control.
    ! success=.false. if any stage would need a non-positive P
    ! (matches this project's earlier pressure-positivity guards).
    !--------------------------------------------------------------
    subroutine rk45_step(r, P, m, h, rho_table, P_table, n, P5, m5, P4, m4, success)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: r, P, m, h
        real(dp), intent(in) :: rho_table(n), P_table(n)
        real(dp), intent(out):: P5, m5, P4, m4
        logical, intent(out) :: success

        real(dp), parameter :: c2=0.2_dp, c3=0.3_dp, c4=0.8_dp, c5=8.0_dp/9.0_dp, c6=1.0_dp, c7=1.0_dp

        real(dp), parameter :: a21=0.2_dp
        real(dp), parameter :: a31=3.0_dp/40.0_dp,     a32=9.0_dp/40.0_dp
        real(dp), parameter :: a41=44.0_dp/45.0_dp,    a42=-56.0_dp/15.0_dp,   a43=32.0_dp/9.0_dp
        real(dp), parameter :: a51=19372.0_dp/6561.0_dp, a52=-25360.0_dp/2187.0_dp, &
                                a53=64448.0_dp/6561.0_dp, a54=-212.0_dp/729.0_dp
        real(dp), parameter :: a61=9017.0_dp/3168.0_dp, a62=-355.0_dp/33.0_dp, &
                                a63=46732.0_dp/5247.0_dp, a64=49.0_dp/176.0_dp, a65=-5103.0_dp/18656.0_dp
        real(dp), parameter :: a71=35.0_dp/384.0_dp, a73=500.0_dp/1113.0_dp, a74=125.0_dp/192.0_dp, &
                                a75=-2187.0_dp/6784.0_dp, a76=11.0_dp/84.0_dp

        real(dp), parameter :: b4_1=5179.0_dp/57600.0_dp, b4_3=7571.0_dp/16695.0_dp, &
                                b4_4=393.0_dp/640.0_dp, b4_5=-92097.0_dp/339200.0_dp, &
                                b4_6=187.0_dp/2100.0_dp, b4_7=1.0_dp/40.0_dp

        real(dp) :: k1P,k1m, k2P,k2m, k3P,k3m, k4P,k4m, k5P,k5m, k6P,k6m, k7P,k7m
        real(dp) :: P_stage

        success = .false.

        if (P <= 0.0_dp) return
        call tov_rhs(r, P, m, rho_table, P_table, n, k1P, k1m)

        P_stage = P + h*a21*k1P
        if (P_stage <= 0.0_dp) return
        call tov_rhs(r+c2*h, P_stage, m+h*a21*k1m, rho_table, P_table, n, k2P, k2m)

        P_stage = P + h*(a31*k1P + a32*k2P)
        if (P_stage <= 0.0_dp) return
        call tov_rhs(r+c3*h, P_stage, m+h*(a31*k1m+a32*k2m), rho_table, P_table, n, k3P, k3m)

        P_stage = P + h*(a41*k1P + a42*k2P + a43*k3P)
        if (P_stage <= 0.0_dp) return
        call tov_rhs(r+c4*h, P_stage, m+h*(a41*k1m+a42*k2m+a43*k3m), rho_table, P_table, n, k4P, k4m)

        P_stage = P + h*(a51*k1P + a52*k2P + a53*k3P + a54*k4P)
        if (P_stage <= 0.0_dp) return
        call tov_rhs(r+c5*h, P_stage, m+h*(a51*k1m+a52*k2m+a53*k3m+a54*k4m), rho_table, P_table, n, k5P, k5m)

        P_stage = P + h*(a61*k1P + a62*k2P + a63*k3P + a64*k4P + a65*k5P)
        if (P_stage <= 0.0_dp) return
        call tov_rhs(r+c6*h, P_stage, m+h*(a61*k1m+a62*k2m+a63*k3m+a64*k4m+a65*k5m), &
                     rho_table, P_table, n, k6P, k6m)

        P5 = P + h*(a71*k1P + a73*k3P + a74*k4P + a75*k5P + a76*k6P)
        if (P5 <= 0.0_dp) return
        m5 = m + h*(a71*k1m + a73*k3m + a74*k4m + a75*k5m + a76*k6m)

        call tov_rhs(r+c7*h, P5, m5, rho_table, P_table, n, k7P, k7m)

        P4 = P + h*(b4_1*k1P + b4_3*k3P + b4_4*k4P + b4_5*k5P + b4_6*k6P + b4_7*k7P)
        m4 = m + h*(b4_1*k1m + b4_3*k3m + b4_4*k4m + b4_5*k5m + b4_6*k6m + b4_7*k7m)

        success = .true.

    end subroutine rk45_step


    !--------------------------------------------------------------
    ! Integrates outward using Dormand-Prince RK45 with embedded
    ! error control (matches scipy's solve_ivp method='RK45').
    ! Step size grows/shrinks based on the actual local truncation
    ! error estimate, not a hand-tuned heuristic -- it naturally
    ! shrinks in the steep crust near the surface. Surface radius
    ! and mass are linearly refined between the last accepted
    ! point and the first sub-surface point.
    !--------------------------------------------------------------
    subroutine integrate_star(rho_table, P_table, d_table, n, rho_c, &
                               R_star, M_star, n_steps_used)

        implicit none

        integer, intent(in)  :: n
        real(dp), intent(in) :: rho_table(n), P_table(n), d_table(n)
        real(dp), intent(in) :: rho_c
        real(dp), intent(out):: R_star, M_star
        integer, intent(out) :: n_steps_used

        real(dp) :: r, P, m
        real(dp) :: r_prev, P_prev, m_prev
        real(dp) :: P_surface
        real(dp) :: h, P5, m5, P4, m4
        real(dp) :: err_P, err_m, err_norm, factor, frac
        real(dp) :: r0

        real(dp), parameter :: rtol = 1.0e-8_dp
        real(dp), parameter :: atol = 1.0_dp
        real(dp), parameter :: h_min = 1.0e-8_dp
        real(dp), parameter :: h_max = 1.0e5_dp

        integer, parameter :: max_steps = 2000000
        integer :: istep
        logical :: success

        P_surface = P_table(1)

        r0 = 100.0_dp
        r  = r0
        P  = eos_P_of_rho(rho_table, P_table, d_table, n, rho_c)
        m  = (4.0_dp/3.0_dp) * pi * r0**3 * rho_c

        h = 1.0_dp   ! initial guess; RK45 error control will adjust it immediately

        r_prev = r
        P_prev = P
        m_prev = m

        istep = 0
        do while (P > P_surface .and. istep < max_steps)

            call rk45_step(r, P, m, h, rho_table, P_table, n, P5, m5, P4, m4, success)

            if (.not. success) then
                h = h * 0.5_dp
                if (h < h_min) then
                    print *, "WARNING: step size underflow, stopping."
                    exit
                end if
                cycle
            end if

            err_P = abs(P5-P4) / (atol + rtol*max(abs(P), abs(P5)))
            err_m = abs(m5-m4) / (atol + rtol*max(abs(m), abs(m5)))
            err_norm = sqrt(0.5_dp*(err_P**2 + err_m**2))

            if (err_norm <= 1.0_dp) then
                r_prev = r
                P_prev = P
                m_prev = m

                r = r + h
                P = P5
                m = m5
                istep = istep + 1

                if (err_norm > 1.0e-12_dp) then
                    factor = 0.9_dp * err_norm**(-0.2_dp)
                else
                    factor = 5.0_dp
                end if
                factor = min(5.0_dp, max(0.2_dp, factor))
                h = min(h_max, h*factor)
            else
                factor = 0.9_dp * err_norm**(-0.2_dp)
                factor = max(0.2_dp, factor)
                h = max(h_min, h*factor)
            end if

        end do

        if (P <= P_surface) then
            frac = (P_prev - P_surface) / (P_prev - P)
            R_star = r_prev + frac*(r - r_prev)
            M_star = m_prev + frac*(m - m_prev)
        else
            R_star = r
            M_star = m
        end if

        n_steps_used = istep

    end subroutine integrate_star



    !--------------------------------------------------------------
    ! Same physics as integrate_star, but also writes the full
    ! (r, rho) radial profile to a file as it integrates -- for
    ! plotting the internal structure of one specific star.
    !--------------------------------------------------------------
    subroutine integrate_star_profile(rho_table, P_table, d_table, n, rho_c, &
                                       R_star, M_star, n_steps_used, profile_file)
        implicit none
        integer, intent(in)  :: n
        real(dp), intent(in) :: rho_table(n), P_table(n), d_table(n)
        real(dp), intent(in) :: rho_c
        real(dp), intent(out):: R_star, M_star
        integer, intent(out) :: n_steps_used
        character(len=*), intent(in) :: profile_file
        real(dp) :: r, P, m
        real(dp) :: r_prev, P_prev, m_prev
        real(dp) :: P_surface
        real(dp) :: h, P5, m5, P4, m4
        real(dp) :: err_P, err_m, err_norm, factor, frac
        real(dp) :: r0, rho_here
        real(dp), parameter :: rtol = 1.0e-8_dp
        real(dp), parameter :: atol = 1.0_dp
        real(dp), parameter :: h_min = 1.0e-8_dp
        real(dp), parameter :: h_max = 1.0e5_dp
        integer, parameter :: max_steps = 2000000
        integer :: istep, unit_out
        logical :: success
        P_surface = P_table(1)
        r0 = 100.0_dp
        r  = r0
        P  = eos_P_of_rho(rho_table, P_table, d_table, n, rho_c)
        m  = (4.0_dp/3.0_dp) * pi * r0**3 * rho_c
        h = 1.0_dp
        r_prev = r
        P_prev = P
        m_prev = m
        open(newunit=unit_out, file=profile_file, status='replace', action='write')
        write(unit_out, '(A)') "# r[km]   rho[g/cm^3]   log10(rho)"
        rho_here = eos_rho_of_P(rho_table, P_table, n, P)
        write(unit_out, '(F14.6,3X,ES14.6,3X,F12.6)') r/1.0e5_dp, rho_here, log10(rho_here)
        istep = 0
        do while (P > P_surface .and. istep < max_steps)
            call rk45_step(r, P, m, h, rho_table, P_table, n, P5, m5, P4, m4, success)
            if (.not. success) then
                h = h * 0.5_dp
                if (h < h_min) exit
                cycle
            end if
            err_P = abs(P5-P4) / (atol + rtol*max(abs(P), abs(P5)))
            err_m = abs(m5-m4) / (atol + rtol*max(abs(m), abs(m5)))
            err_norm = sqrt(0.5_dp*(err_P**2 + err_m**2))
            if (err_norm <= 1.0_dp) then
                r_prev = r
                P_prev = P
                m_prev = m
                r = r + h
                P = P5
                m = m5
                istep = istep + 1
                rho_here = eos_rho_of_P(rho_table, P_table, n, P)
                write(unit_out, '(F14.6,3X,ES14.6,3X,F12.6)') r/1.0e5_dp, rho_here, log10(rho_here)
                if (err_norm > 1.0e-12_dp) then
                    factor = 0.9_dp * err_norm**(-0.2_dp)
                else
                    factor = 5.0_dp
                end if
                factor = min(5.0_dp, max(0.2_dp, factor))
                h = min(h_max, h*factor)
            else
                factor = max(0.2_dp, 0.9_dp*err_norm**(-0.2_dp))
                h = max(h_min, h*factor)
            end if
        end do
        close(unit_out)
        if (P <= P_surface) then
            frac = (P_prev - P_surface) / (P_prev - P)
            R_star = r_prev + frac*(r - r_prev)
            M_star = m_prev + frac*(m - m_prev)
        else
            R_star = r
            M_star = m
        end if
        n_steps_used = istep
    end subroutine integrate_star_profile

end module tov_solver
