program Ns2dHall
    ! ======================================================================
    ! 2D (r,theta) crustal field evolution, Hall + Ohmic,
    ! Nx = Ntheta = 100.  Written from the equations of Skiathas &
    ! Gourgouliatos (2024):
    !
    !   db/dt = curl A,   A = -D curl b  -  G (curl b) x b            (eqs 9,17,19)
    !   b = grad P x grad phi + T grad phi ,  x = r/R
    !   dP/dt = r sin(th) A_phi                                        (eq 20)
    !   dT/dt = sin(th) A_th + r sin(th) dA_th/dr - sin(th) dA_r/dth  (eq 21)
    !
    ! The scalar equations below were obtained by expanding these with
    ! computer algebra (sympy) and are exact, with S(x)=D(x) (Ohmic) and
    ! H(x)=G(x) (Hall), L*P = P_xx + (P_thth - cot(th) P_th)/x^2 :
    !
    !   P_t = S L*P + (H/(x^2 sin)) ( T_th P_x - T_x P_th )
    !   T_t = d/dx( S T_x ) + (S/x^2)( T_thth - cot(th) T_th )         (Ohmic, conservative)
    !         + [ sin A_th + x sin dA_th/dx - sin dA_r/dth ]_(P only)
    !         + (T/(x^2 sin)) [ 2 H cot(th) T_x + (H' - 2H/x) T_th ]
    !   with  A_r  = (H/(x^2 sin^2)) L*P P_x ,   A_th = (H/(x^3 sin^2)) L*P P_th
    !
    ! There is NO T*T_{x,theta} term (it cancels exactly in the curl).
    !
    ! Numerics:
    !   * Ohmic terms: implicit ADI (Peaceman-Rachford type, theta_adi=0.55);
    !     P with the validated routine adi_step_2d_hall_omp (zero velocities),
    !     T with the conservative routine adi_step_2d_T_flux_omp.
    !   * Hall terms: explicit, variable-step Adams-Bashforth 3 (Lagrange
    !     extrapolant integrals, so changing steps are exact), restarted
    !     (Euler, AB2, AB3) after every injection.  Advective derivatives are
    !     upwinded (UPW_ORDER=1 first order, 2 = minmod-limited second order).
    !   * dt <= DT_SAFETY * 0.5 dx^2 / (H_max max(|P|,|T|)); a step-rejection
    !     net (NaN, >1e30, field more than doubling) only guards against
    !     accidents.
    !
    ! Physics inputs: BPS/Potekhin crust model,
    ! conductivity with impurity Q, pasta conductivity / 1000,
    ! repeated flux injection at the pasta midpoint x0, optional crustal field.
    !
    ! Domain: x = 0.923897 (pasta inner edge = crust-core boundary) to 1.0.
    ! Inner boundary P = T = 0 (perfectly conducting core wall), as in
    ! Gourgouliatos & Cumming (2014).  To use the older domain starting at
    ! x = 0.90, set x_c = 0.90_dp in run_one.
    !
    ! Environment (all optional):
    !   RUN_IDX     run only the n-th combination (0 = all)       run = 2*(iQ-1)+iinit
    !   LOGT_MAX    evolve to log10(t/yr)                          (6.0)
    !   INJ_PERIOD  injection period in yr                         (10)
    !   AMP_LOG     log10 g_max of one injection                   (0.8)
    !   DT_SAFETY   multiplier on the whistler step estimate       (1.6)
    !   RHO_CUT     log10 density below which Hall is off          (10.0)
    !   UPW_ORDER   1 or 2                                         (1)
    !   OUTDIR      output directory
    ! ======================================================================
    use utilities
    use tov_constants, only: c_light, pi
    use crust_model
    use cooling_curve
    use g_profile_library
    use psi_hall_omp
    use psi_flux_omp
    use t_of_rho, only: temperature_of_rho, T_interior
    use omp_lib
    implicit none

    integer, parameter :: Nx = 100, Ntheta = 100
    real(dp), parameter :: x0 = 0.926484_dp
    real(dp), parameter :: sigma_right = 0.006_dp/2.3548_dp
    real(dp), parameter :: theta_adi_val = 0.55_dp
    real(dp), parameter :: B0_gauss = 1.0e12_dp / (1.0_dp+9.5_dp)   ! g=1 <-> B0_gauss
    real(dp), parameter :: e_esu = 4.80320425e-10_dp
    real(dp), parameter :: pasta_x_lo = 0.923897_dp, pasta_x_hi = 0.929071_dp
    real(dp), parameter :: reduction_factor = 1000.0_dp
    real(dp), parameter :: x_crust_mid = 0.934_dp, w_crust = 0.006_dp
    real(dp), parameter :: yr_s = 365.25_dp*86400.0_dp
    real(dp), parameter :: t_inj_stop = 1.0e4_dp
    integer,  parameter :: n_Q = 5, n_init = 4

    real(dp) :: Q_vals(n_Q)
    integer  :: init_vals(n_init)
    integer  :: iQ, iinit, n_runs, run_pick, upw_order
    real(dp) :: rho_cut_log, log_t_max, dt_safety, period_in, amp_in, tmp
    character(len=256) :: outdir
    character(len=512) :: key

    call init_crust_model("../../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)
    call init_cooling_curve("../../1d/crust_conductivity/FP_SF_B0.sparse")

    Q_vals    = (/ 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp /)
    init_vals = (/ 0, 1, 2, 3 /)          ! 0 = no pre-existing field, 1 = crustal field

    call env_real("INJ_PERIOD", 10.0_dp, period_in)
    call env_real("AMP_LOG",    0.8_dp,  amp_in)
    call env_real("RHO_CUT",    10.0_dp, rho_cut_log)
    call env_real("LOGT_MAX",   6.0_dp,  log_t_max)
    call env_real("DT_SAFETY",  1.6_dp,  dt_safety)
    call env_real("UPW_ORDER",  1.0_dp,  tmp)
    upw_order = nint(tmp)
    call env_str ("OUTDIR", "output", outdir)
    call env_real("RUN_IDX", 0.0_dp, tmp)
    run_pick = nint(tmp)

    print '(A,I0,A,F6.2,A,F5.2,A,F7.3,A,I0)', " OpenMP threads per run: ", omp_get_max_threads(), &
          "  Hall cutoff log10(rho)=", rho_cut_log, "  log t_max=", log_t_max, &
          "  DT_SAFETY=", dt_safety, "  UPW_ORDER=", upw_order

    n_runs = 0
    do iQ = 1, n_Q
        do iinit = 1, n_init
            n_runs = n_runs + 1
            if (run_pick == 0 .or. n_runs == run_pick) then
                call make_key(Q_vals(iQ), period_in, amp_in, init_vals(iinit), rho_cut_log, key)
                print '(A,I0,2A)', " === run ", n_runs, " : ", trim(key)
                call run_one(Q_vals(iQ), period_in, amp_in, init_vals(iinit), rho_cut_log, &
                             log_t_max, dt_safety, upw_order, trim(outdir), trim(key))
            end if
        end do
    end do

contains

    subroutine env_real(name, default, val)
        character(len=*), intent(in) :: name
        real(dp), intent(in)  :: default
        real(dp), intent(out) :: val
        character(len=64) :: buf
        integer :: st, ios
        val = default
        call get_environment_variable(name, buf, status=st)
        if (st == 0) then
            read(buf, *, iostat=ios) val
            if (ios /= 0) val = default
        end if
    end subroutine env_real

    subroutine env_str(name, default, val)
        character(len=*), intent(in) :: name, default
        character(len=*), intent(out) :: val
        integer :: st
        call get_environment_variable(name, val, status=st)
        if (st /= 0 .or. len_trim(val) == 0) val = default
    end subroutine env_str

    subroutine make_key(Q, freq, logg, init, rcut, key)
        real(dp), intent(in) :: Q, freq, logg, rcut
        integer, intent(in) :: init
        character(len=*), intent(out) :: key
        character(len=8)  :: qs, as, rs
        character(len=12) :: ps
        character(len=40) :: bs
        write(qs, '(F5.3)') Q
        write(as, '(F5.2)') logg
        if (as(1:1) == ' ') as(1:1) = '0'
        if (abs(freq-nint(freq)) < 1.0e-9_dp) then
            write(ps, '(I0)') nint(freq)
        else
            write(ps, '(F0.2)') freq
        end if
        write(rs, '(F4.1)') rcut
        if (init == 0) then
            bs = "_B_none"
        else if (init == 2) then
            bs = "_Datta"
        else if (init == 3) then
            bs = "_Ch7"
        else
            write(bs, '(A,F5.3)') "_B00.00_Kc", x_crust_mid
        end if
        write(key, '(A,A,A,A,A,A,A,A)') "ns2d_H1_S1_N100x100_Q", trim(qs), "_I1_P", trim(ps), &
              "_A", trim(as), trim(bs), "_Rc"//trim(adjustl(rs))//"_V2"
    end subroutine make_key

  
    pure function minmod(a, b) result(m)
        real(dp), intent(in) :: a, b
        real(dp) :: m
        if (a*b <= 0.0_dp) then
            m = 0.0_dp
        else
            m = sign(min(abs(a), abs(b)), a)
        end if
    end function minmod

    ! upwind derivative d/dx of F at (i,j) for advection speed u (x direction)
    function updx(F, i, j, u, h, ord) result(d)
        real(dp), intent(in) :: F(Nx,Ntheta), u, h
        integer, intent(in)  :: i, j, ord
        real(dp) :: d
        if (u >= 0.0_dp) then
            d = (F(i,j)-F(i-1,j))/h
            if (ord >= 2 .and. i >= 3) &
                d = d + 0.5_dp*minmod(F(i,j)-2.0_dp*F(i-1,j)+F(i-2,j), &
                                      F(i+1,j)-2.0_dp*F(i,j)+F(i-1,j))/h
        else
            d = (F(i+1,j)-F(i,j))/h
            if (ord >= 2 .and. i <= Nx-2) &
                d = d - 0.5_dp*minmod(F(i,j)-2.0_dp*F(i+1,j)+F(i+2,j), &
                                      F(i-1,j)-2.0_dp*F(i,j)+F(i+1,j))/h
        end if
    end function updx

    ! upwind derivative d/dtheta
    function updt(F, i, j, u, h, ord) result(d)
        real(dp), intent(in) :: F(Nx,Ntheta), u, h
        integer, intent(in)  :: i, j, ord
        real(dp) :: d
        if (u >= 0.0_dp) then
            d = (F(i,j)-F(i,j-1))/h
            if (ord >= 2 .and. j >= 3) &
                d = d + 0.5_dp*minmod(F(i,j)-2.0_dp*F(i,j-1)+F(i,j-2), &
                                      F(i,j+1)-2.0_dp*F(i,j)+F(i,j-1))/h
        else
            d = (F(i,j+1)-F(i,j))/h
            if (ord >= 2 .and. j <= Ntheta-2) &
                d = d - 0.5_dp*minmod(F(i,j)-2.0_dp*F(i,j+1)+F(i,j+2), &
                                      F(i,j-1)-2.0_dp*F(i,j)+F(i,j+1))/h
        end if
    end function updt

   
    ! Hall right-hand sides  fP = (dP/dt)_Hall , fT = (dT/dt)_Hall
    subroutine hall_rhs(P, T, fP, fT, x, sth, cth, H, Hp, dx, dth, ord)
        real(dp), intent(in)  :: P(Nx,Ntheta), T(Nx,Ntheta)
        real(dp), intent(in)  :: x(Nx), sth(Ntheta), cth(Ntheta), H(Nx), Hp(Nx), dx, dth
        integer,  intent(in)  :: ord
        real(dp), intent(out) :: fP(Nx,Ntheta), fT(Nx,Ntheta)
        real(dp), allocatable :: Ar(:,:), Ath(:,:)
        real(dp) :: lap, dPx, dPt, ps, vx, vt, ux, ut, tsrc
        integer :: i, j

        allocate(Ar(Nx,Ntheta), Ath(Nx,Ntheta))
        Ar = 0.0_dp; Ath = 0.0_dp; fP = 0.0_dp; fT = 0.0_dp

        ! poloidal-current terms  A_r , A_theta  (from P only)
        !$omp parallel do default(shared) private(i,j,lap,dPx,dPt) schedule(static)
        do j = 2, Ntheta-1
            do i = 2, Nx-1
                lap = (P(i+1,j)-2.0_dp*P(i,j)+P(i-1,j))/dx**2 &
                    + (1.0_dp/x(i)**2)*( (P(i,j+1)-2.0_dp*P(i,j)+P(i,j-1))/dth**2 &
                      - (cth(j)/sth(j))*(P(i,j+1)-P(i,j-1))/(2.0_dp*dth) )
                dPx = (P(i+1,j)-P(i-1,j))/(2.0_dp*dx)
                dPt = (P(i,j+1)-P(i,j-1))/(2.0_dp*dth)
                Ar(i,j)  = (H(i)/(x(i)**2*sth(j)**2))*lap*dPx
                Ath(i,j) = (H(i)/(x(i)**3*sth(j)**2))*lap*dPt
            end do
        end do
        !$omp end parallel do

        ! drift of P by the toroidal field
        !$omp parallel do default(shared) private(i,j,ps,vx,vt,ux,ut) schedule(static)
        do j = 2, Ntheta-1
            do i = 2, Nx-1
                ps = H(i)/(x(i)**2*sth(j))
                vx = -ps*(T(i,j+1)-T(i,j-1))/(2.0_dp*dth)
                vt =  ps*(T(i+1,j)-T(i-1,j))/(2.0_dp*dx)
                ux = -T(i,j)*2.0_dp*H(i)*cth(j)/(x(i)**2*sth(j)**2)
                ut = -T(i,j)*(Hp(i)-2.0_dp*H(i)/x(i))/(x(i)**2*sth(j))
                fP(i,j) = -(vx*updx(P,i,j,vx,dx,ord) + vt*updt(P,i,j,vt,dth,ord))
                fT(i,j) = -(ux*updx(T,i,j,ux,dx,ord) + ut*updt(T,i,j,ut,dth,ord))
            end do
        end do
        !$omp end parallel do

        ! T source from P
        !$omp parallel do default(shared) private(i,j,tsrc) schedule(static)
        do j = 3, Ntheta-2
            do i = 3, Nx-2
                tsrc = sth(j)*Ath(i,j) &
                     + x(i)*sth(j)*(Ath(i+1,j)-Ath(i-1,j))/(2.0_dp*dx) &
                     - sth(j)*(Ar(i,j+1)-Ar(i,j-1))/(2.0_dp*dth)
                fT(i,j) = fT(i,j) + tsrc
            end do
        end do
        !$omp end parallel do
        deallocate(Ar, Ath)
    end subroutine hall_rhs

    ! relative magnetic energies (poloidal, toroidal) ~ int B^2 dV
    subroutine energies(P, T, x, sth, dx, dth, epol, etor)
        real(dp), intent(in)  :: P(Nx,Ntheta), T(Nx,Ntheta), x(Nx), sth(Ntheta), dx, dth
        real(dp), intent(out) :: epol, etor
        real(dp) :: dPx, dPt
        integer :: i, j
        epol = 0.0_dp; etor = 0.0_dp
        do j = 2, Ntheta-1
            do i = 2, Nx-1
                dPx = (P(i+1,j)-P(i-1,j))/(2.0_dp*dx)
                dPt = (P(i,j+1)-P(i,j-1))/(2.0_dp*dth)
                epol = epol + (dPx**2 + (dPt/x(i))**2)/sth(j)
                etor = etor + T(i,j)**2/sth(j)
            end do
        end do
        epol = epol*dx*dth; etor = etor*dx*dth
    end subroutine energies

   
    subroutine run_one(Q_imp, period, logg, init_field, rcut_log, logt_max, dt_safety, ord, outdir, key)
        real(dp), intent(in) :: Q_imp, period, logg, rcut_log, logt_max, dt_safety
        integer, intent(in)  :: init_field, ord
        character(len=*), intent(in) :: outdir, key

        real(dp) :: x(Nx), theta(Ntheta), sth(Ntheta), cth(Ntheta), rho_k(Nx)
        real(dp) :: bump(Nx), gcr(Nx), H_of_x(Nx), Hprime(Nx), S_local(Nx)
        real(dp), allocatable :: P(:,:), Tpot(:,:), Psave(:,:), Tsave(:,:), Zv(:,:)
        real(dp), allocatable :: fP(:,:), fT(:,:), fP1(:,:), fT1(:,:), fP2(:,:), fT2(:,:)
        real(dp) :: x_c, x_s, dx, dtheta, R_star_cm, amplitude, sigma_left, gx
        real(dp) :: H_max, dx_min2, n_e_here, rho_cut
        real(dp) :: t_yr, next_inj, t_target, log_t_min
        real(dp) :: t_end, max_frac, field_scale, dt_hall_yr, dt_sec
        real(dp) :: Ts_now, Tb_now, T_here, sig, Ts0, Ts1, t_mid
        real(dp) :: ab0, ab1, ab2, hh1, hh2, ha, hb, hs
        real(dp) :: cap_fac, pm_old, tm_old, pm_new, tm_new, psurf, ssum, epol, etor
        integer  :: i, j, k, unit_out, step_count, n_inj_done, n_out, shrink
        integer  :: nhist, n_rej, n_rej_consec, mt(2), mf(2)
        integer(8) :: clk0, clk1, crate
        logical  :: blown, bad

        allocate(P(Nx,Ntheta), Tpot(Nx,Ntheta), Psave(Nx,Ntheta), Tsave(Nx,Ntheta), Zv(Nx,Ntheta), &
                 fP(Nx,Ntheta), fT(Nx,Ntheta), fP1(Nx,Ntheta), fT1(Nx,Ntheta), fP2(Nx,Ntheta), fT2(Nx,Ntheta))
        Zv = 0; fP = 0; fT = 0; fP1 = 0; fT1 = 0; fP2 = 0; fT2 = 0

        x_c = pasta_x_lo; x_s = 1.0_dp
        if (init_field == 2) x_c = get_x_crust_core()
        if (init_field == 3) x_c = 0.90_dp   ! crust-core boundary 
        dx = (x_s - x_c)/real(Nx-1, dp)
        do i = 1, Nx
            x(i) = x_c + real(i-1,dp)*dx
        end do
        dtheta = pi/real(Ntheta-1, dp)
        do j = 1, Ntheta
            theta(j) = real(j-1,dp)*dtheta
            sth(j) = sin(theta(j)); cth(j) = cos(theta(j))
        end do
        R_star_cm = crust_R_star_km*1.0e5_dp

        ! ---- Hall coefficient H(x) = c B0 / (4 pi e n_e R^2)
        rho_cut = 10.0_dp**rcut_log
        do i = 1, Nx
            rho_k(i) = get_rho_at_r(x(i)*crust_R_star_km)
            if (rho_k(i) >= rho_cut) then
                n_e_here = get_electron_density_at_r(rho_k(i))
                H_of_x(i) = c_light*B0_gauss/(4.0_dp*pi*e_esu*n_e_here*R_star_cm**2)
            else
                H_of_x(i) = 0.0_dp
            end if
        end do
        H_max = maxval(H_of_x)
        Hprime(1) = (H_of_x(2)-H_of_x(1))/dx
        do i = 2, Nx-1
            Hprime(i) = (H_of_x(i+1)-H_of_x(i-1))/(2.0_dp*dx)
        end do
        Hprime(Nx) = (H_of_x(Nx)-H_of_x(Nx-1))/dx
        dx_min2 = min(dx**2, (x_c*dtheta)**2)

        ! ---- injection 
        amplitude = 10.0_dp**logg - 1.0_dp
        sigma_left = sigma_right*sqrt((1.0_dp+amplitude)/amplitude)
        do i = 1, Nx
            gx = x(i) - x0
            if (x(i) >= x0) then
                bump(i) = amplitude*exp(-gx*gx/(2.0_dp*sigma_right**2))
            else
                bump(i) = amplitude*exp(-gx*gx/(2.0_dp*sigma_left**2))
            end if
        end do

        ! ---- pre-existing crustal field 
        gcr = 0.0_dp
        if (init_field == 1) then
            call g_profile_crustal(x, gcr, Nx, dx, x_crust_mid, w_crust, x_c)
            gcr = gcr/gcr(Nx)
        end if
        if (init_field == 2) then
            call g_profile_expelled_datta(x, gcr, Nx, dx, x_c+2.0_dp*0.01_dp, 0.01_dp, x_c)
            gcr = gcr/gcr(Nx)
            bump = 0.0_dp
        end if
        if (init_field == 3) then
            call g_profile_expelled_ch7(x, gcr, Nx, 0.925_dp, sigma_right, 9.5_dp)
            bump = 0.0_dp
        end if

        ! ---- initial state (first injection at t=0)
        do j = 1, Ntheta
            do i = 1, Nx
                P(i,j) = (gcr(i) + bump(i))*sth(j)**2
            end do
        end do
        Tpot = 0.0_dp

        t_yr = 1.0e-3_dp
        n_inj_done = 1
        step_count = 0
        blown = .false.
        nhist = 0
        hh1 = 0.0_dp; hh2 = 0.0_dp
        cap_fac = 1.0_dp
        n_rej = 0; n_rej_consec = 0
        next_inj = period

        open(newunit=unit_out, file=trim(outdir)//"/"//trim(key)//".dat", status="replace", action="write")
        write(unit_out, '(A)') "# "//trim(key)
        write(unit_out, '(A,F6.3,A,F6.2,A,F5.2)') "# Q=", Q_imp, "  period_yr=", period, "  logg=", logg
        write(unit_out, '(A,F6.2,A,ES10.3)') "# Hall cutoff log10(rho)=", rcut_log, "  H_max(1/s)=", H_max
        write(unit_out, '(A,F7.3,A,I0)') "# Ohmic implicit ADI + Hall explicit variable-step AB3;  DT_SAFETY=", &
              dt_safety, "  UPW_ORDER=", ord
        write(unit_out, '(A)') "# equations: Skiathas & Gourgouliatos 2024 eqs 9,17,19-21, expanded exactly (see source header)"
        write(unit_out, '(A)') "# units: g=1 (crustal plateau) <-> B0 = 1e12/10.5 G"
        write(unit_out, '(A)') "# log10(t/yr) g_surface_eq max|P| max|T| nsteps n_inj n_rej E_pol E_tor"
        flush(unit_out)

        call system_clock(clk0, crate)
        log_t_min = -1.0_dp
        do n_out = 1, 400
            t_target = 10.0_dp**(log_t_min + (logt_max-log_t_min)*real(n_out,dp)/400.0_dp)
            if (t_target <= t_yr) cycle

            do while (t_yr < t_target)
                psurf = P(Nx, Ntheta/2+1)
                if (psurf > 1.0e-3_dp) then
                    max_frac = 0.02_dp
                else
                    max_frac = 0.002_dp
                end if
                pm_old = maxval(abs(P))
                tm_old = maxval(abs(Tpot))
                field_scale = max(pm_old, tm_old, 1.0e-6_dp)
                if (H_max > 0.0_dp) then
                    dt_hall_yr = dt_safety*cap_fac*0.5_dp*dx_min2/(H_max*field_scale)/yr_s
                else
                    dt_hall_yr = huge(1.0_dp)
                end if

                t_end = min(t_target, next_inj, t_yr*(1.0_dp+max_frac), t_yr+dt_hall_yr)
                if (t_end <= t_yr) t_end = t_yr + 1.0e-9_dp

                ! surface temperature may not change by more than 5% per step
                Ts0 = Ts_at_age(max(t_yr,1.0e-3_dp))
                do shrink = 1, 40
                    Ts1 = Ts_at_age(t_end)
                    if (abs(Ts1-Ts0) <= 0.05_dp*Ts0) exit
                    t_end = t_yr + 0.5_dp*(t_end-t_yr)
                end do

                dt_sec = (t_end - t_yr)*yr_s
                t_mid = 0.5_dp*(t_yr + t_end)
                Ts_now = Ts_at_age(max(t_mid,1.0e-3_dp))
                Tb_now = T_interior(Ts_now, crust_gs)

                !$omp parallel do default(shared) private(k,T_here,sig) schedule(static)
                do k = 1, Nx
                    T_here = temperature_of_rho(rho_k(k), Tb_now, Ts_now, crust_gs)
                    sig = get_sigma_given_T(rho_k(k), T_here, Q_imp)
                    if (init_field /= 2 .and. init_field /= 3 .and. x(k) >= pasta_x_lo .and. x(k) <= pasta_x_hi) sig = sig/reduction_factor
                    S_local(k) = c_light**2/(4.0_dp*pi*sig*R_star_cm**2)
                end do
                !$omp end parallel do

                Psave = P
                Tsave = Tpot

                ! ---- Hall right-hand sides from the current state ----
                call hall_rhs(P, Tpot, fP, fT, x, sth, cth, H_of_x, Hprime, dx, dtheta, ord)

                ! ---- variable-step Adams-Bashforth weights (seconds) ----
                hs = dt_sec
                ab1 = 0.0_dp; ab2 = 0.0_dp
                if (nhist == 0) then
                    ab0 = hs
                else if (nhist == 1) then
                    ha  = hh1
                    ab0 = (hs*hs/2.0_dp + ha*hs)/ha
                    ab1 = -hs*hs/(2.0_dp*ha)
                else
                    ha  = hh1
                    hb  = hh1 + hh2
                    ab0 = (hs**3/3.0_dp + (ha+hb)*hs**2/2.0_dp + ha*hb*hs)/(ha*hb)
                    ab1 = -(hs**3/3.0_dp + hb*hs**2/2.0_dp)/(ha*(hb-ha))
                    ab2 =  (hs**3/3.0_dp + ha*hs**2/2.0_dp)/(hb*(hb-ha))
                end if

                ! ---- implicit Ohmic step 
                call adi_step_2d_hall_omp(P, x, theta, Nx, Ntheta, dx, dtheta, x_s, &
                                          S_local, Zv, Zv, dt_sec, theta_adi_val)
                call adi_step_2d_T_flux_omp(Tpot, x, theta, Nx, Ntheta, dx, dtheta, &
                                            S_local, dt_sec, theta_adi_val)

                ! ---- explicit Hall increment 
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 2, Ntheta-1
                    do i = 2, Nx-1
                        P(i,j)    = P(i,j)    + ab0*fP(i,j) + ab1*fP1(i,j) + ab2*fP2(i,j)
                        Tpot(i,j) = Tpot(i,j) + ab0*fT(i,j) + ab1*fT1(i,j) + ab2*fT2(i,j)
                    end do
                end do
                !$omp end parallel do

                ! ---- accept / reject (guard against accidents) 
                pm_new = maxval(abs(P))
                tm_new = maxval(abs(Tpot))
                ssum = sum(P) + sum(Tpot)
                bad = (ssum /= ssum) .or. any(P /= P) .or. any(Tpot /= Tpot) &
                      .or. (pm_new > 1.0e30_dp) .or. (tm_new > 1.0e30_dp) &
                      .or. (pm_new > 2.0_dp*pm_old + 1.0_dp) &
                      .or. (tm_new > 2.0_dp*tm_old + 1.0e-3_dp)
                if (bad) then
                    P = Psave
                    Tpot = Tsave
                    cap_fac = max(0.5_dp*cap_fac, 1.0e-8_dp)
                    n_rej = n_rej + 1
                    n_rej_consec = n_rej_consec + 1
                    if (n_rej_consec > 200) then
                        blown = .true.
                        mt = maxloc(abs(Tpot)); mf = maxloc(abs(fT))
                        write(unit_out,'(A,2I5,ES12.4)') "# DIAG max|T| at (i,j) and value:", mt, Tpot(mt(1),mt(2))
                        write(unit_out,'(A,2I5,ES12.4)') "# DIAG max|fT| at (i,j) and value:", mf, fT(mf(1),mf(2))
                        write(unit_out,'(A,3ES12.4)') "# DIAG dt_sec ab0 cap_fac:", dt_sec, ab0, cap_fac
                        exit
                    end if
                    cycle
                end if
                n_rej_consec = 0
                cap_fac = min(1.0_dp, cap_fac*1.005_dp)

                ! shift the Hall history
                fP2 = fP1; fT2 = fT1
                fP1 = fP;  fT1 = fT
                hh2 = hh1; hh1 = dt_sec
                nhist = min(nhist+1, 2)

                t_yr = t_end
                step_count = step_count + 1

                if (t_yr >= t_inj_stop - 1.0e-6_dp .and. next_inj < 1.0e10_dp) next_inj = huge(1.0_dp)
                if (next_inj < 1.0e10_dp) then
                    if (abs(t_yr-next_inj) < 1.0e-6_dp*next_inj .and. t_yr < t_inj_stop-1.0e-6_dp) then
                        do j = 1, Ntheta
                            do i = 1, Nx
                                P(i,j) = P(i,j) + bump(i)*sth(j)**2
                            end do
                        end do
                        n_inj_done = n_inj_done + 1
                        next_inj = next_inj + period
                        nhist = 0          ! injection is a jump: restart the AB histry
                    end if
                end if

                if (mod(step_count, 200000) == 0) then
                    print '(A,I0,A,ES11.4,A,ES11.4,A,ES11.4,A,I0,A,I0,A,ES9.2)', "  step=", step_count, &
                          "  t=", t_yr, "  max|P|=", maxval(abs(P)), "  max|T|=", maxval(abs(Tpot)), &
                          "  inj=", n_inj_done, "  rej=", n_rej, "  cap=", cap_fac
                end if
            end do
            if (blown) exit
            call energies(P, Tpot, x, sth, dx, dtheta, epol, etor)
            write(unit_out, '(F10.5,3X,ES16.8,3X,ES14.6,3X,ES14.6,3X,I10,3X,I7,3X,I8,2(3X,ES13.5))') log10(t_yr), &
                  P(Nx,Ntheta/2+1), maxval(abs(P)), maxval(abs(Tpot)), step_count, n_inj_done, n_rej, epol, etor
            flush(unit_out)
        end do

        call system_clock(clk1)
        if (blown) then
            write(unit_out, '(A,ES12.4)') "# FAILED (200 consecutive rejected steps) at t_yr=", t_yr
            print '(A,ES12.4)', " !!! run failed at t_yr=", t_yr
        end if
        write(unit_out, '(A,I0,A,I0,A,F10.2,A)') "# total steps ", step_count, "  rejected ", n_rej, &
              "  wall ", real(clk1-clk0,dp)/real(crate,dp), " s"
        close(unit_out)
        open(newunit=unit_out, file=trim(outdir)//"/"//trim(key)//"_profiles.dat", status="replace", action="write")
        write(unit_out,'(A)') "# x  P(x,eq)  T(x,theta=pi/4)  T(x,theta=pi/2)  P(x,pi/4)"
        do i = 1, Nx
            write(unit_out,'(F9.5,4ES16.7)') x(i), P(i,Ntheta/2+1), Tpot(i,Ntheta/4+1), Tpot(i,Ntheta/2+1), P(i,Ntheta/4+1)
        end do
        close(unit_out)
        print '(A,I0,A,I0,A,F10.2,A,I0,A,ES11.4)', " finished: steps=", step_count, "  rejected=", n_rej, &
              "  wall=", real(clk1-clk0,dp)/real(crate,dp), " s  injections=", n_inj_done, "  g_surf=", P(Nx,Ntheta/2+1)
        deallocate(P, Tpot, Psave, Tsave, Zv, fP, fT, fP1, fT1, fP2, fT2)
    end subroutine run_one

end program Ns2dHall
