module crust_model

    use utilities
    use pchip
    use density_profile_reader
    use t_of_rho
    use conductivity
    use tov_constants, only: G_grav, c_light, Msun

    implicit none
    private
    public :: init_crust_model
    public :: get_rho_at_r, get_T_at_rho, get_sigma_at_rho, get_sigma_at_r
    public :: get_sigma_given_T, get_electron_density_at_r
    public :: crust_M_star, crust_R_star_km, crust_gs, crust_Tb, crust_Ts
    public :: get_crust_core_radius, get_x_crust_core
    public :: get_x_at_density

    real(dp), allocatable, save :: r_prof(:), rho_prof(:), d_prof(:)
    integer, save :: n_prof
    real(dp), save :: crust_M_star, crust_R_star_km, crust_gs, crust_Tb, crust_Ts
    logical, save :: crust_initialized = .false.

    ! Crust-core boundary, by convention taken at half of nuclear
    ! saturation density (rho_0 ~ 2.8e14 g/cm^3). NOT the neutron
    ! drip density (~4.3e11) 
    real(dp), parameter :: rho_saturation = 2.8e14_dp
    
    ! Originally set to 0.5*rho_saturation; revised
    ! to 3.0e14 g/cm^3 after finding it better matches the crust
    ! width implied by Bhattacharya & Datta (1996) Fig. 1,
    ! AkmalPR star. 
    public :: rho_crust_core
    real(dp), parameter :: rho_crust_core = 1.290e14_dp

contains

   
    subroutine init_crust_model(density_profile_file, Ts_input, M_star_Msun)

        implicit none

        character(len=*), intent(in) :: density_profile_file
        real(dp), intent(in) :: Ts_input       ! surface temp, K
        real(dp), intent(in) :: M_star_Msun    ! stellar mass, solar masses

        real(dp), allocatable :: log10rho_tmp(:)
        real(dp) :: R_star_cm

      
        if (crust_initialized) return

        call read_density_profile(density_profile_file, r_prof, rho_prof, log10rho_tmp, n_prof)
        deallocate(log10rho_tmp)

        allocate(d_prof(n_prof))
        call pchip_coeffs(r_prof, rho_prof, n_prof, d_prof)

        crust_M_star    = M_star_Msun * Msun
        crust_R_star_km = r_prof(n_prof)
        R_star_cm       = crust_R_star_km * 1.0e5_dp

        crust_gs = surface_gravity(crust_M_star, R_star_cm, G_grav, c_light)
        crust_Ts = Ts_input
        crust_Tb = T_interior(crust_Ts, crust_gs)

        crust_initialized = .true.

        print *
        print '(A)', " Crust model initialized:"
        print '(A,F10.4,A)',  "   M_star = ", crust_M_star/Msun, " Msun"
        print '(A,F10.4,A)',  "   R_star = ", crust_R_star_km, " km"
        print '(A,ES14.6,A)', "   g_s    = ", crust_gs, " cm/s^2"
        print '(A,ES14.6,A)', "   T_s    = ", crust_Ts, " K"
        print '(A,ES14.6,A)', "   T_b    = ", crust_Tb, " K"
        print *

    end subroutine init_crust_model


    function get_rho_at_r(r) result(rho)
        implicit none
        real(dp), intent(in) :: r   ! km
        real(dp) :: rho
        if (.not. crust_initialized) then
            print *, "ERROR: crust_model not crust_initialized -- call init_crust_model first."
            stop 1
        end if
        call pchip_eval(r_prof, rho_prof, d_prof, n_prof, r, rho)
    end function get_rho_at_r


    function get_T_at_rho(rho) result(T)
        implicit none
        real(dp), intent(in) :: rho
        real(dp) :: T
        T = temperature_of_rho(rho, crust_Tb, crust_Ts, crust_gs)
    end function get_T_at_rho


  
    function get_sigma_given_T(rho, T, Q) result(sigma)
        implicit none
        real(dp), intent(in) :: rho, T, Q
        real(dp) :: sigma
        real(dp) :: T8, rholog
        real(dp) :: rholog_min, rholog_max

        rholog_min = 3.0_dp + u0(1)     / uc
        rholog_max = 3.0_dp + u0(n_tab) / uc

        rholog = log10(rho)
        if (rholog < rholog_min) rholog = rholog_min
        if (rholog > rholog_max) rholog = rholog_max

        T8 = T/1.0e8_dp
        call sigma_ns(rholog, T8, Q, sigma)
    end function get_sigma_given_T

    function get_electron_density_at_r(rho) result(n_e)
        implicit none
        real(dp), intent(in) :: rho
        real(dp) :: n_e
        real(dp) :: rholog
        real(dp) :: rholog_min, rholog_max
        rholog_min = 3.0_dp + u0(1)     / uc
        rholog_max = 3.0_dp + u0(n_tab) / uc
        rholog = log10(rho)
        if (rholog < rholog_min) rholog = rholog_min
        if (rholog > rholog_max) rholog = rholog_max
        call get_electron_density(rholog, n_e)
    end function get_electron_density_at_r


    function get_sigma_at_rho(rho, Q) result(sigma)
        implicit none
        real(dp), intent(in) :: rho, Q
        real(dp) :: sigma
        real(dp) :: T

        T = get_T_at_rho(rho)
        sigma = get_sigma_given_T(rho, T, Q)
    end function get_sigma_at_rho


    function get_sigma_at_r(r, Q) result(sigma)
        implicit none
        real(dp), intent(in) :: r, Q
        real(dp) :: sigma
        real(dp) :: rho
        rho = get_rho_at_r(r)
        sigma = get_sigma_at_rho(rho, Q)
    end function get_sigma_at_r


    !--------------------------------------------------------------
    ! Radius (km) at which rho(r) = rho_crust_core, found by
    ! bisection on the loaded profile. Not a stored constant
    ! -- recomputed from whatever star is currently initialized, so
    ! it stays correct if a different EOS is loaded later.
    !--------------------------------------------------------------
    function get_crust_core_radius() result(r_cc)
        implicit none
        real(dp) :: r_cc
        real(dp) :: r_lo, r_hi, r_mid, rho_mid
        integer :: iter

        if (.not. crust_initialized) then
            print *, "ERROR: crust_model not initialized -- call init_crust_model first."
            stop 1
        end if

        r_lo = r_prof(1)
        r_hi = r_prof(n_prof)

        if (rho_crust_core > get_rho_at_r(r_lo) .or. rho_crust_core < get_rho_at_r(r_hi)) then
            print *, "ERROR: rho_crust_core is outside the loaded profile's density range."
            stop 1
        end if

        ! rho(r) decreases monotonically with r 
        do iter = 1, 100
            r_mid = 0.5_dp*(r_lo+r_hi)
            rho_mid = get_rho_at_r(r_mid)
            if (rho_mid > rho_crust_core) then
                r_lo = r_mid
            else
                r_hi = r_mid
            end if
            if (abs(r_hi-r_lo) < 1.0e-8_dp) exit
        end do

        r_cc = 0.5_dp*(r_lo+r_hi)

    end function get_crust_core_radius


   
    function get_x_crust_core() result(x_cc)
        implicit none
        real(dp) :: x_cc
        x_cc = get_crust_core_radius() / crust_R_star_km
    end function get_x_crust_core


   
    ! Fractional radius x=r/R at which rho(r) = target_rho, found
    ! by bisection on the actual profile.
  
    function get_x_at_density(target_rho) result(x_out)
        implicit none
        real(dp), intent(in) :: target_rho
        real(dp) :: x_out
        real(dp) :: r_lo, r_hi, r_mid, rho_mid
        integer :: iter

        if (.not. crust_initialized) then
            print *, "ERROR: crust_model not initialized -- call init_crust_model first."
            stop 1
        end if

        r_lo = r_prof(1)
        r_hi = r_prof(n_prof)

        do iter = 1, 100
            r_mid = 0.5_dp*(r_lo+r_hi)
            rho_mid = get_rho_at_r(r_mid)
            if (rho_mid > target_rho) then
                r_lo = r_mid
            else
                r_hi = r_mid
            end if
            if (abs(r_hi-r_lo) < 1.0e-8_dp) exit
        end do

        x_out = 0.5_dp*(r_lo+r_hi) / crust_R_star_km

    end function get_x_at_density

end module crust_model
