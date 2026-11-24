module cooling_curve

    use utilities
    use cubic_spline
    implicit none

    integer, parameter :: n_cool = 44
    real(dp), save :: y2_cool(n_cool)
    real(dp), save :: log_t_tab(n_cool), log_Ts_tab(n_cool)
    logical, save :: cool_initialized = .false.

contains

    
    subroutine init_cooling_curve(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: unit, i

        open(newunit=unit, file=filename, status='old', action='read')
        do i = 1, n_cool
            read(unit,*) log_t_tab(i), log_Ts_tab(i)
        end do
        close(unit)

        call spline_coeffs(log_t_tab, log_Ts_tab, n_cool, y2_cool)
        cool_initialized = .true.

    end subroutine init_cooling_curve


   
    function Ts_at_age(t_yr) result(Ts)
        implicit none
        real(dp), intent(in) :: t_yr
        real(dp) :: Ts
        real(dp) :: age, log_ts_val

        if (.not. cool_initialized) then
            print *, "ERROR: cooling_curve not initialized."
            stop 1
        end if

        if (t_yr <= 0.0_dp) then
            age = log_t_tab(1)
        else
            age = log10(t_yr)
        end if

        if (age <= log_t_tab(1)) then
            log_ts_val = log_Ts_tab(1)
        else if (age >= log_t_tab(n_cool)) then
            ! Extrapolate past the tabulated range, matching BOTH the
            ! value AND slope of the spline at the boundary point
           
            block
                real(dp) :: boundary_val, slope_at_boundary
                real(dp) :: h_local
                h_local = 1.0e-4_dp
                call spline_eval(log_t_tab, log_Ts_tab, y2_cool, n_cool, &
                                  log_t_tab(n_cool), boundary_val)
                call spline_eval(log_t_tab, log_Ts_tab, y2_cool, n_cool, &
                                  log_t_tab(n_cool)-h_local, slope_at_boundary)
                slope_at_boundary = (boundary_val - slope_at_boundary) / h_local
                log_ts_val = boundary_val + slope_at_boundary*(age - log_t_tab(n_cool))
            end block
        else
            call spline_eval(log_t_tab, log_Ts_tab, y2_cool, n_cool, age, log_ts_val)
        end if

        Ts = 10.0_dp**log_ts_val

    end function Ts_at_age

end module cooling_curve
