module density_profile_reader

    use utilities
    implicit none

contains


    ! Reads a 3-column (r[km], rho[g/cm^3], log10(rho)) profile
    ! file, as produced by TOV_Solver.
  
    subroutine read_density_profile(filename, r, rho, log10rho, n)

        implicit none

        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: r(:), rho(:), log10rho(:)
        integer, intent(out) :: n

        integer :: unit, ios, count
        character(len=256) :: line
        real(dp) :: rv, rhov, lv
        logical :: file_exists

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "ERROR: density profile file not found: ", trim(filename)
            stop 1
        end if

        open(newunit=unit, file=filename, status='old', action='read')

        count = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            count = count + 1
        end do

        n = count
        if (n < 2) then
            print *, "ERROR: fewer than 2 usable rows in ", trim(filename)
            stop 1
        end if

        block
            real(dp), allocatable :: r_raw(:), rho_raw(:), log10rho_raw(:)
            integer :: n_dupes

            allocate(r_raw(n), rho_raw(n), log10rho_raw(n))

            rewind(unit)
            count = 0
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle

                read(line, *, iostat=ios) rv, rhov, lv
                if (ios /= 0) cycle

                count = count + 1
                r_raw(count) = rv
                rho_raw(count) = rhov
                log10rho_raw(count) = lv
            end do

            close(unit)

            ! Deduplicate: keep only strictly increasing r. The
            ! written profile can have consecutive rows that round
            ! to the same printed r (adaptive RK45 steps shrink to
            ! sub-file-precision near a stiff surface) -- these
            ! give PCHIP a zero interval width and corrupt every
            ! derivative from there onward. Drop the repeats.
            allocate(r(n), rho(n), log10rho(n))
            r(1) = r_raw(1)
            rho(1) = rho_raw(1)
            log10rho(1) = log10rho_raw(1)
            n_dupes = 0
            n = 1

            block
                integer :: k
                do k = 2, size(r_raw)
                    if (r_raw(k) > r(n)) then
                        n = n + 1
                        r(n) = r_raw(k)
                        rho(n) = rho_raw(k)
                        log10rho(n) = log10rho_raw(k)
                    else
                        n_dupes = n_dupes + 1
                    end if
                end do
            end block

            if (n_dupes > 0) then
                print '(A,I0,A)', " NOTE: dropped ", n_dupes, &
                      " duplicate/non-increasing r rows (sub-precision near surface)."
            end if

            deallocate(r_raw, rho_raw, log10rho_raw)
        end block

        print '(A,I0,A)', " Loaded ", n, " rows from density profile."
        print '(A,F10.4,A,F10.4,A)', "   r range: ", r(1), " to ", r(n), " km"
        print '(A,F8.3,A,F8.3)', "   log10(rho) range: ", log10rho(1), " to ", log10rho(n)

    end subroutine read_density_profile

end module density_profile_reader
