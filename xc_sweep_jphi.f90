program XCSweepJphi
    use utilities
    implicit none

    integer, parameter :: npts = 2000
    real(dp) :: x(npts), g(npts), Jphi(npts)
    real(dp), parameter :: x0 = 0.925_dp, xw = 0.006_dp, amplitude = 10.5_dp
    real(dp) :: dx, gx, x_min
    real(dp) :: xc_candidates(5)
    integer :: i, k, unit_out
    character(len=64) :: fname

    xc_candidates = (/ 0.905_dp, 0.910_dp, 0.915_dp, 0.918_dp, 0.920_dp /)

    do k = 1, 5
        x_min = xc_candidates(k)
        dx = (1.0_dp - x_min)/real(npts-1,dp)
        do i = 1, npts
            x(i) = x_min + real(i-1,dp)*dx
            gx = x(i)-x0
            g(i) = amplitude*exp(-gx*gx/(2.0_dp*xw*xw))
        end do

        do i = 2, npts-1
            Jphi(i) = -(( (g(i+1)-2.0_dp*g(i)+g(i-1))/dx**2 ) - 2.0_dp*g(i)/x(i)**2)
        end do
        Jphi(1) = 0.0_dp
        Jphi(npts) = 0.0_dp

        write(fname, '(A,I0,A)') "output/xc_sweep_", k, ".dat"
        open(newunit=unit_out, file=trim(fname), status="replace", action="write")
        write(unit_out, '(A,F8.4)') "# x_c=", x_min
        do i = 1, npts
            write(unit_out, '(F12.8,3X,ES16.8,3X,ES16.8)') x(i), g(i), Jphi(i)
        end do
        close(unit_out)

        block
            real(dp) :: dip1_min, dip2_min
            integer :: idx0
            idx0 = minloc(abs(x-x0), dim=1)
            dip1_min = minval(Jphi(1:idx0))     ! surface-side dip (x > x0)
            dip2_min = minval(Jphi(idx0:npts))  ! core-side dip (x < x0)
            print '(A,F6.3,A,ES12.4,A,ES12.4,A,F8.4)', &
                  " x_c=", x_min, "  surface-side dip=", dip1_min, &
                  "  core-side dip=", dip2_min, "  ratio(core/surf)=", dip2_min/dip1_min
        end block
    end do

end program XCSweepJphi
