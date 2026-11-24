program PredictJphi
    use utilities
    implicit none

    integer, parameter :: npts = 2000
    real(dp) :: x(npts), g(npts), Jphi(npts)
    real(dp), parameter :: x0 = 0.925_dp
    real(dp), parameter :: fwhm = 0.006_dp
    real(dp), parameter :: xw_right = fwhm / 2.3548_dp   ! surface side: matches Fig 7.1's stated width
    real(dp), parameter :: xw_left = xw_right   ! same width both sides -- fixes skew AND the curvature mismatch causing the mangled Jphi
    real(dp), parameter :: amplitude = 9.5_dp

    real(dp) :: dx, gx
    integer :: i, unit_out

    dx = (1.0_dp - 0.90_dp)/real(npts-1,dp)
    do i = 1, npts
        x(i) = 0.90_dp + real(i-1,dp)*dx
        gx = x(i)-x0
        if (x(i) >= x0) then
            g(i) = 1.0_dp + amplitude*exp(-gx*gx/(2.0_dp*xw_right*xw_right))
        else
            g(i) = (1.0_dp+amplitude)*exp(-gx*gx/(2.0_dp*xw_left*xw_left))
        end if
    end do

    do i = 2, npts-1
        Jphi(i) = -(( (g(i+1)-2.0_dp*g(i)+g(i-1))/dx**2 ) - 2.0_dp*g(i)/x(i)**2)
    end do
    Jphi(1) = 0.0_dp
    Jphi(npts) = 0.0_dp

    open(newunit=unit_out, file="output/predicted_jphi.dat", status="replace", action="write")
    write(unit_out,'(A)') "# x   g(x)   Jphi(x)"
    do i = 1, npts
        write(unit_out, '(F12.8,3X,ES16.8,3X,ES16.8)') x(i), g(i), Jphi(i)
    end do
    close(unit_out)
    print *, "Written to output/predicted_jphi.dat"

end program PredictJphi
