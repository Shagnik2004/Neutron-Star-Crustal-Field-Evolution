program SolveGFromJphi

    use utilities
    use crust_model
    use tridiagonal_solver

    implicit none

    integer, parameter :: N = 2000
    real(dp) :: x(N), Jphi(N), g(N)
    real(dp) :: a(N), b(N), c(N), d(N)
    real(dp) :: x_c, x_s, dx
    real(dp), parameter :: x0 = 0.925_dp
    real(dp), parameter :: fwhm_central = 0.006_dp
    real(dp), parameter :: xw_central = fwhm_central / 2.3548_dp
    real(dp), parameter :: A_central = 1.0_dp
    real(dp), parameter :: A_side = 0.45_dp     ! side-lobe amplitude, relative
    real(dp), parameter :: xw_side = 0.0035_dp  ! side-lobe width
    real(dp), parameter :: delta_side = 0.008_dp ! offset of each side lobe from x0

    real(dp) :: gx
    integer :: i, unit_out

    call init_crust_model("../../TOV_Solver/output/density_profile_1p4Msun.dat", &
                           10.0_dp**4.5_dp, 1.4_dp)

    x_c = 0.90_dp   ! pragmatic domain, wide enough to contain the full assumed Jphi structure
    x_s = 1.0_dp
    print '(A,F10.6,A,F10.6)', " x_c = ", x_c, "   x0 = ", x0

    dx = (x_s - x_c) / real(N-1, dp)
    do i = 1, N
        x(i) = x_c + real(i-1,dp)*dx
    end do

    ! Assumed J_phi: central positive lobe + two negative side lobes,
    ! flanking it symmetrically -- matches the real Fig 7.2 shape.
    do i = 1, N
        gx = x(i)-x0
        Jphi(i) = A_central*exp(-gx*gx/(2.0_dp*xw_central*xw_central)) &
                - A_side*exp(-(x(i)-(x0-delta_side))**2/(2.0_dp*xw_side*xw_side)) &
                - A_side*exp(-(x(i)-(x0+delta_side))**2/(2.0_dp*xw_side*xw_side))
    end do

    ! Solve g_xx - 2g/x^2 = -Jphi, with g=0 at x_c (Dirichlet),
    ! Robin (dg/dx + g/x_s = 0) at x_s, via ghost-point elimination
    ! (same technique used for the diffusion equation's outer BC).
    a(1) = 0.0_dp; b(1) = 1.0_dp; c(1) = 0.0_dp; d(1) = 0.0_dp

    do i = 2, N-1
        a(i) = 1.0_dp/dx**2
        b(i) = -2.0_dp/dx**2 - 2.0_dp/x(i)**2
        c(i) = 1.0_dp/dx**2
        d(i) = -Jphi(i)
    end do

    a(N) = 2.0_dp/dx**2
    b(N) = -2.0_dp/(dx*x_s) - 2.0_dp/dx**2 - 2.0_dp/x_s**2
    c(N) = 0.0_dp
    d(N) = -Jphi(N)

    call thomas_solve(a, b, c, d, g, N)

    print '(A,ES14.6)', " g at surface (unnormalized): ", g(N)
    print '(A,ES14.6)', " g at x_c: ", g(1)

    open(newunit=unit_out, file="output/solve_g_from_jphi.dat", status="replace", action="write")
    write(unit_out, '(A)') "# x   g(x)   Jphi(x)"
    do i = 1, N
        write(unit_out, '(F12.8,3X,ES16.8,3X,ES16.8)') x(i), g(i), Jphi(i)
    end do
    close(unit_out)
    print *, "Written to output/solve_g_from_jphi.dat"

end program SolveGFromJphi
