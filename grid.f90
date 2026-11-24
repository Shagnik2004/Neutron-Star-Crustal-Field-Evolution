module grid

    use utilities
    implicit none

contains

    !--------------------------------------------------------------
    ! Builds a uniform 1D grid x(0:N) from x_min to x_max with
    ! spacing dx (caller-specified, not derived -- keeps the
    ! stability condition r = D*dt/dx^2 easy to reason about
    ! explicitly in main.f90).
    !--------------------------------------------------------------
    subroutine build_grid_fixed_dx(x_min, x_max, dx, x, N)

        implicit none

        real(dp), intent(in)  :: x_min, x_max, dx
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out)  :: N

        integer :: j

        N = nint((x_max - x_min) / dx)

        allocate(x(0:N))
        do j = 0, N
            x(j) = x_min + real(j, dp) * dx
        end do

    end subroutine build_grid_fixed_dx

end module grid
