module psi_hall_omp
    ! OpenMP versions of adi_step_2d_hall (poloidal P) and
    ! adi_step_2d_T_full (toroidal T).  The arithmetic is IDENTICAL to
    ! the validated serial routines in psi_diffusion_solver_2d.f90; only
    ! the independent line solves (one tridiagonal system per theta-line
    ! in the x-sweep, one per x-line in the theta-sweep) are distributed
    ! over threads.  Every line solve reads only the old field and writes
    ! only its own line, so the result does not depend on the thread
    ! count (identical to a 1-thread run).
    use utilities
    use tridiagonal_solver
    implicit none
contains

    subroutine adi_step_2d_hall_omp(psi, x, theta, Nx, Ntheta, dx, dtheta, x_s, S, Vx, Vtheta, dt, theta_adi_in)
        implicit none
        integer, intent(in)     :: Nx, Ntheta
        real(dp), intent(inout) :: psi(Nx, Ntheta)
        real(dp), intent(in)    :: x(Nx), theta(Ntheta), dx, dtheta, x_s, dt
        real(dp), intent(in)    :: S(Nx), Vx(Nx,Ntheta), Vtheta(Nx,Ntheta)
        real(dp), intent(in)    :: theta_adi_in

        real(dp), allocatable :: psi_half(:,:)
        real(dp) :: a(Nx), b(Nx), c(Nx), d(Nx), sol(Nx)
        real(dp) :: at(Ntheta), bt(Ntheta), ct(Ntheta), dt_rhs(Ntheta), solt(Ntheta)
        real(dp) :: rx, rtheta, cot_j, theta_adi, anti_adi, adv_x, adv_theta
        integer :: i, j

        allocate(psi_half(Nx, Ntheta))
        theta_adi = theta_adi_in
        anti_adi = 1.0_dp - theta_adi

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(j,i,a,b,c,d,sol,rx,rtheta,cot_j,adv_x)
        do j = 2, Ntheta-1
            cot_j = cos(theta(j)) / sin(theta(j))
            a(1) = 0.0_dp; b(1) = 1.0_dp; c(1) = 0.0_dp; d(1) = 0.0_dp
            do i = 2, Nx-1
                rx = theta_adi * dt*S(i) / (dx*dx)
                adv_x = theta_adi * dt*Vx(i,j) / dx
                a(i) = -rx
                b(i) = 1.0_dp + 2.0_dp*rx
                c(i) = -rx
                if (Vx(i,j) >= 0.0_dp) then
                    a(i) = a(i) - adv_x
                    b(i) = b(i) + adv_x
                else
                    b(i) = b(i) - adv_x
                    c(i) = c(i) + adv_x
                end if
                rtheta = anti_adi * dt*S(i) / (x(i)**2 * dtheta**2)
                d(i) = psi(i,j) &
                     + rtheta*( (psi(i,j+1) - 2.0_dp*psi(i,j) + psi(i,j-1)) &
                                - cot_j*dtheta*(psi(i,j+1)-psi(i,j-1))/2.0_dp )
            end do
            rx = theta_adi * dt*S(Nx) / (dx*dx)
            a(Nx) = -2.0_dp*rx
            b(Nx) = 1.0_dp + 2.0_dp*rx*(1.0_dp + dx/x_s)
            c(Nx) = 0.0_dp
            rtheta = anti_adi * dt*S(Nx) / (x_s**2 * dtheta**2)
            d(Nx) = psi(Nx,j) &
                  + rtheta*( (psi(Nx,j+1) - 2.0_dp*psi(Nx,j) + psi(Nx,j-1)) &
                             - cot_j*dtheta*(psi(Nx,j+1)-psi(Nx,j-1))/2.0_dp )
            call thomas_solve(a, b, c, d, sol, Nx)
            psi_half(:,j) = sol
        end do
        !$omp end parallel do
        psi_half(:,1) = 0.0_dp
        psi_half(:,Ntheta) = 0.0_dp

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(i,j,at,bt,ct,dt_rhs,solt,rx,rtheta,cot_j,adv_theta)
        do i = 2, Nx-1
            at(1) = 0.0_dp; bt(1) = 1.0_dp; ct(1) = 0.0_dp; dt_rhs(1) = 0.0_dp
            do j = 2, Ntheta-1
                cot_j = cos(theta(j)) / sin(theta(j))
                rtheta = theta_adi * dt*S(i) / (x(i)**2 * dtheta**2)
                adv_theta = theta_adi * dt*Vtheta(i,j) / dtheta
                at(j) = -rtheta * (1.0_dp + cot_j*dtheta/2.0_dp)
                bt(j) = 1.0_dp + 2.0_dp*rtheta
                ct(j) = -rtheta * (1.0_dp - cot_j*dtheta/2.0_dp)
                if (Vtheta(i,j) >= 0.0_dp) then
                    at(j) = at(j) - adv_theta
                    bt(j) = bt(j) + adv_theta
                else
                    bt(j) = bt(j) - adv_theta
                    ct(j) = ct(j) + adv_theta
                end if
                rx = anti_adi * dt*S(i) / (dx*dx)
                dt_rhs(j) = rx*(psi_half(i-1,j) - 2.0_dp*psi_half(i,j) + psi_half(i+1,j)) &
                          + psi_half(i,j)
            end do
            at(Ntheta) = 0.0_dp; bt(Ntheta) = 1.0_dp; ct(Ntheta) = 0.0_dp; dt_rhs(Ntheta) = 0.0_dp
            call thomas_solve(at, bt, ct, dt_rhs, solt, Ntheta)
            psi(i,:) = solt
        end do
        !$omp end parallel do
        psi(1,:) = psi_half(1,:)
        psi(Nx,:) = psi_half(Nx,:)
        psi(:,1) = 0.0_dp
        psi(:,Ntheta) = 0.0_dp
        deallocate(psi_half)
    end subroutine adi_step_2d_hall_omp


    subroutine adi_step_2d_T_full_omp(Tpot, x, theta, Nx, Ntheta, dx, dtheta, S, Wx, Wtheta, dt, theta_adi_in)
        implicit none
        integer, intent(in)     :: Nx, Ntheta
        real(dp), intent(inout) :: Tpot(Nx, Ntheta)
        real(dp), intent(in)    :: x(Nx), theta(Ntheta), dx, dtheta, dt
        real(dp), intent(in)    :: S(Nx), Wx(Nx,Ntheta), Wtheta(Nx,Ntheta)
        real(dp), intent(in)    :: theta_adi_in

        real(dp), allocatable :: T_half(:,:)
        real(dp) :: a(Nx), b(Nx), c(Nx), d(Nx), sol(Nx)
        real(dp) :: at(Ntheta), bt(Ntheta), ct(Ntheta), dt_rhs(Ntheta), solt(Ntheta)
        real(dp) :: rx, rtheta, cot_j, theta_adi, anti_adi, Sprime(Nx), adv_x, adv_theta
        integer :: i, j

        allocate(T_half(Nx, Ntheta))
        theta_adi = theta_adi_in
        anti_adi = 1.0_dp - theta_adi

        Sprime(1) = (S(2)-S(1))/dx
        do i = 2, Nx-1
            Sprime(i) = (S(i+1)-S(i-1))/(2.0_dp*dx)
        end do
        Sprime(Nx) = (S(Nx)-S(Nx-1))/dx

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(j,i,a,b,c,d,sol,rx,rtheta,cot_j,adv_x)
        do j = 2, Ntheta-1
            cot_j = cos(theta(j)) / sin(theta(j))
            a(1) = 0.0_dp; b(1) = 1.0_dp; c(1) = 0.0_dp; d(1) = 0.0_dp
            do i = 2, Nx-1
                rx = theta_adi * dt*S(i) / (dx*dx)
                adv_x = theta_adi * dt*(Sprime(i)+Wx(i,j)) / dx
                a(i) = -rx
                b(i) = 1.0_dp + 2.0_dp*rx
                c(i) = -rx
                if ((Sprime(i)+Wx(i,j)) >= 0.0_dp) then
                    a(i) = a(i) - adv_x
                    b(i) = b(i) + adv_x
                else
                    b(i) = b(i) - adv_x
                    c(i) = c(i) + adv_x
                end if
                rtheta = anti_adi * dt*S(i) / (x(i)**2 * dtheta**2)
                d(i) = Tpot(i,j) &
                     + rtheta*( (Tpot(i,j+1) - 2.0_dp*Tpot(i,j) + Tpot(i,j-1)) &
                                - cot_j*dtheta*(Tpot(i,j+1)-Tpot(i,j-1))/2.0_dp )
            end do
            a(Nx) = 0.0_dp; b(Nx) = 1.0_dp; c(Nx) = 0.0_dp; d(Nx) = 0.0_dp
            call thomas_solve(a, b, c, d, sol, Nx)
            T_half(:,j) = sol
        end do
        !$omp end parallel do
        T_half(:,1) = 0.0_dp
        T_half(:,Ntheta) = 0.0_dp

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(i,j,at,bt,ct,dt_rhs,solt,rx,rtheta,cot_j,adv_theta)
        do i = 2, Nx-1
            at(1) = 0.0_dp; bt(1) = 1.0_dp; ct(1) = 0.0_dp; dt_rhs(1) = 0.0_dp
            do j = 2, Ntheta-1
                cot_j = cos(theta(j)) / sin(theta(j))
                rtheta = theta_adi * dt*S(i) / (x(i)**2 * dtheta**2)
                adv_theta = theta_adi * dt*Wtheta(i,j) / dtheta
                at(j) = -rtheta * (1.0_dp + cot_j*dtheta/2.0_dp)
                bt(j) = 1.0_dp + 2.0_dp*rtheta
                ct(j) = -rtheta * (1.0_dp - cot_j*dtheta/2.0_dp)
                if (Wtheta(i,j) >= 0.0_dp) then
                    at(j) = at(j) - adv_theta
                    bt(j) = bt(j) + adv_theta
                else
                    bt(j) = bt(j) - adv_theta
                    ct(j) = ct(j) + adv_theta
                end if
                rx = anti_adi * dt*S(i) / (dx*dx)
                dt_rhs(j) = rx*(T_half(i-1,j) - 2.0_dp*T_half(i,j) + T_half(i+1,j)) &
                          + T_half(i,j)
            end do
            at(Ntheta) = 0.0_dp; bt(Ntheta) = 1.0_dp; ct(Ntheta) = 0.0_dp; dt_rhs(Ntheta) = 0.0_dp
            call thomas_solve(at, bt, ct, dt_rhs, solt, Ntheta)
            Tpot(i,:) = solt
        end do
        !$omp end parallel do
        Tpot(1,:) = T_half(1,:)
        Tpot(Nx,:) = T_half(Nx,:)
        Tpot(:,1) = 0.0_dp
        Tpot(:,Ntheta) = 0.0_dp
        deallocate(T_half)
    end subroutine adi_step_2d_T_full_omp

end module psi_hall_omp
