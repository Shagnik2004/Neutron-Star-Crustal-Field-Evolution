module psi_flux_omp
    ! Toroidal (T) Ohmic ADI step in CONSERVATIVE form.
    ! Exact Ohmic equation for T = r sin(theta) B_phi with eta(r):
    !     dT/dt = d/dx( S dT/dx ) + (S/x^2) ( T_thth - cot(theta) T_th )
    ! x-fluxes S_{i+-1/2} (T_{i+-1}-T_i)/dx are used in BOTH half steps, so there
    ! is no S'(x) "advection" term, no upwinding and the operator is an M-matrix.
    use utilities
    use tridiagonal_solver
    implicit none
contains

    subroutine adi_step_2d_T_flux_omp(Tpot, x, theta, Nx, Ntheta, dx, dtheta, S, dt, theta_adi_in)
        implicit none
        integer, intent(in)     :: Nx, Ntheta
        real(dp), intent(inout) :: Tpot(Nx, Ntheta)
        real(dp), intent(in)    :: x(Nx), theta(Ntheta), dx, dtheta, dt
        real(dp), intent(in)    :: S(Nx)
        real(dp), intent(in)    :: theta_adi_in

        real(dp), allocatable :: T_half(:,:)
        real(dp) :: a(Nx), b(Nx), c(Nx), d(Nx), sol(Nx)
        real(dp) :: at(Ntheta), bt(Ntheta), ct(Ntheta), dt_rhs(Ntheta), solt(Ntheta)
        real(dp) :: rxm, rxp, rtheta, cot_j, theta_adi, anti_adi, Sm, Sp
        integer :: i, j

        allocate(T_half(Nx, Ntheta))
        theta_adi = theta_adi_in
        anti_adi = 1.0_dp - theta_adi

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(j,i,a,b,c,d,sol,rxm,rxp,rtheta,cot_j,Sm,Sp)
        do j = 2, Ntheta-1
            cot_j = cos(theta(j)) / sin(theta(j))
            a(1) = 0.0_dp; b(1) = 1.0_dp; c(1) = 0.0_dp; d(1) = 0.0_dp
            do i = 2, Nx-1
                Sm = 0.5_dp*(S(i-1)+S(i))
                Sp = 0.5_dp*(S(i)+S(i+1))
                rxm = theta_adi * dt*Sm / (dx*dx)
                rxp = theta_adi * dt*Sp / (dx*dx)
                a(i) = -rxm
                b(i) = 1.0_dp + rxm + rxp
                c(i) = -rxp
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
        !$omp private(i,j,at,bt,ct,dt_rhs,solt,rxm,rxp,rtheta,cot_j,Sm,Sp)
        do i = 2, Nx-1
            at(1) = 0.0_dp; bt(1) = 1.0_dp; ct(1) = 0.0_dp; dt_rhs(1) = 0.0_dp
            Sm = 0.5_dp*(S(i-1)+S(i))
            Sp = 0.5_dp*(S(i)+S(i+1))
            rxm = anti_adi * dt*Sm / (dx*dx)
            rxp = anti_adi * dt*Sp / (dx*dx)
            do j = 2, Ntheta-1
                cot_j = cos(theta(j)) / sin(theta(j))
                rtheta = theta_adi * dt*S(i) / (x(i)**2 * dtheta**2)
                at(j) = -rtheta * (1.0_dp + cot_j*dtheta/2.0_dp)
                bt(j) = 1.0_dp + 2.0_dp*rtheta
                ct(j) = -rtheta * (1.0_dp - cot_j*dtheta/2.0_dp)
                dt_rhs(j) = T_half(i,j) + rxm*(T_half(i-1,j) - T_half(i,j)) &
                                        + rxp*(T_half(i+1,j) - T_half(i,j))
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
    end subroutine adi_step_2d_T_flux_omp

end module psi_flux_omp
