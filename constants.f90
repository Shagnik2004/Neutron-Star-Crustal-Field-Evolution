module tov_constants

    use utilities
    implicit none

    real(dp), parameter :: G_grav  = 6.6743e-8_dp        ! cm^3 g^-1 s^-2
    real(dp), parameter :: c_light = 2.99792458e10_dp    ! cm/s
    real(dp), parameter :: Msun    = 1.98892e33_dp       ! g
    real(dp), parameter :: pi      = 3.14159265358979323846_dp

end module tov_constants