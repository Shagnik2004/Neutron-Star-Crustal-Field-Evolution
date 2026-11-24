program konar_ch7_gprofile

   implicit none

   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: n = 500

   integer :: i
   real(dp) :: x
   real(dp) :: g
   real(dp) :: logg

   open(unit=10,file="output/konar_ch7_gprofile.dat",status="replace")

   write(10,'(A)') "# x   g(x)   log10(g)"

   do i=0,n

      x = real(i,dp)/real(n,dp)

      !---------------------------------------------------------
      ! Chapter 7 trial profile (quintic smoothstep)
      !
      ! x = 0 : crust-core boundary
      ! x = 1 : stellar surface
      !---------------------------------------------------------

      g = 10.d0*x**3 - 15.d0*x**4 + 6.d0*x**5

      if (g <= 1.d-15) then
         logg = -15.d0
      else
         logg = log10(g)
      end if

      write(10,'(3ES20.10)') x,g,logg

   end do

   close(10)

   print *,'Wrote output/konar_ch7_gprofile.dat'

end program konar_ch7_gprofile
