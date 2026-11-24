module binary_search

    use utilities
    implicit none

contains

 

    subroutine find_interval(x, n, x_query, idx)

        implicit none

        integer, intent(in) :: n
        integer, intent(out) :: idx

        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: x_query

        integer :: left
        integer :: right
        integer :: mid

      
        ! boundaries


        if (x_query <= x(1)) then
            idx = 1
            return
        end if

        if (x_query >= x(n)) then
            idx = n-1
            return
        end if

        left = 1
        right = n

        do while (right-left > 1)

            mid = (left+right)/2

            if (x_query >= x(mid)) then
                left = mid
            else
                right = mid
            end if

        end do

        idx = left

    end subroutine find_interval

end module binary_search
