program loop1 

     implicit none 
     integer, parameter :: n = 1000 ! Parameter no change(like const)
      
    real(kind=8), dimension(n) :: x, y ! Two arrays of n-dimension
      integer :: i

      ! loop 
    do i = 1,n
        x(i) = 3.d0 * i
        enddo

    do i = 1,n
        y(n) = 2.d0 * x(i)
        enddo

        print *, "Last y computed:", y(n)

end program loop1
