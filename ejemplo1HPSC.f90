program example1
      implicit none ! all vars are declarate 
      real (kind=8) :: x,y,z ! Doble precision 

      x = 3.d0 ! Doble precision 
      y = 1.d-1 ! is like .1  
      z = x + y
      print *, "z = ", z
      end program example1
