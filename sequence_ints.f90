program arraycons
  implicit none
  integer :: i
  real :: a(10) = (/(i, i=0,9, 1)/)
  print *, a
end program arraycons
