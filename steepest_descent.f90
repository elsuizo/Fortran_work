!*******************************************************
!* Program to demonstrate the use of multi-dimensional *
!*     Steepest Descent Optimization subroutine        *
!* --------------------------------------------------- *
!*  Reference: BASIC Scientific Subroutines, Vol. II   *
!*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris.   *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Example:   Find a local maximum of                  *
!*            F(x,y,z) = sin(x)+2*cos(y)-sin(z)        *
!*                                                     *
!* SAMPLE RUN:                                         *
!*                                                     *
!* How many dimensions: 3                              *
!*                                                     *
!* Convergence criterion: .000000001                   *
!*                                                     *
!* Maximum number of iterations: 50                    *
!*                                                     *
!* Starting constant: 1                                *
!*                                                     *
!* Input the starting point:                           *
!*                                                     *
!*     X[1] = 1                                        *
!*     X[2] = 1                                        *
!*     X[3] = 1                                        *
!*                                                     *
!* The results are:                                    *
!*                                                     *
!*     X(1) = 1.5708084                                *
!*     X(2) = -0.0000139                               *
!*     X(3) = -4.7123879                               *
!*                                                     *
!* Local maximum = 4.0000000                           *
!*                                                     *
!* The number of steps was: 33                         *
!*                                                     *
!*******************************************************
PROGRAM DEMO_STEEPDA

real*8  D(3), Y(3)
real*8  X(10),X1(10)
real*8  e,xk,Eval
integer i,l,m,n

  write(*,"(/' How many dimensions: ')",advance='no')
  read *, l
  write(*,"(' Convergence criterion: ')",advance='no'); 
  read *, e
  write(*,"(' Maximum number of iterations: ')",advance='no')
  read *, m
  write(*,"(' Starting constant: ')",advance='no')
  read *, xk
  print *,' '
  print *,'Input the starting point:'
  print *,' '
  do i=1, l
    write(*,40,advance='no') i
    read *, X(i)
  end do

  call Steepda(l,e,m,xk,X,X1,Y,D,n)   !Call steepest descent optimization subroutine

  print *,' '
  print *,'The results are:'
  print *,' '
  do i=1, l
    write(*,50) i,X(i)
  end do
  write(*,60) Eval(X)
  write(*,70) n

  stop

40 format('    X(',i2,') = ')
50 format('    X(',i2,') = ',F10.7)
60 format(/' Local maximum = ',F10.7)
70 format(' The number of iterations was ',i2//)

end


!********************************************
! Function subroutine              
real*8 Function Eval(X)
  real*8 X(10)
  Eval = dsin(X(1))+2.0*dcos(X(2))-dsin(X(3))
end
!********************************************



  ! Functions called by Steepda()
  Subroutine Utilit1(l,dd,D)
    integer i,l
    real*8 dd,D(3)
    ! Find the magnitude of the gradient
    dd=0.d0
    do i=1, l 
      dd = dd + D(i)*D(i)
    end do
    dd=dsqrt(dd)
    return
  end

  Subroutine Utilit2(l,xk,dd,D,Y,X,X1)
    integer i,l
    real*8 xk,dd,D(3),Y(3),X(10),X1(10),Eval
    ! Update the X(i) 
    do i=1, l
      ! Save old values
      X1(i)=X(i)
      X(i) = X(i) + xk*D(i)/dd
    end do
    Y(3)=Eval(X)
    return
  end

  ! Find approximations of partial derivatives D(i)
  ! by finite differences
  Subroutine Derivatives(l,dd,xk,D,X,Y)
    real*8 a,b,dd,xk,yy,D(3),X(10),Y(3),Eval
    integer i,l
    do i=1, l
      ! Save X(i)
      a=X(i)
      ! Find increment
      b=D(i)*xk/(2.d0*dd)
      ! Move increment in X(i)
      X(i)=X(i)+b
      ! Obtain yy
      yy=Eval(X)
      ! Guard against divide by zero near maximum
      if (dabs(b) < 1.d-12) b=1.d-12
      ! Update D(i)
      D(i)=(yy-Y(3))/b
      ! Guard against locked up derivative
      if (dabs(D(i)) < 1.d-5) D(i)=1.d-5
      ! Restore X(i) and yy
      X(i)=a; yy=Y(3)
    end do
    ! Obtain dd
    call Utilit1(l,dd,D)
    return
  end


!***************************************************
!*    Steepest descent optimization subroutine     *
!* ----------------------------------------------- *
!* This routine finds the local maximum or minimum *
!* of an L-dimensional function using the method   *
!* of steepest decent, or the gradient.            *
!* The function must be available in Function      *
!* Eval(X).In this version, finite differences are *
!* used to calculate the L partial derivatives.    *
!* ----------------------------------------------- *
!* INPUTS:                                         *
!*   l - The dimension of function to study        *
!*   e - The convergence criteria                  *
!*   m - The maximum number of iterations          *
!*   xk - A starting constant                      *
!*   X(i) - Initial values of variables            *
!* OUTPUTS:                                        *
!*   X(i) - The locally optimum set                *
!*   Eval - The local maximum found                *
!*   n - The number of iterations performed,       *
!***************************************************
  Subroutine Steepda(l,e,m,xk,X,X1,Y,D,n)
  ! Labels: 50,51,100,200
  real*8 MACHEPS
  integer i,l,m,n
  real*8 dd,e,xk,X(10),X1(10),D(3),Y(3),Eval
  n=0
  MACHEPS=1.d-15
  !The routine needs three values of Y to get started
  !Generate starting D(i) values
  !These are not even good guesses and slow the program a little
  dd=1.d0
  D(1)=1.d0/dsqrt(dfloat(l))
  do i=2, l  
    D(i)=D(i-1)
  end do
  ! Start initial probe
  do i=1, l
    ! Obtain yy and D[i]
    Y(i)=Eval(X)
    ! Update X[i]
    call Utilit1(l,dd,D)
    call Utilit2(l,xk,dd,D,Y,X,X1)
  end do
  ! We now have a history to base the subsequent search on
  ! Accelerate search if approach is monotonic 
50 if (dabs(Y(2)-Y(1))<MACHEPS) goto 51
  if ((Y(3)-Y(2))/(Y(2)-Y(1))>0.d0) xk=xk*1.2d0
  ! Decelerate if heading the wrong way
51 if (Y(3)<Y(2)) xk=xk/2.d0
  ! Update the Y[i] if value has decreased
  if (Y(3)>Y(2)) goto 100
  ! Restore the X[i]
  do i=1, l
    X(i)=X1(i)
  end do
  goto 200
100 Y(1)=Y(2); Y(2)=Y(3)
  ! Obtain new values
200 Y(3)=Eval(X)
  call Derivatives(l,dd,xk,D,X,Y) ! Get D(i)
  !if dd=0 then the precision limit of the computer has been reached
  if (dd < MACHEPS) return
  ! Update X[i]
  call Utilit2(l,xk,dd,D,Y,X,X1)
  ! Check for maximum iterations and convergence
  n=n+1
  if (n>=m) return
  if (dabs(Y(3)-Y(2))<e) return
  ! Try another iteration
  goto 50
end ! Steepda()


! End of file steepda.f90
