!*******************************************************************************
module kinds 
implicit none
integer, parameter :: dp=kind(0.d0)  ! double precision
end module kinds
!*******************************************************************************
module root
use kinds
implicit none
private
public :: func, bisect, picard, newton_raphson, secant
contains
!*******************************************************************************
subroutine bisect(a,b,x,fx,iter)
!  Bolzano bisection method
!  f(x)=0 continuous in [a,b]
!  if f(a)f(b)<0 then f(x) has at least one root in [a,b]
!  c is the midpoint between a and b
implicit none
integer, parameter :: imax=50
integer :: i
integer, intent(out) :: iter
real(dp), parameter :: tol=1.e-8_dp
real(dp), intent(inout) :: a, b
real(dp), intent(out) :: x, fx
real(dp) :: c, dx
write(*,101) 
101 format (T10,'Bisection method intermediate output',/,   &
            10x,'iter',8x,'x',25x,'f(x)')
do i=1,imax
  c=(a+b)/2
  if(func(c)==0) then 
    exit
  else if(func(a)*func(c)<0.0_dp) then
    b=c
  else 
    a=c
  end if
  dx=abs(a-b)
  if(dx<tol) exit
!-----intermidiate output  
  write(*,*) i, c, func(c)
end do
x=c
fx=func(c)
iter=i
end subroutine bisect
!*******************************************************************************
subroutine picard(x0,x,fx,iter)
!  f(x)=0 is equivalent to x=g(x)
!  x(i)=g(x(k-1)), i=1,2,....
!  Convergence is not good
implicit none
integer, parameter :: imax=50
integer :: i
integer, intent(out) :: iter
real(dp), parameter :: tol=1.e-5_dp
real(dp), intent(in) :: x0
real(dp), intent(out) :: x, fx
real(dp) :: x1, x2, dx
write(*,102) 
102 format (/,T10,'Picard method intermediate output',/,   &
              10x,'iter',8x,'x',25x,'f(x)')
x1=x0
do i=1,imax
  x2=func(x1)+x1
  dx=abs(x2-x1)
  if(dx<tol) exit
  x1=x2
!-----intermediate output  
  write(*,*) i, x2, func(x2) 
end do
x=x2
fx=func(x2)
iter=i
end subroutine picard
!*******************************************************************************
subroutine newton_raphson(x0,x,fx,iter)
!  Newton-Raphson method
!  using intersection of a tangent line (first derivative)
implicit none
integer, parameter :: imax=50
integer :: i
integer, intent(out) :: iter
real(dp), parameter :: tol=1.e-5_dp
real(dp), intent(in) :: x0
real(dp), intent(out) :: x, fx
real(dp) :: x1, x2, dx
write(*,103) 
103 format (/,T10,'Newton-Raphson method intermediate output',/,   &
              10x,'iter',8x,'x',25x,'f(x)')
x1=x0
do i=1,imax
  x2=x1-func(x1)/dfunc(x1)
  dx=abs(x2-x1)
  if(dx<tol) exit
  x1=x2
!-----intermediate output  
  write(*,*) i, x2, func(x2) 
end do
x=x2
fx=func(x2)
iter=i
end subroutine newton_raphson
!*******************************************************************************
subroutine secant(a,b,x,fx,iter)
!  Secant method based on linear interpolation
!  No need for first derivative
implicit none
integer, parameter :: imax=50
integer :: i
integer, intent(out) :: iter
real(dp), parameter :: tol=1.e-5_dp
real(dp), intent(in) :: a, b
real(dp), intent(out) :: x, fx
real(dp) :: x0, x1, x2, dx
write(*,104) 
104 format (/,T10,'Secant method intermediate output',/,   &
              10x,'iter',8x,'x',25x,'f(x)')
x0=a
x1=b
do i=1,imax
  x2=x1-func(x1)*(x1-x0)/(func(x1)-func(x0))
  dx=abs(x2-x1)
  if(dx<tol) exit
  x1=x2
  x0=x1
!-----intermediate output  
  write(*,*) i, x2, func(x2) 
end do
x=x2
fx=func(x2)  
iter=i
end subroutine secant
!*******************************************************************************
function func(x)
!  input function
implicit none 
real(dp), intent(in) :: x
real(dp) :: func
!func = (x-1)**3._dp - 3._dp*x + 2._dp
!func=1._dp/8 * (35._dp * x**4._dp - 30._dp * x**2._dp + 3._dp)
func = x**3._dp + 2._dp*x**2._dp +10._dp*x -20._dp
end function func
!*******************************************************************************
function dfunc(x)
!  first derivative of the function
implicit none
real(dp), intent(in) :: x
real(dp) :: dfunc
dfunc = 2._dp*x**2 + 4._dp*x + 10._dp 
end function dfunc
!*******************************************************************************
end module root
!*******************************************************************************
program main
use kinds
use root
implicit none
integer :: iter
real(dp) :: a, b, x0, x, fx
a=1.0_dp
b=2.0_dp
x0=0.5_dp
call bisect(a,b,x,fx,iter)
write(*,201) x,fx,iter
call picard(x0,x,fx,iter)
write(*,202) x,fx,iter
call newton_raphson(x0,x,fx,iter)
write(*,203) x,fx,iter
call secant(a,b,x,fx,iter)
write(*,204) x,fx,iter
201 format(/,T4,'Result by bisection:',/, &
             3x,'x=',3x,F15.10,/,                &
             3x,'f(x)',1x,F15.10,/,              &
             3x,'iter=',I5)
202 format(/,T4,'Result by Picard:',/, &
             3x,'x=',3x,F15.10,/,                &
             3x,'f(x)',1x,F15.10,/,              &
             3x,'iter=',I5)
203 format(/,T4,'Result by Newton-Raphson :',/, &
             3x,'x=',3x,F15.10,/,                &
             3x,'f(x)',1x,F15.10,/,              &
             3x,'iter=',I5)
204 format(/,T4,'Result by secant :',/, &
             3x,'x=',3x,F15.10,/,                &
             3x,'f(x)',1x,F15.10,/,              &
             3x,'iter=',I5)
end program main
!*******************************************************************************
