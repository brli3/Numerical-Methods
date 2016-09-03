!*******************************************************************************
module kinds
implicit none 
integer, parameter :: dp = kind(0.d0)  !  double precision 
end module kinds
!*******************************************************************************
module adams_bashforth
!  Initial value problem.
!  Adams-Bashforth is a linear multi-step method.
!  Use information of previous steps.
!  A linear combination of previous points and derivatives.
!  Need a single-step method (Runge-Kutta) to drive initially.
!  Because first few points are not available
use kinds
implicit none
public :: ab3, ab3_pece, func, rk45_adaptive

contains
!*******************************************************************************
subroutine ab3_pece(x,y,xstop,m,n)
!  Predict-Evaluate-Correct-Evaluate mode of
!  3rd order Adams-Bashforth method
!  calculate an initial guess with a 3-step method 
!  improve it using a 4-step method
!  also evaluate the function func = y'
implicit none 
integer :: i
integer, intent(in) :: m, n
real(dp), intent(inout) :: x
real(dp), intent(in) :: xstop
real(dp), dimension(m), intent(inout) :: y
real(dp), dimension(m) :: y0, y1, y2, y3, f0, f1, f2, f3
real(dp), dimension(m) :: y3_p, f3_e
real(dp) :: x0, x1, x2, x3, h
!-----step size
h = (xstop-x)/n
!-----initial values
x0 = x
y0 = y
f0 = func(x0,y0,m)
!-----compute first point
call rk45_adaptive(x0,y0,x0+h,m)
x1 = x0
y1 = y0
f1 = func(x1,y1,m)
!-----compute second point
call rk45_adaptive(x1,y1,x1+h,m)
x2 = x1
y2 = y1
f2 = func(x2,y2,m)
do i=3,n
!-----predict
  y3_p = y2 + h*(23._dp*f2 - 16._dp*f1 + 5._dp*f0)/12._dp
  x3 = x2 + h
!-----evaluate
  f3_e = func(x3,y3_p,m)
!-----correct
  y3 = y2 + h*(9._dp*f3_e + 19._dp*f2 - 5._dp*f1 + f0)/24._dp
!-----evaluate
  f3 = func(x3,y3,m)
  write(*,'(<m+1>f10.5)') x3, y3
!-----update
  x2 = x2 + h
  y2 = y3
  f0 = f1
  f1 = f2
  f2 = f3
end do

end subroutine ab3_pece
!*******************************************************************************
subroutine ab3(x,y,xstop,m,n)
!  3 step 3rd order Adams-Bashforth method
!  m-----dimension of odes
!  n-----number of integral division 
!  xstop-----end of integration
implicit none
integer :: i
integer, intent(in) :: m, n
real(dp), intent(in) :: xstop
real(dp), intent(inout) :: x
real(dp), dimension(m), intent(inout) :: y
real(dp), dimension(m) :: y0, y1, y2, y3, f0, f1, f2, f3
real(dp) :: x0, x1, x2, x3, h
!-----step size
h = (xstop-x)/n
!-----initial values
x0 = x
y0 = y
f0 = func(x0,y0,m)
!-----compute first point
call rk45_adaptive(x0,y0,x0+h,m)
x1 = x0
y1 = y0
f1 = func(x1,y1,m)
!-----compute second point
call rk45_adaptive(x1,y1,x1+h,m)
x2 = x1
y2 = y1
f2 = func(x2,y2,m)
do i=3,n
  y3 = y2 + h*(23._dp*f2 - 16._dp*f1 + 5._dp*f0)/12._dp
  x3 = x2 + h
  write(*,'(<m+1>f10.5)') x3, y3
!-----update
  x2 = x2 + h
  y2 = y3
  f0 = f1
  f1 = f2
  f2 = func(x2,y2,m)
end do
end subroutine ab3

!*******************************************************************************
function func(x,y,m)
!  system of ordinary differential equations
!  y' = func = f(x)
implicit none
integer, intent(in) :: m
real(dp), intent(in) :: x
real(dp), dimension(m), intent(in) :: y 
real(dp), dimension(m) :: func
func(1) = y(2)*y(3)
func(2) = -y(1)*y(3)
func(3) = -0.4_dp*y(1)*y(2)
end function func
!*******************************************************************************
subroutine rk45_adaptive(x,y,xstop,m,hstep)
!  Adaptive Runge-Kutta method with embedded Runge-Kutta-Fehlberg formula
implicit none
integer :: i
integer, intent(in) :: m
real(dp), intent(in) :: xstop
real(dp), intent(in), optional :: hstep
real(dp), intent(inout) :: x
real(dp), dimension(m), intent(inout) :: y
real(dp), dimension(m) :: dy
real(dp) :: h, delta, error
real(dp) :: hnext
real(dp), parameter :: tol = 1.0e-6_dp
!-----initial step size
if(present(hstep)) then
  h = min(hstep, xstop-x)
else
  h = min(1.0e-1_dp, xstop-x)
end if
!-----integration
do while(x<xstop)
!-----get increment, error and step correction
  call rkf(x,y,h,m,dy,error)
  delta = 0.9_dp*(tol/error)**0.2_dp
!-----accept current step if error is within tolerance
  if(error>tol) then
!-----recalculate current step with reduced step size  
    h = h * delta
    call rkf(x,y,h,m,dy,error)
  end if
  y = y + dy
  x = x + h
!-----next step size with a margin of safety
  delta = 0.9_dp*(tol/error)**0.2_dp
  if(error/=0.0_dp) then
    hnext = h * delta
  else
    hnext = h
  end if
!-----check if next step is the last one
  if((h>0._dp) == ((x+hnext)>=xstop)) then
    hnext = xstop - x
  end if
  h = hnext
!-----print
  write(*,'(<m+1>E14.5)') x, y
end do

contains

subroutine rkf(x,y,h,m,dy,error)
!  5th order Runge-Kutta-Fehlberg formulas by Cash-Karp
implicit none
integer :: i
integer, intent(in) :: m
real(dp), intent(in) :: x, h
real(dp), dimension(m), intent(in) :: y
real(dp), dimension(m), intent(out):: dy
real(dp), intent(out) :: error
real(dp), dimension(m) :: e
real(dp), dimension(6,m) :: k
real(dp), dimension(2:6) :: a
real(dp), dimension(1) :: b2
real(dp), dimension(2) :: b3
real(dp), dimension(3) :: b4
real(dp), dimension(4) :: b5
real(dp), dimension(5) :: b6
real(dp), dimension(6) :: c, d
real(dp) :: s
!-----list of coefficients
a = ([1._dp/5, 3._dp/10, 3._dp/5, 1._dp, 7._dp/8])
b2 = ([1._dp/5])
b3 = ([3._dp/40, 9._dp/40])
b4 = ([3._dp/10, -9._dp/10, 6._dp/5])
b5 = ([-11._dp/54, 5._dp/2, -70._dp/27, 35._dp/27])
b6 = ([1631._dp/55296, 175._dp/512, 575._dp/13824, 44275._dp/110592, 253._dp/4096 ])
c = ([37._dp/378, 0._dp, 250._dp/621, 125._dp/594, 0._dp, 512._dp/1771])
d = ([2825._dp/27648, 0._dp, 18575._dp/48384, 13525._dp/55296, 277._dp/14336, 1._dp/4])
!-----formulas
k(1,:) = h*func(x, y, m)
k(2,:) = h*func(x+a(2)*h, y+b2(1)*k(1,:), m)
k(3,:) = h*func(x+a(3)*h, y+(b3(1)*k(1,:)+b3(2)*k(2,:)), m)
k(4,:) = h*func(x+a(4)*h, y+(b4(1)*k(1,:)+b4(2)*k(2,:)+b4(3)*k(3,:)), m)
k(5,:) = h*func(x+a(5)*h, y+(b5(1)*k(1,:)+b5(2)*k(2,:)+b5(3)*k(3,:)+b5(4)*k(4,:)), m)
k(6,:) = h*func(x+a(6)*h, y+(b6(1)*k(1,:)+b6(2)*k(2,:)+b6(3)*k(3,:)+b6(4)*k(4,:)+b6(5)*k(5,:)), m)
!-----initialize arrays {dy} and {e} 
dy = 0._dp
e = 0._dp
!-----compute solution increment {dy} and per-step error {e}
do i=1,6
!-----5th order formula
  dy = dy + c(i)*k(i,:)
!-----use 4th order formula implicitly to estimate the truncation error  
  e = e + (c(i)-d(i))*k(i,:)
end do  
!-----rms error
error = sqrt(sum(e**2)/m)
end subroutine rkf
end subroutine rk45_adaptive
!*******************************************************************************
end module adams_bashforth
!*******************************************************************************
program main 
use kinds 
use adams_bashforth
implicit none 
integer :: n = 100
real(dp) :: x, xstop, y(3)
x = 1._dp
xstop = 6.0_dp
y = [1._dp, 1._dp, 1._dp]
write(*,*) 'Result by Adams-Bashforth method: '
call ab3_pece(x,y,xstop,3,n)
end program main
!*******************************************************************************
