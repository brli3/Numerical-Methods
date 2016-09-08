!*******************************************************************************
module kinds
implicit none 
integer, parameter :: dp = kind(0.d0)  !  double precision 
end module kinds
!*******************************************************************************
module runge_kutta
!  Initial value problem.
!  Runge-Kutta is a explicit single step method.
!  Only one previous point and its derivative 
!  are used to determine current value.
!  Both classical and adaptive versions are included.
use kinds 
implicit none 
public :: func, rk4, rk45_adaptive

contains
!*******************************************************************************
subroutine rk4(x,y,xstop,m,n)
!  First order odes with 4th order Runge-Kutta method 
!  Initial value problem {y}' = {f(x,{y})}
!  x-----independent variable
!  y-----dependent variable y(1), y(2), ..., y(m)
!  xstop-----end of integration
!  m-----dimension of odes
!  n-----number of integral steps
implicit none 
integer :: i
integer, intent(in) :: m, n
real(dp), intent(in) :: xstop
real(dp), intent(inout) :: x
real(dp), dimension(m), intent(inout) :: y
real(dp) :: h
h = (xstop-x)/n
do i=1,n
  y = y + rk4_coef(x,y,h,m)
  x = x + h
  write(*,'(<m+1>E14.5)')  x, y
end do

contains

function rk4_coef(x,y,h,m)
!  Coefficients and increament with classical 4th order Runge-Kutta method
implicit none 
integer :: i
integer, intent(in) :: m
real(dp), intent(in) :: x, h
real(dp), dimension(m), intent(in) :: y
real(dp), dimension(m) :: rk4_coef
real(dp), dimension(4,m) :: k
k(1,:) = h*func(x, y, m)
k(2,:) = h*func(x+0.5_dp*h, y+0.5_dp*k(1,:), m)
k(3,:) = h*func(x+0.5_dp*h, y+0.5_dp*k(2,:), m)
k(4,:) = h*func(x+h, y+k(3,:), m)
rk4_coef = (k(1,:) + 2.0_dp*k(2,:) + 2.0_dp*k(3,:) + k(4,:))/6.0_dp
end function rk4_coef
end subroutine rk4
!*******************************************************************************
subroutine rk45_adaptive(xstart,xstop,ystart,ystop,m,hstep)
!  Adaptive Runge-Kutta method with embedded Runge-Kutta-Fehlberg formula
implicit none
integer :: i
integer, intent(in) :: m
real(dp), intent(in) :: xstart
real(dp), intent(in) :: xstop
real(dp), dimension(m), intent(in) :: ystart
real(dp), dimension(m), intent(out) :: ystop
real(dp), intent(in), optional :: hstep
real(dp), dimension(m) :: y, dy
real(dp) :: x, h, delta, error
real(dp) :: hnext
real(dp), parameter :: tol = 1.0e-6_dp
x = xstart
y = ystart
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
ystop = y

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
function func(x,y,m)
!  m-dimension 1st order odes
!  Any higher nth order ode can always be converted to n 1st order odes
implicit none
integer, intent(in) :: m
real(dp), intent(in) :: x
real(dp), dimension(m), intent(in) :: y
real(dp), dimension(m) :: func
!func(1) = y(2)*y(3)
!func(2) = -y(1)*y(3)
!func(3) = -0.51*y(1)*y(2)
!func(1) = y(2)
!func(2) = -9.80665_dp + (65.351e-3_dp)*y(2)**2._dp*exp(-10.53e-5_dp*y(1))
func(1) = y(2)
func(2) = -4.75_dp*y(1) - 10.0_dp*y(2)
end function func
!*******************************************************************************
end module runge_kutta
!*******************************************************************************
program main 
use kinds
use runge_kutta
implicit none
integer :: n
real(dp) :: xstart, h, xstop, tol
real(dp), dimension(2) :: ystart, ysol
xstart = 0.0_dp
xstop = 11.1_dp
ystart = ([-9.0_dp, 0.0_dp])
h = 0.1_dp
write(*,*) 'Result by Runge-Kutta method: '
call rk45_adaptive(xstart,xstop,ystart,ysol,2)
end program main
!*******************************************************************************
