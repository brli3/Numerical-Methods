!*******************************************************************************
module kinds 
implicit none
integer, parameter :: dp=kind(0.d0)  ! double precision
end module kinds
!*******************************************************************************
module shoot_method
!  Boundary value problems
!  Shooting method
!  Guess unknown derivatives at one boundary and 
!  solve an initial value problem.
!  March to the other boundary.
!  Adjust and try until the prescribed boundary value.
use kinds
implicit none
private
public :: shoot_integrate, rk45_adaptive, func 
contains
!*******************************************************************************
subroutine shoot_integrate(xstart,xstop,ystop,yb,u1,u2,m)
!  Shooting integration for 2nd order ode (can be extended for higher order.)
!  with classical Runge-Kutta method. 
!  xstart-----start of integration
!  xstop-----end of integration
!  u1-----1st trial value of unknown init. cond.
!  u2-----2st trial value of unknown init. cond.
!  m-----dimension of odes
!  x,y-----integration output
implicit none
integer, intent(in) :: m
!real(dp), intent(inout) :: u1, u2, xstop 
!real(dp), intent(inout) :: xstart
real(dp), intent(in) :: xstart, xstop
!real(dp), intent(in) :: ystart 
real(dp), intent(in) :: yb
real(dp), dimension(m), intent(out) :: ystop
real(dp), intent(in) :: u1, u2
real(dp) :: usol, xcopy, ycorr(m)
!-----compute correct initial condition
call ridder_res(u1,u2,yb,usol)
write(*,'(t4,a31,/,3(f12.5))') 'Two guesses and computed root: ', u1, u2, usol
call init_cond(usol,ycorr)
call rk45_adaptive(xstart,xstop,ycorr,ystop,m)
contains
!*******************************************************************************
subroutine init_cond(gs,yinit)
!  gs-----a guess of y'
!  yinit-----inital values of the system
real(dp), intent(in) :: gs  ! unknown
real(dp), dimension(m), intent(out) :: yinit
yinit = [0._dp,gs]
end subroutine init_cond
!*******************************************************************************
function resdiff(gs,ybc)
!  y'' = f(x,y,y'); y(a) = alpha; y(b) = beta.
!  boundary residual (difference between the computed and specified
!  boundary value at x = b.
!  computed y(b) = theta(u) 
!  res(u) = theta(u) - beta = 0
implicit none 
real(dp), intent(in) :: gs 
real(dp), intent(in) :: ybc
real(dp), dimension(m) :: y1, y2
real(dp) :: resdiff
call init_cond(gs,y1)
call rk45_adaptive(xstart,xstop,y1,y2,m)
resdiff = y2(m-1) - ybc
end function resdiff
!*******************************************************************************
subroutine ridder_res(a,b,ybc,xsol)
!  Ridder's method using 3 points
!  calculate residual difference function 
implicit none
integer, parameter :: imax=50
integer :: i
integer :: iter
real(dp), parameter :: tol=1.e-5_dp
real(dp), intent(in) :: a, b, ybc
real(dp), intent(out) :: xsol
real(dp) :: fx
real(dp) :: x0, x1, x2, dx, fa, fb, fc, c, sq
!write(*,105)
!105 format (/,T10,'Ridder method intermediate output',/,   &
!              10x,'iter',8x,'x',25x,'f(x)')
x0 = a
x1 = b
fa = resdiff(x0,ybc)
if(fa==0._dp) then
  xsol = x0
  fx = fa
  return
end if
fb = resdiff(x1,ybc)
if(fb==0._dp) then
  xsol = x1
  fx = fb
  return
end if
if(fa*fb>0._dp) stop 'Ridder returns; Root is not bracketed'
do i=1,imax
  c = 0.5_dp*(x0+x1) 
  fc = resdiff(c,ybc)
  sq = sqrt(fc**2-fa*fb) 
  if(sq==0._dp) return
  dx = (c-x0)*fc/sq
  if((fa-fb)<0._dp) dx = -dx
  xsol = c + dx
  fx = resdiff(xsol,ybc)
!-----test for convergence
  if(i>0) then
    if(abs(xsol-c)<tol*max(abs(xsol),1.0_dp)) return
    c = xsol
  end if
!-----Re-bracket the root as tightly as possible
  if(fc*fx>0._dp) then
    if(fa*fx<0._dp) then
      x1 = xsol
      fb = fx
    else
      x0 = xsol
      fa = fx
    end if
  else
    x0 = c
    x1 = xsol
    fa = fc
    fb = fx
  end if
end do
iter = i
end subroutine ridder_res
end subroutine shoot_integrate
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
!-----intermediate result print
!  write(*,'(<m+1>f12.5)') x, y
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
func(1) = y(2)
func(2) = -3._dp*y(1)*y(2)
!func(1) = y(2)
!func(2) = -4.75_dp*y(1) - 10.0_dp*y(2)
end function func
end module shoot_method
!*******************************************************************************
program main
use kinds
use shoot_method
implicit none 
integer :: m 
real(dp) :: xstart, xstop, u1, u2, x, yb, res
real(dp) :: ystart(2), ystop(2)
xstart = 0.0_dp
xstop = 2.00_dp
ystart = [0._dp, 1.5145_dp]   ! with root
yb = 1._dp
u1 = 1._dp
u2 = 2._dp 
call shoot_integrate(xstart,xstop,ystop,yb,u1,u2,2)
!call rk45_adaptive(xstart,xstop,ystart,ystop,2)
write(*,101) xstop, ystop
101 format(T4, 'Result from shooting method, ',"x,y,y': "/,<m+1>(f13.6))
end program main
