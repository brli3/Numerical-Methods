!*******************************************************************************
module kinds 
implicit none
integer, parameter :: dp=kind(0.d0)  ! double precision
end module kinds
!*******************************************************************************
module diff
!  Numerical differentiation by central difference method and five
!  point polynomial fitting. 
!  Richardson extrapolation eliminates truncation errors.
!  Analytical differentiation is often not available 
!  for complex functions and discrete points.
!  e.g. compute velocity and acceleration from discrete distance
!  data points.
use kinds
implicit none 
private 
public :: central_diff, five_diff, richardson
contains
!*******************************************************************************
subroutine central_diff()
!  Central difference method for fist derivative of a function
!  A simple demenstration to show that accuracy decreases with 
!  sufficiently small h.
!  due to subtraction when f(x+h) is too close to f(x-h).
implicit none
integer :: i
real(dp) :: x, h, f1, f2, df 
write(*,*) 'Result by mid-point formula:'
h=1.0_dp
x=sqrt(2._dp)
do i=1,50
  f2=atan(x+h)
  f1=atan(x-h)
  df=(f2-f1)/(2*h)
  h=h/2.0_dp
  write(*,*) i,h,df
end do
write(*,*)
end subroutine central_diff
!*******************************************************************************
subroutine five_diff()
!  Five point polynomial fitting on a uniform grid.
!  Compute f'(0.5) of function f(x)=x^2e^x.
!  f'(i) = (-f(i+2)+8f(i+1)-8f(i-1)+f(i-2))/12h + O(h^4)
!  Details of derivation of eqs are not shown here.
implicit none 
integer :: i, m  ! Midpoint
integer, parameter :: n = 5
real(dp) :: x(n), f(n), h, df
h = 1.0e-2_dp
m = int(n/2) + 1
x(m) = 0.5_dp
!-----set grid around the midpoint
do i=m,n
  x(i) = x(m) + (i-m)*h 
end do
do i=m,1,-1
  x(i) = x(m) - (m-i)*h
end do 
!-----function values
do i=1,n
!  f(i) = sin(x(i)) + cos(x(i)) 
  f(i) = x(i)**2*exp(-x(i))
end do
!-----five point formula
df = -f(m+2)+8._dp*f(m+1)-8._dp*f(m-1)+f(m-2)
df = df/(12._dp*h)
write(*,*) 'Five point poly: '
write(*,'(t8,a2,24x,a2)') 'x', 'y'
write(*,*) x(m), df
write(*,*)
end subroutine five_diff
!*******************************************************************************
subroutine richardson()
!  Richardson extrapolation
!  Eliminate truncation error
implicit none
real(dp) :: x
real(dp) :: h
real(dp) :: f, f_exact
integer :: p
x=0.5
h=0.05
p=2

f=(2.0_dp**p*g(x,h/2.0_dp)-g(x,h))/(2.0_dp**p-1.0_dp)
f_exact=g_exact(x)
write(*,*)
write(*,*) 'Richardson:'
write(*,'(t8,a2,24x,a2)') 'x', 'y'
write(*,*) x, f
write(*,*) 'Analytical result:'
write(*,*) f_exact

contains 
!************************************
function func(x)
implicit none
real(dp), intent(in) :: x
real(dp) :: func
func = x**2*exp(-x)
end function func
!************************************
function g(x,h)
!  Central difference first derivative
implicit none 
real(dp), intent(in) :: x
real(dp), intent(in) :: h
real(dp) :: g
!g = (func(x+h)+func(x-h)-2*func(x))/h**2 
g = (func(x+h)-func(x-h))/(2.0_dp*h)
end function g
!***********************************
function g_exact(x)
implicit none 
real(dp), intent(in) :: x
real(dp) :: g_exact
g_exact = 2*x*exp(-x)-x**2*exp(-x)
end function g_exact
!************************************
end subroutine richardson
!*******************************************************************************
end module diff
!*******************************************************************************
program main
use kinds
use diff
implicit none 
call central_diff()
call five_diff()
call richardson()
end program main
!*******************************************************************************
