!*******************************************************************************
module kinds 
implicit none
integer, parameter :: dp=kind(0.d0)  ! double precision
end module kinds
!*******************************************************************************
module diff
use kinds
implicit none 
private 
public :: central_diff, richardson
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
end subroutine central_diff
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
write(*,*) 'Richardson:'
write(*,*) f
write(*,*)
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
call central_diff
call richardson
end program main
!*******************************************************************************
