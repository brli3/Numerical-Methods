!*******************************************************************************
module kinds 
implicit none 
integer, parameter :: dp=kind(0.d0)  !  double precision
end module kinds
!*******************************************************************************
module intgrl
!  Numerical integration methods including Newton-Cotes (trapezoid and Simpson)
!  and Gauss quadrature (Gauss-Legendre). 
!  Multiple integrals of simple geometrical shapes can be evaluated by Gauss 
!  quadrature in a similar way and is not illustrated here.
use kinds
implicit none
private 
public :: func, intf_exact, trapezoid, auto_simpson, romberg, gauss_legendre, &
          complex_gauss, auto_gauss
contains
!*******************************************************************************
subroutine trapezoid(intf,a,b,n)
!  Trapezoidal rule
!  intf----definite integral of a function
!  a----lower bound
!  b----upper bound
!  n----number of panels
implicit none
integer :: i
integer, intent(in) :: n
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
real(dp) :: h, s

h=(b-a)/n

s=0.0_dp
do i=1,n-1
  s=s+func(a+i*h)
end do
intf = (func(a) + func(b) + 2.0_dp*s)*h/2.0_dp
end subroutine trapezoid
!*******************************************************************************
subroutine auto_simpson(intf,a,b,n)
!  adaptive simpson rule
!  automatic adjustment of panel size
!  n-----initial number of panels
implicit none
integer :: i, imax=20
integer, intent(inout) :: n
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
real(dp) :: intf1, intf2, tol=1.0e-7_dp
do i=1,imax
  call simpson(intf1,a,b,n)
  n=2.0_dp*n
  call simpson(intf2,a,b,n)
  if(abs(intf1-intf2)<tol) exit
end do
intf=intf2

contains
subroutine simpson(intf,a,b,n)
!  integration with Simpson rule
!  higher accuracy than trapezoid
implicit none
integer :: k
integer, intent(in) :: n
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
real(dp) :: h, s1, s2
h=(b-a)/(2.0_dp*n)
!  k=0
s1=func(a+h)
s2=0.0_dp
do k=1,n-1
  s1=s1+func(a+(2.0_dp*k+1.0_dp)*h)
  s2=s2+func(a+2.0_dp*k*h)
end do
intf = (func(a) + func(b) + 4.0_dp*s1 + 2.0_dp*s2)*h/3.0_dp
end subroutine simpson
end subroutine auto_simpson
!*******************************************************************************
subroutine romberg(intf,a,b,tol)
implicit none
integer :: i, j, k, m, imax=40
real(dp), intent(in) :: a, b
real(dp), intent(in) :: tol
real(dp), intent(out) :: intf
real(dp) :: t(1:50,0:49)
real(dp) :: s, temp
t(1,0)=(b-a)/2.0_dp*(func(a)+func(b))
do i=1,imax
  s=0.0_dp
  do j=1,2**(i-1)
    s=s+func(a+(2.0_dp*j-1.0_dp)*(b-a)/2.0_dp**i)
  end do
  t(1,i)=1.0_dp/2*t(1,i-1)+(b-a)/2**i*s
!-----extrapolation
  do m=1,i
    do k=i,1,-1
      temp=4.0_dp**m*t(m,k)-t(m,k-1)
      t(m+1,k-1)=temp/(4.0_dp**m-1.0_dp)
    end do
  end do
  if(abs(t(m,0)-t(m-1,0))<tol) exit
end do
!-----last diagonal term
intf=t(m,0)
end subroutine romberg
!*******************************************************************************
subroutine auto_gauss(intf,a,b,m,n)
!  Recursive Gauss integration 
implicit none 
integer, intent(in) :: m   !  number of maximum redividing times
integer, intent(inout) :: n ! number of panels
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
real(dp), parameter :: tol = 1.0e-8_dp
real(dp) :: s1, s2
integer :: i
n = 2
do i=1,m
  call complex_gauss(s1,a,b,n)
!  Redivision
  n = n*2
  call complex_gauss(s2,a,b,n)
  if(abs(s2-s1)<tol) exit
end do
intf = s2
end subroutine auto_gauss
!*******************************************************************************
subroutine complex_gauss(intf, a, b, n)
!  Complex Gass-Legendre integration
implicit none
integer, intent(in) :: n   !  number of panels
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
real(dp) :: c, d, h, s
integer :: i
h = (b-a)/n
intf = 0.0_dp
do i=1,n
  c = a + (i-1)*h
  d = a + i*h
  call gauss_legendre(s,c,d)
  intf = intf + s
end do
end subroutine complex_gauss
!*******************************************************************************
subroutine gauss_legendre(intf,a,b)
!  Gauss-Legendre integration
implicit none 
real(dp), intent(in) :: a, b
real(dp), intent(out) :: intf
integer, parameter :: n = 5  ! order of accuracy
real(dp), dimension(n) :: t, x, w  
integer :: i
!  standard coeffients from tabel
!  can also be computed 
t = ([-0.9061798459_dp, -0.5384693101_dp, 0.0_dp, 0.538493101_dp, 0.9061798459_dp])
w = ([0.2369268851_dp, 0.4786286705_dp, 0.5688888889_dp, 0.478628705_dp, 0.2369268851_dp])
!  integral bounds transform from [-1,1] to [a.b]
x = (b+a)/2.0_dp + (b-a)*t/2.0_dp
intf = 0.0_dp
do i=1,n
  intf = intf + w(i) * func(x(i))
end do 
intf = intf * (b-a)/2.0_dp
end subroutine gauss_legendre
!*******************************************************************************

!*******************************************************************************
function func(x)
!  function to be evaluated
implicit none
real(dp), intent(in) :: x
real(dp) :: func
func = x**2 + sin(x)
end function func
!************************
function intf_exact(a,b)
!  exact integral value
implicit none
real(dp), intent(in) :: a, b
real(dp) :: intf_exact
intf_exact = (1.0_dp/3*b**3.0_dp - cos(b)) - (1.0_dp/3*a**3.0_dp - cos(a))
end function intf_exact
!*******************************************************************************
end module intgrl

!*******************************************************************************
program main
use kinds
use intgrl
implicit none
integer :: n, m, p
real(dp) :: a, b, intf
n=40
a=-2.0_dp
b=2.0_dp
write(*,101) n
101 format(T4,'Number of panels:',I3)
call trapezoid(intf,a,b,n)
write(*,102) intf
102 format(/,T4,'Result by Trapezoidal rule:',F15.8)
call auto_simpson(intf,a,b,n)
write(*,103) intf, n
103 format(/,T4,'Result by Simpson rule:',F15.8, &
           /,T4,'With adjusted number of panels: ',I4)
call romberg(intf,a,b,1.0e-7_dp)
write(*,104) intf
104 format(/,T4,'Result by Romberg rule:',F15.8)

m=20  !  maximum redivison times 
call auto_gauss(intf,a,b,m,p)
write(*,105) intf, p
105 format(/,T4,'Result by Gauss-Legendre:',F15.8, &
           /,T4,'With actual number of panels:', I4)
          
write(*,106) intf_exact(a,b)
106 format(/,T4,'Exact result:',F15.8)

end program main
!*******************************************************************************
