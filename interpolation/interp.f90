!*******************************************************************************
module kinds
implicit none
integer, parameter :: dp=kind(0.d0)  !double precision 
end module kinds
!*******************************************************************************
module interpolation
use kinds
implicit none
private
public :: lagrange, newton
contains
!*******************************************************************************
subroutine lagrange(n,xdata,ydata,m,x,p)
!  Lagrange's interpolation method
!  xdata,ydata ----- data points (x0,y0), (x1,y1),...,(xn,yn)
!  n ----- number of data points minus 1
!  m ----- number of interpolated points
!  x ----- interpolating points
!  p ----- interpolation result
implicit none
integer :: i
integer, intent(in) :: n, m
real(dp), intent(in) :: xdata(0:n), ydata(0:n), x(0:m)
real(dp), dimension(n), intent(out) :: p(0:m)
do i=0,m
  call lagrange_interp(n,xdata,ydata,x(i),p(i))
end do
contains 
!******************************************************
subroutine lagrange_interp(n,xdata,ydata,x,p)
!  interpolation of a single point
implicit none
integer :: i, j 
integer, intent(in) :: n 
real(dp), intent(in) :: xdata(0:n), ydata(0:n)
real(dp), intent(in) :: x
real(dp), intent(out) :: p
real(dp) :: l(0:n)
!-----weighting function
l(:)=1.0_dp
do i=0,n
  do j=0,n
    if(j==i) cycle
!-----product term    
    l(i)=l(i)*(x-xdata(j))/(xdata(i)-xdata(j))
  end do
end do
!-----sum up for polynomial
p=0.0_dp
do i=0,n
  p=p+ydata(i)*l(i)
end do
end subroutine lagrange_interp
!******************************************************
end subroutine lagrange


!*******************************************************************************
subroutine newton(n,xdata,ydata,m,x,p)
implicit none
integer :: i, j, k
integer, intent(in) :: n, m
real(dp), intent(in) :: xdata(0:n), ydata(0:n), x(0:m)
real(dp), dimension(n), intent(out) :: p(0:m)
do i=0,m
  call newton_interp(n,xdata,ydata,x(i),p(i))
end do
contains
!******************************************************
subroutine newton_interp(n,xdata,ydata,x,y)
implicit none 
integer :: i, j, k
integer, intent(in) :: n
real(dp), intent(in) :: xdata(0:n), ydata(0:n)
real(dp), intent(in) :: x
real(dp), intent(out) :: y
real(dp) :: q(0:n,0:n), b(n+1)
do i=0,n
  q(i,0)=ydata(i)
end do
do i=1,n
  do j=1,i
    q(i,j)=(q(i,j-1)-q(i-1,j-1))/(xdata(i)-xdata(i-j))
  end do
end do
b(n+1)=q(n,n)
do k=n,1,-1
  b(k)=q(k-1,k-1)+b(k+1)*(x-xdata(k-1))
end do
y=b(1)
end subroutine newton_interp
!******************************************************
end subroutine newton
!*******************************************************************************
end module interpolation


!*******************************************************************************
program main
use kinds
use interpolation
implicit none
real(dp), parameter :: pi = 3.141592653589793_dp
real(dp) :: adata(0:3), xdata(0:3), ydata(0:3)
real(dp) :: xdeg(0:2), x(0:2), y(0:2), yreal(0:2)
!-----data points from degree to radi
adata = ([30._dp, 45._dp, 60._dp, 90._dp])
xdata = adata*pi/180._dp
ydata = ([sqrt(3._dp)/2, sqrt(2._dp)/2, 1._dp/2, 0._dp])
!-----points to be interpolated
xdeg = ([47._dp, 53._dp, 79._dp])
x = xdeg*pi/180._dp
write(*,*) '  Computed cos(x) with various interpolation methods'
write(*,101) xdeg
call lagrange(3,xdata,ydata,2,x,y)
write(*,102) y
call newton(3,xdata,ydata,2,x,y)
write(*,103) y
yreal=cos(x)
write(*,104) yreal
101 format(/,T3,'Angle(deg):',T16,3F12.5)
102 format(T3,'Langrange:',T16,3F12.7)
103 format(T3,'Newton:',T16,3F12.7)
104 format(T3,'Real:',T16,3F12.7)
end program main
!*******************************************************************************
