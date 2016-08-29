!*******************************************************************************
module kinds
implicit none
integer, parameter :: dp=kind(0.d0)  !double precision 
end module kinds
!*******************************************************************************
module curvfit
use kinds
private
public :: uptri, gramschmidt, leastsquares, lsqcurvefit
contains
!*******************************************************************************
subroutine lsqcurvefit(x,y,n,c,m)
!  Polynomial curve fitting with least squares
!  x,y ---- input data
!  n ---- dimension of x and y
!  m ---- order of polynomials
!  c ---- coefficients of polynomials
!  f(x) = c1 + c2*x + c3*x^2 + c4*x^3 + .... 
implicit none
integer, intent(in) :: m, n
integer :: i
real(dp), intent(in) :: x(n), y(n)
real(dp), intent(out) ::  c(m)
real(dp) :: a(n,m)
!-----error 
if(m>n) stop 'Warning: order of polynomial is greater than the number of data points'
do i=1,n
  a(i,:)=bv(x(i),m)
end do
call leastsquares(a,y,c,n,m)
contains
!-------------------------------------------------------------------------------
function bv(x,m)
!  base function
implicit none
integer, intent(in) :: m
real(dp), intent(in) :: x
integer :: i
real(dp) :: bv(m)
bv(1)=1.0_dp
do i=2,m
  bv(i)=bv(i-1)*x
end do  
end function bv
!-------------------------------------------------------------------------------
end subroutine lsqcurvefit
!*******************************************************************************
subroutine leastsquares(a,b,x,m,n)
!  Least squares solution for Ax=b where A(m,n)
!  Ax=b, QRx=b, Rx=QTb
!  Back substitution 
implicit none
integer, intent(in) :: m, n
real(dp), dimension(m,n), intent(in) :: a
real(dp), dimension(m), intent(in) :: b
real(dp), dimension(n), intent(out) :: x
real(dp), dimension(m,n) :: q
real(dp), dimension(n,m) :: qt
real(dp), dimension(n,n) :: r
real(dp), dimension(n) ::  qtb
call gramschmidt(a,q,r,m,n)
qt=transpose(q)
!-----Rx=Q'b
qtb=matmul(qt,b)
!-----back substitution
call uptri(r,qtb,x,n)
end subroutine leastsquares
!*******************************************************************************
subroutine gramschmidt(a,q,r,m,n)
!  Corrected Gram-Schmidt method for QR decomposition
!  A=QR  is a rectangular matrix A(m,n)
!  Q: orthogonal matrix, R: upper triagular matrix
!  Q(m,n), R(n,n) 
implicit none
integer, intent(in) :: m, n
integer :: i, j, k
real(dp), dimension(m,n), intent(in) :: a
real(dp), dimension(m,n), intent(out) :: q
real(dp), dimension(n,n), intent(out) :: r
real(dp), dimension(m) :: temp
r(1,1)=sqrt(dot_product(a(:,1),a(:,1)))
q(:,1)=a(:,1)/r(1,1)
do k=2,n
  do j=1,k-1
    r(j,k)=dot_product(q(:,j),a(:,k))
  end do
  temp=a(:,k)
  do j=1,k-1
    temp=temp-q(:,j)*r(j,k)
  end do
  r(k,k)=sqrt(dot_product(temp,temp))
  q(:,k)=temp/r(k,k)
end do
end subroutine gramschmidt
!*******************************************************************************
subroutine uptri(a,b,x,n)
implicit none 
integer :: i,k
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp) :: s
x(n)=b(n)/a(n,n)
do i=n-1,1,-1
  s=0.0_dp
  do k=i+1,n
    s=s+a(i,k)*x(k)
  end do
  x(i)=(b(i)-s)/a(i,i)
end do 
end subroutine uptri
!*******************************************************************************
end module curvfit
!*******************************************************************************
program main
use kinds
use curvfit
implicit none
integer, parameter :: n=7
real(dp), dimension(n) :: x, y
real(dp) :: c1(3), c2(4), c3(9)
x(:)=([-3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp])
y(:)=([4.0_dp, 2.0_dp, 3.0_dp, -0.0_dp, -1.0_dp, -2.0_dp, -5.0_dp])
call lsqcurvefit(x,y,n,c1,3)
call lsqcurvefit(x,y,n,c2,4)
write(*,101) c1, c2
101 format(/T4,'Polynomial curve fitting',//, &
            T3,'2nd order polynomial with coefficients: ', 3(/,F16.11),//, &
            T3,'3rd order polynomial with coefficients: ', 4(/,F16.11),/)
write(*,*)
call lsqcurvefit(x,y,n,c2,9)
end program main
!*******************************************************************************
