!*******************************************************************************
module kinds
implicit none
integer, parameter :: dp=kind(0.d0)  !  double precision
end module kinds
!*******************************************************************************
module nonlinear
use kinds
implicit none
private 
public :: newton 
contains
!*******************************************************************************
subroutine newton()
!  Newton-Raphson iterative methods for systems of equations
!  Solve J(x)dx=-f(x) by Gauss elimination
!  Jacobian matrix is evaluated analytically
!  or by finite difference approximation
integer :: i, imax=50
integer, parameter :: n=2
real(dp) :: x(n), f(n), dx(n), df(n,n)
real(dp) :: tol=1.0e-8_dp
!-----initial guess
x=([2.0_dp, 2.0_dp])
write(*,*) "Nonlinear equations by Newton's method"
do i=1,imax
  f=func(x)
  df=jacob(x)
  call gauss(df,-f,dx,n)
  x=x+dx
  write(*,*) i,x
  if(sqrt(dx(1)**2+dx(2)**2)<tol) exit
end do
end subroutine newton
!******************************************************
function func(x)
!  input function
implicit none
integer, parameter :: n=2
real(dp), intent(in) :: x(n)
real(dp) :: func(n)
func(1) = 6*x(1)**3 + x(1)*x(2) - 3*x(2)**3 -4
func(2) = x(1)**2 - 18*x(1)*x(2)**2 + 16*x(2)**3 + 1
end function func
!******************************************************
function jac(x)
!  Analytically obtained first-derivative Jacobian matrix
implicit none 
integer, parameter :: n=2
real(dp), intent(in) :: x(n) 
real(dp) :: jac(n,n)
jac(1,1) = 18*x(1)**2 + x(2) 
jac(1,2) = x(1) - 9*x(2)**2
jac(2,1) = 2*x(1) - 18*x(2)**2 
jac(2,2) = -36*x(1)*x(2) + 48*x(2)**2 
end function jac
!******************************************************
function jacob(x)
!  Jacobian matrix by finite difference approximation
implicit none 
integer :: i
integer, parameter :: n=2
real(dp) :: x(n)
real(dp) :: f1(n), f2(n), h, temp
real(dp) :: jacob(n,n)
h=1.0e-4_dp
f1=func(x)
do i=1,n
  temp=x(i)
!-----x+h  
  x(i)=temp+h
  f2=func(x)
!-----change x(i) back to obtain the derivative
!     with respect to next variable
  x(i)=temp
  jacob(:,i)=(f2-f1)/h
end do
end function jacob
!*******************************************************************************
subroutine gauss(a,b,x,n)
!  Gauss elimination
implicit none 
integer :: i,k,idmax
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n,1:n+1) :: ab
real(dp), dimension(1:n,1:n) :: aup
real(dp), dimension(1:n) :: bup
real(dp), dimension(1:n+1) :: temp1, temp2
real(dp) :: lam, elmax
ab(1:n,1:n)=a
ab(:,n+1)=b
do k=1,n-1
  elmax=abs(ab(k,k))
  idmax=k
! maximum element in the column
  do i=k+1,n
    if(abs(ab(i,k))>elmax) then
      elmax=ab(i,k)
      idmax=i
    end if
  end do
!  swap two lines
  temp1=ab(k,:)
  temp2=ab(idmax,:)
  ab(k,:)=temp2
  ab(idmax,:)=temp1
  do i=k+1,n  
    lam=ab(i,k)/ab(k,k)
    ab(i,:)=ab(i,:)-lam*ab(k,:)
  end do
end do
aup(:,:)=ab(1:n,1:n)
bup(:)=ab(:,n+1)
call uptri(aup,bup,x,n)
end subroutine gauss
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
end module nonlinear
!*******************************************************************************
program main
use kinds
use nonlinear
implicit none 
call newton()
end program main
