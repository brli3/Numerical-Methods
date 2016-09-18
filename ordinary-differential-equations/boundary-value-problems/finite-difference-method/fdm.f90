!*******************************************************************************
module kinds
implicit none
integer, parameter :: dp=kind(0.d0)  !double precision 
end module kinds
!*******************************************************************************
module finite_diff
!  Finite difference method for B.V.P
!  For linear and nonlinear ODEs
use kinds
implicit none
private
public :: fd_lin, fd_nonlin

contains
!*******************************************************************************
subroutine fd_lin()
!  Sample finite difference method routine to solve a 2nd order linear ODE
!  y'' = -4y + 4x; y(0)=0, y'(pi/2)=0
implicit none
integer :: i
integer, parameter :: n = 11  ! Number of grid points
real(dp), parameter :: pi = 3.14159265359_dp
real(dp) :: xstart, xstop 
real(dp) :: ystop
real(dp), dimension(1:n,1:n) :: a  ! Coefficient matrix
real(dp) :: b(1:n), x(1:n), y(1:n), h, h2
xstop = pi/2
xstart = 0.0_dp
h = (xstop - xstart)/(n-1)  ! Grid space 
h2 = h*h
!-----Initialization
a = 0.0_dp
b = 0.0_dp
x(1) = xstart
b(1) = 0._dp
do i=2,n
  x(i) = x(i-1) + h
  b(i) = 4._dp*h2*x(i)
end do
!-----Compute coefficient matrix 
do i=1,n
  a(i,i) = -2._dp + 4._dp*h2
end do
do i=2,n
  a(i,i-1) = 1.0_dp
end do
do i=1,n-1
  a(i,i+1) = 1.0_dp
end do
!-----Boundary values
a(1,1) = 1.0_dp
a(1,2) = 0.0_dp
a(n,n-1) = 2.0_dp

write(*,'(/,t6,a3,9x,a3,/)') 'x','y'
!-----Solve the tri-diagonal matrix with Thomas algorithm
call gauss(a,b,y,n)
do i=1,n
  write(*, '(2(f12.6))') x(i), y(i)
end do
end subroutine fd_lin
!*******************************************************************************
subroutine fd_nonlin()
!  Solve a non-linear ODE
!  y'' = -3yy'; y(0)=0, y(2)=1
!  Newton-Raphson method is used for the 
!  system of non-linear eqs.
implicit none
integer :: i
integer, parameter :: n = 11  ! Number of grid points
real(dp) :: xstart, xstop, h, x(1:n), y(1:n)
xstart = 0.0_dp
xstop = 2.0_dp
x(1) = xstart
h = (xstop - xstart)/(n-1)  ! Grid space 
do i=2,n
  x(i) = x(i-1) + h
end do
call newton(y)
write(*,'(/,t6,a3,9x,a3,/)') 'x','y'
do i=1,n
  write(*, '(2(f12.6))') x(i), y(i)
end do
end subroutine fd_nonlin
!*******************************************************************************
function func_diff(x,y,yprime)
!  Differential equation 
!  y'' = f(x,y,y') = -3yy'
implicit none
real(dp), intent(in) :: x, y, yprime
real(dp) :: func_diff
func_diff = -3.0_dp*y*yprime
end function func_diff
!*******************************************************************************
function func(y)
!  input function
!  residual function 
implicit none
integer :: i
integer, parameter :: n=11
real(dp), intent(in) :: y(n)
real(dp) :: x(1:n)
real(dp) :: xstart, xstop, h
real(dp) :: func(n)
xstart = 0.0_dp
xstop = 2.0_dp
h = (xstop - xstart)/(n-1)  ! Grid space 

func = 0.0_dp
func(1) = y(1) 
func(n) = y(n) - 1.0 

do i=2,n-1
  func(i) = y(i-1) -2._dp*y(i) + y(i+1) - h*h* &
            func_diff(x(i),y(i),(y(i+1)-y(i-1))/(2._dp*h))
end do
end function func
!******************************************************
subroutine newton(x)
!  Newton-Raphson iterative methods for systems of equations
!  Solve J(x)dx=-f(x) by Gauss elimination
!  Jacobian matrix is evaluated analytically
!  or by finite difference approximation
integer :: i, j
integer, parameter :: imax=50
integer, parameter :: n=11
real(dp), intent(out) :: x(n) 
real(dp) :: f(n), dx(n), df(n,n)
real(dp) :: tol=1.0e-8_dp
!-----initial guess
x = 1.0_dp
write(*,*) "Intermediate results from Newton's method:"
do j=1,n
  write(*,'(f12.6)') x(j)
end do
write(*,*)
do i=1,imax
  f=func(x)
  df=jacob(x)
  call gauss(df,-f,dx,n)
  x=x+dx
!  write(*,'(a6,i2,/,<n>(f12.6),/)') 'Iter: ', i,x
  write(*,'(a6,i3)') 'Iter', i
  do j=1,n
    write(*,'(f12.6)') x(j)
  end do
  write(*,*)
  if(sqrt(dx(1)**2+dx(2)**2)<tol) exit
end do
end subroutine newton
!******************************************************
function jacob(x)
!  Jacobian matrix by finite difference approximation
implicit none
integer :: i
integer, parameter :: n=11
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
subroutine thomas(a,f,x,n)
!  Thomas algorithm
!  only for tri-diagonal matrices Ax=f
!  may include a matrix check procedure
implicit none 
integer :: i
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: f
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n-1) :: c,d
real(dp), dimension(1:n) :: b,u,y
real(dp), dimension(2:n) :: e,l
!  copy A to vectors
do i=1,n
  b(i)=a(i,i)
end do 
do i=1,n-1
  c(i)=a(i,i+1)
end do
do i=2,n
  e(i)=a(i,i-1)
end do
!  decompose 
do i=1,n-1
  d(i)=c(i)
end do
u(1)=b(1)
do i=2,n
  l(i)=e(i)/u(i-1)
  u(i)=b(i)-l(i)*c(i-1)
end do
!  back substitution
y(1)=f(1)
do i=2,n
  y(i)=f(i)-l(i)*y(i-1)
end do
!  back substitution
x(n)=y(n)/u(n)
do i=n-1,1,-1
  x(i)=(y(i)-c(i)*x(i+1))/u(i)
end do
end subroutine thomas
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
end module finite_diff
!******************************************************************************
program main
use kinds
use finite_diff
implicit none
write(*,*) 'Solution for linear example: '
call fd_lin()
write(*,*)
write(*,*) 'Solution for non-linear example: '
call fd_nonlin()
end program main
