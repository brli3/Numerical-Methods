!*******************************************************************************
module kinds
integer, parameter :: dp=kind(0.d0)  !double precision 
end module kinds
!*******************************************************************************
module linsolv
use kinds
private  
public :: uptri, downtri, gauss, doolittle, lusolv, crout, det, cholesky,      &
          mat_eq, mat_inv, gauss_iter, thomas, jacobi, gs, richardson,          &
          conjgrad, tdma
contains
!*******************************************************************************
subroutine mat_inv(a,inv_a,n)
!  calculate the inverse of matrix A
!  use AX=I
implicit none
integer :: i
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n,1:n), intent(out) :: inv_a
real(dp), dimension(1:n,1:n) :: e
do i=1,n
  e(i,i)=1.0_dp
end do
call mat_eq(a,e,inv_a,n,n)
end subroutine mat_inv
!*******************************************************************************
subroutine gauss_iter(a,b,x,n)
!  an simple iterative approach to correct the solution
implicit none 
integer :: i
integer, parameter :: itmax=5
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2, db, dx
!  initial value
call gauss(a,b,x1,n)
!x1=1.0_dp !  test a bad x1 value
do i=1,itmax
  db=matmul(a,x1)-b
!  correction dx
  call gauss(a,db,dx,n)
!  corrected solution x2
  x2=x1-dx
!  update
  x1=x2
end do
x=x1
end subroutine gauss_iter
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
subroutine mat_eq(a,b,x,n,m)
!  Gauss elimination for multiple sets of equations
implicit none 
integer :: i, k, idmax
integer, intent(in) :: n, m
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n,1:m), intent(in) :: b
real(dp), dimension(1:n,1:m), intent(out) :: x
real(dp), dimension(1:n,1:n+m) :: ab
real(dp), dimension(1:n,1:n) :: aup
real(dp), dimension(1:n+m) :: temp1, temp2
real(dp), dimension(1:n) :: temp3, temp4
real(dp) :: lam, elmax
ab(1:n,1:n)=a
ab(:,n+1:n+m)=b
do k=1,n-1
  elmax=abs(ab(k,k))
  idmax=k
! maximum element in a column
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
do i=1,m
  temp3=ab(:,n+i)
  call uptri(aup,temp3,temp4,n)
  x(:,i)=temp4
end do
end subroutine mat_eq
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
subroutine lusolv(a,b,x,n)
!  LU decomposition
implicit none 
integer :: i,k
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n,1:n) :: l, u
real(dp), dimension(1:n) :: y
call doolittle(a,l,u,n)
call downtri(l,b,y,n)
call uptri(u,y,x,n)
end subroutine lusolv
!*******************************************************************************
subroutine crout(a,l,u,n)
!  Crout decomposition
implicit none 
integer :: i,j,k,m
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n,1:n), intent(out) :: l, u
real(dp) :: s
l(:,1)=a(:,1)
u(1,:)=a(1,:)/l(1,1)
do k=2,n
  do i=k,n
    s=0.0_dp
    do m=1,k-1
      s=s+l(i,m)*u(m,k)
    end do
    l(i,k)=a(i,k)-s
  end do
  do j=k+1,n
    s=0.0_dp
    do m=1,k-1
      s=s+l(k,m)*u(m,j)
    end do
    u(k,j)=(a(k,j)-s)/l(k,k)
  end do
  u(k,k)=1
end do
end subroutine crout
!*******************************************************************************
subroutine doolittle(a,l,u,n)
! Doolittle decomposition
implicit none 
integer :: i,j,k,m
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n,1:n), intent(out) :: l, u
real(dp) :: s
u(1,:)=a(1,:)
l(:,1)=a(:,1)/u(1,1)
do k=2,n
  l(k,k)=1.0_dp
  do j=k,n
    s=0.0_dp
    do m=1,k-1
      s=s+l(k,m)*u(m,j)
    end do
    u(k,j)=a(k,j)-s
  end do
  do i=k+1,n
    s=0.0_dp
    do m=1,k-1
      s=s+l(i,m)*u(m,k)
    end do
    l(i,k)=(a(i,k)-s)/u(k,k)
  end do
end do
end subroutine doolittle
!*******************************************************************************
subroutine cholesky(a,b,x,n)
!  Cholesky decomposition for symmetric matrices
implicit none
integer :: i,j,k
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n,1:n) :: l, lt
real(dp), dimension(1:n) :: y
real(dp) :: s
logical :: sym
!  check symetricity
do i=1,n
  do j=1,n
    sym=a(i,j)==a(j,i)
    if (sym==0) stop "Cholesky returns. Matrix is not symmetric."
  end do
end do
!  first column
if(a(1,1)<0) stop "Cholesky returns. Matrix is not positive definite."
l(1,1)=sqrt(a(1,1))
l(2:,1)=a(2:,1)/l(1,1)
!  diagonal term
do j=2,n
  s=0.0_dp
  do k=1,j-1
    s=s+l(j,k)**2
  end do
  if((a(j,j)-s)<0) stop "Cholesky returns. Matrix is not positive definite."
  l(j,j)=sqrt(a(j,j)-s)
! nondiagonal term
  do i=j+1,n
    s=0.0_dp
    do k=1,j-1
      s=s+l(i,k)*l(j,k)
    end do
  l(i,j)=(a(i,j)-s)/l(j,j)
  end do
end do
!  end of decomposition
call downtri(l,b,y,n)
do i=1,n
  do j=1,n
    lt(i,j)=l(j,i)
  end do
end do
call uptri(lt,y,x,n)
end subroutine cholesky
!*******************************************************************************
function det(a,n)
!  calculate determinant of matrix A
implicit none 
integer :: i
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp) :: det
real(dp), dimension(1:n,1:n) :: l, u
call crout(a,l,u,n)
det=1.0_dp
do i=1,n
  det=det*l(i,i)
end do
end function det
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
subroutine downtri(a,b,x,n)
implicit none 
integer :: i,k
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: b
real(dp), dimension(1:n), intent(out) :: x
real(dp) :: s
x(1)=b(1)/a(1,1)
do i=2,n
  s=0.0_dp
  do k=1,i-1
    s=s+a(i,k)*x(k)
  end do
  x(i)=(b(i)-s)/a(i,i)
end do 
end subroutine downtri
!******************************************************************************
subroutine jacobi(a,b,x,x0,n)
!  Jacobi iterative method
implicit none 
integer :: i,j,k
integer, parameter :: imax=500
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: x0, b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2
real(dp) :: s, dx2, tol=1.0e-7_dp
do i=1,n
  if(a(i,i)==0.0_dp) stop "Jacobi returns. Matrix contains zero diagonal element."
end do
x1=x0
do k=1,imax
  do i=1,n
    s=0.0_dp
    do j=1,n
      if(j==i) cycle
      s=s+a(i,j)*x1(j)
    end do
    x2(i)=(b(i)-s)/a(i,i)
  end do
!  accuracy 
  dx2=0.0_dp
  do i=1,n
    dx2=dx2+(x1(i)-x2(i))**2
  end do
  dx2=sqrt(dx2)
  if(dx2<tol) exit
  x1=x2
end do
x=x2
end subroutine jacobi
!*******************************************************************************
subroutine gs(a,b,x,x0,n)
!  Gauss-Seidel iterative method
!  with relaxation 
implicit none 
integer :: i,j,k
integer, parameter :: imax=500
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: x0, b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2
real(dp) :: s, dx2, tol=1.0e-7_dp
real(dp), parameter :: rf=0.5
!  initial values
x1=x0
x2=x1
do k=1,imax
  do i=1,n
    s=0
    do j=1,n
!  update using known values    
      if(j<i) then  
        s=s+a(i,j)*x2(j)
      else if(j>i) then
        s=s+a(i,j)*x1(j)
      end if
    end do
    x2(i)=(b(i)-s)*rf/a(i,i)+(1-rf)*x1(i)
  end do
  dx2=0.0_dp
  do i=1,n
    dx2=dx2+(x1(i)-x2(i))**2
  end do
  dx2=sqrt(dx2)
  if(dx2<tol) exit
  x1=x2
end do
x=x2
end subroutine gs
!*******************************************************************************
subroutine richardson(a,b,x,x0,n)
!  richardson iterative method
implicit none 
integer :: i,j,k
integer, parameter :: imax=500
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: x0, b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2
real(dp) :: s, dx2, tol=1.0e-7_dp
real(dp), parameter :: rf=0.1
!  initial value
x1=x0
do k=1,imax
  do i=1,n
    s=0
    do j=1,n
!  update
      s=s+a(i,j)*x1(j)
    end do
    x2(i)=x1(i)-rf*(s-b(i))
  end do
  dx2=0.0_dp
  do i=1,n
    dx2=dx2+(x1(i)-x2(i))**2
  end do
  dx2=sqrt(dx2)
  if(dx2<tol) exit
  x1=x2
end do
x=x2
end subroutine richardson
!*******************************************************************************
subroutine conjgrad(a,b,x,x0,n)
!  Conjugate gradient method
!  for symmetric matrices only
implicit none 
integer :: i,j,k
integer, parameter :: imax=200
integer, intent(in) :: n
real(dp), dimension(1:n,1:n), intent(in) :: a
real(dp), dimension(1:n), intent(in) :: x0, b
real(dp), dimension(1:n), intent(out) :: x
real(dp), dimension(1:n) :: x1, x2, r0, r1, p0, p1
real(dp) :: alpha, beta, temp1, temp2, temp3, drs
real(dp) :: tol=1.0e-7_dp
real(dp), parameter :: rf=0.6
logical :: sym
do i=1,n
  do j=1,n
    sym=a(i,j)==a(j,i)
    if (sym==0) stop "Conjugate gradient returns. Matrix is not symmetric."
  end do
end do
x1=x0
r0=b-ar(a,x1,n)
p0=r0
do k=1,imax
  temp1=dr(r0,n)
  temp2=rar(a,p0,n)
  alpha=temp1/temp2
  x2=x1+alpha*p0
  drs=sqrt(dr(r0,n))
  if(drs<tol) exit
  r1=r0-alpha*ar(a,p0,n)
  temp3=dr(r1,n)
  beta=temp3/temp1
  p1=r1+beta*p0
!  update all
  r0=r1
  p0=p1
  x1=x2
end do
x=x2
contains
!******************************
function ar(a,r,n)
!  calculate A*r
implicit none 
integer i,n
real(dp), dimension(n,n), intent(in) :: a
real(dp), dimension(n), intent(in) :: r
real(dp), dimension(n) :: temp, ar
temp=0.0_dp
do i=1,n
  do j=1,n
    temp(i)=temp(i)+a(i,j)*r(j)
  end do
end do
ar=temp
end function ar
!******************************
function dr(r,n)
!  calculate square of a vector length
implicit none 
integer, intent(in) :: n
integer :: i
real(dp), dimension(n), intent(in) :: r
real(dp) :: s, dr
s=0.0_dp
do i=1,n
  s=s+r(i)**2
end do
dr=s
end function dr
!******************************
function v1v2(v1,v2,n)
!  vector dot product
implicit none 
integer, intent(in) :: n
integer :: i
real(dp), dimension(n), intent(in) :: v1, v2
real(dp) :: v1v2
v1v2=0.0_dp
do i=1,n
  v1v2=v1v2+v1(i)*v2(i)
end do
end function v1v2
!******************************
function rar(a,r,n)
!  calculate r*Ar
implicit none 
integer, intent(in) :: n
integer :: i
real(dp), dimension(n,n), intent(in) :: a
real(dp), dimension(n), intent(in) :: r
real(dp), dimension(n) :: temp
real(dp) :: rar
temp=ar(a,r,n)
rar=v1v2(r,temp,n)
end function rar
end subroutine conjgrad
!*******************************************************************************
subroutine tdma(istart,jstart,ni,nj,it,jt,phi)
!  Tridiagonal matrix or Thomas algorithm 
!  Iterative and line by line
!  A very simple and general method for CFD solution
implicit none 
integer, intent(in) :: istart, jstart, ni, nj, it, jt
integer :: i, j, jj, nim1, njm1, jsrm1
real, dimension(it,jt), intent(out) :: phi
real, dimension(it,jt) :: ap, an, as, ae, aw, su
real, dimension(2*it) :: a, b, c, d
!-----commence w-e sweep
do i=istart,nim1
  c(jsrm1)=phi(i,jsrm1)
!-----commence s-n transverse
  do j=jstart,njm1
!-----assemble tdma coefficients
    a(j)=an(i,j)     
    b(j)=as(i,j)
    c(j)=aw(i,j)*phi(i-1,j)+ae(i,j)*phi(i+1,j)+su(i,j)
    d(j)=ap(i,j)
!-----calculate coefficients of recurrence formula
    a(j)=a(j)/(d(j)-b(j)*a(j-1))
    c(j)=(b(j)*c(j-1)+c(j))/(d(j)-b(j)*a(j-1))
  end do
!-----obtain new phi
  do jj=jstart,njm1
    j=nj+jsrm1-jj
    phi(i,j)=a(j)*phi(i,j+1)+c(j)
  end do
end do
end subroutine tdma
!*******************************************************************************
end module linsolv 
!*******************************************************************************
program main 
use kinds
use linsolv 
implicit none 
integer, parameter :: n=5
integer :: i, j
real(dp), allocatable :: a(:,:), b(:), x(:), x0(:), inv_a(:,:)
real(dp) :: determinant
allocate(a(n,n), b(n), x(n), x0(n), inv_a(n,n))
a(1,:) = [6.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
a(2,:) = [-1.0_dp, 7.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
a(3,:) = [0.0_dp, -2.0_dp, 8.0_dp, 2.0_dp, 0.0_dp]
a(4,:) = [0.0_dp, 0.0_dp, 3.0_dp, 7.0_dp, -2.0_dp]
a(5,:) = [0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 5.0_dp]
b(:) = [2.0_dp, -3.0_dp, 4.0_dp, -3.0_dp, 1.0_dp]
determinant=det(a,n)
write(*,*) '  Determinant: '
write(*,'(F16.8)') determinant 
write(*,*)
call mat_inv(a,inv_a,n)
write(*,100) ((inv_a(i,j), j=1,n), i=1,n)
100 format(/,'Inverse of A:',/,5(F12.8),/)
call gauss(a,b,x,n)
write(*,101) x
101 format(T5, 'Result by gauss elimination:',/,T4,'x=',5(/7XF12.8),/)
call lusolv(a,b,x,n)
write(*,102) x
102 format(T5, 'Result by LU decomposition:',/,T4,'x=',5(/7XF12.8),/)
call gauss_iter(a,b,x,n)
write(*,103) x
103 format(T5, 'Result by Gauss with iterative correction:',/,T4,'x=',5(/7XF16.10),/)
call thomas(a,b,x,n)
write(*,104) x
104 format(T5, 'Result by thomas algorithm:',/,T4,'x=',5(/7XF12.8),/)
!call cholesky(a,b,x,n)
!write(*,105) x
!105 format(T5, 'Result by Cholesky decomposition:',/,T4,'x=',5(/7XF12.8))
x0(:) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
call jacobi(a,b,x,x0,n)
write(*,106) x
106 format(T5, 'Result by Jacobi:',/,T4,'x=',5(/7xF12.8),/)
call gs(a,b,x,x0,n)
write(*,107) x
107 format(T5, 'Result by Gauss-Seidel:',/,T4,'x=',5(/7xF12.8),/)
call richardson(a,b,x,x0,n)
write(*,108) x
108 format(T5, 'Result by Richardson:',/,T4,'x=',5(/7xF12.8),/)
call conjgrad(a,b,x,x0,n)
write(*,109) x
109 format(T5, 'Result by Conjugate gradient:',/,T4,'x=',5(/7xF12.8),/)
deallocate(a,b,x,x0)
end program main
!*******************************************************************************
