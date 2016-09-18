program main
!----------------------------------------------------------------------
!   Example of an initial value problem;
!	solves 1d transient heat conduction equation
!	dT/dt=alpha*d2T/dx2
!	explicit method
!	iterative mehtod(explicit)
!----------------------------------------------------------------------
implicit none 
integer :: i,ni,nt,i_read,file1,file2
real :: t,alpha,tmax,xmax,delx,delt,s
real, dimension(:), allocatable :: x,temp,temp_cp
file1=1
file2=2
open(unit=file1,file='diff.in')
open(unit=file2,file='diff.dat')
do i_read=1,2
  read(file1,*)
end do
read(file1,*) ni
read(file1,*) delt
read(file1,*) tmax
read(file1,*) xmax
read(file1,*) alpha
!-----initiate
allocate(x(ni),temp(ni),temp_cp(ni))
!-----grid
delx=xmax/ni
x(1)=0.
x(ni)=xmax
do i=2,ni-1
  x(i)=x(i-1)+delx
end do
!-----initial value
nt=int(tmax/delt)
s=alpha*delt/(delx**2)
do i=1,ni
  temp(i)=0.0
end do
t=0.
!-----boundary condition
temp(1)=100.
temp(ni)=100.
!-----time integration
do while (t<=tmax)
  do i=1,ni
    temp_cp(i)=temp(i)
  end do
  do i=2,ni-1
    temp(i)=s*temp_cp(i-1)+(1-2*s)*temp_cp(i)+s*temp_cp(i+1)
  end do
  write(file2,102) t
  write(file2,103) (temp(i),i=1,ni)
  t=t+delt
end do
deallocate(x,temp,temp_cp)
102 format('time:',1X,7F6.1)
103 format(10F8.2)
end program main
