module filter1d

implicit none

contains

subroutine gaussian_filter_1d0(x_in, y_in, x, ynew, delta)
! Filtering 1-D progress variable from laminar flame calculation 
! using a 1-D Gaussian filter.
! Spatial profile of 1-D progress variable needs to be 
! provided in a file for the corresponding fuel type.
! File name should be formated like 'progvar1d_C3H8_phi_1_0.txt'.

! Ruipengyu Li 20/03/2017
! Adapted from Bo
integer, parameter :: dp = selected_real_kind(15)
integer :: nl
integer :: i, j
integer :: iex
integer :: ierr
real(dp), intent(in) :: delta
real(dp), dimension(:), intent(in) :: x_in, y_in
real(dp), dimension(:), intent(out) :: x, ynew
real(dp), allocatable, dimension(:) :: dx, y, x_ex
real(dp) :: pi, gw1, gw2, dx_min, dx_max, gauss, gauss_ex

nl = size(x_in)
allocate(dx(nl-1), y(nl-1), stat=ierr)
x = 0.0_dp
dx = 0.0_dp
y = 0.0_dp
ynew = 0.0_dp
do i=1,nl-1
  x(i) = 0.5_dp * (x_in(i) + x_in(i+1))
  y(i) = 0.5_dp * (y_in(i) + y_in(i+1))
  dx(i) = x_in(i+1) - x_in(i)
end do
x(nl) = x_in(nl)
dx_min = minval(dx)
dx_max = maxval(dx)
iex = int((delta+0.01*dx_min) / dx_min)
allocate(x_ex(-iex:nl+iex), stat=ierr)
pi = 4.0_dp*atan(1.0_dp)
do i=1,nl-1
  x_ex(i) = x(i)
end do
do i=0,-iex,-1
  x_ex(i) = x_ex(i+1) - dx(1)
end do
do i=nl,nl+iex
  x_ex(i) = x_ex(i-1) + dx(nl-1)
end do
do i=1,nl-1
  gauss = 0.0_dp
  gauss_ex = 0.0_dp
  do j=i-iex,i+iex-1
    if(j>=1 .and. j<(nl-1)) then
      gw1 = sqrt(6.0/pi/delta**2) * exp(-6.0*(x(j)-x(i))**2/delta**2)
      gw2 = sqrt(6.0/pi/delta**2) * exp(-6.0*(x(j+1)-x(i))**2/delta**2)
      ynew(i) = ynew(i) + 0.5*(gw1*y(j)+gw2*y(j+1)) * (x(j+1)-x(j))
      gauss = gauss + 0.5*(gw1+gw2) * (x(j+1)-x(j))
    else
      gw1 = sqrt(6.0/pi/delta**2) * exp(-6.0*(x_ex(j)-x_ex(i))**2/delta**2)
      gw2 = sqrt(6.0/pi/delta**2) * exp(-6.0*(x_ex(j+1)-x_ex(i))**2/delta**2)
      gauss_ex = gauss_ex + 0.5*(gw1+gw2) * (x_ex(j+1)-x_ex(j))
    end if
  end do
  ynew(i) = ynew(i) * (gauss+gauss_ex)/gauss
  ynew(nl) = ynew(nl-1)
end do
deallocate(dx, y, x_ex, stat=ierr)
end subroutine gaussian_filter_1d0

end module filter1d

program main

use filter1d

implicit none

integer, parameter :: dp = selected_real_kind(15)
integer :: i, n
integer :: ierr
real(dp), allocatable, dimension(:) :: x, y
real(dp), allocatable, dimension(:) :: xnew, ynew
real(dp) :: dx, delta
n = 100
allocate(x(n), y(n), stat=ierr)
allocate(xnew(n), ynew(n), stat=ierr)
x = 0.0_dp
dx = 1.0_dp/n
do i=2,n
  x(i) = x(1) + (i-1)*dx
end do
y = 0.0_dp
where(x > 0.5_dp) y = 1.0_dp
delta = 1.0e-1_dp
call gaussian_filter_1d0(x, y, xnew, ynew, delta)
write(*,'(2f12.5)') (xnew(i), ynew(i), i=1,n) 
end program main
