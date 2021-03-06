
module constants
implicit none
real,parameter    :: gamma    = 5./3.,&
					 Courant  = .5  , &
					 length   = 1.  ,&
					 tmax     = 10. ,&
					 qo       = 5.  ,&  
					 delta    = 1.e-30 ,&
					 uL       = 1.e-30 ,&      !(uL=.5 for piston , uL=.0 for shock tube)
					 uR       = 1.e-30  	 !(uR=.0 for piston , uR=.0 for shock tube)
integer,parameter :: jmax 	  = 100  ,& ! number of cells
					 itterMax = 1000  ,& ! Maximum itterations
					 writeStep= 10    ,& ! Itter step to write to file
					 IC       = 2       ! Initial condiiton (1-pisition,2 Shock)
end module

program main 
use constants
implicit none 
real :: dt=0.00005 , time=.0  , dmi
real,dimension(0:jmax+1) :: rho,p,e,q,c,dx,dm,r
real,dimension(-1:jmax+1) :: u,x,uT
integer :: i=0 , j

call OpenFiles()
call InitDomain(x,dx,r)
call InitCond(rho,e,p,c,q,u,dm,dx,x)

do while (time.LT.tmax.AND.i.LT.itterMax)
	call Boundary(u) 			! Impose BC on u
	x(-1)  = x(-1) + dt*u(-1)	! Move left most ghost cell
	!q = (/(qo*(dm(j)/(x(j)-x(j-1)))*(MIN(0.0,u(j)-u(j-1)))**2,j=0,jmax+1,1)/) ! Artificial viscious pressure
	q = (/(qo*rho(j)*ABS(MIN(0.0,u(j)-u(j-1)))*(.1*sqrt(gamma*p(j)/rho(j))+ABS(u(j)-u(j-1))),j=0,jmax+1,1)/) ! Artificial viscious pressure
	do j=0,jmax  
		x(j)   = x(j) + u(j)*dt
		dmi    = .5*(rho(j+1)+rho(j))*(r(j+1)-r(j))    ! Mass in "momentum zone" of cell interface of j and j+1 cell
		u(j)  = u(j) + (dt/dmi)*(p(j)+q(j)-p(j+1)-q(j+1)) ! Update the velocity in a temporary array
		rho(j) = dm(j)/(x(j)-x(j-1))
		e(j)   = e(j) + (dt/dm(j))*(p(j)+q(j))*(u(j-1)-u(j))
		p(j)   = (gamma-1.)*rho(j)*e(j)

	end do 
	r   = .5*(x(0:jmax+1)+x(-1:jmax))! Find new cell centers
	dx  =  x(0:jmax+1)-x(-1:jmax)   ! Find new cell widths
	i = i + 1						! Update itteration counter
	time = time+ dt					! update timer
	c = SQRT(gamma*p/rho)
	call Output(rho,e,p,u,dm,time,x,r,i)
	call Timestep(u,dx,dt,c,i)
	if (MINVAL(dx).LT..0) then
		print*,'Negative volume'
		exit 
	end if 
end do 
print*,'Terminated at:',i
call CloseFiles()
end program

subroutine InitDomain(x,dx,r)
use Constants
implicit none
real,intent(inout),dimension(-1:jmax+1) :: x
real,intent(inout),dimension(0:jmax+1)  :: dx,r
integer :: i
x   = (/((i*length/(SIZE(x)-1)),i=0,SIZE(x)-1,1)/)  ! Make domain uniformly spaced
dx  =  x(0:jmax+1)-x(-1:jmax)   					! Spacing of cells
r   = .5*(x(0:jmax+1)+x(-1:jmax))         			! Cell centers
end subroutine 

subroutine Boundary(u)
use Constants
implicit none
real,intent(inout),dimension(-1:jmax+1) :: u
if (IC.EQ.1) then
	u(-1)     = uL ! Left cell must move with piston
	u(jmax+1) = uR ! Right cell must move with right edge
else if (IC.EQ.2) then
	u(-1) = u(0)
	u(jmax+1) = u(jmax)
end if 
end subroutine

subroutine InitCond(rho,e,p,c,q,u,dm,dx,x)
use Constants
implicit none 
real,intent(inout),dimension(0:jmax+1) :: rho,e,p,c,q,dm,dx
real,intent(inout),dimension(-1:jmax+1):: u,x
integer :: i

if (IC.EQ.1) then
	rho = 1.
	p   = 1.
	e   = p/((gamma-1.)*rho)
	u(-1)        = uL     ! Required to be equal to the piston speed on left (moving)
	u(0:jmax)    = 0.
	u(jmax+1)    = uR     ! Required to be equal to the piston speed on right (stationary)
else if (IC.EQ.2) then
	do i=0,jmax+1
		if (x(i).LT..5) then
			rho(i) = 1. 
			p(i)   = 1.
			u(i)   = delta
		else 
			rho(i) = .125
			p(i)   = .1
			u(i)   = delta
		end if 
	end do 
	u(-1) = delta
	e = p/((gamma-1.)*rho)
end if 
dm  = rho*dx
end subroutine 

subroutine Timestep(u,dx,dt,c,i)
use Constants 
implicit none
real,intent(inout) :: dt 
real,intent(in),dimension(-1:jmax+1) :: u
real,intent(in),dimension(0:jmax+1)  :: dx,c
integer,intent(in) :: i
real :: dt_temp 
dt_temp = dt                            ! Previous time step
dt = COURANT*MINVAL(dx/(ABS(.5*(u(-1:jmax)+u(0:jmax+1)))+c))
if (i.NE.0) then
	dt = MIN(dt,1.1*dt_temp)  				   ! Doesnt allow time step to grow too fast
end if 
end subroutine 

subroutine Output(d,e,p,u,dm,time,x,r,itter)
use Constants
implicit none
real,intent(in),dimension(0:jmax+1)   :: d,e,p,dm,r
real,intent(in),dimension(-1:jmax+1) :: x,u
integer,intent(in) :: itter
real,intent(in) :: time
integer :: st,en,choice ,j
choice = 0 ! Set choice to 0 to prin out real domain, 1 to print out entire domain (including ghost cells)

st = 1
en = jmax


if (MOD(itter,writeStep)==0) then ! Only write out itter modulus writeStep (speeds up plotting)
    !write(1,"(1D12.5)") d(1:jmax) ! Rho 
    write(1,*) d(1:jmax)
	write(2,*) u(0:jmax) ! Velocity
    write(3,*) e(1:jmax) ! Internal Energy
    write(4,*) p(1:jmax) ! Pressure
	write(7,*) x(0:jmax)
	!write(7,"(6(1PD12.5))") (u(j),(dt),(p(j)+q(j)-p(j+1)-q(j+1)),q(j),dm(j),x(j) ,j=-1,jmax)  ! interface
	write(8,*) time      ! Current Time
    write(9,*) r(1:jmax) ! Lagrange Cell Ceters
    write(10,*) SUM(dm)               ! Total mass
	write(11,*) itter
end if 
end subroutine

subroutine OpenFiles()
implicit none
! Open Files to write out to 
open(unit=1,file='Output/Density.txt')
open(unit=2,file='Output/Velocity.txt')
open(unit=3,file='Output/InternalEnergy.txt')
open(unit=4,file='Output/Pressure.txt')
open(unit=7,file='Output/LagnCellEdge.txt')
open(unit=8,file='Output/CurrentTime.txt')
open(unit=9,file='Output/LagnCellCenter.txt')
open(unit=10,file='Output/TotalMass.txt')
open(unit=11,file='Output/Step.txt')
end subroutine

subroutine CloseFiles()
implicit none
! Open Files to write out to 
close(1) ! Close density
close(2) ! Close velocity
close(3) ! Close Internal e
close(4) ! Close pressure
close(7) ! Close cell interface
close(8) ! Close dt for my solution
close(9) ! Close Lagrage cell ceter
close(10)! Close total mass
end subroutine
