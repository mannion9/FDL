! 3) Try the different characteristic speeds, those used in q, and dt
! 5) Check shock tube initial condtion

module constants
implicit none
real,parameter    :: gamma    = 1.4 ,&
					 Courant  = .9  ,&
					 length   = 1.  ,&
					 tmax     = 2. ,&
					 qo       = 2.  ,&
					 delta    = 1.e-30 ,&
					 uL       = .5 ,&      !(uL=.5 for piston , uL=.0 for shock tube)
					 uR       = 1.e-30  	 !(uR=.0 for piston , uR=.0 for shock tube)
integer,parameter :: jmax 	  = 100  ,& ! number of cells
					 itterMax = 200  ,& ! Maximum itterations
					 writeStep= 10    ,& ! Itter step to write to file
					 IC       = 1       ! Initial condiiton (1-pisition,2 Shock)
end module

program main 
use constants
implicit none 
real :: dt=0.05 , time=.0  , dmi
real,dimension(0:jmax+1) :: rho,p,e,q,c,dx,dm,r
real,dimension(-1:jmax+1) :: u,x,uT
integer :: i , j

call OpenFiles()

!Create Domain
x   = (/((i*length/(SIZE(x)-1)),i=0,SIZE(x)-1,1)/)  ! Make domain uniformly spaced
dx  =  x(0:jmax+1)-x(-1:jmax)   					! Spacing of cells
r   = .5*(x(0:jmax+1)+x(-1:jmax))         			! Cell centers
! Initial Conditions
call Initialcondition(rho,e,p,c,q,u,dm,dx,x)
i = 0
do while (time.LT.tmax.AND.i.LT.itterMax)
	call Boundary(u)
	call Timestep(u,dx,dt,c)
	call Output(rho,e,p,u,dm,time,x,r,i)
	x(-1)  = x(-1) + dt*u(-1)
	do j=0,jmax  
		dmi    = .5*(rho(j+1)+rho(j))*(r(j+1)-r(j))    ! Mass in "momentum zone" of cell interface of j and j+1 cell
		x(j)   = x(j) + u(j)*dt
		uT(j)  = u(j) + (dt/dmi)*(p(j)+q(j)-p(j+1)-q(j+1)) ! Update the velocity in a temporary array
		rho(j) = dm(j)/(x(j)-x(j-1))
		e(j)   = e(j) + (dt/dm(j))*(q(j)+p(j))*(u(j-1)-u(j))
		p(j)   = (gamma-1.)*rho(j)*e(j)
			
		if (uT(j-1)-uT(j).GT..0) then  ! If updated cell is compressed add viscous pressure
			! Linear
			!c(j)   = sqrt(gamma*p(j)/rho(j))
			!q(j) = qo*c(j)*rho(j)*(uT(j-1)-uT(j))
			! NonLinear
			q(j) = qo*rho(j)*(uT(j-1)-uT(j))**(2)
			! Combo
			!c(j)   = sqrt(gamma*p(j)/rho(j))
			!q(j) = (qo/10)*c(j)*rho(j)*(uT(j-1)-uT(j))
			!q(j) = q(j)+qo*rho(j)*(uT(j-1)-uT(j))**(2)
		else 
			q(j) = .0
		end if 
		c(j)   = sqrt(gamma*(p(j)+q(j))/rho(j))
	end do 
	u   = uT						 ! Update the velocity 
	r   = .5*(x(0:jmax+1)+x(-1:jmax))! Find new cell centers
	dx  =  x(0:jmax+1)-x(-1:jmax)   ! Find new cell widths
	i = i + 1						! Update itteration counter
	time = time+ dt					! update timer
	if (MINVAL(dx).LT..0) then
		print*,'Negative volume'
	end if 
end do 
print*,i
call CloseFiles()
end program

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

subroutine InitialCondition(rho,e,p,c,q,u,dm,dx,x)
use Constants
implicit none 
real,intent(inout),dimension(0:jmax+1) :: rho,e,p,c,q,dm,dx
real,intent(inout),dimension(-1:jmax+1):: u,x
integer :: i
if (IC.EQ.1) then
	rho = 1.
	p   = delta
	e   = p/((gamma-1.)*rho)
	u(-1)        = uL     ! Required to be equal to the piston speed on left (moving)
	u(0:jmax)    = delta
	u(jmax+1)    = uR ! Required to be equal to the piston speed on right (stationary)
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


do i=0,jmax+1
	if (u(i-1)-u(i).GT..0) then
		! Linear
		!c   = sqrt(gamma*p/rho)
		!q(i) = qo*c(i)*rho(i)*(u(i-1)-u(i))
		! Nonlinear
		q(i) = qo*rho(i)*(u(i-1)-u(i))**(2)
		! Combo
		!c   = sqrt(gamma*p/rho)
		!q(i) = (qo/10)*c(i)*rho(i)*ABS(u(i-1)-u(i))
		!q(i) = q(i) + qo*rho(i)*(u(i-1)-u(i))**(2)
	else 
		q(i) = .0
	end if 
end do 
c   = sqrt(gamma*(p+q)/rho)
dm  = rho*dx
end subroutine 

subroutine Timestep(u,dx,dt,c)
use Constants 
implicit none
real,intent(inout) :: dt 
real,intent(in),dimension(-1:jmax+1) :: u
real,intent(in),dimension(0:jmax+1)  :: dx,c
real :: dt_temp 
dt_temp = dt                            ! Previous time step
dt = COURANT*MINVAL(dx)/(MAXVAL(ABS(u(-1:jmax))+c))
dt = MIN(dt,2.*dt_temp)  				   ! Doesnt allow time step to grow too fast
end subroutine 

subroutine Output(d,e,p,u,dm,time,x,r,itter)
use Constants
implicit none
real,intent(in),dimension(0:jmax+1)   :: d,e,p,dm,r
real,intent(in),dimension(-1:jmax+1) :: x,u
integer,intent(in) :: itter
real,intent(in) :: time
integer :: st,en,choice 
choice = 0 ! Set choice to 0 to prin out real domain, 1 to print out entire domain (including ghost cells)

st = 1
en = jmax


if (MOD(itter,writeStep)==0) then ! Only write out itter modulus writeStep (speeds up plotting)
    write(1,*) d(1:jmax) ! Rho 
    write(2,*) u(0:jmax) ! Velocity
    write(3,*) e(1:jmax) ! Internal Energy
    write(4,*) p(1:jmax) ! Pressure
	write(7,*) x(0:jmax) ! interface
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