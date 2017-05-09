! The length of the time array is "unkown" for this program because this depends on the MAIN programs result and t_max,
! but it is limited by itterMAX from the MAIN program. We deal with this after loading and assining values to time
! we find the MAXLOC(time); this is the final time.
! Primative W = (rho,u,p)^(T)
module constants1
use constants
implicit none
real,parameter :: gm1   = gamma-1.  ,&
                  gp1   = gamma+1.
end module

program RiemannSolver
use Constants
use constants1
implicit none
integer,parameter :: res = 1.
real :: p_star,u_star,f,fp,fl,fr,guess,S,center,newton,a_r,a_l,dv
real,dimension(jmin:jmax)    :: r,d,u,p,e
real,dimension(0:itterMax)     :: time = 0.d0
real :: dl,pl,ul,dr,pr,ur
integer :: i,j,Ntmax

call OpenFiles()
call ReadInputs(dl,pl,ul,dr,pr,ur,center,time,Ntmax,r(jmin:jmax))
a_l = SQRT(gamma*pl/dl)
a_r = SQRT(gamma*pr/dR)

! Find p_star and u_star
guess = .5*(pl+pr)
!guess = .5*(pl+pr)-.125*(ur-ul)*(dl+dr)*(a_l+a_r)
p_star = newton(guess,dl,pl,ul,dr,pr,ur)
CALL f_func(p_star,dl,pl,ul,dr,pr,ur,f,fp,fl,fr)
u_star = .5*(ul+ur)+.5*(fr-fl)

dv = (r(jmax)-r(jmin))/REAL(res*jmax-jmin)
r  = (/(real(i)*dv,i=jmin,res*jmax,1)/)
write(10,*) r
! Time equal zero
do i=jmin,res*jmax
	if (r(i).LE.center) then
		d(i) = dl
		u(i) = ul
		p(i) = pl
		e(i) = p(i)/(gm1*d(i))
	else
		d(i) = dr
		u(i) = ur
		p(i) = pr
		e(i) = p(i)/(gm1*d(i))
	end if
end do
write(4,*) d
write(7,*) u
write(8,*) p
write(9,*) e

! Future time
do j=0,Ntmax
	do i=jmin,res*jmax
        S = (r(i)-center)/time(j+1)  ! The j+1 is required because t(0)=0.00, but we have already written out t=0.0 data
		CALL sample(S,u_star,p_star,dl,pl,ul,dr,pr,ur,d(i),u(i),p(i))
		e(i) = p(i)/(gm1*d(i))
	end do
	write(4,*) d
	write(7,*) u
	write(8,*) p
	write(9,*) e
end do
call CloseFiles()
end program

subroutine sample(S,u_star,p_star,dl,pl,ul,dr,pr,ur,d,u,p)
!-------------------------------------------------
! Sample determines for the value of the primates
! in each at a specifc location
! INPUTS:
!   S     : r/t is the location in r-t space
!   u_star: U in starred region
!   p_star: P in starred region
!   w_L   : primative variables left
!   w_R   : primative variables right
! OUTPUTS:
!   w     : primative variables at r/t
!-------------------------------------------------
use Constants
use constants1
implicit none
real,intent(in) :: S,u_star,p_star
real,intent(in) :: dl,pl,ul,dr,pr,ur
real,intent(inout) :: d,u,p
real :: a_l,a_r, a_l_star,a_r_star, &
        S_l,S_l_tl,S_l_hd, &
        S_r,S_r_tl,S_r_hd, &
        rho_l_star_shock,rho_r_star_shock, &
        rho_l_star_rare,rho_r_star_rare, &
        rho_l_fan,rho_r_fan, &
        u_l_fan,u_r_fan,p_l_fan,p_r_fan
real,dimension(0:2) :: w_l_star,w_r_star,w_l_fan_star,w_r_fan_star, w_r_fan,w_l_fan,w

a_l = SQRT(gamma*pl/dl)   ! Toro 3.6
a_r = SQRT(gamma*pr/dr)
a_l_star = a_l*(p_star/pl)**(gm1/(2.*gamma))  ! Toro 4.54/4.61
a_r_star = a_r*(p_star/pr)**(gm1/(2.*gamma))
rho_l_star_shock = dl*(((p_star/pl)+gm1/gp1)/(gm1*p_star/(gp1*pl)+1)) ! Toro 4.50/4.57
rho_r_star_shock = dr*(((p_star/pr)+gm1/gp1)/(gm1*p_star/(gp1*pr)+1))
rho_l_star_rare  = dl*(p_star/pl)**(1./gamma) ! Toro 4.53/4.60
rho_r_star_rare  = dr*(p_star/pr)**(1./gamma)
rho_l_fan = dl*(2./gp1+gm1/(a_l*gp1)*(ul-S))**(2./gm1) ! Toro 4.56/4.63
rho_r_fan = dr*(2./gp1-gm1/(a_r*gp1)*(ur-S))**(2./gm1)
u_l_fan = (2./gp1)*(a_l+(gm1/2.)*ul+S)   ! Toro 4.56/4.63
u_r_fan = (2./gp1)*(-a_r+(gm1/2.)*ur+S)
p_l_fan = pl*(2./gp1+gm1/(gp1*a_l)*(ul-S) )**(2.*gamma/gm1)
p_r_fan = pr*(2./gp1-gm1/(gp1*a_r)*(ur-S) )**(2.*gamma/gm1)
S_l = ul - a_l*SQRT(gp1*p_star/(2.*gamma*pl)+(gm1/(2.*gamma))) ! Toro Eq. 4.52
S_r = ur + a_r*SQRT(gp1*p_star/(2.*gamma*pr)+(gm1/(2.*gamma))) ! Toro Eq. 4.52
S_l_hd = ul - a_l       ! Toro 4.55
S_l_tl = u_star - a_l_star
S_r_hd = ur + a_r       ! Toro 4.62
S_r_tl = u_star + a_r_star
! Shock
w_l_star = (/rho_l_star_shock,u_star,p_star/)
! Rarefaction (in starred region)
w_l_fan_star = (/rho_l_star_rare,u_star,p_star/)
! Rarefaction (within fan)
w_l_fan = (/rho_l_fan,u_l_fan,p_l_fan/)
! Shock
w_r_star = (/rho_r_star_shock,u_star,p_star/)
! Rarefaction (in starred region)
w_r_fan_star = (/rho_r_star_rare,u_star,p_star/)
! Rarefaction (in fan)
w_r_fan = (/rho_r_fan,u_r_fan,p_r_fan/)

if (S.LE.u_star) then
  ! Left of discontinuity
  if (p_star.GT.pl) then
    ! Shock
    if (S.LE.S_l) then
      ! To the left of the shock
      d = dl
      u = ul
      p = pl
    else
      ! To the right of the shock
      d = rho_l_star_shock
      u = u_star
      p = p_star
      !w= w_l_star  ! Will have conditions in star region
    end if
  else
    ! Rarefaction
    if (S.LE.S_l_hd) then
      ! To the right of head
      d = dl
      u = ul
      p = pl
    else if (S.GE.S_l_hd .AND. S.LE.S_l_tl) then
     ! Inside the fan
      d = rho_l_fan
      u = u_l_fan
      p = p_l_fan
    else
     ! In the star region
     d = rho_l_star_rare
     u = u_star
     p = p_star
    end if
  end if
else if (S.GE.u_star) then
  ! Right of disconinuity
  if (p_star>pr) then
    ! Shock
    if (S.GE.S_r) then
      ! To the right of the shock
      d = dr
      u = ur
      p = pr
    else
      ! To the left of the shock
      d = rho_r_star_shock
      u = u_star
      p = p_star
      !w = w_r_star
    end if
  else
    ! Rarefaction
    if (S.GE.S_r_tl .AND. S.LE.S_r_hd) then
      ! Inside the fan
      d = rho_r_fan
      u = u_r_fan
      p = p_r_fan
      !w = w_r_fan
    else if (S.GE.S_r_hd) then
      ! To the right of head
      d = dr
      u = ur
      p = pr
    else
      ! In the star region
     d = rho_r_star_rare
     u = u_star
     p = p_star
      !w = w_r_fan_star
    end if
  end if
end if
end subroutine

function newton(guess,dl,pl,ul,dr,pr,ur) result(p_star)
!---------------------------------------------
! Newton determines the root of an equation
! using newtons method.
!
! INPUTS:
!   guess : intial guess of value of root
!   w_l   : primative variables left     (0)
!   w_l   : primative variables right    (0)
!   f_func: the function we seek root of
! OUTPUT
!   p_star: value of p in stared region
!--------------------------------------------
implicit none
real :: guess,p_star
real :: dl,pl,ul,dr,pr,ur
real :: f,fp,fl,fr,it, itold, TOL=10.**(-6) , CHA
integer :: i, itter_max=20,break=0
it = guess
CHA =  guess
do i=1,itter_max
  itold = it
  if (CHA > TOL) then !.AND. i<itter_max) then
    CALL f_func(it,dl,pl,ul,dr,pr,ur,f,fp,fl,fr)
    it = itold - f/fp
    if (it < 0) then  ! Requires that the solution be positive for physical
      it = itold
    end if
    CHA = 2.*abs(it-itold)/(it+itold)
  else
    exit
  end if
end do
p_star = it
end  function newton

subroutine f_func(p,dl,pl,ul,dr,pr,ur,f,fp,f_l,f_r)
!---------------------------------------
! The function, whose root returns p*
! in the star region.
!
! INPUTS:
!   p  : pressure in sarted region (0)
!   w_l: primative variables left  (0)
!   w_r: primative variables right (0:2)
!   recall - W = (rho,u,p)^(T)
! OUTPUT:
!   f  : scalar (0)
!   fp : derivative of f
!   fl : function fl
!   fr : function fr
!---------------------------------------
use Constants
use constants1
implicit none
real,intent(in):: dl,pl,ul,dr,pr,ur
real,intent(in) :: p
real,intent(inout) :: f,fp,f_l,f_r
real :: f_lp,f_rp,del_u,A_l,A_r,B_l,B_r,al,ar
real :: A,B,ak,rho,pk,fk,fkp
integer:: i
del_u  = ur-ul
A_l    = (2./gp1)/dl   ! Toro 4.8
A_r    = (2./gp1)/dr
B_l    = (gm1/gp1)*pl
B_r    = (gm1/gp1)*pr
al     = SQRT(gamma*pl/dl)   ! Toro 3.6
ar     = SQRT(gamma*pr/dr)
do i=1,2
  if (i==1) then
    A    = A_l
    B    = B_l
    ak    = al
    rho  = dl
    pk   = pl
  else
    A    = A_r
    B    = B_r
    ak    = a_r
    rho = dr
    pk   = pr
  end if
  if (p>pk) then
    ! shock
      fk  = (p-pk)*SQRT(A/(p+B))
      fkp = SQRT(A/(B+p))*(1-(p-pk)/(2.*(B+p)))
  else
    ! rarefaction
      fk  = (2.*ak)/gm1*((p/pk)**(gm1/(2.*gamma))-1.)
      fkp = (1./(rho*ak))*(p/pk)**(-1.*gm1/(2.*gamma))
  end if
  if (i==1) then
    f_l  = fk
    f_lp = fkp
  else
    f_r  = fk
    f_rp = fkp
  end if
end do
f   = f_l + f_r + del_u
fp  = f_lp + f_rp
end subroutine

subroutine ReadInputs(dl,pl,ul,dr,pr,ur,center,time,Ntmax,r)
use Constants
implicit none
real,intent(inout),dimension(jmin:jmax) 	 :: r
real,intent(inout),dimension(0:itterMax) :: time
real :: dl,pl,ul,dr,pr,ur
integer,intent(inout) :: Ntmax
real,intent(inout) 	  :: center
integer :: i,check
real    :: reader

! Read in left and right state
read(1,*) dl
read(1,*) ul
reaD(1,*) pl
read(1,*) dr
read(1,*) ur
read(1,*) pr


! Read in discontinuity location
read(1,*) center

! Read in times
do i=0,itterMAX
	read(2,*,IOSTAT=check) reader
	if (check>0 .OR. check<0) then ! This statment will stop the loop if we have reached the end of the file
		exit
	else
		time(i) = reader
	end if
end do
Ntmax = MAXLOC(time,1)! This is the index of the last time value
! Read in cell centers
read(3,*) r
end subroutine

subroutine OpenFiles()
implicit none
! Open input files
open(unit=1,file='Inputs/InitialState.txt')
open(unit=2,file='Inputs/CurrentTime.txt')
open(unit=3,file='Inputs/LagnCellCenter.txt')
open(unit=4,file='Output/rho.txt')
open(unit=7,file='Output/velocity.txt')
open(unit=8,file='Output/pressure.txt')
open(unit=9,file='Output/energy.txt')
open(unit=10,file='Output/Position.txt')
end subroutine

subroutine CloseFiles()
implicit none
close(1)
close(2)
close(3)
close(4)
close(7)
close(8)
close(9)
close(10)
end subroutine
