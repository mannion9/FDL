subroutine InitDomain(x,dv,r)
use Constants
implicit none
real,intent(inout),dimension(imin:imax+1) :: x
real,intent(inout),dimension(imin:imax  ) :: dv,r
integer :: i
dv = (xmax-xmin)/REAL(N-1)
r  =(/(real(i-jmin)*dv(i),i=imin,imax,1)/)										! The i-jmin would be just i, but we need to be able to start our counting from any number for LILAC
x(imin:imax+1) = (/(dv(jmin)*(.5+real(i-jmin)),i=imin-1,imax,1)/)
end subroutine


subroutine InitCond(d,e,h,p,u,m,dmi,dv,x,qa,r)
use Constants
implicit none
real,intent(inout),dimension(imin:imax  ) :: d,e,h,m,dv,qa,r,p
real,intent(inout),dimension(imin:imax+1) :: u,x,dmi
real :: dL,dR,uL,uR,pL,pR,x0=(xmax-xmin)/2.
integer :: i

! Piston problem
if (IC.EQ.0) then
	dL = 1.
	dR = dL
	pL = 1.
	pR = pL
	uL = delta
	uR = uL
else if (IC.EQ.1) then
	! Modified sod test
	dL = 1.
	dR = dL/8.
	pL = 1.
	pR = pL/10.
	uL = delta !0.75
	uR = delta
else if (IC.EQ.2) then
	dL = 1.
	dR = dL
	pL = 0.4
	pR = pL
	uL = -2.
	uR = -uL
else if (IC.EQ.3) then
	dL = 1.
	dR = dL
	pL = 1000.
	pR = 0.01
	uL = delta
	uR = uL
else if (IC.EQ.4) then
	dL = 5.99924
	dR = 5.99242
	pL = 460.894
	pR = 46.0950
	uL = 19.5975
	uR = -6.19633
else if (IC.EQ.5) then
	dL = 1.
	dR = dL
	pL = 1000.
	pR = 0.01
	uL = -19.59745
	uR = uL
else if (IC.EQ.6) then
! Strong compression
	dL = 1.4
	dR = dL
	pL = .7143
	pR = pL
	uL = 5.0
	uR = -uL
end if

do i = imin,imax
	if (r(i).LT.x0) then
		d(i) = dL
		p(i) = pL
		u(i) = uL
	else
		d(i) = dR
		p(i) = pR
		u(i) = uR
	end if
end do
u(imax+1)=uR

! Sedov Blast
if (IC.EQ.7) then
	d = 1.
	u   = delta
	p(imin:INT(N/2)-1) = delta
	p(INT(N/2))        = 3*(1.4-1.)*d(3)  ! Give internal energy of .3 to orgin cell (this just tells us the pressure, given .3 internal energy)
	p(INT(N/2)+1:imax) = delta
end if

! Collela-Woodward Blast
if (IC.EQ.8) then
	d = 1.
	u   = .0
	do i=imin,imax
		if (r(i).LT.0.1*(xmax-xmin)) then
			p(i)   = 1000.
		else if (r(i).GE.0.1*(xmax-xmin).AND.r(i).LT.0.9*(xmax-xmin)) then   	   ! Left state
			p(i) = 0.01
		else  ! Right state
			p(i) = 100.0
		end if
	end do
end if

if (IC.EQ.9) then
	d = 1.
	p = 1.
	u = 1.
end if

! Assumes ideal gas at initial state (makes intial conditions easier to match with others)
e  = p/((gamma-1.)*d)
h  = e + .5*(.5*(u(imin+1:imax+1)+u(imin:imax)))**2
m  = d*dv
dmi(jmin:jmax+1) = .5*(m(jmin-1:jmax)+m(jmin:jmax+1))
qa = (/(qo*d(i)*MIN(0.0,(u(i+1)-u(i)))**2.,i=imin,imax,1)/)

! Write this file for Riemann Solver
write(13,*) dL
write(13,*) uL
write(13,*) pL
write(13,*) dR
write(13,*) uR
write(13,*) pR
write(13,*) x0
end subroutine
