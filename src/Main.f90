program main
use constants
implicit none
real                        :: dtf=0.0001,dth=0.0001, time=.0,etemp,t,fpdv,cv=3./2.
real,dimension(imin:imax  ) :: da,db,pa,pb,e,qa,qb,h,dva,dvb,dm,r    ! Located at cell centers
real,dimension(imin:imax+1) :: u,x,dmi                               ! Located at cell edges
integer                     :: i=0 , j
call OpenFiles()
call InitDomain(x,dva,r)
call InitCond(da,e,h,pa,u,dm,dmi,dva,x,qa,r)
call Output(da,e,h,pa,u,dm,time,x,r,i)

do i=0,itterMax-1
	call Boundary(da,u,pa,dva,dm)
	call Timestep(u,dva,dth,dtf,pa,da)

	! Move ghost cells to prevent negative volumes in ghost cells
		!(imin:jmin-1)   = x(imin  :jmin-1) + dth*u(jmin)    ! original
	! x(jmax+2:imax+1) = x(jmax+2:imax+1) + dth*u(jmax+1)
	x(imin:jmin-1)   = x(imin  :jmin-1) + dtf*u(jmin)
	x(jmax+2:imax+1) = x(jmax+2:imax+1) + dtf*u(jmax+1)
	dmi(jmin:jmax+1) = (/(.5*(da(j-1)+da(j))*(r(j)-r(j-1)),j=jmin,jmax+1,1)/)

	do j=jmin,jmax+1
	  ! u(j) = u(j) - (dtf/dmi(j))*(pa(j)+qa(j)-pa(j-1)-qa(j-1)) ! original
			u(j) = u(j) - (dth/dmi(j))*(pa(j)+qa(j)-pa(j-1)-qa(j-1))
	end do

	do j=jmin,jmax+1
	  ! x(j) = x(j) + dth*u(j)  ! orginal
			x(j) = x(j) + dth*u(j)
	end do

	dvb =  x(imin+1:imax+1)-x(imin:imax)    ! Find new cell widths

	do j=jmin,jmax
	  db(j) = dm(j)/dvb(j)
	end do

	qb = (/(qo*.5*(db(j)+da(j))*MIN(0.0,(u(j+1)-u(j)))**2.,j=imin,imax,1)/)

	do j= jmin,jmax
	  etemp = e(j) - pa(j)*(1./db(j)-1./da(j))
	  pb(j) = (gamma-1.)*db(j)*etemp
	  e(j)  = e(j) - .5*(pa(j)+qa(j)+pb(j)+qb(j))*(1./db(j)-1./da(j))
	end do

	! do j= jmin,jmax
	!   etemp = e(j) + .5*(u(j)+u(j+1))**2.
	!   etemp = etemp + dth*(u(j)*.5*(pb(j)+pb(j+1))-u(j+1)*.5*(pb(j)+pb(j-1)))/dm(j)
	!   e(j)  = etemp - .5/dm(j)*(u(j)+u(j+1))**2.
	! end do

	! Solving with temperature
	! do j = jmin,jmax
	! 	t = e(j)/cv
	! 	fpdv = (gamma-1.)*(dvb(j)-dva(j))/(dvb(j)+dva(j))
	! 	t = t*(1.-fpdv)/(1.+fpdv) - (dvb(j)+dva(j))*qb(j)/(cv*(1.+fpdv))
	! 	e(j) = t*cv
	! end do

	! Solving using just volumes
	! do j= jmin,jmax
	!   etemp = e(j) - pa(j)*(2.*(dvb(j)-dva(j))/(dvb(j)+dva(j)))
	!   pb(j) = (gamma-1.)*db(j)*etemp
	!   e(j)  = e(j) - .5*(pa(j)+qa(j)+pb(j)+qb(j))*(2.*(dvb(j)-dva(j))/(dvb(j)+dva(j)))
	! end do

	  pb  = (/(MAX((gamma-1.)*db(j)*e(j),delta),j=imin,imax,1)/)
	  r  = .5*(x(imin+1:imax+1)+x(imin:imax))! Find new cell centers
	  ! time = time + dth					             ! update timer  original
		time = time + dth					             ! update timer


	  if (MOD(i,writeStep).EQ.0) call Output(da,e,h,pb,u,dm,time,x,r,i)
	if (MINVAL(dvb).LT..0)  print*,'Negative volume in cell number',MINLOC(dvb)
	if (dtf.LT.1.e-6)       print*,'Time step has become smaller than 10^(-6)'
	if (time.GT.tmax.OR.MINVAL(dvb).LT.delta.OR.dtf.LT.1.e-6) exit
	da = db
	pa = pb
	dva= dvb
		qa = qb

end do
print*,'Terminated at step:',i
call CloseFiles()
end program
