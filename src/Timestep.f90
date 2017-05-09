subroutine Timestep(u,dv,dth,dtf,p,d)
use Constants
implicit none
real,intent(inout) :: dth,dtf
real,intent(in),dimension(imin:imax+1) :: u
real,intent(in),dimension(imin:imax  ) :: dv,p,d
real,dimension(imin:imax) :: a
real :: dtha
dtha = dth
a = SQRT(gamma*p/d)
dth = COURANT*MINVAL(dv/(a+abs(.5*(u(imin:imax)+u(imin+1:imax+1)))))
! dth = COURANT*MINVAL(dv/a)
dth = MIN(dth,1.1*dtha)
dtf = .5*(dth+dtha)
end subroutine
