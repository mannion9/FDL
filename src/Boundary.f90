subroutine Boundary(d,u,p,dv,dm)
use Constants
implicit none
real,intent(inout),dimension(imin:imax) :: d,p,dv,dm
real,intent(inout),dimension(imin:imax+1) :: u

if (BC.EQ.'Free') then
    u(imin:jmin)     = u(jmin)
    u(jmax+1:imax+1) = u(jmax)
else if (BC.EQ.'Reflective') then
    u(imin:jmin)     = -u(jmin)
    u(jmax+1:imax+1) = -u(jmax)
end if

if (IC.EQ.0) then
    u(imin:jmin)     = 100.
    u(jmax+1:imax+1) = u(jmax)
end if
p(imin:jmin-1) = p(jmin)
p(jmax+1:imax) = p(jmax)
d(imin:jmax-1) = dm(imin:jmax-1)/dv(imin:jmax-1)
d(jmax+1:imax) = dm(jmax+1:imax)/dv(jmax+1:imax)
end subroutine
