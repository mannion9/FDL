subroutine Output(d,e,h,p,u,dm,time,x,r,itter)
use Constants
implicit none
real,intent(in),dimension(imin:imax  ) :: d,e,h,p,dm,r
real,intent(in),dimension(imin:imax+1) :: x,u
integer,intent(in) :: itter
real,intent(in) :: time
write(1,*) d(jmin:jmax)
write(2,*) .5*(u(jmin+1:jmax+1)+u(jmin:jmax)) ! Average velocity at cell center
write(3,*) e(jmin:jmax) ! Internal Energy
write(4,*) p(jmin:jmax) ! Pressure
write(7,*) x(jmin:jmax+1)
write(8,*) time      ! Current Time
write(9,*) r(jmin:jmax) ! Lagrange Cell Ceters
write(10,*) SUM(dm)               ! Total mass
write(11,*) itter
write(12,*) h(jmin:jmax)
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
open(unit=12,file='Output/TotalEnergy.txt')
open(unit=13,file='Output/InitialState.txt')
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
close(11)! Close steps
close(12)! Close total energy
end subroutine
