module constants
implicit none
!----------------------- USER INPUTS ----------------------------------------- !
integer,parameter :: N         = 100    ,& ! Number of cells+1
					 itterMax  = 4000   ,& ! Maximum itterations
					 writeStep = 10    ,& ! Step interval between writing out data
					 IC        =  1        ! Initial condition  (not working in 4)
real,parameter    :: Courant   = .4     ,&
					 qo        = 4.0    ,& ! Viscious pressure constant (Emperically found)
					 xmin      = 0.0    ,&
					 xmax	   = 1.0    ,&
					 tmax      = .3        ! Maximum time
character (LEN=*),parameter :: energy_solver = 'internal' ,& ! ('total','internal')
							   BC            = 'free'     !('free','refl')
!----------------------- Program Constants ----------------------------------------- !
real,parameter    :: gamma = 5./3.        ,& ! Assuming ideal gas EOS
					 delta = 1.e-30          ! A small number
integer,parameter :: jmin  = 0            ,& ! Index of left  most real  cell center
					 jmax  = (N-1)+jmin   ,& ! Index of right most real  cell center
					 imin  = jmin-1       ,& ! Index of left  most ghost cell center
					 imax  = jmax+1          ! Index of right most ghost cell center
end module
