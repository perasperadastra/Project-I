
!	Initialization parameters
!---------------------------
integer::nc=6 !Number of particles per side
double precision::density=1.d0 !Lomgitude of the size box

integer::structure=1 !For SC:1, for fcc:2, for diamond 3.

double precision::temp= 300.0d0  !Temperature (K)
integer::vel_opt=1 ! if 1, implemented bimodal. Else, set up to 0

!Simulation set up
!-----------------

! Total simulation steps
integer::ntimes=100000

! Multiple of steps at which the thermodynamic properties are saved
integer::everyt=100

! Time step at which melting is from initial structure is done
! and rdf can be computed
integer::tmelt=1000 

!Time step (ps)
double precision:: dt= 0.1d0 

integer::thermo=1 !For no thermostat:0, for active thermostat:1

!Lennard-Jones parameters
!----------------
double precision:: epsilon = 10.22d0, sigma = 2.556d0

!Andersen thermostat random number generator seed
!----------------
integer :: rng_seed = 14
