!-----------------------------------------------------------------------------------------------------------------------------------
! variable declaration module
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
module varinit

!	variables with values to be loaded from input files

!	system property values
!	fluid/solvent
!											density, dynamic viscosity, bulk number density, hard-core radius
	double precision :: 					rho, mu, rho_inf, sigma
!	electrolyte
!											number of species
	integer :: 								nspecies
!											electroneutral number density, relative permittivity
	double precision :: 					n0, relperm
!											diffusion coefficients, valencies
	double precision, allocatable :: 		diffcoeff(:), valence(:)
!	solid
!											Young's modulus, Poisson's ratio, Hamaker's constant &
!											solvation pressure tuning parameters
	double precision :: 					capE, nu, capAsfw, &
											capLambda1, capGamma11, capGamma12, phi1, capLambda2, capGamma21, capGamma22, phi2
!											zeta potential (is dependent on the electrolyte properties as well)
	double precision, allocatable ::		zeta(:)
!	universal constants
!											universal charge, Boltzmann constant, permittivity of free space
	double precision :: 					q, kcapB, permvac
	double precision, parameter :: 			pi=3.141592653589793
!	geometry parameters
!											probe oscillation amplitude, mean separation of probe from origin, 
!											undeformed substrate thickness, probe radius, probe oscillation frequency, temperature
	double precision :: 					h0, capD, capL, capR, omega, capT
!	simulation parameters
!											consider adaptive grid for r, y, and t (0 = no; 1 = yes)
	integer :: 								adr, ady, adt
!											m and n (skew parameters of	original and uniform grids) for r, y, and t
	double precision :: 					madr, nadr, mady, nady, madt, nadt
!											switch to not iterate over hd pressure solutions (0 = false; 1 = true), neglect
!											effect of substrate speed in hd pressure solution (0 = false; 1 = true) 
!											applicable only for semi-analytical methodology; patch the solution at r=1 with 
!											solution obtained in the last time-step (0 = false, 1 = true)
	integer :: 								nophd, phditer, phdpatch
!											obtain centerline solution from axisymmetry for hd pressure, substrate deflection, 
!											molecular pressures
	integer :: 								phdaxi, profaxi, molaxi
!											scheme to be used for numerical solution (TBD), scheme to be used for analytical 
!											solution (0 = NR; 1 = NR-direction; 2 = brute-force; 3 = marching in the direction 
!											determined by the pressure at the guess considered; 4 = bisection method by taking 
!											equal steps in each direction starting at the guess and executing bisection-method on
!											the bracket identified first; 5 = same as 4 but with steps normalized according to 
!											energy, i.e. integration of absolution pressure times deflection being), 
!											guess to be taken for the numerical solution (0 = zero; 1 = analytical solution of 
!											current time-step; 2 = numerical solution of last time step), guess to be taken for the 
!											analytical solution (0 = zero; 1 = guess taken as solution of last time step with 
!											re-scaled; 2 = guess taken as solution of last time step without rescaling), start 
!											oscillations from the bottom (0 = false; 1 = true), adapt node-wise error tolerance with
!											the amplification 
	integer :: 								schemenum, schemeana, numguess, anaguess, revtime, toladapt
!											which solution to obtain (0 = analytical; 1 = numerical; 2 = both), reduce output 
!											written to file (0 = none; 1 = retain centerline, quarterly, and, dimensional r, l, 
!											and p values; 2 = retain only centerline, and quarterly values), show progress of
!											simulation on the terminal (0 = false; 1 = true), perform a sanity check on smallness of 
!											appropriate parameters (0 = false; 1 = true)
	integer :: 								soln, minout, termout, sanity
!											error tolerance for node-wise iterations, error tolerance for hd pressure 
!											convergence, threshold for switching to semi-analytical methodology based on 
!											amplification eta*M, threshold on eta*M for switching to semi-analytical methodology 
!											based on divergence N, threshold on N for switching to semi-analytical methodology 
!											based on divergence N, step size for marching or bisection method of root-finding for l,
!											energy-based normalizer exponent for step-size for the bisection method
	double precision :: 					tolnod, tolphd, thresamp, thresdiv1, thresdiv2, dl, bsecshrink
!											relaxation parameter for root-finding scheme, for nodes and for hd pressure 
	double precision :: 					relaxnod, relaxphd
!											number of grid points in r-coordinate, y-coordinate, z-coordinate, and t-coordinate
	integer :: 								nr, ny, nz, nt
!											minimum and maximum values of r, y (reversed), z, and t (values should be restricted
!											in [0,1], [0,-1] for r and t
	double precision :: 					rmin, rmax, ymin, ymax, tmin, tmax		
	
!	variables to be used during simulation and written to output

!	field variables (analytical)
!											analytical solution of total pressure, hd pressure, and deflection
	double precision, allocatable :: 		pana(:,:), pstana(:,:), lana(:,:)
!											analytical solution of EDL disjoining pressure, van der Waals disjoining pressure and
!											solvation pressure
	double precision, allocatable ::		pDLana(:,:), pvdWana(:,:), pSana(:,:)
!	field variables (analytical, dimensional)
!											dimensional value of analytical solution of total pressure, hd pressure, and 
!											deflection dimensional value of analytical solution of EDL disjoining pressure, van 
!											der Waals disjoining pressure and solvation pressure
	double precision, allocatable :: 		panad(:,:), pstanad(:,:), lanad(:,:)
	double precision, allocatable ::		pDLanad(:,:), pvdWanad(:,:), pSanad(:,:)
!	solution grid
!											time, radial coordinate, and probe profile grid
	double precision, allocatable :: 		t(:), r(:), capH(:)
!											dimensional value of time, radial coordinate, and probe profile grid
	double precision, allocatable :: 		td(:), rd(:), capHd(:)
!	miscellaneous
!	system parameters (to be calculated)
!	non-dimensional parameters from the mathematical formulation (refer to the manuscript/paper)
	double precision :: 					epsi, alpha, delta, capK, kappa, eta, gamm
	double precision, allocatable ::		epsi1(:), kappa1(:), eta1(:), gamma1(:)
	double precision, allocatable ::		coeffmom(:,:), coeffnernst(:,:), coeffpoisson(:), coeffre(:), coeffsolidbc(:,:)
!	index iterators
	integer :: 								itemp, itemp1
	integer :: 								izeta, ispecies
	integer :: 								iterphd, iternod
	integer :: 								ir, iy, iz, it, jr, jy, jz, jt, ir1, iy1, iz1, it1, ir2, iy2, iz2, it2
!	temporary value holders
	double precision :: 					containers(20)
	double precision :: 					atemp(1:20,-2:2,-2:2)
!	temporary value holders (variables)
	double precision, allocatable :: 		lana1(:,:), pDLana1(:,:), pvdWana1(:,:), pSana1(:,:)
	double precision, allocatable :: 		pSanaa(:,:), pSanab(:,:)
	integer ::								nt1
	double precision ::						tmax1
	double precision :: 					lanakcapL, lanakcapR, lanacapLl, lanacapLr, lanacapRl, lanacapRr
	double precision :: 					lanakcapLprev, lanakcapRprev, lanacapLlprev, lanacapLrprev, lanacapRlprev, lanacapRrprev
	double precision ::						lanak, lanakprev, lanakmin
	double precision :: 					panakcapL, panakcapR, panacapLl, panacapLr, panacapRl, panacapRr
	double precision :: 					panakcapLprev, panakcapRprev, panacapLlprev, panacapLrprev, panacapRlprev, panacapRrprev
	double precision ::						panak, panakprev
!	environment generation variables
	integer*4 :: 							today(3), now(3)
	character (len=8) :: 					datestamp
	character (len=15) :: 					siminstance
	character (len=15) :: 					probtype
	character (len=2) :: 					dir
	logical :: 								exists(20)
	character :: 							chartemp
	integer,allocatable :: 					regime(:)
!	last time step and last iteration containers (system variables)
	double precision, allocatable :: 		lanalast(:,:), lanalast2(:,:)
	double precision :: 					lanaprev
!	derivative value holders
	double precision :: 					dlanadt
	double precision :: 					dlana0dt, dlana1dt, dlana3dt, dlana0dtprev, dlana1dtprev, dlana3dtprev
!	error and residual variables
	double precision :: 					Fnod, Jnod, Fnodprev, Jnodprev, Fphd, Jphd, Fphdprev, Jphdprev
	double precision :: 					FnodkcapL, FnodkcapR, FnodcapLl, FnodcapLr, FnodcapRl, FnodcapRr
	double precision :: 					FnodkcapLprev, FnodkcapRprev, FnodcapLlprev, FnodcapLrprev, FnodcapRlprev, FnodcapRrprev
	double precision ::						Fnodk, Fnodkprev
	double precision :: 					EcapL, EcapR, dlanak, dlanakcapL, dlanakcapR
	double precision, allocatable :: 		capMamp(:), capNdiv(:)
!	support variables
	integer, allocatable :: 				ipivtest(:)
	integer :: 								iline
end module varinit
!-----------------------------------------------------------------------------------------------------------------------------------
