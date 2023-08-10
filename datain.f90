!-----------------------------------------------------------------------------------------------------------------------------------
! data importing subroutine
!
! imports system material properties, geometry parameters, grid parameters, and simulation parameteres
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine datain

	use varinit
	implicit none

!	importing solvent properties
	open(unit=3,file='../../../../in/dlss/solvent.csv',status='old',action='read')
	read(3,*)
	read(3,*)	rho, mu, rho_inf, sigma
	close(3)

!	importing electrolyte properties
	open(unit=1,file='../../../../in/dlss/electrolyte.csv',status='old',action='read')
	read(1,*)
	read(1,*)	nspecies
	close(1)
	allocate(diffcoeff(1:nspecies))
	allocate(valence(1:nspecies))
	open(unit=1,file='../../../../in/dlss/electrolyte.csv',status='old',action='read')
	read(1,*)
	read(1,*)	nspecies, n0, relperm, (diffcoeff(ispecies), ispecies=1,nspecies), &
				(valence(ispecies), ispecies=1,nspecies)
	close(1)
	
!	importing solid properties
	open(unit=2,file='../../../../in/dlss/solid.csv',status='old',action='read')
	allocate(zeta(2))
	read(2,*)
	read(2,*)	capE, nu, capAsfw, capLambda1, capGamma11, capGamma12, phi1, capLambda2, capGamma21, capGamma22, phi2, &
				zeta(1), zeta(2)
	close(2)

!	importing universal constants
	open(unit=7,file='../../../../in/dlss/univ.csv',status='old',action='read')
	read(7,*)
	read(7,*)	q, kcapB, permvac

!	importing geometry parameteres
	open(unit=4,file='../../../../in/dlss/geometry.csv',status='old',action='read')
	read(4,*)
	read(4,*)	h0, capD, capL, capR, omega, capT
	close(4)


!	importing grid parameteres
	open(unit=5,file='../../../../in/dlss/grid.csv',status='old',action='read')
	read(5,*)
	read(5,*)	adr, ady, adt, madr, nadr, mady, nady, madt, nadt
	read(5,*)
	read(5,*)	nophd, phditer, phdpatch
	read(5,*)
	read(5,*)	phdaxi, profaxi, molaxi
	read(5,*)
	read(5,*)	schemenum, schemeana, numguess, anaguess, revtime, toladapt
	read(5,*)
	read(5,*)	soln, minout, termout, sanity
	read(5,*)
	read(5,*)	tolnod, tolphd, thresamp, thresdiv1, thresdiv2, dl, bsecshrink
	read(5,*)
	read(5,*)	relaxnod, relaxphd
	read(5,*)
	read(5,*)	nr, ny, nz, nt
	read(5,*)
	read(5,*)	rmin, rmax, ymin, ymax, tmin, tmax	
	close(5)

!	importing problem prob
	open(unit=60,file='../../../../in/dlss/type.csv',status='old',action='read')
	read(60,'(A)')	probtype
	close(60)

!	getting timestamp
	call idate(today)
	call itime(now)
	write(datestamp,10)			today(3),today(2),today(1)
	write(siminstance,20)		now(1),now(2),now(3)
	siminstance = trim(probtype)//"_"//siminstance

	call execute_command_line(	"mkdir -p ../../../../out/dlss_num/"//datestamp//"/"//trim(siminstance),exitstat=itemp)
	call execute_command_line(	"mkdir -p ../../../../out/dlss_num/"//datestamp//"/"//trim(siminstance)// &
								"/jacobian",exitstat=itemp)

	10 format (I4,I2.2,I2.2)
	20 format (I2.2,I2.2,I2.2)

end subroutine datain
!-----------------------------------------------------------------------------------------------------------------------------------
