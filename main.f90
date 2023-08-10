!-----------------------------------------------------------------------------------------------------------------------------------
! dynamic loading on soft-substrate solver for x-z coordinate system (transformed from r-z co-ordinate system) (wrapper)
!
! author -              pratyaksh karan (iit kharagpur; pratyakshkaran @ gmail.com)
! date -         		2019-01-17
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
program main

!	variable declaration
	use varinit
	use derivative
	implicit none

!	calling subroutines
!	importing input data from in-files
	call datain
!	allocating size for matrices and initializing
	call allocinit
!	calculating pertinent parameter values for the solution mechanism
	call paramcalc
!	generating mesh
	call meshgen
!	time-marching
	call march
!	finish with writing out any relevent data
	call finout
!	calling plot data generators
	call plotana
	call plotanad
!	purging heavy output files
	if (minout .ne. 0) then
		call purgemin
	end if

end program main
!-----------------------------------------------------------------------------------------------------------------------------------
