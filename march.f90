!-----------------------------------------------------------------------------------------------------------------------------------
! time marching subroutine
!
! marches across time steps and writes time step output to file (calls to the respective subroutines), and temporarily stores last
! time-step data for non-steady terms
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine march

    use varinit
    implicit none
	
	do ir = 1,nr
		capH(ir) =			1.0 + 0.5*r(ir)**2.0	
	end do	
	do it = 1,nt
		call paramcalc1
		if (termout .eq. 1) then
			write(*,*) "t = ",t(it)
		end if
		call anainst
		call anaout
		lanalast2 = lanalast
		lanalast = lana
		containers(20) = pstana(3,nr)
		containers(19) = pstana(1,nr)
	end do

end subroutine march
!-----------------------------------------------------------------------------------------------------------------------------------
