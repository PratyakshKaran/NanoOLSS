!-----------------------------------------------------------------------------------------------------------------------------------
! mesh generation subroutine
!
! generates the mesh for x-y coordinate system
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine meshgen

	use varinit
	implicit none

!	generating t-grid
	if (revtime .eq. 0) then
		if (adt .eq. 0) then
			do it = 1,nt
				t(it) = 		tmin+((tmax-tmin)/(nt-1))*(it-1)
			end do
		else
			nt1 = (nt-1)/2+1
			tmax1 = tmin+((tmax-tmin)/2.0)
			do it1 = 1,nt1
				t(it1) = 		(1.0/(nt1-1))*(it1-1)
			end do
			containers(2) =		(madt*nadt*(tmax1-tmin)-nadt*tmax1)/ &
								((tmin-madt*(tmax1-tmin))*tmax1+madt*nadt*(tmax1-tmin)**2.0)
			containers(3) =		((madt-nadt)*(tmax1-tmin)-tmin)/ &
								((tmin-madt*(tmax1-tmin))*tmax1+madt*nadt*(tmax1-tmin)**2.0)
			containers(1) =		-containers(2)*tmin
			do it1 = 1,nt1
				t(it1) = 		(containers(1)-t(it1))/(t(it1)*containers(3)-containers(2))
				t(nt-it1+1) = 	tmax-t(it1)
			end do	
		end if
	else
		if (adt .eq. 0) then
			do it = 1,nt
				t(it) = 		omega*pi+tmin+((tmax-tmin)/(nt-1))*(it-1)
			end do
		else
			nt1 = (nt-1)/2+1
			tmax1 = tmin+((tmax-tmin)/2.0)
			do it1 = 1,nt1
				t(it1) = 		(1.0/(nt1-1))*(it1-1)
			end do
			containers(2) =		((1.0-madt)*(1.0-nadt)*(tmax1-tmin)-(1.0-nadt)*tmax1)/ &
								((tmin-(1.0-madt)*(tmax1-tmin))*tmax1+(1.0-madt)*(1.0-nadt)*(tmax1-tmin)**2.0)
			containers(3) =		(((1.0-madt)-(1.0-nadt))*(tmax1-tmin)-tmin)/ &
								((tmin-(1.0-madt)*(tmax1-tmin))*tmax1+(1.0-madt)*(1.0-nadt)*(tmax1-tmin)**2.0)
			containers(1) =		-containers(2)*tmin
			do it1 = 1,nt1
				t(it1) = 		omega*pi+((containers(1)-t(it1))/(t(it1)*containers(3)-containers(2)))
				t(nt-it1+1) = 	omega*2.0*pi+tmax-t(it1)				
			end do
		end if
	end if
	
!	generating r-grid
	if (adr .eq. 0) then
		do ir = 1,nr
			r(ir) = 			rmin+((rmax-rmin)/(nr-1))*(ir-1)
		end do
	else
		do ir = 1,nr
			r(ir) = 			(1.0/(nr-1))*(ir-1)
		end do
		containers(2) =			(madr*nadr*(rmax-rmin)-nadr*rmax)/((rmin-madr*(rmax-rmin))*rmax+madr*nadr*(rmax-rmin)**2.0)
		containers(3) =			((madr-nadr)*(rmax-rmin)-rmin)/((rmin-madr*(rmax-rmin))*rmax+madr*nadr*(rmax-rmin)**2.0)
		containers(1) =			-containers(2)*rmin
		do ir = 1,nr
			r(ir) = 			(containers(1)-r(ir))/(r(ir)*containers(3)-containers(2))
		end do
	end if
	
end subroutine meshgen
!-----------------------------------------------------------------------------------------------------------------------------------
