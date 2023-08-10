!-----------------------------------------------------------------------------------------------------------------------------------
! output write subroutine analytical
!
! writes individual time step solution to output files for analytical solution
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine anaout

    use varinit
    implicit none

	open(unit=11,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lana.dat', &
                status='old',position="append",action='write')
	open(unit=12,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pana.dat', &
                status='old',position="append",action='write')
	open(unit=13,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstana.dat', &
                status='old',position="append",action='write')
	open(unit=14,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLana.dat', &
                status='old',position="append",action='write')
	open(unit=15,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWana.dat', &
                status='old',position="append",action='write')
	open(unit=16,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSana.dat', &
                status='old',position="append",action='write')
	open(unit=41,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat', &
                status='old',position="append",action='write')
	open(unit=42,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat', &
                status='old',position="append",action='write')
	open(unit=43,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstanad.dat', &
                status='old',position="append",action='write')
	open(unit=44,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLanad.dat', &
                status='old',position="append",action='write')
	open(unit=45,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWanad.dat', &
                status='old',position="append",action='write')
	open(unit=46,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSanad.dat', &
                status='old',position="append",action='write')
	open(unit=72,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.dat', &
							status='old',position="append",action='write')
    open(unit=73,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/Hd.dat', &
                            status='old',position="append",action='write')

    write(11,*)  	(lana(3,ir), ir=1,nr)
    write(11,*)  	(lana(0,ir), ir=1,nr)
    write(11,*)  	(lana(1,ir), ir=1,nr)
    write(11,*)  	(lana(2,ir), ir=1,nr)
    write(12,*)  	(pana(3,ir), ir=1,nr)
    write(12,*)  	(pana(0,ir), ir=1,nr)
    write(12,*)  	(pana(1,ir), ir=1,nr)
    write(12,*)  	(pana(2,ir), ir=1,nr)
    write(13,*)  	(pstana(3,ir), ir=1,nr)
    write(13,*)  	(pstana(0,ir), ir=1,nr)
    write(13,*)  	(pstana(1,ir), ir=1,nr)
    write(13,*)  	(pstana(2,ir), ir=1,nr)
    write(14,*) 	(pDLana(3,ir), ir=1,nr)
    write(14,*) 	(pDLana(0,ir), ir=1,nr)
    write(14,*) 	(pDLana(1,ir), ir=1,nr)
    write(14,*) 	(pDLana(2,ir), ir=1,nr)
    write(15,*) 	(pvdWana(3,ir), ir=1,nr)
    write(15,*) 	(pvdWana(0,ir), ir=1,nr)
    write(15,*) 	(pvdWana(1,ir), ir=1,nr)
    write(15,*) 	(pvdWana(2,ir), ir=1,nr)
    write(16,*) 	(pSana(3,ir), ir=1,nr)
    write(16,*) 	(pSana(0,ir), ir=1,nr)
    write(16,*) 	(pSana(1,ir), ir=1,nr)
    write(16,*) 	(pSana(2,ir), ir=1,nr)
    write(41,*)  	(lanad(3,ir), ir=1,nr)
    write(41,*)  	(lanad(0,ir), ir=1,nr)
    write(41,*)  	(lanad(1,ir), ir=1,nr)
    write(41,*)  	(lanad(2,ir), ir=1,nr)
    write(42,*)  	(panad(3,ir), ir=1,nr)
    write(42,*)  	(panad(0,ir), ir=1,nr)
    write(42,*)  	(panad(1,ir), ir=1,nr)
    write(42,*)  	(panad(2,ir), ir=1,nr)
    write(43,*)  	(pstanad(3,ir), ir=1,nr)
    write(43,*)  	(pstanad(0,ir), ir=1,nr)
    write(43,*)  	(pstanad(1,ir), ir=1,nr)
    write(43,*)  	(pstanad(2,ir), ir=1,nr)
    write(44,*) 	(pDLanad(3,ir), ir=1,nr)
    write(44,*) 	(pDLanad(0,ir), ir=1,nr)
    write(44,*) 	(pDLanad(1,ir), ir=1,nr)
    write(44,*) 	(pDLanad(2,ir), ir=1,nr)
    write(45,*) 	(pvdWanad(3,ir), ir=1,nr)
    write(45,*) 	(pvdWanad(0,ir), ir=1,nr)
    write(45,*) 	(pvdWanad(1,ir), ir=1,nr)
    write(45,*) 	(pvdWanad(2,ir), ir=1,nr)
    write(46,*) 	(pSanad(3,ir), ir=1,nr)
    write(46,*) 	(pSanad(0,ir), ir=1,nr)
    write(46,*) 	(pSanad(1,ir), ir=1,nr)
    write(46,*) 	(pSanad(2,ir), ir=1,nr)
	write(72,*)		(rd(ir), ir=1,nr)
    write(73,*) 	(capHd(ir), ir=1,nr)
	
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(41)
    close(42)
    close(43)
    close(44)
    close(45)
    close(46)
    close(72)
    close(73)
    close(92)
	close(95)

end subroutine anaout
!-----------------------------------------------------------------------------------------------------------------------------------
