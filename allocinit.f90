!-----------------------------------------------------------------------------------------------------------------------------------
! matrix space allocation and initialization subroutine
!
! allocates the dimension to  matrices, sets values to required variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine allocinit

	use varinit
	implicit none

!	Allocating dimension to matrices
	allocate(pana				(0:3,1:nr))
	allocate(pstana				(0:3,1:nr))
	allocate(lana				(0:3,1:nr))
	allocate(pDLana				(0:3,1:nr))
	allocate(pvdWana			(0:3,1:nr))
	allocate(pSana				(0:3,1:nr))
	allocate(panad				(0:3,1:nr))
	allocate(pstanad			(0:3,1:nr))
	allocate(lanad				(0:3,1:nr))
	allocate(pDLanad			(0:3,1:nr))
	allocate(pvdWanad			(0:3,1:nr))
	allocate(pSanad				(0:3,1:nr))
	allocate(t					(1:nt))
	allocate(r					(1:nr))
	allocate(capH				(1:nr))
	allocate(td					(1:nt))
	allocate(rd					(1:nr))
	allocate(capHd				(1:nr))
	allocate(epsi1				(1:nt))
	allocate(kappa1				(1:nt))
	allocate(gamma1				(1:nt))
	allocate(eta1				(1:nt))
	allocate(coeffmom			(1:3,1:13))
	allocate(coeffnernst		(1:nspecies,1:10))
	allocate(coeffpoisson		(1:4))
	allocate(coeffre			(1:11))
	allocate(coeffsolidbc		(1:3,1:10))
	allocate(lana1				(0:3,1:nr))
	allocate(pDLana1			(0:3,1:nr))
	allocate(pvdWana1			(0:3,1:nr))
	allocate(pSana1				(0:3,1:nr))
	allocate(pSanaa				(0:3,1:nr))
	allocate(pSanab				(0:3,1:nr))
	allocate(lanalast			(0:3,1:nr))
	allocate(lanalast2			(0:3,1:nr))
	allocate(capMamp			(1:nt))
	allocate(capNdiv			(1:nt))
	allocate(regime				(1:nt))

!	Initialization of Matrices
	pana =			0.0
	pstana =		0.0
	lana =			0.0
	pDLana =		0.0
	pvdWana =		0.0
	pSana = 		0.0
	panad =			0.0
	pstanad =		0.0
	lanad =			0.0
	pDLanad =		0.0
	pvdWanad =		0.0
	pSanad = 		0.0
	t =				0.0
	r =				0.0
	capH = 			0.0
	td =			0.0
	rd =			0.0
	capHd = 		0.0
	epsi1 =			0.0
	gamma1 = 		0.0
	kappa1 = 		0.0
	eta1 = 			0.0
	coeffmom = 		0.0
	coeffnernst = 	0.0
	coeffpoisson = 	0.0
	coeffre = 		0.0
	coeffsolidbc = 	0.0
	lana1 = 		0.0
	pDLana1 = 		0.0
	pvdWana1 = 		0.0
	pSana1 = 		0.0
	pSanaa = 		0.0
	pSanab = 		0.0
	lanalast = 		0.0
	lanalast2 = 	0.0
	capMamp =		0.0
	capNdiv =		0.0
	regime =		0

! 	creating files for dumping output
	open(	unit=11,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lana.dat', &
			status='new',action='write')
	open(	unit=12,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pana.dat', &
			status='new',action='write')
	open(	unit=13,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstana.dat', &
			status='new',action='write')
	open(	unit=14,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLana.dat', &
			status='new',action='write')
	open(	unit=15,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWana.dat', &
			status='new',action='write')
	open(	unit=16,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSana.dat', &
			status='new',action='write')
	open(	unit=41,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat', &
			status='new',action='write')
	open(	unit=42,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat', &
			status='new',action='write')
	open(	unit=43,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstanad.dat', &
			status='new',action='write')
	open(	unit=44,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLanad.dat', &
			status='new',action='write')
	open(	unit=45,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWanad.dat', &
			status='new',action='write')
	open(	unit=46,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSanad.dat', &
			status='new',action='write')
	open(	unit=71,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/t.dat', &
			status='new',action='write')
	open(	unit=72,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/r.dat', &
			status='new',action='write')
	open(	unit=73,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/H.dat', &
			status='new',action='write')
	open(	unit=81,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/td.dat', &
			status='new',action='write')
	open(	unit=82,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.dat', &
			status='new',action='write')
	open(	unit=83,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/Hd.dat', &
			status='new',action='write')
	open(	unit=91,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.dat', &
			status='new',action='write')
	open(	unit=92,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/mumps.log', &
			status='new',action='write')
	open(	unit=93,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.dat', &
			status='new',action='write')
	open(	unit=94,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/main.log', &
			status='new',action='write')
	open(	unit=95,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/anawarn.log', &
			status='new',action='write')
	
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
	close(71)
	close(72)
	close(73)
	close(81)
	close(82)
	close(83)
	close(91)
	close(92)
	close(93)
	close(94)
	close(95)

end subroutine allocinit
!-----------------------------------------------------------------------------------------------------------------------------------
