!-----------------------------------------------------------------------------------------------------------------------------------
! plot data generator subroutine from analytical solution
!
! assorts data to be plotted
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine plotana

	use varinit
	implicit none
	
	integer :: iph
	double precision :: dtmpar(1:20)
	integer :: itmpar(1:20)
	integer :: tph(1:4)
	double precision, allocatable ::	capF(:)
	double precision, allocatable :: 	lcl(:), 		pcl(:), 	pstcl(:), 		pDLcl(:), 		pvdWcl(:), 		pScl(:)	
	double precision, allocatable :: 	lcl0(:),		pcl0(:), 	pstcl0(:), 		pDLcl0(:), 		pvdWcl0(:), 	pScl0(:)
	double precision, allocatable ::	pinst(:)
	double precision, allocatable :: 	lph(:,:), 		pph(:,:), 	pstph(:,:), 	pDLph(:,:), 	pvdWph(:,:), 	pSph(:,:)
	double precision, allocatable :: 	lph0(:,:), 		pph0(:,:), 	pstph0(:,:), 	pDLph0(:,:),	pvdWph0(:,:), 	pSph0(:,:)
	
	allocate(capF(1:nt))
	allocate(lcl(1:nt))
	allocate(pcl(1:nt))
	allocate(pstcl(1:nt))
	allocate(pDLcl(1:nt))
	allocate(pvdWcl(1:nt))
	allocate(pScl(1:nt))
	allocate(lcl0(1:nt))
	allocate(pcl0(1:nt))
	allocate(pstcl0(1:nt))
	allocate(pDLcl0(1:nt))
	allocate(pvdWcl0(1:nt))
	allocate(pScl0(1:nt))
	allocate(pinst(1:nr))
	allocate(lph(1:4,1:nr))
	allocate(pph(1:4,1:nr))
	allocate(pstph(1:4,1:nr))
	allocate(pDLph(1:4,1:nr))
	allocate(pvdWph(1:4,1:nr))
	allocate(pSph(1:4,1:nr))
	allocate(lph0(1:4,1:nr))
	allocate(pph0(1:4,1:nr))
	allocate(pstph0(1:4,1:nr))
	allocate(pDLph0(1:4,1:nr))
	allocate(pvdWph0(1:4,1:nr))
	allocate(pSph0(1:4,1:nr))
	
	tph(1) = 1
	tph(2) = tph(1)
	do while(t(tph(2)) .le. (t(nt)/4.0))
		tph(2) = tph(2)+1
	end do
	tph(2) = tph(2)-1
	tph(3) = (nt-1)/2+1
	tph(4) = tph(3)
	do while(t(tph(4)) .le. ((3.0*t(nt))/4.0))
		tph(4) = tph(4)+1
	end do
	tph(4) = tph(4)-1
	open(unit=11,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lana.dat',status='old',action='read')
	do it = 1,nt
		read(11,*)	lcl(it)
		read(11,*)	lcl0(it)
		read(11,*)
		read(11,*)
	end do
	close(11)
	open(unit=11,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(11,*)	(lph(iph,ir),ir=1,nr)
			read(11,*)	(lph0(iph,ir),ir=1,nr)
			read(11,*)
			read(11,*)
			iph = iph + 1
		else
			read(11,*)
			read(11,*)
			read(11,*)
			read(11,*)
		end if
	end do
	close(11)
	open(unit=13,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstana.dat',status='old',action='read')
	do it = 1,nt
		read(13,*)	pstcl(it)
		read(13,*)	pstcl0(it)
		read(13,*)
		read(13,*)
	end do
	close(13)
	open(unit=13,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(13,*)	(pstph(iph,ir),ir=1,nr)
			read(13,*)	(pstph0(iph,ir),ir=1,nr)
			read(13,*)
			read(13,*)
			iph = iph + 1
		else
			read(13,*)
			read(13,*)
			read(13,*)
			read(13,*)
		end if
	end do
	close(13)
	open(unit=14,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLana.dat',status='old',action='read')
	do it = 1,nt
		read(14,*)	pDLcl(it)
		read(14,*)	pDLcl0(it)
		read(14,*)
		read(14,*)
	end do
	close(14)
	open(unit=14,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(14,*)	(pDLph(iph,ir),ir=1,nr)
			read(14,*)	(pDLph0(iph,ir),ir=1,nr)
			read(14,*)
			read(14,*)
			iph = iph + 1
		else
			read(14,*)
			read(14,*)
			read(14,*)
			read(14,*)
		end if
	end do
	close(14)
	open(unit=15,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWana.dat',status='old',action='read')
	do it = 1,nt
		read(15,*)	pvdWcl(it)
		read(15,*)	pvdWcl0(it)
		read(15,*)
		read(15,*)
	end do
	close(15)
	open(unit=15,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(15,*)	(pvdWph(iph,ir),ir=1,nr)
			read(15,*)	(pvdWph0(iph,ir),ir=1,nr)
			read(15,*)
			read(15,*)
			iph = iph + 1
		else
			read(15,*)
			read(15,*)
			read(15,*)
			read(15,*)
		end if
	end do
	close(15)
	open(unit=16,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSana.dat',status='old',action='read')
	do it = 1,nt
		read(16,*)	pScl(it)
		read(16,*)	pScl0(it)
		read(16,*)
		read(16,*)
	end do
	close(16)
	open(unit=16,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(16,*)	(pSph(iph,ir),ir=1,nr)
			read(16,*)	(pSph0(iph,ir),ir=1,nr)
			read(16,*)
			read(16,*)
			iph = iph + 1
		else
			read(16,*)
			read(16,*)
			read(16,*)
			read(16,*)
		end if
	end do
	close(16)
	open(unit=12,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pana.dat',status='old',action='read')
	do it = 1,nt
		read(12,*)	pcl(it)
		read(12,*)	pcl0(it)
		read(12,*)
		read(12,*)
	end do
	close(12)
	open(unit=12,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pana.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(12,*)	(pph(iph,ir),ir=1,nr)
			read(12,*)	(pph0(iph,ir),ir=1,nr)
			read(12,*)
			read(12,*)
			iph = iph + 1
		else
			read(12,*)
			read(12,*)
			read(12,*)
			read(12,*)
		end if
	end do
	close(12)
	open(unit=12,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pana.dat',status='old',action='read')
	do it = 1,nt
		read(12,*)	(pinst(ir),ir=1,nr)
		do iline = 2,4
			read(12,*)
		end do
		capF(it) = 0.0
		do ir = 1,nr-1
			capF(it) = capF(it) + 0.25*(pinst(ir)+pinst(ir+1))*(r(ir)+r(ir+1))*(r(ir+1)-r(ir))
		end do
	end do
	close(12)
	capF = capF*2.0*pi
	
	call execute_command_line('rm -r '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot',exitstat=iline)
	call execute_command_line('mkdir -p '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot',exitstat=iline)
	call execute_command_line('cp '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/Hana.dat '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot/Hana.dat',exitstat=iline)
	call execute_command_line('cp '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/r.dat '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot/r.dat',exitstat=iline)
	call execute_command_line('cp '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/t.dat '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot/t.dat',exitstat=iline)
										
	open(unit=111,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lcl.dat', &
											status='new',action='write')
	open(unit=112,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pcl.dat', &
											status='new',action='write')
	open(unit=113,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstcl.dat', &
											status='new',action='write')
	open(unit=114,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLcl.dat', &
											status='new',action='write')
	open(unit=115,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWcl.dat', &
											status='new',action='write')
	open(unit=116,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pScl.dat', &
											status='new',action='write')
	open(unit=121,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lcl0.dat', &
											status='new',action='write')
	open(unit=122,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pcl0.dat', &
											status='new',action='write')
	open(unit=123,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstcl0.dat', &
											status='new',action='write')
	open(unit=124,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLcl0.dat', &
											status='new',action='write')
	open(unit=125,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWcl0.dat', &
											status='new',action='write')
	open(unit=126,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pScl0.dat', &
											status='new',action='write')
	open(unit=161,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lph.dat', &
											status='new',action='write')
	open(unit=162,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pph.dat', &
											status='new',action='write')
	open(unit=163,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstph.dat', &
											status='new',action='write')
	open(unit=164,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLph.dat', &
											status='new',action='write')
	open(unit=165,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWph.dat', &
											status='new',action='write')
	open(unit=166,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pSph.dat', &
											status='new',action='write')
	open(unit=171,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lph0.dat', &
											status='new',action='write')
	open(unit=172,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pph0.dat', &
											status='new',action='write')
	open(unit=173,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstph0.dat', &
											status='new',action='write')
	open(unit=174,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLph0.dat', &
											status='new',action='write')
	open(unit=175,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWph0.dat', &
											status='new',action='write')
	open(unit=176,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pSph0.dat', &
											status='new',action='write')
	open(unit=211,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/tph.dat', &
											status='new',action='write')
	open(unit=221,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/F.dat', &
											status='new',action='write')

	write(111,*)		(lcl(it),it=1,nt)
	write(112,*)		(pcl(it),it=1,nt)
	write(113,*)		(pstcl(it),it=1,nt)
	write(114,*)		(pDLcl(it),it=1,nt)
	write(115,*)		(pvdWcl(it),it=1,nt)
	write(116,*)		(pScl(it),it=1,nt)
	write(121,*)		(lcl0(it),it=1,nt)
	write(122,*)		(pcl0(it),it=1,nt)
	write(123,*)		(pstcl0(it),it=1,nt)
	write(124,*)		(pDLcl0(it),it=1,nt)
	write(125,*)		(pvdWcl0(it),it=1,nt)
	write(126,*)		(pScl0(it),it=1,nt)
	do iph = 1,4
		write(161,*)	(lph(iph,ir),ir=1,nr)
		write(162,*)	(pph(iph,ir),ir=1,nr)
		write(163,*)	(pstph(iph,ir),ir=1,nr)
		write(164,*)	(pDLph(iph,ir),ir=1,nr)
		write(165,*)	(pvdWph(iph,ir),ir=1,nr)
		write(166,*)	(pSph(iph,ir),ir=1,nr)
		write(171,*)	(lph0(iph,ir),ir=1,nr)
		write(172,*)	(pph0(iph,ir),ir=1,nr)
		write(173,*)	(pstph0(iph,ir),ir=1,nr)
		write(174,*)	(pDlph0(iph,ir),ir=1,nr)
		write(175,*)	(pvdWph0(iph,ir),ir=1,nr)
		write(176,*)	(pSph0(iph,ir),ir=1,nr)
	end do
	write(211,*) 		(tph(iph), iph=1,4)
	write(221,*)		(capF(it),it=1,nt)

	close(111)
	close(112)
	close(113)
	close(114)
	close(115)
	close(116)
	close(121)
	close(122)
	close(123)
	close(124)
	close(125)
	close(126)
	close(161)
	close(162)
	close(163)
	close(164)
	close(165)
	close(166)
	close(171)
	close(172)
	close(173)
	close(174)
	close(175)
	close(176)
	close(211)
	close(221)

end subroutine plotana
