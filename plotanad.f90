!-----------------------------------------------------------------------------------------------------------------------------------
! plot data generator subroutine from analytical solution (dimensional values)
!
! assorts data to be plotted
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine plotanad

	use varinit
	implicit none
	
	integer :: iph
	double precision :: dtmpar(1:20)
	integer :: itmpar(1:20)
	integer :: tph(1:4)
	double precision, allocatable ::	capHcl(:), 		capF(:)
	double precision, allocatable :: 	lcl(:), 		pcl(:), 	pstcl(:), 		pDLcl(:), 		pvdWcl(:), 		pScl(:)	
	double precision, allocatable :: 	lcl0(:),		pcl0(:), 	pstcl0(:), 		pDLcl0(:), 		pvdWcl0(:), 	pScl0(:)
	double precision, allocatable ::	capHph(:,:),	pinst(:)
	double precision, allocatable :: 	lph(:,:), 		pph(:,:), 	pstph(:,:), 	pDLph(:,:), 	pvdWph(:,:), 	pSph(:,:)
	double precision, allocatable :: 	lph0(:,:), 		pph0(:,:), 	pstph0(:,:), 	pDLph0(:,:),	pvdWph0(:,:), 	pSph0(:,:)
	
	allocate(capF(1:nt))
	allocate(capHcl(1:nt))
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
	allocate(capHph(1:4,1:nr))
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

	open(unit=73,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/Hd.dat',status='old',action='read')
	do it = 1,nt
		read(73,*)		capHcl(it)
	end do
	close(73)
	open(unit=73,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/Hd.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(73,*)	(capHph(iph,ir),ir=1,nr)
			iph = iph + 1
		else
			read(73,*)
		end if
	end do
	close(73)
	open(unit=41,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat',status='old',action='read')
	do it = 1,nt
		read(41,*)		lcl(it)
		read(41,*)		lcl0(it)
		read(41,*)
		read(41,*)
	end do
	close(41)
	open(unit=41,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(41,*)	(lph(iph,ir),ir=1,nr)
			read(41,*)	(lph0(iph,ir),ir=1,nr)
			read(41,*)
			read(41,*)
			iph = iph + 1
		else
			read(41,*)
			read(41,*)
			read(41,*)
			read(41,*)
		end if
	end do
	close(41)
	open(unit=43,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstanad.dat',status='old',action='read')
	do it = 1,nt
		read(43,*)		pstcl(it)
		read(43,*)		pstcl0(it)
		read(43,*)
		read(43,*)
	end do
	close(43)
	open(unit=43,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pstanad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(43,*)	(pstph(iph,ir),ir=1,nr)
			read(43,*)	(pstph0(iph,ir),ir=1,nr)
			read(43,*)
			read(43,*)
			iph = iph + 1
		else
			read(43,*)
			read(43,*)
			read(43,*)
			read(43,*)
		end if
	end do
	close(43)
	open(unit=44,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLanad.dat',status='old',action='read')
	do it = 1,nt
		read(44,*)		pDLcl(it)
		read(44,*)		pDLcl0(it)
		read(44,*)
		read(44,*)
	end do
	close(44)
	open(unit=44,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pDLanad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(44,*)	(pDLph(iph,ir),ir=1,nr)
			read(44,*)	(pDLph0(iph,ir),ir=1,nr)
			read(44,*)
			read(44,*)
			iph = iph + 1
		else
			read(44,*)
			read(44,*)
			read(44,*)
			read(44,*)
		end if
	end do
	close(44)
	open(unit=45,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWanad.dat',status='old',action='read')
	do it = 1,nt
		read(45,*)		pvdWcl(it)
		read(45,*)		pvdWcl0(it)
		read(45,*)
		read(45,*)
	end do
	close(45)
	open(unit=45,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pvdWanad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(45,*)	(pvdWph(iph,ir),ir=1,nr)
			read(45,*)	(pvdWph0(iph,ir),ir=1,nr)
			read(45,*)
			read(45,*)
			iph = iph + 1
		else
			read(45,*)
			read(45,*)
			read(45,*)
			read(45,*)
		end if
	end do
	close(45)
	open(unit=46,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSanad.dat',status='old',action='read')
	do it = 1,nt
		read(46,*)		pScl(it)
		read(46,*)		pScl0(it)
		read(46,*)
		read(46,*)
	end do
	close(46)
	open(unit=46,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/pSanad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(46,*)	(pSph(iph,ir),ir=1,nr)
			read(46,*)	(pSph0(iph,ir),ir=1,nr)
			read(46,*)
			read(46,*)
			iph = iph + 1
		else
			read(46,*)
			read(46,*)
			read(46,*)
			read(46,*)
		end if
	end do
	close(46)
	open(unit=42,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat',status='old',action='read')
	do it = 1,nt
		read(42,*)		pcl(it)
		read(42,*)		pcl0(it)
		read(42,*)
		read(42,*)
	end do
	close(42)
	open(unit=42,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat',status='old',action='read')
	iph = 1
	do it = 1,nt
		if (it .eq. tph(iph)) then
			read(42,*)	(pph(iph,ir),ir=1,nr)
			read(42,*)	(pph0(iph,ir),ir=1,nr)
			read(42,*)
			read(42,*)
			iph = iph + 1
		else
			read(42,*)
			read(42,*)
			read(42,*)
			read(42,*)
		end if
	end do
	close(42)
	open(unit=42,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat',status='old',action='read')
	open(unit=72,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.dat',status='old',action='read')
	do it = 1,nt
		read(72,*)		(rd(ir),ir=1,nr)
		read(42,*)		(pinst(ir),ir=1,nr)
		do iline = 2,4
			read(42,*)
		end do
		capF(it) = 0.0
		do ir = 1,nr-1
			capF(it) = 	capF(it) + 0.25*(pinst(ir)+pinst(ir+1))*(rd(ir)+rd(ir+1))*(rd(ir+1)-rd(ir))
		end do
	end do
	close(42)
	close(72)
	capF = capF*2.0*pi
		
	call execute_command_line('mkdir -p '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot',exitstat=iline)
	call execute_command_line('cp '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/rd.dat '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot/rd.dat',exitstat=iline)
	call execute_command_line('cp '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/td.dat '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
										'/forplot/td.dat',exitstat=iline)

	open(unit=611,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lcld.dat', &
											status='new',action='write')
	open(unit=612,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pcld.dat', &
											status='new',action='write')
	open(unit=613,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstcld.dat', &
											status='new',action='write')
	open(unit=614,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLcld.dat', &
											status='new',action='write')
	open(unit=615,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWcld.dat', &
											status='new',action='write')
	open(unit=616,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pScld.dat', &
											status='new',action='write')
	open(unit=621,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lcl0d.dat', &
											status='new',action='write')
	open(unit=622,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pcl0d.dat', &
											status='new',action='write')
	open(unit=623,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstcl0d.dat', &
											status='new',action='write')
	open(unit=624,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLcl0d.dat', &
											status='new',action='write')
	open(unit=625,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWcl0d.dat', &
											status='new',action='write')
	open(unit=626,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pScl0d.dat', &
											status='new',action='write')
	open(unit=661,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lphd.dat', &
											status='new',action='write')
	open(unit=662,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pphd.dat', &
											status='new',action='write')
	open(unit=663,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstphd.dat', &
											status='new',action='write')
	open(unit=664,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLphd.dat', &
											status='new',action='write')
	open(unit=665,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWphd.dat', &
											status='new',action='write')
	open(unit=666,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pSphd.dat', &
											status='new',action='write')
	open(unit=671,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/lph0d.dat', &
											status='new',action='write')
	open(unit=672,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pph0d.dat', &
											status='new',action='write')
	open(unit=673,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pstph0d.dat', &
											status='new',action='write')
	open(unit=674,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pDLph0d.dat', &
											status='new',action='write')
	open(unit=675,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pvdWph0d.dat', &
											status='new',action='write')
	open(unit=676,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/pSph0d.dat', &
											status='new',action='write')
	open(unit=721,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/Fd.dat', &
											status='new',action='write')
	open(unit=761,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/Hcld.dat', &
											status='new',action='write')
	open(unit=771,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/forplot/Hphd.dat', &
											status='new',action='write')



	write(611,*)		(lcl(it),it=1,nt)
	write(612,*)		(pcl(it),it=1,nt)
	write(613,*)		(pstcl(it),it=1,nt)
	write(614,*)		(pDLcl(it),it=1,nt)
	write(615,*)		(pvdWcl(it),it=1,nt)
	write(616,*)		(pScl(it),it=1,nt)
	write(621,*)		(lcl0(it),it=1,nt)
	write(622,*)		(pcl0(it),it=1,nt)
	write(623,*)		(pstcl0(it),it=1,nt)
	write(624,*)		(pDLcl0(it),it=1,nt)
	write(625,*)		(pvdWcl0(it),it=1,nt)
	write(626,*)		(pScl0(it),it=1,nt)	
	do iph = 1,4
		write(661,*)	(lph(iph,ir),ir=1,nr)
		write(662,*)	(pph(iph,ir),ir=1,nr)
		write(663,*)	(pstph(iph,ir),ir=1,nr)
		write(664,*)	(pDLph(iph,ir),ir=1,nr)
		write(665,*)	(pvdWph(iph,ir),ir=1,nr)
		write(666,*)	(pSph(iph,ir),ir=1,nr)
		write(671,*)	(lph0(iph,ir),ir=1,nr)
		write(672,*)	(pph0(iph,ir),ir=1,nr)
		write(673,*)	(pstph0(iph,ir),ir=1,nr)
		write(674,*)	(pDlph0(iph,ir),ir=1,nr)
		write(675,*)	(pvdWph0(iph,ir),ir=1,nr)
		write(676,*)	(pSph0(iph,ir),ir=1,nr)
		write(771,*)	(capHph(iph,ir),ir=1,nr)
	end do
	write(721,*)		(capF(it),it=1,nt)
	write(761,*)		(capHcl(it),it=1,nt)
	
	close(611)
	close(612)
	close(613)
	close(614)
	close(615)
	close(616)
	close(621)
	close(622)
	close(623)
	close(624)
	close(625)
	close(626)
	close(661)
	close(662)
	close(663)
	close(664)
	close(665)
	close(666)
	close(671)
	close(672)
	close(673)
	close(674)
	close(675)
	close(676)
	close(771)
	close(721)
	close(761)

end subroutine plotanad
