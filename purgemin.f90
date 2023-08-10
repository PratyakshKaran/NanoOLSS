!-----------------------------------------------------------------------------------------------------------------------------------
! purging excess output subroutine
!
! purges non-required output data if needed for space
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine purgemin

	use varinit
    implicit none
!	integer :: iline
	
	call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.dat '// &
								'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.txt',exitstat=iline)
	call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.dat '// &
								'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.txt',exitstat=iline)
	if (minout .eq. 1) then
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.dat '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.txt',exitstat=iline)
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.txt',exitstat=iline)
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.txt',exitstat=iline)
	end if
	call execute_command_line(		'rm '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/*.dat',exitstat=iline)
	call execute_command_line(		'rm -r '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)// &
									'/jacobian',exitstat=iline)
	call execute_command_line(		'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.txt '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.dat',exitstat=iline)
	call execute_command_line(		'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.txt '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.dat',exitstat=iline)
	if (minout .eq. 1) then
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.txt '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/rd.dat',exitstat=iline)
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.txt '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/panad.dat',exitstat=iline)
		call execute_command_line(	'mv '//'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.txt '// &
									'../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/lanad.dat',exitstat=iline)
	end if

end subroutine purgemin
