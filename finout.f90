!-----------------------------------------------------------------------------------------------------------------------------------
! final subroutine
!
! writes all remaining pertinent variables to output files (i.e. those not written for each time-step) and final subroutine for the
! program
!-----------------------------------------------------------------------------------------------------------------------------------!-----------------------------------------------------------------------------------------------------------------------------------
subroutine finout

    use varinit
    implicit none

    open(unit=1,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param.dat',status='old', &
			                        position='append', action='write')
    open(unit=5001,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/param1.dat',status='old', &
			                        position='append', action='write')
    open(unit=1004,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/r.dat',status='old', &
                                    position='append', action='write')
    open(unit=50,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/t.dat',status='old', &
                                    position='append', action='write')
    open(unit=51,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/td.dat',status='old', &
                                    position='append', action='write')
    open(unit=1015,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/H.dat', status='old', &
    								position="append",action='write')

    write(1,*)  	'solvent: '
    write(1,*)  	'rho, mu, rho_inf, sigma'
    write(1,*)  	rho, mu, rho_inf, sigma
    write(1,*)  	'electrolyte: '
    write(1,*)  	'nspecies, n0, (diffcoeff(ispecies), ispecies=1,nspecies), ', &
					'(valence(ispecies), ispecies=1,nspecies), relperm'
    write(1,*)  	nspecies, n0, (diffcoeff(ispecies), ispecies=1,nspecies), &
            		(valence(ispecies), ispecies=1,nspecies), relperm
    write(1,*)  	'solid: '
    write(1,*)  	'capE, nu, capAsfw, capLambda1, capGamma11, capGamma12, phi1, capLambda2, capGamma21, capGamma22, phi2,'//&
					'zeta(1), zeta(2)'
    write(1,*)  	capE, nu, capAsfw, capLambda1, capGamma11, capGamma12, phi1, capLambda2, capGamma21, capGamma22, phi2, &
					zeta(1), zeta(2)
	write(1,*) 		'univ:'
	write(1,*)		'q, kcapB, permvac'
	write(1,*)		q, kcapB, permvac
    write(1,*)  	'geometry: '
    write(1,*)  	'h0, capD, capL, capR, omega, capT'
    write(1,*)  	h0, capD, capL, capR, omega, capT
    write(1,*)  	'grid: '
	write(1,*)		'adr, ady, adt, madr, nadr, mady, nady, madt, nadt'
	write(1,*)		adr, ady, adt, madr, nadr, mady, nady, madt, nadt
	write(1,*)		'nophd, phditer, phdpatch'
	write(1,*)		nophd, phditer, phdpatch
	write(1,*)		'phdaxi, profaxi, molaxi'
	write(1,*)		phdaxi, profaxi, molaxi
	write(1,*)		'schemenum, schemeana, numguess, anaguess, revtime, toladapt'
	write(1,*)		schemenum, schemeana, numguess, anaguess, revtime, toladapt
	write(1,*)		'soln, minout, termout, sanity'
	write(1,*)		soln, minout, termout, sanity
	write(1,*)		'tolnod, tolphd, thresamp, thresdiv1, thresdiv2, dl, bsecshrink'
	write(1,*)		tolnod, tolphd, thresamp, thresdiv1, thresdiv2, dl, bsecshrink
	write(1,*)		'relaxnod, relaxphd'
	write(1,*)		relaxnod, relaxphd
	write(1,*)		'nr, ny, nz, nt'
	write(1,*)		nr, ny, nz, nt
	write(1,*)		'rmin, rmax, ymin, ymax, tmin, tmax	'
	write(1,*)		rmin, rmax, ymin, ymax, tmin, tmax
	write(1,*)		'================================'
    write(1,*)  	'probtype: '
    write(1,*)  	probtype
    write(1,*)  	'parameters: '
    write(1,*)  	'epsi, alpha, delta, capK, kappa, eta, gamm'
    write(1,*)  	epsi, alpha, delta, capK, kappa, eta, gamm
    write(1,*)  	'coeffmom'
    do itemp = 1,3
        write(1,*)  (coeffmom(itemp,itemp1), itemp1=1,13)
    end do
    write(1,*)  'nernst'
    do itemp = 1,nspecies
        write(1,*)  (coeffnernst(itemp,itemp1), itemp1=1,10)
    end do
    write(1,*)  'poisson'
    write(1,*)  	(coeffpoisson(itemp), itemp=1,4)
    write(1,*)  're'
    write(1,*)  	(coeffre(itemp), itemp=1,11)
    do itemp = 1,3
        write(1,*)  (coeffsolidbc(itemp,itemp1), itemp1=1,10)
    end do
	
	write(5001,*)	(epsi1(it), it=1,nt)
	write(5001,*)	(kappa1(it), it=1,nt)
	write(5001,*)	(eta1(it), it=1,nt)
	write(5001,*)	(gamma1(it), it=1,nt)
	write(5001,*)	(capMamp(it), it=1,nt)
	write(5001,*)	(capNdiv(it), it=1,nt)
	write(5001,*)	(regime(it), it=1,nt)
	
	write(1004,*)  	(r(ir), ir=1,nr)
    write(50,*)  	(t(it), it=1,nt)
    write(51,*)  	(td(it), it=1,nt)
    write(1015,*) 	(capH(ir), ir=1,nr)

    close(1)
    close(5001)
    close(1004)
    close(50)
    close(51)
	close(1015)
	
end subroutine finout
!-----------------------------------------------------------------------------------------------------------------------------------
