!-----------------------------------------------------------------------------------------------------------------------------------
! parameter calculation subroutine
!
! calculates pertinent parameter values as required by the solution mechanism
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine paramcalc

	use varinit
	implicit none

!	calculating aspect ratios, elasticity co-efficients and Debye length
	epsi =		capD/capR
	alpha =		h0/capD
	delta =		capL/capR
	capK = 		sqrt((q**2.0*n0*(abs(valence(1)*valence(2)**2.0)+abs(valence(1)**2.0*valence(2))))/(relperm*permvac*kcapB*capT))
	kappa = 	(epsi/sqrt(epsi))
	eta = 		(kappa*delta)/epsi
	gamm = 		delta/sqrt(epsi)

!	calculating coefficient values of momentum equations
	do itemp = 1,2
		coeffmom(itemp,1) = 	(sqrt(epsi)**2.0*epsi**2.0*rho*omega*capR**2.0)/mu
		coeffmom(itemp,2) = 	(alpha*sqrt(epsi)**2.0*epsi**2.0*rho*omega*capR**2.0)/mu
		coeffmom(itemp,3) = 	(alpha*sqrt(epsi)**2.0*epsi**2.0*rho*omega*capR**2.0)/mu
		coeffmom(itemp,4) = 	(alpha*sqrt(epsi)**2.0*epsi**2.0*rho*omega*capR**2.0)/mu
		coeffmom(itemp,5) = 	(alpha*sqrt(epsi)**2.0*epsi**2.0*rho*omega*capR**2.0)/mu
		coeffmom(itemp,6) = 	1.0
		coeffmom(itemp,7) = 	epsi**2.0
		coeffmom(itemp,8) = 	epsi**2.0
		coeffmom(itemp,9) = 	sqrt(epsi)**2.0
		coeffmom(itemp,10) = 	epsi**2.0
		coeffmom(itemp,11) = 	epsi**2.0
		coeffmom(itemp,12) = 	(sqrt(epsi)*epsi*n0*zeta(1))/(mu*omega*alpha)
		coeffmom(itemp,13) = 	((sqrt(epsi)*epsi**2.0*capR*rho*0.0)/(mu*omega*alpha))
	end do
	coeffmom(3,1) = 	(sqrt(epsi)**2.0*epsi**4.0*rho*omega*capR**2.0)/mu
	coeffmom(3,2) = 	(alpha*sqrt(epsi)**2.0*epsi**4.0*rho*omega*capR**2.0)/mu
	coeffmom(3,3) = 	(alpha*sqrt(epsi)**2.0*epsi**4.0*rho*omega*capR**2.0)/mu
	coeffmom(3,4) = 	(alpha*sqrt(epsi)**2.0*epsi**4.0*rho*omega*capR**2.0)/mu
	coeffmom(3,5) = 	0.0
	coeffmom(3,6) = 	sqrt(epsi)**2.0
	coeffmom(3,7) = 	epsi**4.0
	coeffmom(3,8) = 	epsi**4.0
	coeffmom(3,9) = 	sqrt(epsi)**2.0*epsi**2.0
	coeffmom(3,10) = 	0.0
	coeffmom(3,11) = 	0.0
	coeffmom(3,12) = 	(sqrt(epsi)**2.0*epsi**2.0*n0*zeta(1))/(mu*omega*alpha)
	coeffmom(3,13) = 	((sqrt(epsi)**2.0*epsi**3.0*capR*rho*0.0)/(mu*omega*alpha))

!	calculating coefficient values of species conservation equation
	do ispecies = 1,nspecies
		coeffnernst(ispecies,1) = 	sqrt(epsi)**2.0*epsi**2.0*omega*capR**2.0
		coeffnernst(ispecies,2) = 	alpha*sqrt(epsi)**2.0*epsi**2.0*omega*capR**2.0
		coeffnernst(ispecies,3) = 	alpha*sqrt(epsi)**2.0*epsi**2.0*omega*capR**2.0
		coeffnernst(ispecies,4) = 	alpha*sqrt(epsi)**2.0*epsi**2.0*omega*capR**2.0
		coeffnernst(ispecies,5) = 	diffcoeff(ispecies)*epsi**2.0
		coeffnernst(ispecies,6) = 	diffcoeff(ispecies)*epsi**2.0
		coeffnernst(ispecies,7) = 	diffcoeff(ispecies)*sqrt(epsi)**2.0
		coeffnernst(ispecies,8) = 	((diffcoeff(ispecies)*q*valence(ispecies)*zeta(1))/(kcapB*capT))*epsi**2.0
		coeffnernst(ispecies,9) = 	((diffcoeff(ispecies)*q*valence(ispecies)*zeta(1))/(kcapB*capT))*epsi**2.0
		coeffnernst(ispecies,10) = 	((diffcoeff(ispecies)*q*valence(ispecies)*zeta(1))/(kcapB*capT))*sqrt(epsi)**2.0
	end do

!	calculating coefficients of poisson equation
 	coeffpoisson(1) = 	epsi**2.0
	coeffpoisson(2) = 	epsi**2.0
	coeffpoisson(3) = 	sqrt(epsi)**2.0
	coeffpoisson(4) = 	(eta**2.0*epsi**2.0*n0*capR**2.0*q)/(permvac*relperm*zeta(1))

!	calculating coefficients of reynolds equation
	coeffre(1) =		1.0
	coeffre(2) = 		eta/alpha
	coeffre(3) = 		(eta*gamm)/alpha
	coeffre(4) = 		(eta*gamm**2.0)/alpha
	coeffre(5) =		1.0
	coeffre(6) =		eta
	coeffre(7) = 		eta*gamm
	coeffre(8) =		eta*gamm**2.0
	coeffre(9) = 		1.0
	coeffre(10) = 		eta
	coeffre(11) = 		eta**2.0

!	calculating coeffients of solid b.c.s
	do itemp = 1,2
		coeffsolidbc(itemp,1) =		1.0
		coeffsolidbc(itemp,2) =		gamm
		coeffsolidbc(itemp,3) = 	gamm*eta*(epsi/delta)
		coeffsolidbc(itemp,4) = 	gamm**2.0*eta*(epsi/delta)
		coeffsolidbc(itemp,5) = 	gamm*eta*(epsi/delta)
		coeffsolidbc(itemp,6) = 	gamm*eta**2.0*(epsi/delta)
		coeffsolidbc(itemp,7) = 	gamm*eta**3.0*(epsi/delta)
		coeffsolidbc(itemp,8) = 	gamm*eta*((epsi**3.0)/delta)
		coeffsolidbc(itemp,9) =		epsi*sqrt(epsi)
		coeffsolidbc(itemp,10) =	(epsi**3.0)/sqrt(epsi)
	end do
	coeffsolidbc(3,1) =		1.0
	coeffsolidbc(3,2) =		gamm
	coeffsolidbc(3,3) = 	gamm*eta*(epsi/delta)
	coeffsolidbc(3,4) = 	gamm**2.0*eta*(epsi/delta)
	coeffsolidbc(3,5) = 	1.0
	coeffsolidbc(3,6) = 	eta
	coeffsolidbc(3,7) = 	eta**2.0
	coeffsolidbc(3,8) = 	gamm*eta*((epsi**2.0*sqrt(epsi))/delta)
	coeffsolidbc(3,9) =	 	gamm*eta*((epsi**4.0)/(sqrt(epsi)*delta))
	coeffsolidbc(3,10) =	2.0*epsi**2.0

end subroutine paramcalc
!-----------------------------------------------------------------------------------------------------------------------------------
