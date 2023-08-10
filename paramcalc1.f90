!-----------------------------------------------------------------------------------------------------------------------------------
! parameter calculation subroutine for rectified subroutine
!
! calculates pertinent parameter values as required by the rectified solution mechanism
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine paramcalc1

	use varinit
	implicit none

!	calculating aspect ratios, elasticity co-efficients and Debye length
	epsi1(it) =		epsi*(1.0+alpha*cos(t(it)))
	kappa1(it) = 	((mu*omega*alpha*epsi)/(epsi1(it)**2.0*capE))*(((1.0+nu)*(1.0-2*nu))/(1.0-nu))
	eta1(it) = 		(kappa1(it)*delta)/epsi1(it)
	gamma1(it) = 	delta/sqrt(epsi1(it))

end subroutine paramcalc1
