!-----------------------------------------------------------------------------------------------------------------------------------
! time step evaluation subroutine for analytical solution
!
! obtains solution for a particular time step as per the scheme specified (Newton-Raphson or Source-term linearization)
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine anainst

!   including required modules
    use varinit
    use derivative
    implicit none

!   including the required header files
    include '/usr/local/openmpi/include/mpif.h'
    include '/usr/local/mumps/include/dmumps_struc.h'
		
!   declaring variables for mumps
	type (dmumps_struc) :: mumps_par
	integer :: ierr

	if (termout .eq. 1) then
		write(*,*) 		"analytical solution, time step ",it
	end if	
	
   open(	unit=92,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/mumps.log', &
         	status='old',position='append',action='write')
   open(	unit=95,file='../../../../out/dlss_num/'//datestamp//'/'//trim(siminstance)//'/anawarn.log', &
        	status='old',position='append',action='write')
	
!	obtaining solution when molecular forces dominate over hydrodynamic forces
!	taken to be separation that is given by the	grid parameter thres times the particle size in solvation force formulation
	capMamp(it) =	1.0+abs((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
					((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
					exp(-epsi1(it)*capK*capR)) + &
					abs(-(capAsfw/(6.0*pi*epsi1(it)*epsi*mu*omega*alpha*capR**3.0))) + &
					abs(-((epsi1(it)**2.0*capLambda1*rho_inf*kcapB*capT)/(mu*omega*alpha*epsi))* &
					exp(-((epsi1(it)*capR)/(capGamma11*sigma)))) + &
					abs(-((epsi1(it)**2.0*capLambda2*rho_inf*kcapB*capT)/(mu*omega*alpha*epsi))* &
					exp(-((epsi1(it)*capR)/(capGamma21*sigma))))				
	capNdiv(it) =	max(max(abs(eta1(it)*epsi1(it)*capK*capR*capMamp(it)),abs(3.0*eta1(it)*capMamp(it))), &
					max(abs(((eta1(it)*epsi1(it)*capR*capMamp(it))/sigma)*((1.0/capGamma11)+((2.0*pi)/capGamma12)* &
					tan(((2.0*pi*epsi1(it)*capR*(1.0+eta1(it)*capMamp(it)))/(capGamma12*sigma))+phi1))), &
					abs(((eta1(it)*epsi1(it)*capR*capMamp(it))/sigma)*((1.0/capGamma21)+((2.0*pi)/capGamma22)* &
					tan(((2.0*pi*epsi1(it)*capR*(1.0+eta1(it)*capMamp(it)))/(capGamma22*sigma))+phi2)))))
	
	if ((eta1(it)*capMamp(it) > thresamp) .or. ((eta1(it)*capMamp(it) > thresdiv1) .and. (capNdiv(it) > thresdiv2))) then
		regime(it) = 	1
	else
		regime(it) = 	0
	end if
	
	if (regime(it) .eq. 1) then
!		setting guess value as per condition
		if (anaguess .eq. 0) then
			lana = 		0.0
			pstana =	0.0
			pDLana = 	0.0
			pSana =		0.0
			pvdWana =	0.0
			pana =		0.0
		elseif (anaguess .eq. 1) then
			if (it .ne. 1) then
				lana = 		(kappa1(it-1)/kappa1(it))*lana
				pstana =	((epsi1(it)**2.0*kappa1(it))/(epsi1(it-1)**2.0*kappa1(it)))*pstana
				pDLana =	((epsi1(it)**2.0*kappa1(it))/(epsi1(it-1)**2.0*kappa1(it)))*pDLana
				pvdWana =	((epsi1(it)**2.0*kappa1(it))/(epsi1(it-1)**2.0*kappa1(it)))*pvdWana
				pSana =		((epsi1(it)**2.0*kappa1(it))/(epsi1(it-1)**2.0*kappa1(it)))*pSana
				pana =		((epsi1(it)**2.0*kappa1(it))/(epsi1(it-1)**2.0*kappa1(it)))*pana
			end if
		end if
		Fphd = 1.0e10
		iterphd = 0
!		iterating for enveloping solution of hd pressure
		do while ((Fphd) .ge. (tolphd))
			iterphd = iterphd+1
			do ir = 1,nr
!				iterating for node-wise solution of deflection
				Fnod = 1.0e10
				iternod = 0
				if (toladapt .eq. 1) then
					containers(1) = tolnod/capMamp(it)
				else
					containers(1) = tolnod
				end if		
				do while ((abs(Fnod)) .ge. (containers(1)))
					Jnodprev = Jnod
					iternod = iternod+1
					pDLana(0,ir) = 		(((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
										((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
										exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lana(0,ir)))
					pvdWana(0,ir) = 	-capAsfw/ &
										(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lana(0,ir))**3.0)
					pSanaa(0,ir) =   	-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
										exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma11*sigma))* &
										cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma12*sigma))+phi1)
					pSanab(0,ir) =   	-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
										exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma21*sigma))* &
										cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma22*sigma))+phi2)
					pSana(0,ir) =		pSanaa(0,ir)+pSanab(0,ir)
!					collating all pressure components to get total for leading order
					pana(0,ir) = 		pstana(0,ir)+pDLana(0,ir)+pvdWana(0,ir)+pSana(0,ir)
					if (((schemeana .eq. 0) .or. (schemeana .eq. 1)) .or. (schemeana .eq. 2)) then
						Fnod =				lana(0,ir)- pana(0,ir)* &
											(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))
!						solving using NR, NR-dir, or BF by putting J (i.e. F_prime) accordingly
						if (schemeana .ne. 2) then
							Jnod =			1.0-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))* &
											( (-(eta1(it)*epsi1(it)*capK*capR)*pDLana(0,ir)) + &
											(-((3.0*eta1(it))/(capH(ir)+eta1(it)*lana(0,ir)))*pvdWana(0,ir)) + &
											(-(((eta1(it)*epsi1(it)*capR)/sigma)*((1.0/capGamma11)+((2.0*pi)/capGamma12)* &
											tan(((2.0*pi*epsi1(it)*capR* &
											(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma12*sigma))+phi1)))* &
											pSanaa(0,ir)) + &
											(-(((eta1(it)*epsi1(it)*capR)/sigma)*((1.0/capGamma21)+((2.0*pi)/capGamma22)* &
											tan(((2.0*pi*epsi1(it)*capR* &
											(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma22*sigma))+phi2)))* &
											pSanab(0,ir)) )
							if ((schemeana .eq. 0) .and. (abs(Jnod) .ge. 1.0e6)) then
								do while (abs(Jnod) .ge. 1.0e4)
									Jnod = 		Jnod/10.0
								end do							
							end if
							if (schemeana .eq. 1) then
								Jnod =		Jnod/abs(Jnod)
							end if					
						else
							Jnod =			1.0
						end if
						lana(0,ir) =	lana(0,ir) - relaxnod*(Fnod/Jnod)
						if (termout .eq. 1) then
							write(*,*)	"time						",			t(it),	"	time step	",		it	
							write(*,*)	"root-find error				",		Fnod,	"	node no		",		ir, &
										"	iteration no	",					iternod
							write(*,*)	"hydrodynamic error				",		Fphd,	"	iteration no	",	iterphd
							write(*,*)	"amplification switch				",	(eta1(it)*capMamp(it) > thresamp)
							write(*,*)	"divergence precheck				",	(capMamp(it) > thresdiv1)
							write(*,*)	"divergence switch				",		(capNdiv(it) > thresdiv2)
						else
							write(95,*)	"time						",			t(it),	"	time step	",		it	
							write(95,*)	"root-find error				",		Fnod,	"	node no		",		ir, &
										"	iteration no	",					iternod
							write(95,*)	"hydrodynamic error				",		Fphd,	"	iteration no	",	iterphd
							write(95,*)	"amplification switch				",	(eta1(it)*capMamp(it) > thresamp)
							write(95,*)	"divergence precheck				",	(capMamp(it) > thresdiv1)
							write(95,*)	"divergence switch				",		(capNdiv(it) > thresdiv2)

						end if
					else
						dlanak = 			dl
						dlanakcapL = 		dl
						dlanakcapR =		dl
						lanak = 			lana(0,ir)
						panak = 			pana(0,ir)
						Fnod =				lanak-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panak		
						Fnodk =				Fnod
!						solving by proceeding in the direction of pressure (very prone to divergence)
						if (schemeana .eq. 3) then
							do while (Fnodk*Fnodkprev .ge. 0.0)
								Fnodkprev =		Fnodk
								panakprev = 	panak
								lanakprev =		lanak
								if (pana(0,ir) .ge. 0) then
									lanak = 	lanak+dlanak
									panak =		((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanak))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lanak)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma12*sigma))+phi1))+&
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma22*sigma))+phi2))+&
												pstana(0,ir)
									if (panak .gt. pana(0,ir)) then
										lanak = lanakprev-dlanak
									end if
								else
									lanak = 	lanak-dlanak
									panak =		((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanak))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lanak)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma12*sigma))+phi1))+&
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/ &
												(capGamma22*sigma))+phi2))+pstana(0,ir)					
								end if
								Fnodk = 		lanak-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
												(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panak
							end do
							Fnod =				Fnodk
							lana(0,ir) = 		lanak
							pana(0,ir) = 		panak
							pDLana(0,ir) = 		(((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lana(0,ir)))
							pvdWana(0,ir) = 	-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lana(0,ir))**3.0)
							pSanaa(0,ir) =   	-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/ &
												(capGamma12*sigma))+phi1)
							pSanab(0,ir) =   	-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/ &
												(capGamma22*sigma))+phi2)
							pSana(0,ir) =		pSanaa(0,ir)+pSanab(0,ir)
							pana(0,ir) = 		pstana(0,ir)+pDLana(0,ir)+pvdWana(0,ir)+pSana(0,ir)
!						solving by scope-incremental bisection
						else
							lanacapLr =			lanak
							lanacapLl =			lanak 
							lanacapRr =			lanak
							lanacapRl =			lanak
							FnodcapLl = 		1.0
							FnodcapLr = 		1.0
							FnodcapRl = 		1.0
							FnodcapRr = 		1.0
							do while (	((FnodcapRl*FnodcapRr .gt. 0.0) .and. (FnodcapLl*FnodcapLr .gt. 0.0)) .or. &
										((FnodcapRl*FnodcapRr .lt. 0.0) .and. (FnodcapLl*FnodcapLr .lt. 0.0)) )
								lanacapLr =		lanacapLl
								lanacapLl =		lanacapLl-dlanakcapL
								lanacapRl =		lanacapRr
								lanacapRr =		lanacapRr+dlanakcapR
								panacapLl = 	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanacapLl))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0* &
												(capH(ir)+eta1(it)*lanacapLl)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLl))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*& 
												(capH(ir)+eta1(it)*lanacapLl))/(capGamma12*sigma))+phi1))+ &
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLl))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLl))/ &
												(capGamma22*sigma))+phi2))+pstana(0,ir)
								panacapLr = 	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanacapLr))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0* &
												(capH(ir)+eta1(it)*lanacapLr)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLr))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*& 
												(capH(ir)+eta1(it)*lanacapLr))/(capGamma12*sigma))+phi1))+ &
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLr))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapLr))/ &
												(capGamma22*sigma))+phi2))+pstana(0,ir)
								panacapRl = 	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanacapRl))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0* &
												(capH(ir)+eta1(it)*lanacapRl)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRl))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*& 
												(capH(ir)+eta1(it)*lanacapRl))/(capGamma12*sigma))+phi1))+ &
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRl))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRl))/ &
												(capGamma22*sigma))+phi2))+pstana(0,ir)
								panacapRr = 	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
												((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
												exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanacapRr))) + &
												(-capAsfw/ &
												(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0* &
												(capH(ir)+eta1(it)*lanacapRr)**3.0)) +&
												(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRr))/(capGamma11*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*& 
												(capH(ir)+eta1(it)*lanacapRr))/(capGamma12*sigma))+phi1))+ &
												(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
												exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRr))/(capGamma21*sigma))* &
												cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanacapRr))/ &
												(capGamma22*sigma))+phi2))+pstana(0,ir)
								FnodcapLl = 	lanacapLl-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
												(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panacapLl
								FnodcapLr = 	lanacapLr-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
												(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panacapLr
								FnodcapRl = 	lanacapRl-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
												(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panacapRl
								FnodcapRr = 	lanacapRr-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
												(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panacapRr
								if (((FnodcapRl*FnodcapRr .lt. 0.0) .and. (FnodcapLl*FnodcapLr .lt. 0.0))) then
									lanacapLl =		lanacapLr
									lanacapLr =		lanacapLr+dlanakcapL
									lanacapRr =		lanacapRl
									lanacapRl =		lanacapRl-dlanakcapR
									dlanakcapL =	bsecshrink*dlanakcapL
									dlanakcapR =	bsecshrink*dlanakcapR
								end if
							end do
							if ((FnodcapRl*FnodcapRr .le. 0.0)) then
								lanakcapL = 	lanacapRl
								lanakcapR =		lanacapRr
							elseif ((FnodcapLl*FnodcapLr .le. 0.0)) then
								lanakcapL = 	lanacapLl
								lanakcapR =		lanacapLr
							else
								if (termout .eq. 0) then
									write(*,*) 'ERROR IN PROGRAM - VIEWSCOPING FOR ROOT IS NOT WORKING CORRECTLY'
									write(95,*) 'ERROR IN PROGRAM - VIEWSCOPING FOR ROOT IS NOT WORKING CORRECTLY'
								else
									write(95,*) 'ERROR IN PROGRAM - VIEWSCOPING FOR ROOT IS NOT WORKING CORRECTLY'
								end if
							end if
							do while (abs(Fnodk) .ge. containers(1))
								lanak = 	0.5*(lanakcapL+lanakcapR)
								panak =		((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
											((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
											exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanak))) + &
											(-capAsfw/ &
											(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lanak)**3.0)) +&
											(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma11*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma12*sigma))+phi1))+&
											(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/(capGamma21*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanak))/ &
											(capGamma22*sigma))+phi2))+pstana(0,ir)
								panakcapL =	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
											((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
											exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanakcapL))) + &
											(-capAsfw/ &
											(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lanakcapL)**3.0)) +&
											(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapL))/(capGamma11*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapL))/(capGamma12*sigma))+phi1))+&
											(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapL))/(capGamma21*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapL))/ &
											(capGamma22*sigma))+phi2))+pstana(0,ir)
								panakcapR =	((((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
											((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
											exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lanakcapR))) + &
											(-capAsfw/ &
											(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*(capH(ir)+eta1(it)*lanakcapR)**3.0)) +&
											(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapR))/(capGamma11*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapR))/(capGamma12*sigma))+phi1))+&
											(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
											exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapR))/(capGamma21*sigma))* &
											cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lanakcapR))/ &
											(capGamma22*sigma))+phi2))+pstana(0,ir)
								Fnodk =		lanak-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panak
								FnodkcapL =	lanakcapL-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panakcapL
								FnodkcapR =	lanakcapR-(((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
											(((1.0+nu)*(1.0-2.0*nu))/(1.0-nu)))*panakcapR
								if (FnodkcapL*Fnodk .lt. 0.0) then
									lanakcapR =		lanak
								else
									lanakcapL =		lanak
								end if
								Fnod =				Fnodk
								lana(0,ir) = 		lanak
								pana(0,ir) = 		panak
								pDLana(0,ir) = 		(((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
													((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
													exp(-epsi1(it)*capK*capR*(capH(ir)+eta1(it)*lana(0,ir)))
								pvdWana(0,ir) = 	-capAsfw/ &
													(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0* &
													(capH(ir)+eta1(it)*lana(0,ir))**3.0)
								pSanaa(0,ir) =   	-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
													exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma11*sigma))* &
													cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/ &
													(capGamma12*sigma))+phi1)
								pSanab(0,ir) =   	-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
													exp(-(epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/(capGamma21*sigma))* &
													cos(((2.0*pi*epsi1(it)*capR*(capH(ir)+eta1(it)*lana(0,ir)))/ &
													(capGamma22*sigma))+phi2)
								pSana(0,ir) =		pSanaa(0,ir)+pSanab(0,ir)
								pana(0,ir) = 		pstana(0,ir)+pDLana(0,ir)+pvdWana(0,ir)+pSana(0,ir)
							end do
						end if
					end if
					if (termout .eq. 1) then
							write(*,*)	"time						",			t(it),	"	time step	",		it	
							write(*,*)	"altr-root error				",		Fnod,	"	node no		",		ir, &
										"	iteration no	",					iternod
							write(*,*)	"hydrodynamic error				",		Fphd,	"	iteration no	",	iterphd
							write(*,*)	"amplification switch				",	(eta1(it)*capMamp(it) > thresamp)
							write(*,*)	"divergence precheck				",	(capMamp(it) > thresdiv1)
							write(*,*)	"divergence switch				",		(capNdiv(it) > thresdiv2)
					else
							write(95,*)	"time						",			t(it),	"	time step	",		it	
							write(95,*)	"altr-root error				",		Fnod,	"	node no		",		ir, &
										"	iteration no	",					iternod
							write(95,*)	"hydrodynamic error				",		Fphd,	"	iteration no	",	iterphd
							write(95,*)	"amplification switch				",	(eta1(it)*capMamp(it) > thresamp)
							write(95,*)	"divergence precheck				",	(capMamp(it) > thresdiv1)
							write(95,*)	"divergence switch				",		(capNdiv(it) > thresdiv2)
					end if
					Fnod = 	0.0
				end do
			end do
			do ir = 1,nr
				lana(1,ir) = 		0.0
				pstana(1,ir) = 		0.0
				pDLana(1,ir) = 		0.0
				pvdWana(1,ir) = 	0.0
				pSana(1,ir) = 		0.0
				pana(1,ir) = 		0.0
				lana(2,ir) = 		0.0
				pstana(2,ir) = 		0.0
				pDLana(2,ir) = 		0.0
				pvdWana(2,ir) = 	0.0
				pSana(2,ir) = 		0.0
				pana(2,ir) = 		0.0
			end do
!			obtaining solution for hydrodynamic pressure
			do ir = 1,nr
				lana(3,ir) = 		lana(0,ir)
			end do
!	      	mumps initialization steps
			call mpi_init(ierr)
			mumps_par%comm = 	mpi_comm_world
			mumps_par%job = 		-1
			mumps_par%sym = 		0
			mumps_par%par = 		1
			call dmumps(mumps_par)
!		  	disabling terminal output from mumps
			mumps_par%ICNTL(1) = 	-1
			mumps_par%ICNTL(2) = 	-1
			mumps_par%ICNTL(3) = 	-1
			mumps_par%ICNTL(4) = 	-1
			mumps_par%ICNTL(5) = 	-1
			if (mumps_par%infog(1).lt.0) then
				write(92,*) "error in mumps initialization"
				write(92,*)  mumps_par%INFOG(1)
				write(92,*)  mumps_par%INFOG(2)
				if (termout .eq. 1) then
					write(*,*) "error in mumps initialization"
					write(*,*)  mumps_par%INFOG(1)
					write(*,*)  mumps_par%INFOG(2)
				end if
				goto 450
			end if
			if (mumps_par%myid .eq. 0) then
				mumps_par%n = nr
				mumps_par%nnz = 3*(nr-1)+1
				allocate(mumps_par%irn(mumps_par%nnz))
				allocate(mumps_par%jcn(mumps_par%nnz))
				allocate(mumps_par%a(mumps_par%nnz))
				allocate(mumps_par%rhs(mumps_par%n))
			end if
			ir = 1
			itemp = 1
			mumps_par%irn(itemp) = 		ir
			mumps_par%jcn(itemp) = 		ir+2
			mumps_par%a(itemp) = 		(-(r(ir+1)-r(ir))**2.0)* &
										(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir+2)-r(ir)))
			itemp = itemp+1
			mumps_par%irn(itemp) = 		ir
			mumps_par%jcn(itemp) = 		ir+1
			mumps_par%a(itemp) = 		((r(ir+2)-r(ir))**2.0)* &
										(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir+2)-r(ir)))
			itemp = itemp+1
			mumps_par%irn(itemp) = 		ir
			mumps_par%jcn(itemp) = 		ir
			mumps_par%a(itemp) = 		(-((r(ir+2)-r(ir))**2.0-(r(ir+1)-r(ir))**2.0))* &
										(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir+2)-r(ir)))
			itemp = itemp+1
			mumps_par%rhs(ir) = 		(((r(ir+1)-r(ir))*(r(ir+2)-r(ir))*(r(ir+2)-r(ir+1)))*0.0)* &
										(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir+2)-r(ir)))
			containers(1) = 0.0
			do ir = 2,nr-1
				mumps_par%irn(itemp) = 	ir
				mumps_par%jcn(itemp) = 	ir+1
				mumps_par%a(itemp) = 	((r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**3.0)*(2.0*(r(ir-1)-r(ir))) + &
										((capH(ir)+eta1(it)*lana(3,ir))**3.0+3.0*r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**2.0* &
										(r(ir)+eta1(it)*differential1d(lana(3,:),'x','cd',ir,0,'x',r,r)))* &
										(-(r(ir-1)-r(ir))**2.0))* &
										((r(ir+1)-r(ir-1))**2.0)
				itemp = itemp+1
				mumps_par%irn(itemp) = 	ir
				mumps_par%jcn(itemp) = 	ir
				mumps_par%a(itemp) = 	((r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**3.0)*(2.0*(r(ir+1)-r(ir-1))) + &
										((capH(ir)+eta1(it)*lana(3,ir))**3.0+3.0*r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**2.0* &
										(r(ir)+eta1(it)*differential1d(lana(3,:),'x','cd',ir,0,'x',r,r)))* &
										((r(ir-1)-r(ir))**2.0-(r(ir+1)-r(ir))**2.0))* &
										((r(ir+1)-r(ir-1))**2.0)
				itemp = itemp+1
				mumps_par%irn(itemp) = 	ir
				mumps_par%jcn(itemp) = 	ir-1
				mumps_par%a(itemp) = 	((r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**3.0)*(-2.0*(r(ir+1)-r(ir))) + &
										((capH(ir)+eta1(it)*lana(3,ir))**3.0+3.0*r(ir)*(capH(ir)+eta1(it)*lana(3,ir))**2.0* &
										(r(ir)+eta1(it)*differential1d(lana(3,:),'x','cd',ir,0,'x',r,r)))* &
										((r(ir+1)-r(ir))**2.0))* &
										((r(ir+1)-r(ir-1))**2.0)
				itemp = itemp+1
				if (it .eq. 1) then
					dlana3dt = 0
				else
					if (phditer .eq. 0) then
						if (it .eq. 2) then
							dlana3dt = 	((lana(3,ir)-lanalast(3,ir))/(t(it)-t(it-1)))
						else
							dlana3dt = 	(-(t(it-1)-t(it))**2.0*lanalast2(3,ir)+(t(it-2)-t(it))**2.0*lanalast(3,ir) - &
										(((t(it-2)-t(it))**2.0)-((t(it-1)-t(it))**2.0))*lana(3,ir))/ &
										((t(it-2)-t(it))*(t(it-1)-t(it))*(t(it-2)-t(it-1)))
						end if
					else 
						if (phditer .eq. 1) then
							dlana3dt = 	0.0
						end if
					end if
				end if
				mumps_par%rhs(ir) = 	(((12.0*r(ir)*(epsi1(it)/epsi)*(eta1(it)/alpha)*dlana3dt)-(12.0*r(ir)*sin(t(it))))* &
										((r(ir+1)-r(ir))*(r(ir-1)-r(ir))*(r(ir+1)-r(ir-1))))* &
										((r(ir+1)-r(ir-1))**2.0)
			end do
			ir = nr
			mumps_par%irn(itemp) = 		ir
			mumps_par%jcn(itemp) = 		ir
			mumps_par%a(itemp) = 		(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir)-r(ir-2))**3.0)
			mumps_par%rhs(ir) = 		0.0*(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir)-r(ir-2))**3.0)
			if (it .ne. 1) then	
				if (phdpatch .eq. 1) then
					mumps_par%rhs(ir) = 	((4.0/3.0)*sin(t(it)))*((epsi1(it)/epsi1(it-1))**2.0)* &
											(((rmax-rmin)/2.0)*(capH(ir)+eta1(it)*lana(3,ir))**3.0*(r(ir)-r(ir-2))**3.0)
				end if
			end if
			mumps_par%JOB = 6
			call dmumps(mumps_par)
			if (mumps_par%infog(1).lt.0) then
				write(92,*) "error in mumps execution"
				write(92,*)  mumps_par%INFOG(1)
				write(92,*)  mumps_par%INFOG(2)
				if (termout .eq. 1) then
					write(*,*) "error in mumps execution"
					write(*,*)  mumps_par%INFOG(1)
					write(*,*)  mumps_par%INFOG(2)
				end if
				goto 450
			end if
			Fphdprev = Fphd
			Fphd = 0.0
			do ir = 1,nr
				Fphd = 	Fphd+abs(pstana(3,ir)-((1.0-relaxphd)*pstana(3,ir)+relaxphd*mumps_par%rhs(ir)))
			end do
			Fphd = 		Fphd/nr
			do ir = 1,nr
				pstana(3,ir) = (1.0-relaxphd)*pstana(3,ir)+relaxphd*mumps_par%rhs(ir)
				pstana(0,ir) = (1.0-relaxphd)*pstana(0,ir)+relaxphd*mumps_par%rhs(ir)
				pstana(1,ir) = 0.0
				pstana(2,ir) = 0.0
			end do
			if (nophd .eq. 1) then
				Fphd = 0.0
			end if
!      	 	mumps finalization steps
			if (mumps_par%myid.eq.0) then
				deallocate(mumps_par%irn)
				deallocate(mumps_par%jcn)
				deallocate(mumps_par%a)
				deallocate(mumps_par%rhs)
			end if
			mumps_par%job = -2
			call dmumps(mumps_par)
			if(mumps_par%INFOG(1).lt.0) then
				write(92,*) "error in mumps finalization"
				write(92,*)  mumps_par%INFOG(1)
				write(92,*)  mumps_par%INFOG(2)
				if (termout .eq. 1) then
					write(*,*) "error in mumps finalization"
					write(*,*)  mumps_par%INFOG(1)
					write(*,*)  mumps_par%INFOG(2)
				end if
				goto 450
			end if
			450 call mpi_finalize(ierr)
		end do
		if (termout .eq. 1) then
			write(*,*)	"time-step 						",it,	"	time		",t(it),	"	hd error		",Fphd
		else
			write(95,*)	"time-step 						",it,	"	time		",t(it),	"	hd error		",Fphd
		end if
!		solution combination
		do ir = 1,nr
			pDLana(3,ir) = 		pDLana(0,ir)
			pvdWana(3,ir) = 	pvdWana(0,ir)
			pSana(3,ir) = 		pSana(0,ir)
			pana(3,ir) = 		pana(0,ir)
		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	obtaining solution for hydrodynamic pressure regime
	else
!		zeroeth order solution
!		obtaining pressure terms
		do ir = 1,nr
			pstana(0,ir) = 	(3.0*sin(t(it)))*((1/(capH(ir)**2.0)))
			pDLana(0,ir) =  (((64.0*epsi1(it)**2.0*n0*kcapB*capT)/(mu*omega*alpha*epsi))* &
							((tanh((q*0.5*(zeta(1)+zeta(2)))/(4.0*kcapB*capT)))**2.0))* &
							exp(-epsi1(it)*capK*capR*capH(ir))
			pvdWana(0,ir) = -capAsfw/(6.0*pi*epsi1(it)*epsi*alpha*mu*omega*capR**3.0*capH(ir)**3.0)
			pSanaa(0,ir) = 	(-((capLambda1*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
							exp(-(epsi1(it)*capR*capH(ir))/(capGamma11*sigma))* &
							cos(((2.0*pi*epsi1(it)*capR*capH(ir))/(capGamma12*sigma))+phi1))
			pSanab(0,ir) = 	(-((capLambda2*rho_inf*kcapB*capT*epsi1(it)**2)/(mu*omega*alpha*epsi))* &
							exp(-(epsi1(it)*capR*capH(ir))/(capGamma21*sigma))* &
							cos(((2.0*pi*epsi1(it)*capR*capH(ir))/(capGamma22*sigma))+phi2))
			pSana(0,ir) = 	pSanaa(0,ir)+pSanab(0,ir)
		end do
		do ir = 1,nr
			pana(0,ir) =	pstana(0,ir)+pDLana(0,ir)+pvdWana(0,ir)+pSana(0,ir)
		end do
		do ir = 1,nr
			lana(0,ir) =	((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
							(((1.0-2.0*nu)*(1.0+nu))/(1.0-nu))*pana(0,ir)
		end do
!		first order
!		obtaining pressure terms
!		iterating for convergence of hydrodynamic pressure
!		obtaining solution for first order hydrodynamic pressure from current iteration
!	    mumps initialization steps
		call mpi_init(ierr)
		mumps_par%comm = mpi_comm_world
		mumps_par%job = -1
		mumps_par%sym = 0
		mumps_par%par = 1
		call dmumps(mumps_par)
!		disabling terminal output from mumps
		mumps_par%ICNTL(1) = -1
		mumps_par%ICNTL(2) = -1
		mumps_par%ICNTL(3) = -1
		mumps_par%ICNTL(4) = -1
		mumps_par%ICNTL(5) = -1
		if (mumps_par%infog(1).lt.0) then
			write(92,*) "error in mumps initialization"
			write(92,*)  mumps_par%INFOG(1)
			write(92,*)  mumps_par%INFOG(2)
			if (termout .eq. 1) then
				write(*,*) "error in mumps initialization"
				write(*,*)  mumps_par%INFOG(1)
				write(*,*)  mumps_par%INFOG(2)
			end if
			goto 600
		end if
		if (mumps_par%myid .eq. 0) then
			mumps_par%n = nr
			mumps_par%nnz = 3*(nr-1)+1
			allocate(mumps_par%irn(mumps_par%nnz))
			allocate(mumps_par%jcn(mumps_par%nnz))
			allocate(mumps_par%a(mumps_par%nnz))
			allocate(mumps_par%rhs(mumps_par%n))
		end if
		ir = 1
		itemp = 1
		mumps_par%irn(itemp) = 		ir
		mumps_par%jcn(itemp) = 		ir+2
		mumps_par%a(itemp) = 		(-(r(ir+1)-r(ir))**2.0)* &
									(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir+2)-r(ir)))
		itemp = itemp+1
		mumps_par%irn(itemp) = 		ir
		mumps_par%jcn(itemp) = 		ir+1
		mumps_par%a(itemp) = 		((r(ir+2)-r(ir))**2.0)* &
									(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir+2)-r(ir)))
		itemp = itemp+1
		mumps_par%irn(itemp) = 		ir
		mumps_par%jcn(itemp) = 		ir
		mumps_par%a(itemp) = 		(-((r(ir+2)-r(ir))**2.0-(r(ir+1)-r(ir))**2.0))* &
									(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir+2)-r(ir)))
		itemp = itemp+1
		mumps_par%rhs(ir) = 		(((r(ir+1)-r(ir))*(r(ir+2)-r(ir))*(r(ir+2)-r(ir+1)))*0.0)* &
									(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir+2)-r(ir)))
		containers(1) = 0.0
		do ir = 2,nr-1
			mumps_par%irn(itemp) = 	ir
			mumps_par%jcn(itemp) = 	ir+1
			mumps_par%a(itemp) = 	((r(ir)*capH(ir)**3.0)*(2.0*(r(ir-1)-r(ir))) + &
									(capH(ir)**3.0+3.0*r(ir)*2.0*capH(ir)**2.0)* &
									(-(r(ir-1)-r(ir))**2.0))* &
									((r(ir+1)-r(ir-1))**2.0)
			itemp = itemp+1
			mumps_par%irn(itemp) = 	ir
			mumps_par%jcn(itemp) = 	ir
			mumps_par%a(itemp) = 	((r(ir)*capH(ir)**3.0)*(2.0*(r(ir+1)-r(ir-1))) + &
									(capH(ir)**3.0+3.0*r(ir)*2.0*capH(ir)**2.0)* &
									((r(ir-1)-r(ir))**2.0-(r(ir+1)-r(ir))**2.0))* &
									((r(ir+1)-r(ir-1))**2.0)
			itemp = itemp+1
			mumps_par%irn(itemp) = 	ir
			mumps_par%jcn(itemp) = 	ir-1
			mumps_par%a(itemp) = 	((r(ir)*capH(ir)**3.0)*(-2.0*(r(ir+1)-r(ir))) + &
									(capH(ir)**3.0+3.0*r(ir)*2.0*capH(ir)**2.0)* &
									((r(ir+1)-r(ir))**2.0))* &
									((r(ir+1)-r(ir-1))**2.0)
			itemp = itemp+1
			if (it .eq. 1) then
				dlana0dt = 0
			else
				if (phditer .eq. 0) then
					if (it .eq. 2) then
						dlana0dt = 	((lana(0,ir)-lanalast(0,ir))/(t(it)-t(it-1)))
					else
						dlana0dt = 	(-(t(it-1)-t(it))**2.0*lanalast2(0,ir)+(t(it-2)-t(it))**2.0*lanalast(0,ir) - &
									(((t(it-2)-t(it))**2.0)-((t(it-1)-t(it))**2.0))*lana(0,ir))/ &
									((t(it-2)-t(it))*(t(it-1)-t(it))*(t(it-2)-t(it-1)))
					end if
				else 
					if (phditer .eq. 1) then
						dlana0dt = 	0.0
					end if
				end if
			end if
			mumps_par%rhs(ir) = 	(((12.0*r(ir)*(epsi1(it)/epsi)*(1.0/alpha)*dlana0dt)- &
									3.0*r(ir)*capH(ir)**3.0*lana(0,ir)*differential1d(pstana(0,:),'xx','cd',ir,0,'x',r,r)- &
									3.0*(capH(ir)**2.0*lana(0,ir)+2.0*r(ir)**2.0*capH(ir)*lana(0,ir)+r(ir)*capH(ir)* &
									differential1d(lana(0,:),'x','cd',ir,0,'x',r,r))* &
									differential1d(pstana(0,:),'x','cd',ir,0,'x',r,r))* &
									((r(ir+1)-r(ir))*(r(ir-1)-r(ir))*(r(ir+1)-r(ir-1))))* &
									((r(ir+1)-r(ir-1))**2.0)
		end do
		ir = nr
		mumps_par%irn(itemp) = 		ir
		mumps_par%jcn(itemp) = 		ir
		mumps_par%a(itemp) = 		(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir)-r(ir-2))**3.0)
		mumps_par%rhs(ir) = 		0.0*(((rmax-rmin)/2.0)*capH(ir)**3.0*(r(ir)-r(ir-2))**3.0)
		mumps_par%JOB = 6
		call dmumps(mumps_par)
		if (mumps_par%infog(1).lt.0) then
			write(92,*) "error in mumps execution"
			write(92,*)  mumps_par%INFOG(1)
			write(92,*)  mumps_par%INFOG(2)
			if (termout .eq. 1) then
				write(*,*) "error in mumps execution"
				write(*,*)  mumps_par%INFOG(1)
				write(*,*)  mumps_par%INFOG(2)
			end if
			goto 600
		end if
		do ir = 1,nr
			pstana(1,ir) = mumps_par%rhs(ir)
		end do
!       mumps finalization steps
		if (mumps_par%myid.eq.0) then
			deallocate(mumps_par%irn)
			deallocate(mumps_par%jcn)
			deallocate(mumps_par%a)
			deallocate(mumps_par%rhs)
		end if
		mumps_par%job = -2
		call dmumps(mumps_par)
		if(mumps_par%INFOG(1).lt.0) then
			write(92,*) "error in mumps finalization"
			write(92,*)  mumps_par%INFOG(1)
			write(92,*)  mumps_par%INFOG(2)
			if (termout .eq. 1) then
				write(*,*) "error in mumps finalization"
				write(*,*)  mumps_par%INFOG(1)
				write(*,*)  mumps_par%INFOG(2)
			end if
			goto 600
		end if				
		600 call mpi_finalize(ierr)
!		obtaining pressure terms
		do ir = 1,nr
			pDLana(1,ir) = 	-capK*epsi1(it)*capR*lana(0,ir)*pDLana(0,ir)
			pvdWana(1,ir) = -((3.0*lana(0,ir))/capH(ir))*pvdWana(0,ir)
			pSanaa(1,ir) = 	-((epsi1(it)*capR*lana(0,ir))/sigma)* &
							((1.0/capGamma11)+((2.0*pi)/capGamma12)* &
							tan(((2.0*pi*epsi1(it)*capR*capH(ir))/(capGamma12*sigma))+phi1))*pSanaa(0,ir)
			pSanab(1,ir) = 	-((epsi1(it)*capR*lana(0,ir))/sigma)* &
							((1.0/capGamma21)+((2.0*pi)/capGamma22)* &
							tan(((2.0*pi*epsi1(it)*capR*capH(ir))/(capGamma22*sigma))+phi2))*pSanab(0,ir)
			pSana(1,ir) =	pSanaa(1,ir)+pSanab(1,ir)
			pana(1,ir) = 	pstana(1,ir) + pDLana(1,ir) + pvdWana(1,ir) + pSana(1,ir)		
		end do
		do ir = 1,nr
			lana(1,ir) = 	((mu*omega*alpha*epsi)/(epsi1(it)**2.0*kappa1(it)*capE))* &
							(((1.0-2.0*nu)*(1.0+nu))/(1.0-nu))*pana(1,ir)
		end do
		do ir = 1,nr
			pstana(2,ir) = 	0.0
			pDLana(2,ir) =  0.0
			pvdWana(2,ir) = 0.0
			pSana(2,ir) =   0.0
			pana(2,ir) =   	0.0
			lana(2,ir) =   	0.0			
		end do
!		solution combination
		do ir = 1,nr
			lana(3,ir) = 		lana(0,ir)+		eta1(it)*lana(1,ir)+		eta1(it)**2.0*lana(2,ir)
			pstana(3,ir) = 		pstana(0,ir)+	eta1(it)*pstana(1,ir)+		eta1(it)**2.0*pstana(2,ir) 
			pDLana(3,ir) = 		pDLana(0,ir)+	eta1(it)*pDLana(1,ir)+		eta1(it)**2.0*pDLana(2,ir)
			pvdWana(3,ir) = 	pvdWana(0,ir)+	eta1(it)*pvdWana(1,ir)+		eta1(it)**2.0*pvdWana(2,ir)
			pSana(3,ir) = 		pSana(0,ir)+	eta1(it)*pSana(1,ir)+		eta1(it)**2.0*pSana(2,ir)
			pana(3,ir) = 		pana(0,ir)+		eta1(it)*pana(1,ir)+		eta1(it)**2.0*pana(2,ir)
		end do
	end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (phdaxi .eq. 1) then
		ir = 1
		pstana(0,ir) = 	(-(r(ir+1)-r(ir))**2.0*pstana(0,ir+2)+(r(ir+2)-r(ir))**2.0*pstana(0,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pstana(1,ir) = 	(-(r(ir+1)-r(ir))**2.0*pstana(1,ir+2)+(r(ir+2)-r(ir))**2.0*pstana(1,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pstana(2,ir) = 	(-(r(ir+1)-r(ir))**2.0*pstana(3,ir+2)+(r(ir+2)-r(ir))**2.0*pstana(2,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pstana(3,ir) = 	(-(r(ir+1)-r(ir))**2.0*pstana(4,ir+2)+(r(ir+2)-r(ir))**2.0*pstana(3,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
	end if
	if (profaxi .eq. 1) then
		ir = 1
		lana(0,ir) = 	(-(r(ir+1)-r(ir))**2.0*lana(0,ir+2)+(r(ir+2)-r(ir))**2.0*lana(0,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		lana(1,ir) = 	(-(r(ir+1)-r(ir))**2.0*lana(1,ir+2)+(r(ir+2)-r(ir))**2.0*lana(1,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		lana(2,ir) = 	(-(r(ir+1)-r(ir))**2.0*lana(3,ir+2)+(r(ir+2)-r(ir))**2.0*lana(2,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		lana(3,ir) = 	(-(r(ir+1)-r(ir))**2.0*lana(4,ir+2)+(r(ir+2)-r(ir))**2.0*lana(3,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
	end if
	if (molaxi .eq. 1) then
		ir = 1
		pDLana(0,ir) = 	(-(r(ir+1)-r(ir))**2.0*pDLana(0,ir+2)+(r(ir+2)-r(ir))**2.0*pDLana(0,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pDLana(1,ir) = 	(-(r(ir+1)-r(ir))**2.0*pDLana(1,ir+2)+(r(ir+2)-r(ir))**2.0*pDLana(1,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pDLana(2,ir) = 	(-(r(ir+1)-r(ir))**2.0*pDLana(3,ir+2)+(r(ir+2)-r(ir))**2.0*pDLana(2,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pDLana(3,ir) = 	(-(r(ir+1)-r(ir))**2.0*pDLana(4,ir+2)+(r(ir+2)-r(ir))**2.0*pDLana(3,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pvdWana(0,ir) = (-(r(ir+1)-r(ir))**2.0*pvdWana(0,ir+2)+(r(ir+2)-r(ir))**2.0*pvdWana(0,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pvdWana(1,ir) = (-(r(ir+1)-r(ir))**2.0*pvdWana(1,ir+2)+(r(ir+2)-r(ir))**2.0*pvdWana(1,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pvdWana(2,ir) = (-(r(ir+1)-r(ir))**2.0*pvdWana(3,ir+2)+(r(ir+2)-r(ir))**2.0*pvdWana(2,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pvdWana(3,ir) = (-(r(ir+1)-r(ir))**2.0*pvdWana(4,ir+2)+(r(ir+2)-r(ir))**2.0*pvdWana(3,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pSana(0,ir) = 	(-(r(ir+1)-r(ir))**2.0*pSana(0,ir+2)+(r(ir+2)-r(ir))**2.0*pSana(0,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pSana(1,ir) = 	(-(r(ir+1)-r(ir))**2.0*pSana(1,ir+2)+(r(ir+2)-r(ir))**2.0*pSana(1,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pSana(2,ir) = 	(-(r(ir+1)-r(ir))**2.0*pSana(3,ir+2)+(r(ir+2)-r(ir))**2.0*pSana(2,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
		pSana(3,ir) = 	(-(r(ir+1)-r(ir))**2.0*pSana(4,ir+2)+(r(ir+2)-r(ir))**2.0*pSana(3,ir+1))/ &
						((r(ir+2)+r(ir+1)-2*r(ir))*(r(ir+2)-r(ir+1)))
	end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
!	obtaining dimensional values
	do ir = 1,nr
		rd(ir) = 					sqrt(epsi1(it))*capR*r(ir)
	end do
	do ir = 1,nr
		do itemp = 0,3
			panad(itemp,ir) = 		pana(itemp,ir)*((mu*omega*alpha*epsi)/(epsi1(it)**2.0))
			pstanad(itemp,ir) = 	pstana(itemp,ir)*((mu*omega*alpha*epsi)/(epsi1(it)**2.0))
			pDLanad(itemp,ir) = 	pDLana(itemp,ir)*((mu*omega*alpha*epsi)/(epsi1(it)**2.0))
			pvdWanad(itemp,ir) = 	pvdWana(itemp,ir)*((mu*omega*alpha*epsi)/(epsi1(it)**2.0))
			pSanad(itemp,ir) = 		pSana(itemp,ir)*((mu*omega*alpha*epsi)/(epsi1(it)**2.0))
			lanad(itemp,ir) = 		lana(itemp,ir)*kappa1(it)*capL			
		end do
		capHd(ir) = 				capH(ir)*epsi1(it)*capR
	end do
	td(it) = (1.0/omega)*t(it)

end subroutine anainst
!-----------------------------------------------------------------------------------------------------------------------------------
