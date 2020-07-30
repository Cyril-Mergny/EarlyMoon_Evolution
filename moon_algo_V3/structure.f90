SUBROUTINE Rheology
    use variables
	implicit none
    integer :: i,ncc,j,idx
    real(kind=dp) ::  dr
	real(kind=dp) :: threshold_eta

	!Assigning values to arrays profiles
	ncc = 0
	do ilay = 1,5
		if (ilay .eq. 1) then ; dr = Rlim(1)/(n_layer(1)-1) ;
		elseif (ilay .eq. 3 .or. ilay .eq. 4 .or. ilay .eq. 5) then ; dr = drmag
		else ; dr = (Rlim(ilay)-Rlim(ilay-1))/n_layer(ilay)
		endif
		do i = 1, n_layer(ilay)
			ncc = ncc+1
			if (ncc .eq. 1) then ; r_in(1) = 0._dp
			else ; r_in(ncc) = r_in(ncc-1)+dr
			endif
			rho_in(ncc) = rho_lay(ilay)
			vp_in(ncc) = vp_lay(ilay)
			vs_in(ncc) = vs_lay(ilay)
			eta_in(ncc) = eta_lay(ilay)
			D_in(ncc) = lambda_lay(ilay)/rho_lay(ilay)/cp_lay(ilay)
			enddo
		enddo
	!Change interface radius to exact rad
	r_in(n_lim(3)) = Rlim(3)
	r_in(n_lim(4)) = Rlim(4)
	
	open(10,file="Outputs/melt.txt",position="append")
	!Crust Melting and Viscosity changes with temperature
	do i = n_lim(4)+1,n_lim(5)
		!Calcul meltfraction
		phi_melt = -22._dp/(T(i)-2134+0.1724*(r_in(i)/1e3)+1.374d-4*(r_in(i)/1.e3)**2) -0.05_dp
		!phi_melt = 0
		if (phi_melt .lt. 0) phi_melt = 0._dp
		if (phi_melt .gt. 1) phi_melt = 1._dp
		!Melt
		if (phi_melt .gt. 1.e-4) then
			if (phi_melt .gt. 0.301) phi_melt = 0.301_dp !tides.f90 can not handle phi>0.3
			CALL Melt_models(i)
		else !No melt: Arrhenius law
			Tsolidus= 2134._dp -0.1724_dp*(Rlim(4)/1e3) - 1.3714d-4*(Rlim(4)/1e3)**2 - 4.4_dp/0.01_dp
			eta_in(i) = etacrust_base*exp(Enj/Rgaz/Tbase*(Tsolidus/T(i) - 1))
			threshold_eta = 1d30
			if (eta_in(i) .gt. threshold_eta) then
				eta_in(i) = threshold_eta
				endif
			endif
		write(10,*) phi_melt
		enddo
	
	write(10,*) "-"
	close(10)
		
END SUBROUTINE Rheology


SUBROUTINE Internal_structure
	!**************************************************************
	!Ouputs the limit Radius (Rlim) of each interface (e.g. crust/mantle)
	!and the density in each layers
	!**************************************************************
	use variables
	implicit none
	!Radius Interface
	Rlim(5)=rtot
	Rlim(4)=Rlim(5)-rcrust
	Rlim(3)=Rlim(4)-rmag
	Rlim(2)=Rlim(3)-rast
	Rlim(1)=rcore
	!nbr of pts in each layer
	n_lim(1) = n_layer(1)
	n_lim(2) = n_lim(1)+n_layer(2)
	n_lim(3) = n_lim(2)+n_layer(3)
	n_lim(4) = n_lim(3)+n_layer(4)
	n_lim(5) = n_lim(4)+n_layer(5)
	
	! Crust
	rhocrust=2900_dp
	rho_lay(5)=rhocrust
	mu_lay(5)=mu_crust
	K_lay(5)=K_crust
	eta_lay(5) = etacrust_base
	cp_lay(5) = cpcrust
	lambda_lay(5) = lambda_mant
	! Magma ocean
	rhomant= 3300_dp
	rho_lay(4)=rhomant
	mu_lay(4)= 0
	K_lay(4)=  K_mant
	eta_lay(4) = 0
	cp_lay(4) = cpcrust
	lambda_lay(4) = lambda_mant
	! Astenophere
	!call Melt_models(i)
	rhomant= 3300_dp
	rho_lay(3)=rhomant
	mu_lay(3)=mu_mant
	K_lay(3)=K_mant
	eta_lay(3)=etamant
	cp_lay(3) = cpcrust
	lambda_lay(3) = lambda_mant
	! Mantle
	rho_lay(2)=rhomant
	mu_lay(2)=mu_mant
	K_lay(2)=K_mant
	eta_lay(2)=etamant
	cp_lay(2) = cpcrust
	lambda_lay(2) = lambda_mant
	! Core
	rhocore= 7800_dp
	rho_lay(1) = rhocore
	mu_lay(1) = mu_core
	K_lay(1) = K_core
	eta_lay(1) = etacore
	cp_lay(1) = cpcrust
	lambda_lay(1) = lambda_mant

	vs_lay=sqrt(mu_lay/rho_lay)
	vp_lay=sqrt((K_lay+4._dp/3._dp*mu_lay)/rho_lay)

END SUBROUTINE Internal_structure

SUBROUTINE Melt_models(i)
	! For partially melt layers
	! Change viscosity and shear modulus
	! with melt fraction
    use variables
    implicit none
    integer :: i
    real(kind=dp) :: alphamelt,mu_ast,gamma,zeta_mu,zeta_eta,delta,mu_melt
    real(kind=dp) :: phi_crit,phi_crit_mu,errorfunc,errorfunc_mu,poisson_ratio

	!"Composite law"
	alphamelt = 30._dp
	poisson_ratio = 0.25_dp
	mu_ast = 1e1
	gamma = 5._dp
	zeta_mu = 1e-6
	zeta_eta = 1e-10
	delta = 13._dp-gamma
	phi_crit = 0.56_dp
	phi_crit_mu = 0.6_dp
	errorfunc=(1-zeta_eta)*erf(((sqrt(pi))/(2*(1-zeta_eta)))*((1-phi_melt)/phi_crit)*(1+((1-phi_melt)/phi_crit)**gamma))
	errorfunc_mu=(1-zeta_mu)*erf(((sqrt(pi))/(2*(1-zeta_mu)))*((1-phi_melt)/phi_crit_mu) *(1+((1-phi_melt)/phi_crit_mu)**gamma))

	if (phi_melt.lt.0.30) then
		eta_in(i) = etacrust_base*exp(-alphamelt*phi_melt)
		mu_melt =((1/mu_lay(5))+((phi_melt/mu_lay(5))*((40-24*poisson_ratio)/15)))**(-1) !rond
	elseif (phi_melt.gt.0.29) then
		eta_in(i) = ((1+((1-phi_melt)/phi_crit)**(delta))/(1-errorfunc)**(2.5*phi_crit))
		mu_melt = mu_ast*(((1+((1-phi_melt)/phi_crit_mu)**(delta))/(1-errorfunc_mu)**(2.5*phi_crit_mu)))
		endif
	vs_in(i)=sqrt(mu_melt/rho_lay(5))
	vp_in(i)=sqrt((K_lay(5)+4._dp/3._dp*mu_melt)/rho_lay(5))
	!higher conductivity in mushy law
	D_in(i) = D_in(i) + D_in(i)*(phi_melt)*50._dp

    
END SUBROUTINE Melt_models




