
SUBROUTINE conduction
    use variables
    use param_y, only : htide
    implicit none
    integer :: i, outloop
    real(kind=dp) :: dt,maxDeltaT, thresh_min,thresh_max
    real(kind=dp), dimension(ncond) :: D,Cp,H,R
    real(kind=dp), dimension(4,ncond) :: Diag
	real(kind=dp), dimension(2,ncond) :: Tcond
	
	!Copy of arrays for conduction
	D = D_in(n_lim(4):n)
	D(1) = D(2)
    Cp = Cpcrust
    H = htide(n_lim(4):n)/rho_in(n_lim(4):n)
    R = r_in(n_lim(4):n)
    Tcond(2,:) = T(n_lim(4):n)
	
	!Adaptative timestep
	thresh_max  = 10_dp
	thresh_min  = 3_dp
	outloop = 0
	
	do while (outloop .eq. 0)
	
		dt = timestep/1.d3
		mode = "base_temp"
		!Main conduction Loop
		do i = 1,int(timestep/dt)
			Tcond(1,:) = Tcond(2,:)
			call Matrice(R,H,dt,Diag,D,Cp,Tcond)
			call tridag(Diag(3,:),Diag(2,:),Diag(1,:),Diag(4,:),Tcond(2,:),ncond)
			enddo
			
		outloop = 1
		!Computing max diff bw prev & current Temp
		maxDeltaT  = maxval(abs(T(n_lim(4):n)-Tcond(2,:)))
		!Reset temperature if timestep too big
		if (maxDeltaT .gt. thresh_max) then
			print *, "Reducing timestep (too big), calculing temperature profile again"
			print *,"maxDT=",maxDeltaT, "timestep=",timestep/1d6/year
			timestep = timestep/2._dp
			Tcond(2,:) = T(n_lim(4):n)
			outloop = 0
			endif
		enddo

	!Increasing timestep when Temp changes negligeable
	if (maxDeltaT .lt. thresh_min .and.  dcrust_rate*timestep/drmag .lt. 0.5 .and. n_layer(4) .gt. 0) then
		timestep = timestep + timestep*(thresh_min-maxDeltaT)*2._dp
		endif
	if (timestep .gt. 10.e6*year) timestep = 10.e6*year
		
		
	!Assigning new values to public array T
	T(n_lim(4):n) = Tcond(2,:)
	
	print *, "phi", -22._dp/(Tbase-2134._dp + 0.1724*(Rlim(4)/1e3)+1.3714d-4*(Rlim(4)/1.e3)**2) -0.05_dp
	print *," Tmagma (K) ", Tbase
	
	print *,"timestep (Myr) =",timestep/1d6/year
	print *,"time (Myr) = ", time/1d6/year
	print *, "a =",a/Rearth,", e=",ecc
	print *, n_layer
	print *,"--------------------------------------------------"

END SUBROUTINE conduction


!**********CRISTALLISATION/FUSION**********

SUBROUTINE cristallisation
    ! Calculate dissipated nrj
    ! by radioactivity Hussmann2010,
    ! abundance C in Laneuville2013
    use variables
    implicit none
	real(kind=dp) :: Delta1,Delta11,Vocean,dVocean
    real(kind=dp) :: L,dVolcrust,dVolast,dast,dmant_rate
    
    !Bilan nrj -> change in lid radius
    L = 5.d5
    f = 0.2_dp
    Vocean = 4._dp/3._dp*pi*((Rtot-rcrust)**3-(Rtot-rcrust-rmag)**3)
    CALL derive_Tbase
	CALL power_crust
	dcrust_rate = -lambda_lay(5)*DTbase/L/rhocrust - (Prad_mant+Prad_core+Ptide_mant)/(4*pi*Rlim(4)**2)/L/rhocrust
	dcrust_rate = dcrust_rate/(1 + rhomant/rhocrust*(1-f)/f + Cpcrust*rhomant*Vocean*DTocean/(4*pi*Rlim(4)**2)/L/rhocrust)
	if (n_layer(4) .le. 0 ) dcrust_rate = 0

	!Changing radius values
	rcrust = rcrust + dcrust_rate*timestep
	
	!Cristallisation mantle
	dast = 0
	if (dcrust_rate .gt. 0) then
		dVolcrust = 4._dp*pi*(Rtot-rcrust)**2._dp*(dcrust_rate*timestep)
		dVolast = (1._dp-f)/f*dVolcrust
		dast = dVolast/(4._dp*pi*(Rtot-rast)**2)
		!dmant_rate = (1._dp-f)/f * (Rlim(4)/Rlim(3))**2 *dcrust_rate
		endif
	!Changing radius ocean
	rmag = rmag - dcrust_rate*timestep -dast
	rast = rast + dast
	
	!Change of ocean Temperature

	if (n_layer(4) - (nint(rcrust/drmag)-n_layer(5)) - (nint(rast/drmag)-n_layer(3)) .le. 0) then
		dVocean = -Vocean/timestep
		Vocean = 0._dp
	else
		dVocean = (4._dp/3._dp*pi*((Rtot-rcrust)**3-(Rtot-rcrust-rmag)**3)-Vocean)/timestep
		Vocean = 4._dp/3._dp*pi*((Rtot-rcrust)**3-(Rtot-rcrust-rmag)**3)
		endif
		
	DTocean = 0.1724_dp*dcrust_rate/1e3 + 1.3714d-4*(Rtot-rcrust)*2*dcrust_rate/1e6
	DTocean = (DTocean + 0.88_dp/Vocean_ini*dVocean/(0.2_dp*Vocean/Vocean_ini +0.01_dp)**2)
	
	!DTocean = -Tsolidus/f * 4._dp*pi*(Rtot-rcrust)**2._dp*dcrust_rate/Vocean_ini * timestep
	Tbase = Tbase + DTocean*timestep
	T(n_lim(4)) = Tbase
	
	
	!Changing Temperature grid when interface moves
	!cristallisation
	CALL derive_Tbase
	if (nint(rcrust/drmag) .gt. n_layer(5)) then
		!new T point + interpolation
		T(n_lim(4)-1) = Tbase
		Delta1 = r_in(n_lim(4))-r_in(n_lim(4)-1)
		Delta11 = r_in(n_lim(4)+1)-r_in(n_lim(4)-1)
		T(n_lim(4)) = DTbase*(Delta11**2*Delta1-Delta1**2*Delta11) + Delta1**2*T(n_lim(4)+1)+(Delta11**2-Delta1**2)*T(n_lim(4)-1)
		T(n_lim(4)) = T(n_lim(4))/Delta11**2
		endif
	!Fusion: deepest temp and assign it to neighb
	if (nint(rcrust/drmag) .lt. n_layer(5)) then
		T(n_lim(4)+1) = Tbase
		endif
		
	!Changing mapping
	!Carefull when ocean dissapears
	if (n_layer(4) - (nint(rcrust/drmag)-n_layer(5)) - (nint(rast/drmag)-n_layer(3)) .lt. 0) then
		n_layer(5) = n_layer(5) + n_layer(4) !The mantle takes it all
		n_layer(4) = 0
	else
		n_layer(4) = n_layer(4) - (nint(rcrust/drmag)-n_layer(5)) - (nint(rast/drmag)-n_layer(3))
		n_layer(5) = nint(rcrust/drmag)
		n_layer(3) = nint(rast/drmag)
		endif
	!New interface limits
	Rlim(5)=rtot
	Rlim(4)=Rlim(5)-rcrust
	Rlim(3)=Rlim(4)-rmag
	Rlim(2)=Rlim(3)-rast
	Rlim(1)=rcore
	!New spatial grid
	n_lim(3) = n_lim(2)+n_layer(3)
	n_lim(4) = n_lim(3)+n_layer(4)
	n_lim(5) = n_lim(4)+n_layer(5)
	ncond = n_layer(5)+1
	
END SUBROUTINE cristallisation

SUBROUTINE derive_Tbase
	! Three points derivate of the
	! Temperature at base of crust
    use variables
    implicit none
    integer :: i
    real(kind=dp) :: Delta1,Delta11
    real(kind=dp), dimension(ncond) :: Tcond,R
    
    Tcond = T(n_lim(4):n)
    R = r_in(n_lim(4):n)
    Delta1 = R(2)-R(1)
    Delta11 = R(3)-R(1)
    DTbase = Delta11**2*Tcond(2)-Delta1**2*Tcond(3)-(Delta11**2-Delta1**2)*Tcond(1)
    DTbase = DTbase/(Delta11**2*Delta1-Delta1**2*Delta11)
END SUBROUTINE derive_Tbase


!**********SUBROUTINES**********

SUBROUTINE Matrice(R,H,dt,Diag,D,Cp,Tcond)
    !Output diagonals and soruce terms
    !required to solve the thermic equations
    use variables, only : mode,ncond,n_layer
    use prec
	implicit none
	integer :: i
	real(kind=dp), dimension(2,ncond) :: Tcond
    real(kind=dp), dimension(ncond) :: D,Cp ; real :: D_p,D_m
    real(kind=dp),intent(in) :: dt
    real(kind=dp), dimension(ncond):: R,H
    real(kind=dp), dimension(4,ncond), intent(out) :: Diag
    
    if (mode == "base_temp") then
    !Boundary term - base
    Diag(1,1) = 0
    Diag(2,1) = 0.5
    Diag(3,1) = 0
    Diag(4,1) = 0.5*Tcond(1,1)
    elseif (mode == "base_flux") then
    Diag(1,1) = -D(1)*dt/(R(2)-R(1))**2
    Diag(2,1) = 0.5 + D(1)*dt/(R(2)-R(1))**2
    Diag(3,1) = 0
    Diag(4,1) = 0.5*Tcond(1,1)
    endif
    !Main Loop
    do i = 2,ncond-1
        D_p = (D(i)+D(i+1))/2.d0
        D_m = (D(i-1)+D(i))/2.d0
		!if (i .eq. n_layer(4)-1 .or. i .eq. n_layer(4)) then ; D_p = D(i) ;D_m = D(i) ; endif ! no propagation of D at interface
        Diag(1,i) = -dt/(R(i+1)-R(i-1))*(D(i)/R(i) + D_p/(R(i+1)-R(i)))
        Diag(2,i) = 0.5 + dt/(R(i+1)-R(i-1))*(D_p/(R(i+1)-R(i)) + D_m/(R(i)-R(i-1)))
        Diag(3,i) = dt/(R(i+1)-R(i-1))*(D(i)/R(i) - D_m/(R(i)-R(i-1)))
        Diag(4,i) = 0.5*Tcond(1,i) + H(i)*dt/2.d0/Cp(i)
        enddo
    !Boundary term -surface
    Diag(1,ncond) = 0
    Diag(2,ncond) = 0.5
    Diag(3,ncond) = 0
    Diag(4,ncond) = 0.5*Tcond(1,ncond)
    return
END SUBROUTINE Matrice


SUBROUTINE tridag(a,b,c,r,u,ncond)
    !Solve tridiagonal system Mu=r
    !Return the solution vector u
    integer :: ncond, j
    double precision, dimension(ncond) :: a, b, c, r, gam,u
    real :: bet
    bet=b(1)
    u(1)=r(1)/bet
    do j = 2,ncond
        gam(j) = c(j-1)/bet
        bet = b(j)-a(j)*gam(j)
        u(j) = (r(j)-a(j)*u(j-1))/bet
        enddo
    do j = ncond-1,1,-1
        u(j) = u(j)-gam(j+1)*u(j+1)
        enddo
    return
END SUBROUTINE tridag



SUBROUTINE initialise_temp
    !Create an initial temperature profile
    use variables
    integer :: i
    real(kind=dp) :: Delta_T
    
	!Linear Temperature
	Tbase = 2134._dp -0.1724_dp*(Rlim(4)/1e3) - 1.3714d-4*(Rlim(4)/1e3)**2 - 4.4_dp/(0.2_dp+0.01_dp)
	Tsolidus= 2134._dp -0.1724_dp*(Rlim(4)/1e3) - 1.3714d-4*(Rlim(4)/1e3)**2 - 4.4_dp/0.01_dp
	T = 0
    do i = n_lim(4),n
		T(i) = (i-n_lim(4))*(Tsurf-Tbase)/(ncond-1) + Tbase
        enddo
END SUBROUTINE initialise_temp





