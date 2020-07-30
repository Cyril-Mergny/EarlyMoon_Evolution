SUBROUTINE call_orbit
	! Calls orbitevol module made by Jonathan
	use orbit
	use variables
	
	kloveSat(:) = real(klove2)
	
	CALL power_crust
	Pdiss = Pcrust+Ptide_mant

!	print *,Rlim(4),n_layer(4),n_layer(5)
!	print *,ecc,a/Rearth,PlanSpin,Pdiss/1e12
!	print *,"------------------------------------"
	
	xincl = 0.001_dp
	beta = 0.001_dp
	CALL calc_orbitalEvolution(ecc, w_in, xincl, beta, PlanSpin, timestep)
	
	! upload period
    period =  2._dp*pi/w_in

	
END SUBROUTINE call_orbit




SUBROUTINE get_htide
    !Calculates htide
    use variables, only : zeta,ecc,w_in,rtot,n,ilay,nlay,n_lim
    use param_y, only : htide,HH,R,MUC
    implicit none
    real(kind(0.d0)) :: i1,beta, f0, f1, f2
    integer :: i

	!Zeta coefficient for htide
    beta=sqrt(1.d0-ecc*ecc)
    f0=1.d0+31.d0/2.d0*ecc**2.d0+255.d0/8.d0*ecc**4.d0+185.d0/16.d0*ecc**6.d0+25.d0/64.d0*ecc**8.d0
    f1=1.d0+15.d0/2.d0*ecc**2.d0+45.d0/8.d0*ecc**4.d0+5.d0/16.d0*ecc**6.d0
    f2=1.d0+3.d0*ecc**2.d0+3.d0/8.d0*ecc**4.d0
    zeta=2.d0/7.d0*f0/beta**15.d0-4.d0/7.d0*f1/beta**12.d0*cos(i1)+1.d0/7.d0*f2/beta**9.d0*(1.d0+(cos(i1))**2.d0)
    htide = 21.d0/10.d0*REAL(HH(2,:))*IMAG(MUC)*(rtot)**4/R**2*(w_in**5)*zeta
    htide(1) = 0.d0
    !Remap htide size to erase
    !double values at interface
	do ilay = 1, nlay-1
		do i = n_lim(ilay)+1, n+nlay
		htide(i) = htide(i+1)
		enddo
		enddo
END SUBROUTINE get_htide



SUBROUTINE power_crust
    ! intern power in the crust
    ! from G. Tobie thesis
    use param_y, only : R,htide
    use variables
    integer :: i

    Pcrust = 0
    do i = n_lim(4), n-1
        Pcrust = Pcrust + 4._dp*pi*((R(i+1)+R(i))/2.D0)**2*(REALPART(htide(i+1))+REALPART(htide(i)))/2.D0*(R(i+1)-R(i))
        enddo
    do i= n_lim(4)-1, n
        Pcrust = Pcrust + 4._dp*pi*((R(i-1)+R(i))/2.D0)**2*(HTIDE(i-1)+HTIDE(i))/2.D0*(R(i)-R(i-1))
        enddo
    Pcrust = Pcrust/2
    
    Ptide_mant = 0
    do i= n_lim(1), n_lim(4)-1
        Ptide_mant = Ptide_mant + 4._dp*pi*((R(i+1)+R(i))/2.D0)**2*(REALPART(htide(i+1))+REALPART(htide(i)))/2.D0*(R(i+1)-R(i))
        enddo
    do i= n_lim(1)-1, n_lim(4)
        Ptide_mant = Ptide_mant + 4._dp*pi*((R(i-1)+R(i))/2.D0)**2*(HTIDE(i-1)+HTIDE(i))/2.D0*(R(i)-R(i-1))
        enddo
    Ptide_mant = Ptide_mant/2
    
END SUBROUTINE power_crust





SUBROUTINE radioactive
    ! Calculate dissipated nrj
    ! by radioactivity Hussmann2010,
    ! abundance C in Laneuville2013
    use variables
    implicit none
    integer :: i
    real(kind=dp), dimension(4) :: C,H,compo,tau
    real(kind=dp) :: lambda

    compo(1) = 0.992745_dp !U238
    compo(2) = 7.2d-3 !U235
    compo(3) = 1._dp !Th
    compo(4) = 1.17d-4 !K

    H(1) = 9.48d-5
    H(2) = 5.69d-4
    H(3) = 2.69d-5
    H(4) = 2.92d-5
    
    !Laneuville2013
	C(1) = 25.1d-9
	C(2) = 25.1d-9
	C(3) = 25.1d-9*3.7
	C(4) = 25.1d-9*2500._dp

    tau(1) = 4.468d9*year
    tau(2) = 0.7038d9*year
    tau(3) = 14.05d9*year
    tau(4) = 1.277d9*year
    
    ! Massique Radioactivity Power
    Prad_mass = 0
    do i = 1,4
        lambda  = log(0.5)/tau(i)
        Prad_mass = Prad_mass + compo(i)*C(i)*H(i)*exp(lambda*(time-tpr))
        enddo
    Prad_mant = Prad_mass * 4/3*pi*(Rlim(4)**3-Rlim(1)**3)*rho_lay(3)
    Prad_crust = Prad_mass * 4/3*pi*(Rlim(5)**3-Rlim(4)**3)*rho_lay(5)
    Prad_core = Prad_mass * 4/3*pi*(Rlim(1)**3)*rho_lay(1)
END SUBROUTINE radioactive



