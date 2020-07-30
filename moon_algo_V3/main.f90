! *************************
! Main program to start the
! moon evolution algorithm
! *************************

MODULE prec
   integer, parameter :: dp = SELECTED_REAL_KIND(8,3)
END MODULE prec

MODULE variables
    use prec
    !Visco
    integer :: visco
    real(kind=dp) :: alpha_gam,drmag,mean_melt
    !Melt
    integer :: melt
    real(kind=dp) :: phi_melt
	!Geometry
    real(kind=dp) :: rtot,rcore,rast,rcrust,rmag
    real(kind=dp) :: ecc,period,w_in
    !Internal structure
    real(kind=dp) :: vsmant, vpmant, etamant, mu_mant, K_mant, rhomant,lambda_mant
    real(kind=dp) :: vscrust, vpcrust, etacrust_base, mu_crust, K_crust, rhocrust
    real(kind=dp) :: vscore, vpcore, etacore, mu_core, K_core, rhocore
    real(kind=dp) :: etamag
    real(kind=dp) :: enj,alpha_mag
    !Layers
    integer :: n,n_in,ilay
    integer, parameter :: nlay=5
    integer, dimension(nlay) :: n_layer(nlay) ,n_lim(nlay)
    real(kind=dp), dimension(nlay) :: Rlim, rho_lay, K_lay, vp_lay
    real(kind=dp), dimension(nlay) :: vs_lay, eta_lay, mu_lay,cp_lay,lambda_lay
    real(kind=dp), allocatable :: r_in(:),rho_in(:),eta_in(:),vs_in(:),vp_in(:),D_in(:)
    !Tides & Power
    integer :: nloops
    real(kind=dp) :: zeta
    complex(kind=dp) :: klove2
    real(kind=dp), dimension(2) :: klovePlan
    real(kind=dp) :: Ptide_tot,Ptide_crust,Ptide_mant, Ptide_ast
    real(kind=dp) :: Prad_mass, Prad_mant,Prad_crust,Prad_core
    !Cristallisation
	real(kind=dp) :: Tf,dcrust,dmag,Vocean_ini,f,dcrust_rate,DTocean
    !Conduction
    integer :: ncond
    character(len=15) :: mode
    real(kind=dp), allocatable :: T(:)
    real(kind=dp) :: timestep,time, Tsurf, Tbase,DTbase,DTsurf,Fconv
    real(kind=dp) :: Pbase,Psurf,Cpcrust,Tsolidus
    !Orbital param
	real(kind=dp) :: planMass,satMass,planRad,satRad,satMoIr,planMoIr,planQ_1,satQ_1,planSpin
	REAL(KIND=dp) :: a,xincl,beta,spin
    !Physic constants
    real(kind=dp), parameter :: pi = acos(-1._dp)
    real(kind=dp), parameter :: Rgaz =  8.31446261815324_dp
    real(kind=dp), parameter :: Mmoon =  7.34767309d22
	real(kind=dp), parameter :: Mearth =  5.972d24
	real(kind=dp), parameter :: Rearth =  6371d3
    real(kind=dp), parameter :: gmoon = 1.625_dp
    real(kind=dp), parameter :: year =  365._dp*86400._dp
    real(kind=dp), parameter :: tpr = 4.6d9*year !radioactiv
    real(kind=dp), parameter :: Grav   = 6.67384e-11

END MODULE variables


!MAIN PROGRAM
PROGRAM planet_profile
	 
	use variables
	implicit none
    integer :: i,j

    !Initialising program
    CALL Import_param
    CALL rewrite_files

    allocate(r_in(n),rho_in(n),eta_in(n),vs_in(n),vp_in(n),D_in(n),T(n))
	CALL internal_structure
    CALL initialise_temp
	CALL write_Tcrustprofile

    !Main Loop
    do i = 1,nloops
		time = time + timestep
        CALL rheology
		CALL write_rlim
		CALL write_radius
        CALL tides(r_in,rho_in,vp_in,vs_in,eta_in,n,w_in)
        CALL get_htide
        CALL conduction
		CALL radioactive
		CALL cristallisation
		CALL call_orbit

        !Writing options
        Call write_structure
        CALL write_Hprofile
        CALL write_Tcrustprofile
        CALL power_crust
        enddo

   ! CALL write_power
    CALL write_radius
    !CALL radioactive
    !CALL cristallisation

    deallocate(r_in,rho_in,eta_in,vs_in,vp_in,D_in,T)

END PROGRAM planet_profile



!****************SUBROUTINES**********************


SUBROUTINE Import_param
    !Import parameters
    use variables
    implicit none
    character(len = 15) :: junk
    integer :: i
    !Import rheological parameters
    open(11,file='input_parameters.in')
	read(11,*) junk !Orbit parameters
    read(11,*) junk,nloops
    read(11,*) junk,timestep
    read(11,*) junk,a
    read(11,*) junk,ecc
	read(11,*) junk !SatGeometry
    read(11,*) junk,rtot
    read(11,*) junk,rcore
    read(11,*) junk,rast
    read(11,*) junk,rmag
    read(11,*) junk,rcrust
	read(11,*) junk,n_layer(1),n_layer(2),n_layer(3),n_layer(4),n_layer(5)
	read(11,*) junk ! Plan/Sat ppt
	read(11,*) junk, planRad
	read(11,*) junk, planMass
	read(11,*) junk, satMass
	read(11,*) junk, planSpin
	read(11,*) junk, satMoIr
	read(11,*) junk, planMoIr
	read(11,*) junk, planQ_1
	read(11,*) junk, klovePlan(2)
	read(11,*) junk ! Sat Rheology
	read(11,*) junk, etacrust_base
	read(11,*) junk, etamant
	read(11,*) junk, etacore
	read(11,*) junk, K_crust
	read(11,*) junk, K_mant
	read(11,*) junk, K_core
	read(11,*) junk, mu_crust
	read(11,*) junk, mu_mant
	read(11,*) junk, mu_core
	read(11,*) junk, lambda_mant
	read(11,*) junk, Cpcrust
	read(11,*) junk !misc
	read(11,*) junk, Tsurf
	read(11,*) junk, Tbase
	read(11,*) junk, Enj
    
    !Assigning other variables
    !Planet
	satRad = Rtot
    ! period from demi grand axe
    a = a * Rearth
    period = 2._dp*pi*SQRT(a**3/(Grav*(planMass+satMass)))
    w_in = 2._dp*pi/period
    !Layers
    !Magma and crust must have same grid
    drmag = rmag/real(n_layer(4))
    n_layer(5) = nint(rcrust/drmag)
    n_layer(3) = nint(rast/drmag)
	n = n_layer(1)+n_layer(2)+n_layer(3)+n_layer(4)+n_layer(5)
	!conduction
	ncond = n_layer(5)+1
	Vocean_ini = 4._dp/3._dp*pi*(Rtot**3-(Rtot-1000.e3)**3)
	DTocean = 0._dp
    !timestep
    time = 0
    timestep = timestep*1.d6*year
    ! /!\ alpha_gam shouldnt be defined here
    alpha_gam = 0.30_dp
   
END SUBROUTINE Import_param



SUBROUTINE rewrite_files
	use variables, only: n
	
	open(10,file='Outputs/melt.txt')
	write(10,*) "melt"
	close(10)
	
	open(11,file='Outputs/structure.txt')
	write(11,*) n,"points, structure"
	close(11)
	
	open(12,file='Outputs/H_profile.txt')
	write(12,*) n
	close(12)

    open(16,file="Outputs/Tcrust_profile.txt")
    write(16,*) "Tcrust_profile"
    close(16)
    
    open(17,file="Outputs/radius.txt")
    write(17,*) "Radius conduction"
    close(17)
    
    open(20,file='Outputs/rlim_evol.txt')
    write(20,*) "Rlim4", " Time", " ecc", " a", " kreal"," kimag", "test"
    close(20)
END SUBROUTINE rewrite_files




