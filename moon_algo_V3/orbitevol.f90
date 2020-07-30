MODULE orbit
!-----------------------------------------------------------------------------------------------------------------------------------
!   Contains all the procedure and data/variables used for the computation of (tidally-driven) orbital evolution
!   BE CAREFUL: all the orbital model is written in a dimensional framework!
!   Note: all (default) quantities are in SI units (i.e. mks), except in 'calc_orbitalEvolution' where time is in Earth year or day
!   (and angles in °) and in 'param.in'.
!-----------------------------------------------------------------------------------------------------------------------------------

!   MODULES USED
use variables

! Variables
IMPLICIT NONE
PUBLIC
INTEGER,       PARAMETER :: lmax=2                         !max. degree in Kaula's development (for effect of dissipation in planet)
REAL(KIND=dp), PARAMETER :: e_lowlimit=1e-6_dp             !perhaps machine-dependent => TO IMPROVE!!
REAL(KIND=dp), DIMENSION(lmax) :: kloveSat
REAL(KIND=dp)   :: Porb, Pspin  ,Pdiss               !semi-maj. axis, orbital period of the sat., spin period of the planet
REAL(KIND=dp) :: planPower=0._dp, satplanPow=0._dp,satplanL=0._dp !dissipated power inside the planet, total power of the planet-sat system


!Private procedures of the module
PRIVATE :: calc_aeiDeriv,F_spinDeriv,F_obliqDeriv, F_Flmp, F_kn, F_Glpq

CONTAINS

SUBROUTINE calc_orbitalEvolution(e, xn, beta, xincl, spin, timestep)
	!   Computes the orbital evolution of the satellite around its parent planet: gives new values of semi-major axis, eccentricity
	!   ('e'), obliquity of the central body ('beta'), inclination of the satellite's orbit ('xincl') and central planet spin ('spin',
	!   in rad/s), after a [long-term] evolution of 'timestep' seconds.
	!   NOTE: only 1:1 spin-orbit resonance (for now)

	! Arguments
	REAL(KIND=dp), INTENT(INOUT) :: e, xn, beta, xincl, spin
	REAL(KIND=dp), INTENT(IN)    :: timestep
	! Local variables
	INTEGER       :: i, nbtimestep
	REAL(KIND=dp) :: aOld,   eOld,      xinclOld,  dt2, finaldt, dadtPlan, dedtPlan, dadt, dedt, dxincldt, dspin_dt, dbeta_dt, &
					 powOrb, satpowRot, planpowRot, xLOrb, planLrot, satLrot,dtorb,eps
	
	dtorb = 40e3
	eps = EPSILON(1._dp)
	! Main block
	finaldt    = MOD(timestep,(dtorb*year))/year                         !'nbtimstep+1'-th (i.e. final) time step value
	nbtimestep = FLOOR(timestep/(dtorb*year)) + COUNT([finaldt>eps])    !total number of (in years-)timesteps

	!dt2 = dtorb
	dt2 = timestep/year
	DO i=1,1 !nbtimestep
	   !IF ((i == nbtimestep).AND.(finaldt>eps))  dt2 = finaldt
	   aOld     = a
	   eOld     = e
	   xinclOld = xincl
	   CALL calc_aeiDeriv(e,a, xn, xincl+beta, spin, dadtPlan, dedtPlan, dadt, dedt, dxincldt, planPower)

	   !evolution of orbital elements+spin:
	   a     = a     + dt2 * year*dadt
	   e     = e     + dt2 * year*dedt
	   xincl = xincl + dt2 * year*dxincldt/pi/180._dp   !!Ross (1989)  inclination with respect to inertial plane (e.g. ecliptic)
	   dspin_dt  = F_spinDeriv(aOld, eOld, xinclOld+beta, dadtPlan, dedtPlan, dxincldt)
	   dbeta_dt  = F_obliqDeriv(aOld, eOld, xinclOld, spin, beta, dadtPlan, dedtPlan, dxincldt, dspin_dt)
	   spin     = spin  + dt2 * year*dspin_dt
	   beta     = beta  + dt2 * year*dbeta_dt

	   !IF (e < e_lowlimit) e = e_lowlimit    !safety: TO IMPROVE!
	   
	   !secondary quantities:
	   xn    = xn*(aOld/a)**1.5_dp    !mean motion (Kepler's 3rd law)
	   Porb  = pi/(xn*43200._dp)      !revolution period (in Earth days) evolution
	   Pspin = pi/(spin*43200._dp)    !spin period (in Earth days) evolution

	   !powers (used to check conservation of energy: satplanPow = -(Pdiss+planPower)
	   powOrb     = .5_dp*Grav*(planMass+satMass)*satMass*dadt/a**2
	   satpowRot  = -1.5_dp*(satMoIr*satMass*satRad**2)*Grav*(satMass+planMass)*dadt/a**4
	   planpowRot = (planMoIr*planMass*planRad**2)*spin*dspin_dt
	   satplanPow = powOrb + satpowRot + planpowRot    !saved every convective timestep (see 'write_results_means' in main program)

	   !angular momenta (projection on the normal to the inertial plane; satplanL must be conserved):
	   xLOrb    = planMass*satMass*SQRT(Grav/(planMass+satMass)*a*(1._dp-e**2))
	   planLrot = spin*planMoIr*planMass*planRad**2
	   satLrot  = xn*satMoIr*satMass*satRad**2    !1:1 spin-orbit resonance
	   satplanL = (xLOrb+satLrot)*COS(xincl*pi/180._dp) + planLrot*COS(beta*pi/180._dp)    !saved every convective timestep
													
	END DO
END SUBROUTINE calc_orbitalEvolution


SUBROUTINE calc_aeiDeriv(e,a, xn, xincl, spin, da_dtPlan, de_dtPlan, da_dt, de_dt, di_dt, planPow)
	!Computes time derivatives (in qtt/second) of the semi-major axis, eccentricty ('e') and orbit inclination of the satellite
	!('xincl') RELATIVE TO THE CENTRAL PLANET'S EQUATOR, as a balance of tidal dissipation in the satellite and in its parent
	!planet. 'xn' is the mean motion, and 'spin' the spin angular frequency. Also yields the dissipated power inside the planet
	!('planPow').
	!RESULTS: 'da_dt', 'de_dt', 'di_dt' and the two components 'da_dtPlan' and 'de_dtPlan' needed to compute the spin evolution.

	! Arguments
	REAL(KIND=dp), INTENT(IN)  :: e,a, xn, xincl, spin
	REAL(KIND=dp), INTENT(OUT) :: da_dtPlan, de_dtPlan, da_dt, de_dt, di_dt, planPow
	! Local variables
	INTEGER            :: k, l, m, p, q
	REAL(KIND=dp)      :: xirad, Eorbit, afact, efact, da_dtSat, de_dtSat, dadtI_lmpq,dedtI_lmpq,dadtII_lmpq, dedtII_lmpq, didt_lmpq
	REAL(KIND=dp)      :: fact, Flmp, Glpq, Klmpq, wlmpq, sinEpslmpq_I,sinEpslmpq_II
	
	! Main block
	xirad  = xincl*pi/180._dp    !inclination in radians
	Eorbit = -.5_dp*Grav*(planMass+satMass)*satMass/a
	
	!Effect of dissipation
	da_dtSat = 0._dp
	de_dtSat = 0._dp
	!IF (e > e_lowlimit) THEN
		!Planete
		de_dtPlan = e*(57._dp/8._dp)*klovePlan(2)*PlanQ_1*(satMass/PlanMass)*(planRad/a)**5*xn
		da_dtPlan = a*(8._dp/19._dp)*(1+51._dp*e**2/4._dp)*(de_dtPlan/e)
		!Satellite
		de_dtSat = -e*(21._dp/2._dp)*kloveSat(2)*satQ_1*(planMass/satMass)*(satRad/a)**5*xn
		da_dtSat = a*2*e**2*de_dtSat/e
		!END IF
		
	!print *,"deplan=", de_dtPlan,"desat=",de_dtSat
	!print *,"daplan=", da_dtPlan,"dasat=",da_dtSat
	da_dt = da_dtSat + da_dtPlan
	de_dt = de_dtSat + de_dtPlan
	di_dt = 0 !didt_lmpq
		
END SUBROUTINE calc_aeiDeriv



FUNCTION F_Flmp(l, m, p, i)
	!Computes Kaula's inclination function (cf. Kaula, 1961) for
	!given degree 'l', angular order 'm' and 'p' and inclination 'i';

	! Arguments and function
	INTEGER,       INTENT(IN) :: l,m,p
	REAL(KIND=dp), INTENT(IN) :: i
	REAL(KIND=dp) :: F_Flmp
	! Local variables
	INTEGER       :: t, s, c, tmax, smax, cmin, cmax, k
	REAL(KIND=dp) :: Flmp_l, Flmp_s

	! Main block
	IF ( 2*p<=l-m ) THEN
	   tmax = p
	ELSE
	   tmax = INT(0.5_dp*REAL(l-m,KIND=dp))
	END IF
	smax = m
	F_Flmp = 0._dp

	DO t=0,tmax
	   DO s=0,smax
		  IF ( m+t >= p+s ) THEN
			 cmin = 0
		  ELSE
			 cmin = p-t-m+s
		  END IF

		  IF ( m+p+t >= l+s ) THEN
			 cmax = l-m-2*t+s
		  ELSE
			 cmax = p-t
		  END IF

		  Flmp_l   = REAL(PRODUCT([(k,k=2,2*l-2*t)]),KIND=dp) / ( REAL(PRODUCT([(k,k=2,t)]),KIND=dp)                        &
				   * REAL(PRODUCT([(k,k=2,l-t)]),KIND=dp) * 2._dp**(2*(l-t))*REAL(PRODUCT([(k,k=2,l-m-2*t)]),KIND=dp) )     &
				   * SIN(i)**(l-m-2*t)
		  Flmp_s   = F_kn(s,m)*COS(i)**s

		  DO c=cmin,cmax
			 F_Flmp = F_Flmp + Flmp_l*Flmp_s*F_kn(c,l-m-2*t+s)*F_kn(p-t-c,m-s)*(-1._dp)**(c-INT(0.5_dp*REAL(l-m,KIND=dp)))
		  END DO
	   END DO
	END DO
END FUNCTION F_Flmp



FUNCTION F_kn(k,n)
	!Gives a real quantity necessary to compute Kaula's
	!inclination function Flmp (see function 'F_Flmp' just above);
	
	! Arguments and function
	INTEGER, INTENT(IN) :: k, n
	REAL(KIND=dp)       :: F_kn
	! Local variables
	INTEGER :: j

	! Main block
	IF (n >= 0) THEN
	   IF (n-k<0) THEN    !caution for n-k<0!!
		  F_kn = 0._dp
	   ELSE
		  F_kn = REAL(PRODUCT([(j,j=2,n)]),KIND=dp) / REAL(PRODUCT([(j,j=2,k)])*PRODUCT([(j,j=2,n-k)]),KIND=dp)
	   END IF
	ELSE
	   F_kn = REAL(PRODUCT([(j,j=2,-n+k-1)]),KIND=dp) / REAL(PRODUCT([(j,j=2,k)])*PRODUCT([(j,j=2,-n-1)]),KIND=dp)*(-1._dp)**k
	END IF
END FUNCTION F_kn


FUNCTION F_Glpq(l,p,q,e)
	! Computes Kaula's eccentricity function (cf. Kaula, 1961;
	! method of Szeto and Lambeck (1980)) for given degree 'l',
	! angular order'p' and 'q' and eccentricity 'e';

	! Arguments and function
	INTEGER,       INTENT(IN) :: l, p, q
	REAL(KIND=dp), INTENT(IN) :: e
	REAL(KIND=dp)             :: F_Glpq,eps
	! Local variables
	INTEGER       :: p2, q2, h, n, k
	REAL(KIND=dp) :: b

	eps = EPSILON(1._dp)
	! Main block
	b = e/(1._dp+SQRT(1._dp-e**2._dp))
	F_Glpq = 0._dp
	IF ( 2*p <= l ) THEN
	   p2 = p
	   q2 = q
	ELSE
	   p2 = l-p
	   q2 = -q
	END IF

	IF ( abs(e)<eps ) THEN    !caution for e=0
	   IF ( q == 0 ) THEN
		  F_Glpq = 1._dp
	   ELSE
		  F_Glpq = 0._dp
	   END IF
	ELSE
	   IF ( l-2*p+q == 0 ) THEN
		  DO k=0,p2-1
			 F_Glpq = F_kn(2*k+l-2*p2,l-1)*F_kn(k,2*k+l-2*p2)*(e/2._dp)**(2*k+l-2*p2)
		  END DO
		  F_Glpq = F_Glpq/((1._dp-e**2)**(REAL(l,KIND=dp)-0.5_dp))
	   ELSE
		  DO k=0,lmax
			 IF ( q2 > 0 ) THEN
				h = k+q2
				n = k
			 ELSE
				h = k
				n = k-q2
			 END IF
			 F_Glpq = F_Glpq + F_Plpqk(l,p2,q2,h,e)*F_Qlpqk(l,p2,q2,n,e)*b**(2*k)
		  END DO
		  F_Glpq = F_Glpq*(-b)**(ABS(q))*(1._dp+b**2)**l
	   END IF
	END IF


	CONTAINS

	FUNCTION F_Plpqk(l, p2, q2, h, e)
		!   Gives a real quantity necessary to compute Kaula's
		!   eccentricity function Glpq (see function 'F_Glpq' just above);
		INTEGER, INTENT(IN) :: l, p2, q2, h
		REAL(KIND=dp), INTENT(IN) :: e
		REAL(KIND=dp) :: F_Plpqk

		! Local variables
		INTEGER       :: i, r
		REAL(KIND=dp) :: b

		! Main block
		b = e/(1._dp+SQRT(1._dp-e**2))
		F_Plpqk = 0._dp
		DO r=0,h
		   F_Plpqk = F_Plpqk + F_kn(h-r,2*p2-2*l)*(-1._dp)**r/REAL(PRODUCT([(i,i=2,r)]),KIND=dp)*(REAL(l-2*p2+q2,KIND=dp)*e/(2._dp*b))**r
		END DO
	END FUNCTION F_Plpqk


	FUNCTION F_Qlpqk(l, p2, q2, n, e)
		! Arguments and function
		INTEGER,       INTENT(IN) :: l, p2, q2, n
		REAL(KIND=dp), INTENT(IN) :: e
		REAL(KIND=dp)             :: F_Qlpqk

		! Local variables
		INTEGER       :: i, r
		REAL(KIND=dp) :: b

		! Main block
		b = e/(1._dp+SQRT(1._dp-e**2))
		F_Qlpqk = 0._dp
		DO r=0,n
		!!O   F_Qlpqk = F_Qlpqk + F_kn(n-r,-2*p2)/PRODUCT([(i,i=2,r)])*(REAL(l-2*p2+q2,KIND=dp)*e/(2._dp*b))**r
		   F_Qlpqk = F_Qlpqk + F_kn(n-r,-2*p2)/REAL(PRODUCT([(i,i=2,r)]),KIND=dp)*(REAL(l-2*p2+q2,KIND=dp)*e/(2._dp*b))**r
		END DO
	END FUNCTION F_Qlpqk

END FUNCTION F_Glpq



FUNCTION F_spinDeriv(a, e, xi, da_dtPlan, de_dtPlan, di_dt)
	!Computes the central planet's spin derivative (conservation of angular momentum) as afunction of current semi-major axis ('a'),
	!eccentricity ('e'), RELATIVE inclination ('xi', i.e. inclination with respect to the equator of the central planet)
	!and their respective time derivatives due to dissipation in the planet ('dadt_Plan', 'dedt_Plan'; 'di_dt' is the derivative of
	!the TRUE inclination).

	! Arguments and function
	REAL(KIND=dp), INTENT(IN) :: a, e, xi, da_dtPlan, de_dtPlan, di_dt
	REAL(KIND=dp)             :: F_spinDeriv
	! Local variables
	REAL(KIND=dp) :: xirad, xn, satMoi, satLspin, satL, fact, dsatL_dt    !'satL' & 'fact'(='satL' time deriv.)  usually negligible!

	! Main block
	xirad  = xi*pi/180._dp    !residual inclination in radians
	xn     = SQRT(Grav*(planMass+satMass)/a**3)
	satMoi = satMoIr*satMass*satRad**2

	satLspin = xn*satMoi    !1:1 spin-orbit resonance
	satL    = planMass*satMass*SQRT(Grav/(planMass+satMass)*a*(1._dp-e**2)) + satLspin

	fact     = 1.5_dp*satMoi*xn*da_dtPlan/a
	dsatL_dt = planMass*satMass*SQRT(Grav/(planMass+satMass))                                  &
			 * ( SQRT((1._dp-e**2)/a)*da_dtPlan*.5_dp - SQRT(a/(1._dp-e**2))*de_dtPlan*e ) -fact
	F_spinDeriv = (satL*di_dt*SIN(xirad)-dsatL_dt*COS(xirad))/(planMoIr*planMass*planRad**2)
END FUNCTION F_spinDeriv


FUNCTION F_obliqDeriv(a, e, xincl, spin, obliq, da_dtPlan, de_dtPlan, di_dt, dspin_dt)
	!   Computes the central planet's obliquity derivative (conservation of angular momentum) as a function of current semi-major axi
	!   ('a'), eccentricity ('e'), TRUE inclination ('xi'; with resepct to inertial plane), spin ('spin') and obliquity ('obliq') of the
	!   central planet and the associated derivatives ('da_dtPlan', 'de_dtPlan', 'di_dt', and 'dspin_dt')

	! Arguments and function
	REAL(KIND=dp), INTENT(IN) :: a, e, xincl, spin, obliq, da_dtPlan, de_dtPlan, di_dt, dspin_dt
	REAL(KIND=dp)             :: F_obliqDeriv
	! Local variables
	REAL(KIND=dp) :: obliqRad, xinclRad, xn, satMoi, satL, fact, dsatL_dt, planLderiv_C

	! Main block
	obliqRad = obliq*pi/180._dp
	xinclRad = xincl*pi/180._dp

	!Same as in function 'F_spinDeriv': => TO IMPROVE: put 'satL' and 'dsatL_dt' as !global variables!! => avoid double computation!!
	xn           = SQRT(Grav*(planMass+satMass)/a**3)
	satMoi       = satMoIr*satMass*satRad**2
	satL         = planMass*satMass*SQRT(Grav/(planMass+satMass)*a*(1._dp-e**2)) + xn*satMoi
	fact         = 1.5_dp*satMoi*xn*da_dtPlan/a
	dsatL_dt     = planMass*satMass*SQRT(Grav/(planMass+satMass))                                  &    !Kaula (1964),
				 * ( SQRT((1._dp-e**2)/a)*da_dtPlan*.5_dp - SQRT(a/(1._dp-e**2))*de_dtPlan*e ) -fact    !eq.(59) [corrected]

	planLderiv_C = (satL*di_dt*SIN(xinclRad)-dsatL_dt*COS(xinclRad))/(planMoIr*planMass*planRad**2)
    F_obliqDeriv = (dspin_dt/TAN(obliqRad) - planLderiv_C/SIN(obliqRad))/spin
END FUNCTION F_obliqDeriv


END MODULE orbit
