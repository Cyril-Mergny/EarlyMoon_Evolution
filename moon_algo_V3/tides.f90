!**********************************************************************
!Calcul des fonctions radiales yi de Takeushi et Saito 1972 
!Inclure les modules prec,param_y et interf dans le programme principal
!Mathieu Le Feuvre/Gabriel Tobie
!**********************************************************************

MODULE prec
   integer, parameter :: dp = SELECTED_REAL_KIND(8,3)
END MODULE prec

MODULE param_Y
   use prec
   integer, parameter :: ncmax=25000,ndmax=20,nrk=1000
   real(kind=dp):: w,w0, wmin,wmax,wprec,Tmin,Tmax,Tprec,Ts, ecc
   integer :: nw, acc, n
   real(kind=dp):: O
   real(kind=dp):: GG
   real(kind=dp):: pi, period
   real(kind=dp), parameter::epsilon=1e-5
   real(kind=dp):: rn,rhon,vn,wn, gcore, alpha_c, rhoatm,  Rsat
   real(kind=dp):: Mtot, Rtot, rho_moy, Rmant, Rcore, Rast, MOI, Itot, Mrock, Irock
   real(kind=dp):: deep, colatt, lont, dcolatt,dlont, surft,  Mcore, Msat 
   real(kind=dp):: Det,Det1,Det2,wdet
   real(kind=dp), dimension(ncmax) :: r,rho,eta,vs,vp,T,P,d,mu,Ka,g,mass,colat,lon
   complex(kind=dp), dimension(ncmax) :: lamc,muc,vsc,vpc
   complex(kind=dp) :: I1,I2,I3,Y1S
   integer :: nc,ndis,visco,ord,ndep,nsl,ordmin,ordmax, iprofile, nconv, melt
   integer, dimension(ndmax) :: pdis,psl
   integer :: forc
   character(100) :: ficmod   
   complex(kind=dp) :: Fh2, mubar, h2ana, k2ana, h2anaG
   complex(kind=dp), dimension(6,3) :: Yij
   complex(kind=dp), dimension(6,3,ncmax) ::  Yijr
   complex(kind=dp), dimension(6,ncmax) ::  Yi
   integer :: ni, ncore,nmant, ncrust, nast
   complex(kind=dp) :: dy1dr(ncmax),dy5dr(ncmax),&
   dI1(ncmax),dI2(ncmax),HH(4,ncmax), &
   strain_energy(ncmax), strain_energy_moy(ncmax), htide(ncmax)
   integer :: maxRK
   real(kind=dp):: precRK
 
   integer :: idum
END MODULE

MODULE param_B
   use prec
   real(kind=dp) :: EdJdT,Ealphq,Etaul,Etauh,Etaum,Ema,Emj,Emv,Ee,Ev 
   real(kind=dp) :: tr,dr,dJdt,mj,alphq,taul,tauh,taum,ma,mv,e,v
   real(kind=dp), parameter :: rgaz=8.314472_dp,trF=950._dp+273._dp,drF=1.e-5_dp,dJdtF=9.1e-4_dp,&
   mjF=0.16_dp,alphqF=0.27_dp,taulF=3.981e-3_dp,tauhF=5.26e6_dp,&
   taumF=4.31e6_dp,maF=1.09_dp,mvF=2.1_dp,eF=5.05e5_dp,vF=1.2e-5_dp,&
   deltaF=1.4_dp
END MODULE param_B


MODULE param_rheo
   use prec
   real(kind=dp), parameter :: Ea=50.e3_dp, vis_ref=1.e15_dp, T_ref=255._dp
   real(kind=dp), parameter :: mu_E=3.3e9_dp, K_E=10.e9_dp
   real(kind=dp) ::  vis_temp, bocean, bice, vsmant, vpmant, etamant, rhomant, bmantsup
   real(kind=dp) :: vscrust, vpcrust, etacrust, rhocrust, bcrust, rhocore, rhoocean, etaice, bast
   real(kind=dp) :: mu_core, K_core, etacore, K_mant, mu_mant, mu_ocean, K_ocean, etaocean
   real(kind=dp) :: K_crust, mu_crust
   real(kind=dp) :: rhoice, mu_ice, K_ice, rho_rock
   real(kind=dp) :: A_0, B_0, A_m, A_c, B_m, B_c
   real(kind=dp), parameter :: Rg=8.314472d0
   complex(kind=dp) :: Jcomp

    !Melt
    integer, parameter :: nmelt=101
    integer :: imelt
    real(kind=dp) :: phi_melt
    real(kind=dp) :: mu_mant_ast, etamant_ast

    !Internal structure
    integer, parameter :: nlay=4
    real(kind=dp), dimension(nlay) :: Rlim, rho_lay, K_lay, vp_lay, vs_lay, eta_lay, mu_lay
END MODULE

MODULE interf
   implicit none
   INTERFACE       
       Function BURGER(ww,TT,PP,dd)
       use prec
       real(kind=dp), intent(in) :: ww,TT,PP,dd
       complex(kind=dp) :: BURGER
       End Function BURGER
       
       Function Phi(m,x)
       use prec
       real(kind=dp), intent(in) :: m
       complex(kind=dp), intent(in) :: x
       complex(kind=dp) :: Phi
       End Function Phi
      
       Function Psi(m,x)
       use prec
       real(kind=dp), intent(in) :: m
       complex(kind=dp), intent(in) :: x
       complex(kind=dp) :: Psi
       End Function Psi
         
!       Subroutine Dy(k,rr,yy,dYdr)
!       use prec
!       integer, intent(in) :: k
!       real(kind=dp), intent(in) :: rr
!       complex(kind=dp), dimension(6,3), intent(inout) :: yy
!       complex(kind=dp), dimension(6,3), intent(out) :: dYdr
!       End subroutine Dy

       Function gr(r,k)
       use prec
       real(kind=dp), intent(in) :: r
       integer, intent(in) :: k
       real(kind=dp) :: gr
       End function gr
                   
      SUBROUTINE trapzd(x,y,ia,ib,s,ordre)
      use prec
      real(kind=dp), dimension(101), intent(in) :: x,y
      integer, intent(in) :: ordre,ia,ib
      real(kind=dp), intent(out) :: s
      END SUBROUTINE trapzd 
      
    SUBROUTINE Integ(k,ystart)
    use prec
    integer, INTENT(IN) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: ystart
	END SUBROUTINE Integ
	
	SUBROUTINE rkqs(k,yy,dydx,x,htry,eps,yscal,hdid,hnext)
	use prec
	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: yy
	complex(kind=dp), DIMENSION(6,3), INTENT(IN) :: dydx,yscal
	real(kind=dp), INTENT(INOUT) :: x
	real(kind=dp), INTENT(IN) :: htry,eps
	real(kind=dp), INTENT(OUT) :: hdid,hnext
	END SUBROUTINE rkqs

   	SUBROUTINE rkck(k,yy,dydx,x,h,yout,yerr)
   	use prec
   	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(IN) :: yy,dydx
	real(kind=dp), INTENT(IN) :: x,h
	complex(kind=dp), DIMENSION(6,3), INTENT(OUT) :: yout,yerr
	END SUBROUTINE rkck
	
	SUBROUTINE RK4(k,yy)
    use prec
	integer, INTENT(IN) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: yy	
	END SUBROUTINE rk4
	
	FUNCTION invM(X1)
    use prec
    complex(kind=dp), intent(in), dimension(:,:) :: X1
    complex(kind=dp), dimension(size(X1,1),size(X1,1)) :: invM
    END FUNCTION invM
    
    function ran(idum)
    INTEGER, INTENT(INOUT) :: idum
    REAL :: ran
    end function ran
              
    FUNCTION gammln(xx) 
    
    use prec 
    REAL(kind=dp), intent(in) :: xx
    REAL (kind=dp) :: gammln

    END FUNCTION gammln   

    FUNCTION J_Maxwell(mut,vis,w_forc) 

     use prec 
     REAL(kind=dp), intent(in) :: mut,vis, w_forc
     complex(kind=dp) :: J_Maxwell
    
    END FUNCTION J_Maxwell


    FUNCTION J_Andrade(mut,vis,w_forc) 

     use prec     
     REAL(kind=dp), intent(in) :: mut,vis, w_forc
     complex(kind=dp) :: J_Andrade
     
    END FUNCTION J_Andrade

   END INTERFACE
END MODULE interf


SUBROUTINE tides(r_in,rho_in,vp_in,vs_in,eta_in,n_in,w_in) !,k2_out,Dt_out,Q_out)
   use param_Y
   use prec
   use interf
   use variables, only : nlay,n_lim,ilay,klove2,satQ_1
   
   implicit none
   integer:: n_in, i,n_
   real(kind=dp),dimension(n_in) :: r_in, rho_in, vp_in, vs_in, eta_in
   real(kind=dp),dimension(n_in+nlay) :: r_, rho_, vp_, vs_, eta_
   real(kind=dp) :: Imk2_calc, Imk2_ast, Imk2_ice, Imk2_mant, Qcalc, Qcalc2, w_in, k2_out, Dt_out, Dt_out2, Q_out
   complex(kind=dp) :: k2calc


do i = 1,n_in
r_(i)=r_in(i)
rho_(i)=rho_in(i)
vp_(i)=vp_in(i)
vs_(i)=vs_in(i)
eta_(i)=eta_in(i)
enddo
n_=n_in

!Makes double value at interface
!Important for good transition bw two layers
do ilay = 1,nlay-1
	n_ = n_+1
	do i = n_, n_lim(ilay)+ilay, -1
		r_(i) = r_(i-1)
		rho_(i) = rho_(i-1)
		vp_(i) = vp_(i-1)
		vs_(i) = vs_(i-1)
		eta_(i) = eta_(i-1)
		enddo
	rho_(n_lim(ilay)+ilay) = rho_(n_lim(ilay)+1+ilay)
	vp_(n_lim(ilay)+ilay) = vp_(n_lim(ilay)+1+ilay)
	vs_(n_lim(ilay)+ilay) = vs_(n_lim(ilay)+1+ilay)
	eta_(n_lim(ilay)+ilay) = eta_(n_lim(ilay)+1+ilay)
	enddo



!-- Input parameters-----------------------
 
   Ord=2
! Degree of the perturbing potential

  forc=1
! free oscillation=0 or forced oscillation=1 (tides)


 visco=2
! viscoelastic model assumption: 0=purely elastic; 1=Maxwell; 2=Andrade

 acc=0
! treatment of the liquid layer: 0=quasi-static; 1=dynamic 

 precRK=1.e-10
 maxRK=10000000
!-precision and maximum steps used in the Runge-Kunta scheme

 pi=acos(-1._dp)
 GG=6.67384e-11_dp


!--Define Structure Model 

  w=w_in
  r(1:n_in) = r_
  rho(1:n_in)=rho_
  vp(1:n_in)=vp_
  vs(1:n_in)=vs_
  eta(1:n_in)=eta_
  n=n_


   Call Model_interface

!--Normalize 
   Call Norm(1)
!--Calculate viscolelastic parameters
   Call Viscel   
!--Calculate gravity for each shell        
   call Gravite
!--Calculate yi   
   Call Y6
!--Denormalize 
   Call Norm(-1)   
!--Calculate energy integrals and dissipation
   Call I123

!--Calculate determinant
   Call Deter   

!---Determine k2 and Dt from complex Y functions ------


Imk2_calc=0._dp
Imk2_mant=0._dp
Imk2_ast=0._dp


!Pice=Ptot-Pcrust-Pmant

Imk2_calc=Imk2_calc/2.d0
Imk2_mant=Imk2_mant/2.d0
Imk2_ast=Imk2_ast/2.d0


k2_out = REAL(Yi(5,nc))-1._dp
k2calc = Yi(5,nc)-1._dp



Qcalc=-sqrt(k2calc*conjg(k2calc))/IMAG(Yi(5,nc))
Qcalc2=sqrt(k2calc*conjg(k2calc))/Imk2_calc

Dt_out=asin(1.d0/Qcalc)/w
Dt_out2=asin(1.d0/Qcalc2)/w

Q_out = Qcalc
satQ_1 = 1/Q_out
klove2 = k2calc

END SUBROUTINE tides


SUBROUTINE Model_interface
   use param_Y
   use prec
   implicit none

   integer :: i

!--Define constants and allocate
   O=real(ord)



! Interface position and nature

	nsl=0
	Do i=2,n
	if(vs(i) == 0. .and. vs(i-1) .ne. 0.) nsl=nsl+1
	Enddo

   pdis(:)=0
   psl(:)=0 ! solide/liquide
   ndis=0 ! discontinuite
   nsl=0  
   nc=1     	 

   do nc=1,n
    If(nc>ncmax .or. ndis>ndmax) then
      write(*,*) 'too many shells', nc, ndis
      stop
      return
    Endif
    If (nc>1 .and. r(nc)==r(nc-1)) then
       ndis=ndis+1 
       pdis(ndis)=nc-1
       If (vs(nc)==0. .and. vs(nc-1) .ne. 0.) then
        nsl=nsl+1
        psl(nsl)=nc-1
       Endif
    Endif
    enddo
    nc=nc-1

!----Normalisation constants
     rn=r(n)
     rhon=3000._dp 
     vn=4000._dp
     wn=vn/rn

!     wn=1._dp
!     rn=1._dp
!     rhon=1._dp
!     vn=1._dp

END SUBROUTINE Model_interface


SUBROUTINE Deter
   use param_Y
   use prec
   implicit none

 If(vs(nc)==0.) then
   Det=yijr(2,2,nc)*yijr(6,1,nc)-yijr(2,1,nc)*yijr(6,2,nc)
 Else   
   Det=yijr(2,1,nc)*(yijr(4,2,nc)*yijr(6,3,nc)-yijr(6,2,nc)*yijr(4,3,nc))&
   -yijr(2,2,nc)*(yijr(4,1,nc)*yijr(6,3,nc)-yijr(6,1,nc)*yijr(4,3,nc))&
   +yijr(2,3,nc)*(yijr(4,1,nc)*yijr(6,2,nc)-yijr(6,1,nc)*yijr(4,2,nc))
 Endif   

END SUBROUTINE Deter
   
SUBROUTINE Norm(dir)
   use param_Y
   use prec
   implicit none
   integer, intent(in) :: dir

   If (dir==1) then 
!----Normalisation 
     w=w/wn
     GG=GG*(rhon/wn**2)
     r(:)=r(:)/rn
  	 vp(:)=vp(:)/vn
   	 vs(:)=vs(:)/vn
     rho(:)=rho(:)/rhon 	
     eta(:)=eta(:)/(rhon*vn*rn) 
   
   Else if (dir==-1) then
!----Denormalize
     Yi(1,:)=Yi(1,:)/(wn**2*rn)
     Yi(2,:)=Yi(2,:)*rhon 
     Yi(3,:)=Yi(3,:)/(wn**2*rn)   
     Yi(4,:)=Yi(4,:)*rhon 
!    Yi(5,:)=Yi(5,:)
     Yi(6,:)=Yi(6,:)/rn  

     lamc(:)=lamc(:)*rhon*vn**2
     muc(:)=muc(:)*rhon*vn**2
     mu(:)=mu(:)*rhon*vn**2 
     Ka(:)=Ka(:)*rhon*vn**2 

     GG=GG/(rhon/wn**2)
     r(:)=r(:)*rn
  	 vp(:)=vp(:)*vn
   	 vs(:)=vs(:)*vn
     rho(:)=rho(:)*rhon 	
     eta(:)=eta(:)*(rhon*vn*rn) 
     w=w*wn
     G(:)=G(:)/(rhon/wn**2)/rn**2*rhon*rn**3
   
   Else
   write(*,*) '_Norm : nothing done! '
   Endif
    
END SUBROUTINE Norm

SUBROUTINE Viscel
   use param_Y
   use param_rheo
   use prec
   use interf
   implicit none
   integer :: i
   
   !----Complex variables
     mu(:)=rho(:)*vs(:)**2
     Ka(:)=rho(:)*(vp(:)**2-4._dp/3._dp*vs(:)**2)
     
     DO i=1,nc

     if ((vs(i)==0.d0).or.(mu(i)==0.d0).or.(eta(i)==0.d0)) then 
       muc(i)=0.d0
     else
      if(visco==0) muc(i)=cmplx(mu(i),0.d0,kind=dp)
      if(visco==1) then 
        Jcomp=J_Maxwell(mu(i),eta(i),w)
        muc(i)=cmplx(1.D0,kind=dp)/Jcomp
      endif
      if(visco==2) then 
        Jcomp=J_Andrade(mu(i),eta(i),w)
        muc(i)=cmplx(1.D0,kind=dp)/Jcomp
      endif
     endif
      lamc(i)=Ka(i)-2._dp/3._dp*muc(i)
      vpc(i)=sqrt((lamc(i)+2._dp*muc(i))/rho(i))
      vsc(i)=sqrt(muc(i)/rho(i))


  

     ENDDO      
END SUBROUTINE Viscel

FUNCTION BURGER(ww,TT,PP,dd)
!*****Ecrit par Antoine Mocquet
!*****Rigidité viscoelastique d'apres Faul & Jackson 2005
      use param_B
      use param_Y
      use interf, only: trapzd
      use prec
      implicit none
      real(kind=dp), intent(in) :: ww,TT,PP,dd 
      complex(kind=dp) :: BURGER

      real(kind=dp), dimension(2) :: som
      real(kind=dp), dimension(101) :: ay,by2,denom,xinte

      real(kind=dp), dimension(101,2) :: byr
      real(kind=dp) :: gu,qm1
      real(kind=dp) :: taur,delta
      real(kind=dp) :: cortau,x,coef1,coef2,grain,y0,dlju,taum2,tauh2,s
      integer ::iy,itest,kk,ii
      
      taur  = exp(-e/(rgaz*tr)) 
      
      delta = deltaF*(tauhF/taulF)*(taul/tauh)
 
 !****************Corrections           
   !   ma=ma*2.
   !   alphq=alphq/2.
   !   v = 1.4e-5_dp + Ev
 !**************** 
 
      cortau=exp((e+PP*v)/(rgaz*TT))
      
      If(cortau>1.e200_dp) then
     ! write(*,*) 'Burger Warning, T trop basse'
      cortau=1.e200_dp
      Endif

!--------Le rapport taul/tauh est indépendant de la taille
!--------de grain et des conditions thermodynamiques

         x=taul/tauh
         coef1=alphq*delta/(1._dp-(x**alphq))
         y0=log(x)

         DO iy=1,101
            ay(iy)    = y0*(1._dp-real(iy-1)*1.e-2_dp)
            byr(iy,1) = exp(alphq*ay(iy))
            byr(iy,2) = exp(ay(iy))*byr(iy,1)
            by2(iy)   = exp(2._dp*ay(iy))
         ENDDO
         
         grain=dd/dr

! Anharmonicité anélastique

         dlju=0._dp

         IF(TT .gt. tr) then
         dlju = dJdt*(grain**(-mj))*(TT-tr)
         ENDIF
         
! Correction de la pression et de la température
! et de la taille de grain pour les temps de relaxation

         tauh2=tauh*taur*cortau*(grain**ma)
         taum2=taum*taur*cortau*(grain**mv)

! Calcul des coefficients pré-intégrales

          coef2 = coef1*tauh2

          DO iy=1,101
             denom(iy) = 1._dp+by2(iy)*((ww*tauh2)**2)
          ENDDO

! Integration sur les temps de relaxation

          DO ii=1,2
             som(ii)=0._dp
             itest=0
             DO iy=1,101
                xinte(iy)=byr(iy,ii)/denom(iy)
             ENDDO
             CALL trapzd(ay,xinte,1,101,s,1)
             kk=2
100          CALL trapzd(ay,xinte,1,101,s,kk)
             kk=kk+1
             itest=itest+2**(kk-2)
             If(itest.le.100) go to 100
             If(ii.eq.1) then
                  som(ii)=1._dp+dlju+coef1*s
             Else
                  som(ii)=ww*coef2*s+(1._dp/(ww*taum2))
             Endif
           ENDDO

! Perturbation de la rigidité en GPa
! et de l'attenuation Qmu-1

            gu  = 1._dp/(sqrt((som(1)**2)+(som(2)**2)))
            qm1 = som(2)/som(1)

! Expression complexe de la variation de la rigidité dmu
! muc(w) = mu*dmu avec BURGER=dmu= ( J1/(J1^2 + J2^2) + i J2/(J1^2 + J2^2) )
!gu=module burger , qm1=Im/Re burger
            BURGER=cmplx(som(1)/(som(1)**2+som(2)**2),som(2)/(som(1)**2+som(2)**2))
      
END FUNCTION BURGER

SUBROUTINE trapzd(x,y,ia,ib,s,ordre)
! Calcul d'intégrale à l'ordre "ordre" par la méthode des trapèzes
      use prec
      implicit none
      real(kind=dp), dimension(101), intent(in) :: x,y
      integer, intent(in) :: ordre,ia,ib
      real(kind=dp), intent(out) :: s
      real(kind=dp) :: pas,tnm,del,sum,rere,rr
      integer :: it,k,j

      if(ordre.eq.1) then
         s=0.5_dp*(x(ib)-x(ia))*(y(ia)+y(ib))
      else
         pas = (x(ib)-x(ia))/float(ib-ia)
         it  = 2**(ordre-2)
         tnm = float(it)
         del = (x(ib)-x(ia))/tnm
         rr   = x(ia)+0.5_dp*del
         sum = 0._dp
         
         do j=1,it
            rere = (rr-x(ia))/pas
            k    = ia  + int(rere)
            sum  = sum + y(k)
            rr    = rr   + del
         enddo
         s=0.5_dp*(s+(x(ib)-x(ia))*sum/tnm)
      endif     
END SUBROUTINE trapzd 

  
SUBROUTINE Gravite
     use param_Y
     use prec
     implicit none 
     integer :: i

     mass(1)=0._dp
     G(1)=0._dp

     DO  i=2,nc
         mass(i)=mass(i-1)+(4._dp/3._dp)*pi*(r(i)**3-r(i-1)**3)*rho(i)
         G(i)=GG*mass(i)/(r(i)**2)
     ENDDO
     
     
END SUBROUTINE Gravite  
    
SUBROUTINE Y6
! Calcul des fonctions y1-6 d'après Takeushi et Saito 1972
   use param_Y
   use interf
   use prec
   implicit none      
   integer :: i,j,k    

!++++ Initialisation de la solution au centre ++++ 
   Yij=0._dp
   Yijr=0._dp
   
   ni=2
   CALL Init(ni)
         
   Yijr(:,:,ni)=Yij(:,:)

!++++ Integration du centre vers la surface ++++
   Call UpY

   Yijr(:,:,nc)=Yij(:,:)


!++++ Fonctions propres de la surface vers le centre
   Call Eigen   
      
END SUBROUTINE Y6

!****************************************************
SUBROUTINE Init(k)
! Solution analytique initiale au sommet de la première couche    
      use param_Y
      use prec
      use interf, only: Phi,Psi
      implicit none 
      integer :: i,j,k
      real(kind=dp):: O1
      complex(kind=dp) :: gam,kk,x,f,h

      gam=4._dp/3._dp*pi*GG*rho(k)
      O1=O+1._dp


      IF ((acc==0).and.(vs(k)==0._dp)) THEN

        Yij(:,:)=0._dp
        Yij(5,1)=1._dp
        Yij(6,1)=2._dp*(O-1._dp)/r(k)
!!!!! Yij(6,1) correspond dans ce cas a la variable y7 telle que definie dans Saito (1974) Eq. (17) !!!!!!!

      ELSE

! Solution de Pekeris et Jarosh (1958) pour une sphere homogene compressible
! Troisieme solution (existe uniquement ds un solide)

     If(vs(k)==0._dp) then  !liquide     
     
      Yij(:,3)=0._dp     
      
     ELSE

      kk=(w**2+4._dp*gam)/vpc(k)**2+w**2/vsc(k)**2&
      -sqrt((w**2/vsc(k)**2-(w**2+4._dp*gam)/vpc(k)**2)**2+4._dp*O*(O+1._dp)*gam**2/(vpc(k)**2*vsc(k)**2)) 
      kk=0.5_dp*kk
      kk=sqrt(kk)
     
      f=vsc(k)**2/gam*(kk**2-w**2/vsc(k)**2)
          
      h=f-(O+1._dp)
      x=kk*r(k)
       
      Yij(1,3)=-r(k)/(2._dp*O+3._dp)*(0.5_dp*O*h*Psi(O,x)+f*Phi(O1,x))
      
      Yij(2,3)=-(lamc(k)+2._dp*muc(k))*f*Phi(O,x)&
      +muc(k)/(2._dp*O+3._dp)*(-O*(O-1._dp)*h*Psi(O,x)+2._dp*(2._dp*f+O*(O+1._dp))*Phi(O1,x))
      
      Yij(3,3)=-r(k)/(2._dp*O+3._dp)*(0.5_dp*h*Psi(O,x)-Phi(O1,x))
      
      Yij(4,3)=muc(k)*(Phi(O,x)-1._dp/(2._dp*O+3._dp)*((O-1._dp)*h*Psi(O,x)+2._dp*(f+1._dp)*Phi(O1,x)))
      
      Yij(5,3)=(vpc(k)**2*f-(O+1._dp)*vsc(k)**2)&
      -r(k)**2*3._dp*gam*f/(2._dp*(2._dp*O+3._dp))*Psi(O,x)
      
      Yij(6,3)=1._dp/r(k)*(2._dp*O+1._dp)*Yij(5,3)&
      +r(k)*3._dp*O*gam*h/(2._dp*(2._dp*O+3._dp))*Psi(O,x)      
     ENDIF
     
!    Deuxieme solution
      Yij(1,2)=O/r(k)
      Yij(2,2)=2._dp*muc(k)*O*(O-1._dp)/r(k)**2
      Yij(3,2)=1._dp/r(k)
      Yij(4,2)=2._dp*muc(k)*(O-1._dp)/r(k)**2
      Yij(5,2)=O*gam-w**2
      Yij(6,2)=((2._dp*O+1._dp)*Yij(5,2)-3._dp*O*gam)/r(k)
       
!    Première Solution
   IF (vsc(k)==0._dp) then  
      kk=1._dp/vpc(k)**2*(w**2+4._dp*gam-(O*(O+1._dp)*gam**2)/w**2)
      kk=sqrt(kk)
      f=-w**2/gam
   ELSE       
      kk=(w**2+4._dp*gam)/vpc(k)**2+w**2/vsc(k)**2&
      +sqrt((w**2/vsc(k)**2-(w**2+4._dp*gam)/vpc(k)**2)**2+4._dp*O*(O+1._dp)*gam**2/(vpc(k)**2*vsc(k)**2)) 
      kk=0.5_dp*kk
      kk=sqrt(kk)
           
      f=vsc(k)**2/gam*(kk**2-w**2/vsc(k)**2)            
   ENDIF  
      
      x=kk*r(k)
      h=f-(O+1._dp)
       
      Yij(1,1)=-r(k)/(2._dp*O+3._dp)*(0.5_dp*O*h*Psi(O,x)+f*Phi(O1,x))
      
      Yij(2,1)=-(lamc(k)+2._dp*muc(k))*f*Phi(O,x)&
      +muc(k)/(2._dp*O+3._dp)*(-O*(O-1._dp)*h*Psi(O,x)+2._dp*(2._dp*f+O*(O+1._dp))*Phi(O1,x))
      
      Yij(3,1)=-r(k)/(2._dp*O+3._dp)*(0.5_dp*h*Psi(O,x)-Phi(O1,x))
      
      Yij(4,1)=muc(k)*(Phi(O,x)-1._dp/(2._dp*O+3._dp)*((O-1._dp)*h*Psi(O,x)+2._dp*(f+1._dp)*Phi(O1,x)))
      
      Yij(5,1)=(vpc(k)**2*f-(O+1._dp)*vsc(k)**2)&
      -r(k)**2*3._dp*gam*f/(2._dp*(2._dp*O+3._dp))*Psi(O,x)
      
      Yij(6,1)=1._dp/r(k)*(2._dp*O+1._dp)*Yij(5,1)&
      +r(k)*3._dp*O*gam*h/(2._dp*(2._dp*O+3._dp))*Psi(O,x) 

    
    ENDIF
!+++++++ Normalisation +++++++   

      Yij(:,:)=Yij(:,:)*(r(k)/r(nc))**Ord
                          
END SUBROUTINE Init

!****************************************************          
      Function Phi(m,x)
      use prec
      implicit none
      real(kind=dp), intent(in) :: m
      complex(kind=dp), intent(in) :: x
      complex(kind=dp) :: Phi
      
      Phi=1._dp-x**2/(2._dp*(2._dp*m+3._dp))+x**4/(8._dp*(2._dp*m+3._dp)*(2._dp*m+5._dp))
           
      End Function Phi
!****************************************************     
      Function Psi(m,x)
      use prec
      implicit none
      real(kind=dp), intent(in) :: m
      complex(kind=dp), intent(in) :: x
      complex(kind=dp) :: Psi
      
      Psi=1._dp-x**2/(4._dp*(2._dp*m+5._dp))+x**4/(12._dp*(2._dp*m+5._dp)*(2._dp*m+7._dp))
           
      End Function Psi    
!****************************************************
      SUBROUTINE UpY
!     INTEGRATION DES FONCTIONS YIJ DU CENTRE VERS LA SURFACE 
      use param_Y  
      use prec 
      use interf
      implicit none       
      integer :: i,j,k
      complex(kind=dp) :: As(3), Q21s, Q31s
               
    DO k=ni+1,nc
            If(vs(k)==0._dp .and. vs(k-1) /= 0._dp) then !Transition Sol/Liq
              IF (acc==1) THEN 
             !!!! Following Takeuchi and Saito (1972)
             Yij(1,1:2)=Yij(1,1:2)-Yij(4,1:2)/Yij(4,3)*Yij(1,3) 
             Yij(2,1:2)=Yij(2,1:2)-Yij(4,1:2)/Yij(4,3)*Yij(2,3) 
             Yij(5,1:2)=Yij(5,1:2)-Yij(4,1:2)/Yij(4,3)*Yij(5,3) 
             Yij(6,1:2)=Yij(6,1:2)-Yij(4,1:2)/Yij(4,3)*Yij(6,3) 
             Yij(:,3)=0._dp            
             Yij(4,1:2)=0._dp
             Yij(3,1:2)=-1._dp/(w**2*rho(k)*r(k))*(Yij(2,1:2)-rho(k)*(g(k)*Yij(1,1:2)-Yij(5,1:2)))
             ELSE
             !!!!! Following Saito (1974), Eq. (21)
             As=Yij(2,1:3)-rho(k)*(g(k)*Yij(1,1:3)-Yij(5,1:3))
             Q21s=(Yij(4,3)*As(1)-Yij(4,1)*As(3))/(Yij(4,2)*As(3)-Yij(4,3)*As(2))
             Q31s=(Yij(4,2)*As(1)-Yij(4,1)*As(2))/(Yij(4,3)*As(2)-Yij(4,2)*As(3))  
             Yij(5,1)=Yij(5,1)+Q21s*Yij(5,2)+Q31s*Yij(5,3)
             Yij(6,1)=Yij(6,1)+Q21s*Yij(6,2)+Q31s*Yij(6,3)+4._dp*GG*pi/g(k)*(Yij(2,1)+Q21s*Yij(2,2)+Q31s*Yij(2,3))
             Yij(:,2)=0._dp
             Yij(:,3)=0._dp
             Yij(:,2)=0._dp
             Yij(1:4,1)=0._dp
              ENDIF

            goto 333
            Endif
            
            IF (vs(k-1)==0._dp .and. vs(k)/=0._dp) then !Transition Liq/Sol
              IF (acc==1) THEN 
                Yij(3,1)=0._dp
                Yij(4,1)=0._dp
                Yij(3,2)=0._dp
                Yij(4,2)=0._dp
                Yij(:,3)=0._dp
                Yij(3,3)=1._dp
	      ELSE 
                Yij(1,1)=0._dp
                Yij(2,1)=-rho(k)*Yij(5,1)
                Yij(6,1)=Yij(6,1)+4._dp*pi*GG*rho(k)/g(k)*Yij(5,1)
                Yij(1,2)=1._dp
                Yij(2,2)=rho(k)*g(k)*Yij(1,2)
                Yij(6,2)=-4._dp*pi*GG*rho(k)*Yij(1,2)
                Yij(3,3)=1._dp
              ENDIF
            goto 333
            Endif

            IF (r(k-1)==r(k)) goto 333     
   
            IF ((vs(k)==0._dp).and.(acc==1)) THEN
		CALL RK4(k,Yij)
	    ELSE
		CALL Integ(k,Yij)
            ENDIF

333   continue          


            Yijr(:,:,k)=Yij(:,:)  
            

    ENDDO
   
    END SUBROUTINE UpY
    
!****************************************************
    SUBROUTINE RK4(k,yy)
    use param_Y
    use prec
!    use interf
	IMPLICIT NONE
	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: yy	
    real(kind=dp) :: r1,r2
	real(kind=dp) :: h,h6,hhh,xh,x
	complex(kind=dp), DIMENSION(6,3) :: dym,dyt,yt,dydx
	integer :: nstep,l,kk
	
	r1=r(k-1)
     r2=r(k)
	
	nstep=nrk
	h=r2-r1 	
	h=h/real(nstep)
	
	x=r1
	Do l=1,nstep
  
    call Dy(k,x,yy,dydx)	
	hhh=h*0.5_dp
	h6=h/6.0_dp
	xh=x+hhh
	yt=yy+hhh*dydx
	call Dy(k,xh,yt,dyt)
	yt=yy+hhh*dyt
    call Dy(k,xh,yt,dym)
	yt=yy+h*dym
	dym=dyt+dym
    call Dy(k,x+h,yt,dyt)
	yy=yy+h6*(dydx+dyt+2.0_dp*dym)
	
	x=x+h
	Enddo
	
	END SUBROUTINE rk4
!************************************************************      
   SUBROUTINE Integ(k,ystart)
    use param_Y
    use prec
!    use interf
	IMPLICIT NONE
	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: ystart
	real(kind=dp) :: x1,x2
	real(kind=dp), PARAMETER :: TINY=1.0e-30_dp
	INTEGER :: MAXSTP
	INTEGER :: nstp
	real(kind=dp) :: h,hdid,hnext,x,xsav
	complex(kind=dp), DIMENSION(6,3) :: dydx,yy,yscal
	INTEGER :: nok,nbad,kount,i,j,kk
	real(kind=dp) :: eps,h1,hmin
 
    x1=r(k-1)
    x2=r(k)
		
	MAXSTP=maxRK
	eps=precRK
	 
	hmin=0._dp
	
	!premier pas choisi, sera reduit si l'erreur est trop grande
	h1=(x2-x1)/1000._dp
    
    x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	yy(:,:)=ystart(:,:)
	
	do nstp=1,MAXSTP
		call Dy(k,x,yy,dydx)		
		!terme a partir duquel l'ereur est evaluée : err propto Delta(y)/yscal
		yscal(:,:)=abs(yy(:,:))+abs(h*dydx(:,:))+TINY

		
		!Si on est dans le liquide, certains y restent nuls, donc yscal=0 et err=inf
		!yscal est pose artificiellement grand pour ces y 
		!pour rendre l'erreur nulle 
		Do i=1,6
		Do j=1,3
		If (yy(i,j) == 0._dp .and. dydx(i,j)==0._dp) yscal(i,j)=1.e100_dp
		Enddo
		Enddo

		
		if ((x+h-x2)*(x+h-x1) > 0.0_dp) h=x2-x !si on depasse l'intervalle, h est reduit
		call rkqs(k,yy,dydx,x,h,eps,yscal,hdid,hnext) !retourne hnext<h si erreur trop grande, >h si erreur trop faible
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		if ((x-x2)*(x2-x1) >= 0.0_dp) then   !si on depasse l'intervalle, solution retenue
			ystart(:,:)=yy(:,:)
			RETURN
		end if
		if (abs(hnext) < hmin)&
			write(*,*) ('stepsize smaller than minimum in odeint')
		h=hnext
	end do
!write(*,*) nok
	write(*,*) 'too many steps in runge-kutta for such a precision'
	stop
	 
	END SUBROUTINE Integ
!****************************************************	
    SUBROUTINE rkqs(k,yy,dydx,x,htry,eps,yscal,hdid,hnext)
    use param_Y
    use prec
!    use interf
	IMPLICIT NONE
	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(INOUT) :: yy
	complex(kind=dp), DIMENSION(6,3), INTENT(IN) :: dydx,yscal
	real(kind=dp), INTENT(INOUT) :: x
	real(kind=dp), INTENT(IN) :: htry,eps
	real(kind=dp), INTENT(OUT) :: hdid,hnext
	INTEGER :: ndum,i,j
	real(kind=dp) :: errmax,h,htemp,xnew
	complex(kind=dp), DIMENSION(6,3) :: yerr,ytemp
	REAL(kind=dp), PARAMETER :: SAFETY=0.9_dp,PGROW=-0.2_dp,PSHRNK=-0.25_dp,ERRCON=1.89e-4_dp
	

	h=htry
	do
		call rkck(k,yy,dydx,x,h,ytemp,yerr)
		
		errmax=maxval(abs(yerr(:,:)/yscal(:,:)))/eps

!                write(*,*) 'errmax',errmax, 'radius =',x*rn,'h',h*rn!!!!!!!!!!
		if (errmax <= 1.0_dp) exit !si err trop grande, on diminue le pas
		htemp=SAFETY*h*(errmax**PSHRNK)
		h=sign(max(abs(htemp),0.1_dp*abs(h)),h)
		xnew=x+h 
		if (xnew == x) then
		write(*,*) 'stepsize underflow in rkqs, radius =',x*rn!!!!!!!!!
		stop
		endif
	end do
	if (errmax > ERRCON) then !si err trop petite, on agrandit le pas
		hnext=SAFETY*h*(errmax**PGROW)


!    if (errmax==maxval(abs(yerr(3,1:2)/yscal(3,1:2)))/eps) &
!                write(*,*) 'errmax',errmax, k
!                write(*,*) 'hnext=', hnext
	else
		hnext=5.0_dp*h
	end if
	hdid=h
	x=x+h
	yy(:,:)=ytemp(:,:)
!        if (vs(k).eq.0._dp)  write(*,*) 'TEST - rkqs', k, hdid, yscal(1,2), ytemp(1,2),yerr(1,2)

	END SUBROUTINE rkqs
!****************************************************
   	SUBROUTINE rkck(k,yy,dydx,x,h,yout,yerr)
    use prec
!    use interf
	IMPLICIT NONE
	integer, intent(in) :: k
	complex(kind=dp), DIMENSION(6,3), INTENT(IN) :: yy,dydx
	real(kind=dp), INTENT(IN) :: x,h
	complex(kind=dp), DIMENSION(6,3), INTENT(OUT) :: yout,yerr
	INTEGER :: ndum
	complex(kind=dp), DIMENSION(6,3) :: ak2,ak3,ak4,ak5,ak6,ytemp
	real(kind=dp), PARAMETER :: A2=0.2_dp ,A3=0.3_dp ,A4=0.6_dp ,A5=1.0_dp ,&
		A6=0.875_dp ,B21=0.2_dp ,B31=3.0_dp /40.0_dp ,B32=9.0_dp /40.0_dp ,&
		B41=0.3_dp ,B42=-0.9_dp ,B43=1.2_dp ,B51=-11.0_dp /54.0_dp ,&
		B52=2.5_dp ,B53=-70.0_dp /27.0_dp ,B54=35.0_dp /27.0_dp ,&
		B61=1631.0_dp /55296.0_dp ,B62=175.0_dp /512.0_dp ,&
		B63=575.0_dp /13824.0_dp,B64=44275.0_dp /110592.0_dp ,&
		B65=253.0_dp /4096.0_dp ,C1=37.0_dp /378.0_dp ,&
		C3=250.0_dp /621.0_dp ,C4=125.0_dp /594.0_dp ,&
		C6=512.0_dp /1771.0_dp ,DC1=C1-2825.0_dp /27648.0_dp ,&
		DC3=C3-18575.0_dp /48384.0_dp ,DC4=C4-13525.0_dp /55296.0_dp ,&
		DC5=-277.0_dp /14336.0_dp ,DC6=C6-0.25_dp 

	ytemp=yy+B21*h*dydx
	call Dy(k,x+A2*h,ytemp,ak2)
	ytemp=yy+h*(B31*dydx+B32*ak2)
	call Dy(k,x+A3*h,ytemp,ak3)
	ytemp=yy+h*(B41*dydx+B42*ak2+B43*ak3)
	call Dy(k,x+A4*h,ytemp,ak4)
	ytemp=yy+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	!	write(*,*)'A5',(x+A5*h)*1740000._dp!!!!!!!!!!
	call Dy(k,x+A5*h,ytemp,ak5)
	ytemp=yy+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call Dy(k,x+A6*h,ytemp,ak6)
	yout=yy+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
 END SUBROUTINE rkck
!****************************************************
 SUBROUTINE Dy(k,rr,yy,dYdr)
    use param_Y
    use prec
    use interf
    implicit none
    integer, intent(in) :: k
    integer :: i
    real(kind=dp), intent(in) :: rr
    complex(kind=dp), dimension(6,3), intent(inout) :: yy
    complex(kind=dp), dimension(6,3), intent(out) :: dYdr
    complex(kind=dp), dimension(6,3) :: soliq, yypred
    complex(kind=dp), dimension(3) :: Y3
    real(kind=dp) :: dr,rhor,grr
    complex(kind=dp) :: lamr,mur
   
    
 dr=(rr-r(k-1))/(r(k)-r(k-1))
 rhor=rho(k-1)+dr*(rho(k)-rho(k-1))
 mur=muc(k-1)+dr*(muc(k)-muc(k-1))
 lamr=lamc(k-1)+dr*(lamc(k)-lamc(k-1))
 grr=gr(rr,k)!g(k-1)+dr*(g(k)-g(k-1))!
         

   IF ((vs(k)==0._dp).and.(acc==0)) THEN
    soliq(:,:)=0._dp
    soliq(5,1)=1._dp
    soliq(6,1)=1._dp
    dYdr(:,:)=0.d0
    dYdr(5,1)=(4._dp*pi*GG*rhor/grr-(O+1._dp)/rr)*yy(5,1)+yy(6,1)
    dYdr(6,1)=2._dp*(O-1._dp)/rr*4._dp*pi*GG*rhor/grr*yy(5,1)+((O-1._dp)/rr-4._dp*pi*GG*rhor/grr)*yy(6,1)
!!!!!!!! dYdr(6,1) et yy(6,1) correspondent a la variable y7 de Saito (1974) et y6 dans Takeushi et Saito (1972)
!   write(*,*) 'TEST', rr, real(dYdr(5,1)), real(yy(5,1)), real(dYdr(6,1)), real(yy(6,1))
   ELSE

   soliq(:,:)=1._dp
   Y3(:)=yy(3,:)



   IF (vs(k)==0._dp) then
      soliq(:,3)=0._dp
      soliq(3,:)=0._dp
      soliq(4,:)=0._dp
      Y3(1:2)=-1._dp/(w**2*rhor*rr)*(yy(2,1:2)-rhor*(grr*yy(1,1:2)-yy(5,1:2)))
      Y3(3)=0._dp
      yy(3,:)=Y3(:)
   ENDIF
            
   dYdr(1,:)=1._dp/(lamr+2._dp*mur)*(yy(2,:)-lamr/rr*(2._dp*yy(1,:)-O*(O+1._dp)*Y3(:)))

   dYdr(2,:)=+2._dp/rr*(lamr*dYdr(1,:)-yy(2,:))&
             +1._dp/rr*(2._dp/rr*(lamr+mur)-rhor*grr)*(2._dp*yy(1,:)-O*(O+1.)*Y3(:))&
             +O*(O+1)/rr*yy(4,:)-rhor*(yy(6,:)-(O+1.)/rr*yy(5,:)+2._dp*grr/rr*yy(1,:))
   dYdr(2,:)=dYdr(2,:)-w**2*rhor*yy(1,:)


   IF (vs(k)==0._dp) THEN
       
        dYdr(3,:)=0._dp
       
        dYdr(4,:)=0._dp

    else

      dYdr(3,:)=1._dp/mur*yy(4,:)+1._dp/rr*(Y3(:)-yy(1,:))

      dYdr(4,:)=-lamr/rr*dYdr(1,:)&
       -(lamr+2._dp*mur)/rr**2*(2._dp*yy(1,:)-O*(O+1._dp)*Y3(:))&
       +2._dp*mur/rr**2*(yy(1,:)-Y3(:))-3._dp/rr*yy(4,:)&
       -rhor/rr*(yy(5,:)-grr*yy(1,:))

      dYdr(4,:)=dYdr(4,:)-w**2*rhor*Y3(:)

     endif
  
     dYdr(5,:)=yy(6,:)+4._dp*pi*GG*rhor*yy(1,:)-(O+1._dp)/rr*yy(5,:)
    
     dYdr(6,:)=(O-1._dp)/rr*(yy(6,:)+4._dp*pi*GG*rhor*yy(1,:))&
     +4._dp*pi*GG*rhor/rr*(2._dp*yy(1,:)-O*(O+1._dp)*Y3(:)) 
   
     ENDIF


     dYdR(:,:)=soliq(:,:)*dYdR(:,:)

  
 END SUBROUTINE Dy
!****************************************************   
 FUNCTION gR(rr,k)
   use param_Y
   use prec
   implicit none
   real(kind=dp), intent(in) :: rr
   integer, intent(in) :: k
   real(kind=dp) :: gR
   real(kind=dp) :: maR
   
    maR=mass(k-1)+4._dp*pi/3._dp*(rr**3-r(k-1)**3)*rho(k) !maK(k-1)=masse au sommet de la couche inf
    gR=GG*MaR/(rr**2)
 
 END FUNCTION gR
!****************************************************   
 SUBROUTINE Eigen
 use param_Y
 use prec
 use interf
 complex(kind=dp) :: Mat(3,3),Minv(3,3),deem,c1,c2,c3,c3new,c2new
 complex(kind=dp) :: As(3), Q21s, Q31s
 real(kind=dp) :: F1,F2,F3
 



   Mat(1,:)=Yijr(2,:,nc)
   Mat(2,:)=Yijr(4,:,nc)
   
    F1=0._dp
    F2=0._dp

!    If (forc==2) F1=-(2._dp*O+1._dp)/R(nc)*g(nc)/(4._dp*pi*GG)
    If (forc==2) F1=-3058.0419077634106_dp

    If(forc==0) then
    Mat(3,:)=Yijr(1,:,nc)
    F3=1._dp ! amplitude indeterminee sans source
    else
    Mat(3,:)=Yijr(6,:,nc)
!    F3=(1._dp+2._dp*O)*(r(nc)**(Ord-1)) 
    F3=(1._dp+2._dp*O)*(r(nc)**(-1)) ! Forcage
    Endif

  if (vs(nc)/=0._dp) then 


   Minv(:,:)=invM(Mat)
          
   ! Y2(R) =  c1*ys(2,1) + c2*ys(2,2) + c3*ys(2,3) = F1
   ! Y4(R) =  c1*ys(4,1) + c2*ys(4,2) + c3*ys(4,3) = F2
   ! Y6(R) =  c1*ys(6,1) + c2*ys(6,2) + c3*ys(6,3) = F3
        
    c1=F1*Minv(1,1)+F2*Minv(1,2)+F3*Minv(1,3)
    c2=F1*Minv(2,1)+F2*Minv(2,2)+F3*Minv(2,3)
    c3=F1*Minv(3,1)+F2*Minv(3,2)+F3*Minv(3,3)
    
   endif 

  if (vs(nc)==0._dp) then 
    if (acc==1) then
    c1=(F1*Yijr(6,2,nc)-F3*Yijr(2,2,nc))/(Yijr(2,1,nc)*Yijr(6,2,nc)-Yijr(2,2,nc)*Yijr(6,1,nc))
    c2=(F1*Yijr(6,1,nc)-F3*Yijr(2,1,nc))/(Yijr(2,2,nc)*Yijr(6,1,nc)-Yijr(2,1,nc)*Yijr(6,2,nc))
    c3=0.d0
    else
    c1=F3/Yijr(6,1,nc)
    c2=0._dp
    c3=0._dp
    endif
 endif

    Yi(:,nc)=c1*Yijr(:,1,nc)+c2*Yijr(:,2,nc)+c3*Yijr(:,3,nc) !Fonctions propres en surface

    
    c3new=c3
    c2new=c2

    DO i=nc-1,1,-1 ! Redefinition de c3 aux interfaces (de la surface vers le centre)    
     If(vs(i) /= 0._dp .and. vs(i+1) == 0._dp)  then 
      if (acc==1) c3new=-(c1*Yijr(4,1,i)+c2*Yijr(4,2,i))/Yijr(4,3,i)  
      if (acc==0) then
         Yij(:,:)=Yijr(:,:,i)
         As=Yij(2,1:3)-rho(i+1)*(g(i+1)*Yij(1,1:3)-Yij(5,1:3))
         Q21s=(Yij(4,3)*As(1)-Yij(4,1)*As(3))/(Yij(4,2)*As(3)-Yij(4,3)*As(2))
         Q31s=(Yij(4,2)*As(1)-Yij(4,1)*As(2))/(Yij(4,3)*As(2)-Yij(4,2)*As(3))        
        c3new=Q31s*c1
        c2new=Q21s*c1 
      endif
     endif

     Yi(:,i)=c1*Yijr(:,1,i)+c2new*Yijr(:,2,i)+c3new*Yijr(:,3,i)   


    Enddo
    
!    Yi(:,:)=Yi(:,:)/r(nc)**Ord
              
 END SUBROUTINE Eigen
 
 	FUNCTION invM(X1)
 	use prec
	      implicit none
          complex(kind=dp), intent(in), dimension(:,:) :: X1
          complex(kind=dp), dimension(size(X1,1),size(X1,1)) :: invM !matrice carree
          !--------------------
          complex(kind=dp), dimension(size(X1,1),size(X1,1)) :: W
          complex(kind=dp), dimension(size(X1,1),size(X1,1)) :: Y
          integer i, j, k, n
        
          W(:,:) = X1(:,:)
          invM(:,:)=0._dp
          Do i=1,size(X1,1)
          invM(i,i)=1._dp !identite
          Enddo
          

        ! perform Gauss elimination on W*I=1/X1
          n = size(X1,2)
          do k = 1,n-1
            do i=k+1,n
              W(i,k) = W(i,k)/W(k,k)
              invM(i,:) = invM(i,:) - W(i,k) * invM(k,:)
            end do
            do j=k+1,n
              do i=k+1,n
                W(i,j) = W(i,j) - W(i,k) * W(k,j)
              end do
            end do
          end do

        ! perform back substitution on invM
          do k = n,1,-1
            invM(k,:) = invM(k,:) / W(k,k)
            do i=1,k-1
              invM(i,:) = invM(i,:) - W(i,k) * invM(k,:)
            end do
          end do
     END FUNCTION invM

SUBROUTINE I123
     use param_Y
     use prec 
     implicit none 
     integer :: j
                
     dI1(:)=rho(:)*(conjg(yi(1,:))*yi(1,:)+O*(O+1._dp)*conjg(yi(3,:))*yi(3,:))*r(:)**2
 
  !   dy1dr(1)=0.
  !   dy5dr(1)=0.
  !   Do j=2,nc
  !   if(r(j)==r(j-1)) then
  !   dy1dr(j)=(yi(1,j)-yi(1,j-2))/(r(j)-r(j-2))
  !   dy5dr(j)=(yi(5,j)-yi(5,j-2))/(r(j)-r(j-2))
  !   else
  !   dy1dr(j)=(yi(1,j)-yi(1,j-1))/(r(j)-r(j-1))
  !   dy5dr(j)=(yi(5,j)-yi(5,j-1))/(r(j)-r(j-1))
  !   endif
  !  Enddo
     
      dy1dr(:)=(yi(2,:)-lamc(:)/r(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))/(lamc(:)+2._dp*muc(:))     
      dy5dr(:)=yi(6,:)+4._dp*pi*GG*rho(:)*yi(1,:)-(O+1._dp)/r(:)*yi(5,:)
 
     HH(1,:)=r(:)**2/(conjg(Ka(:)+4._dp/3._dp*muc(:))*(Ka(:)+4._dp/3._dp*muc(:)))&
     *conjg(yi(2,:)-(Ka(:)-2._dp/3._dp*muc(:))/r(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))&
     *(yi(2,:)-(Ka(:)-2._dp/3._dp*muc(:))/r(:)*(2._dp*yi(1,:)-O*(O+1.)*yi(3,:)))&
     +2._dp*r(:)*real(dy1dr(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))&
     +conjg(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:))*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:))
 
     HH(2,:)=4._dp/3._dp*r(:)**2/(conjg(Ka(:)+4._dp/3._dp*muc(:))*(Ka(:)+4._dp/3._dp*muc(:)))&
     *conjg(yi(2,:)-(Ka(:)-2._dp/3._dp*muc(:))/r(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))&
     *(yi(2,:)-(Ka(:)-2._dp/3._dp*muc(:))/r(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))&
     -4._dp/3._dp*r(:)*real(dy1dr(:)*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))&
     +1._dp/3._dp*conjg(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:))*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:))&
     +O*(O+1._dp)*r(:)**2*conjg(yi(4,:))*yi(4,:)/(conjg(muc(:))*muc(:))&
     +O*(O**2-1._dp)*(O+2._dp)*conjg(yi(3,:))*yi(3,:)

     


     HH(3,:)=2._dp*(O+1._dp)*r(:)*real(conjg(yi(1,:))*yi(5,:))&
     -2._dp*O*(O+1._dp)*r(:)*real(conjg(yi(3,:))*yi(5,:))&
     -2._dp*G(:)*r(:)*real(conjg(yi(1,:))*(2._dp*yi(1,:)-O*(O+1._dp)*yi(3,:)))

     HH(4,:)=r(:)*conjg(yi(6,:))*yi(6,:)/(4._dp*pi*GG)
 
     dI2(:)=Ka(:)*HH(1,:)+muc(:)*HH(2,:)+rho(:)*HH(3,:)+HH(4,:)
 
     I1=dI1(2)*r(2)
     I2=dI2(2)*r(2)
     Y1S=yi(1,2) !dy1dr(2)*r(2)
 
     Do j=3,nc
     I1=I1+dI1(j)*(r(j)-r(j-1))
     I2=I2+dI2(j)*(r(j)-r(j-1))

     Y1S=Y1S+dy1dr(j)*(r(j)-r(j-1))
     if(r(j)==r(j-1)) Y1S=Y1S+(yi(1,j)-yi(1,j-1))
     Enddo
     
     I3=r(nc)**2*(conjg(yi(1,nc))*yi(2,nc)+O*(O+1._dp)*conjg(yi(3,nc))*yi(4,nc)+conjg(yi(5,nc))*yi(6,nc)/(4._dp*pi*GG))
   
      HH(2,1)=0.D0

!***Write Hmu*mu in file*********
open(102,file='Outputs/Hmu_mu.dat')
do j=2,nc
	if (muc(j)==0.d0) HH(2,j)=0.D0
	write(102,*) r(j), real(HH(2,j))*imag(muc(j)), period/(24._dp*3600._dp), nc
enddo


END SUBROUTINE I123

 !*********************    
     FUNCTION ran(idum)
     IMPLICIT NONE
     INTEGER, INTENT(INOUT) :: idum
     REAL :: ran
     INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
     REAL, SAVE :: am
     INTEGER, SAVE :: ix=-1,iy=-1,k
     if (idum <= 0 .or. iy < 0) then
          am=nearest(1.0D0,-1.0D0)/IM
          iy=ior(ieor(888889999,abs(idum)),1)
          ix=ieor(777755555,abs(idum))
          idum=abs(idum)+1
     end if
     ix=ieor(ix,ishft(ix,13))
     ix=ieor(ix,ishft(ix,-17))
     ix=ieor(ix,ishft(ix,5))
     k=iy/IQ
     iy=IA*(iy-k*IQ)-IR*k
     if (iy < 0) iy=iy+IM
     ran=am*ior(iand(IM,ieor(ix,iy)),1)
     END FUNCTION ran
!*********************
   !*********************

     FUNCTION gammln(xx) 
     
     
     use prec 
    IMPLICIT NONE
    INTEGER :: j
    REAL(kind=dp), intent(in) :: xx
    REAL(kind=dp) :: gammln
    REAL (kind=dp) :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp 

    DATA cof,stp/76.18009172947146_dp,-86.50532032941677_dp, &
     24.01409824083091_dp,-1.231739572450155_dp,.1208650973866179e-2_dp,& 
    -.5395239384953e-5_dp,2.5066282746310005_dp/ 

    x=xx
    y=x 
    tmp=x+5.5_dp
    tmp=(x+0.5_dp)*log(tmp)-tmp 
    ser=1.000000000190015_dp 
    do j=1,6
     y=y+1._dp
    ser=ser+cof(j)/y 
    enddo

    gammln=tmp+log(stp*ser/x) 

    END FUNCTION gammln   


!*******************
     FUNCTION J_Maxwell(mut,vis,w_forc) 

     use prec 
     use param_rheo     

     IMPLICIT NONE
 
         
     REAL(kind=dp), intent(in) :: mut, vis, w_forc
     complex(kind=dp) :: J_Maxwell
     
     J_Maxwell=cmplx(1._dp/mut,-1._dp/(vis*w_forc),kind=dp)

    END FUNCTION J_Maxwell

!*********************

     FUNCTION J_Andrade(mut,vis,w_forc) 

! following Castillo-Rogez et al. (2011)

     use prec 
     use param_rheo
     use variables, only : alpha_gam
  

     IMPLICIT NONE

     REAL(kind=dp), intent(in) :: mut,vis, w_forc
     REAL(kind=dp) :: gammln
     complex(kind=dp) :: J_Andrade

     J_Andrade=cmplx(1._dp/mut,-1._dp/(vis*w_forc), kind=dp)
     J_Andrade=J_Andrade+cmplx(1._dp/mut,kind=dp)&
                        *(cmplx(0._dp,w_forc*vis/mut,kind=dp))**(-alpha_gam)&
                        *cmplx(exp(gammln(1._dp+alpha_gam)),kind=dp)

    END FUNCTION J_Andrade


!************************
