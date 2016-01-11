
MODULE commondata
!AY wrote
!<Energy>
! eki     : kinetic energy
! eppp    : SCF( =0.5*P(H+F) )  + CI energy                          [kJ/mol]
! ebetag  : nuclear core-core interaction = sum_uv[Z_u*Z_v*Gamma_uv] [kJ/mol]
! eqcff   : potential energy of MM classical part
! eground :
! ecore
! 
  !USE class_atom

  IMPLICIT NONE
  SAVE

  ! Constants and Unit Conversions
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323d0
  DOUBLE PRECISION, PARAMETER :: boltz=1.380662d-23   !(J/K)
  DOUBLE PRECISION, PARAMETER :: avogad=6.022045d023     !Avogadro's number
  DOUBLE PRECISION, PARAMETER :: amu=1.6605655d-27      !mass units in kg
  DOUBLE PRECISION, PARAMETER :: a0=.529177209   !bohr radius/angstrom
  DOUBLE PRECISION, PARAMETER :: hbar=0.0635078    !hbar ps kJ/mol
  
  !APW new variables meant to collect and streamline this mess
  !TYPE(ATOM), ALLOCATABLE, DIMENSION(:) :: atoms

  LOGICAL :: log_cicalc,log_exdyn,usenose,log_solv,log_scalevel,&
       log_dumpvel,log_confine,log_rescale_i,log_dderiv,log_delre,&
       log_ppv,log_restart,log_addcos,log_C60shift
  INTEGER :: nnuc,nelec,norb,npi,istep,nsteps,nsolve,iseed,iscale,istat,irestart,navg,nexstate,idsx
    
  DOUBLE PRECISION :: deltat,ttol,Tmb,unitl,unitm,unite,d_screen,eki,eppp,ebetag,eqcff,eground,ecore
  
END MODULE commondata

!----------------------------------------!

MODULE CLASSICAL
!AY wrote
!<id number>
! pi_iden(i) : conversion from a pi atom id number to a atom id number for all
! atype(i)   : atom type of i-th atom
!<Basic Variables>
! amass(i)  :  mass of i-th atom       []
! vel(3,i)  :  velocity of i-th atom   [A/ps]
! pos(3,i)  :  position of i-th atom   [A]
! fra(3,i)  :  force of i-th atom      [kJ/mol/A]
  IMPLICIT NONE
  SAVE
  
      
  CHARACTER*1, ALLOCATABLE, DIMENSION(:) :: pispec,species
  
  INTEGER :: nphis,nbonds,nthetas,nphis_in,npispec,ngamma,nphi_in,nspec
  INTEGER :: ntheta_type,nbond_type,nphi_type,nnonb_type,nbeta_type
  INTEGER :: b_sw,p_sw,t_sw,tp_sw,nb_sw
  !INTEGER, DIMENSION(nspec) :: pispec_code,theta_equiv,icvn
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atype,pi_iden,ncos,phi_in,iatg,pispec_code,&
       theta_equiv,icvn
  !APW delre
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ilink,bond_code,theta_code,phi_code,nonb_code,beta_code
  !INTEGER, DIMENSION(nspec,nspec) :: pair_code
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bonds,thetas,phis,nonbs,links,pair_code
  DOUBLE PRECISION :: dfield,unitm,mass_tot
  !DOUBLE PRECISION, DIMENSION(nspec) :: mass
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: amass,qatom,new_do,mass
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pos,deriv,vel,vel_l,pos_m,pos_l,fra
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta_pms,theta_cub,bond_pms,phi_pms,nonb_pms,alpha,beta,gamma
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: deriv_bond,deriv_theta,deriv_phi
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: dderiv
  
  !APW for uniquely defined torsional potentials
  INTEGER :: nphib
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: phib_pairs
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: phib_pms
  
  !c:::  force switch
  data b_sw /1/
  data t_sw /1/
  data tp_sw/1/
  data p_sw /1/
  data nb_sw/1/
  
CONTAINS
  FUNCTION length(dx)
    !APW needs to be modified for PBC (periodic boundary conditions)
    
    IMPLICIT NONE
    
    DOUBLE PRECISION :: length
    DOUBLE PRECISION, DIMENSION(3) :: dx
    
    length=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
    
    return
    
  END FUNCTION length

  !c     costhe calculates the  cosine of theta given the vector component
  !c     and lengths
  FUNCTION costhe(da,db,ra,rb)
    !APW needs to be modified for PBC (periodic boundary conditions)
    
    IMPLICIT NONE
    
    DOUBLE PRECISION :: costhe,ra,rb
    DOUBLE PRECISION, DIMENSION(3) :: da,db
    
    costhe=(da(1)*db(1)+da(2)*db(2)+da(3)*db(3))/(ra*rb)
    
  END FUNCTION costhe

  FUNCTION crosprod(a,b)

    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: a,b
    DOUBLE PRECISION, DIMENSION(3) :: crosprod
    crosprod(1)=a(2)*b(3)-a(3)*b(2)
    crosprod(2)=-a(1)*b(3)+a(3)*b(1)
    crosprod(3)=a(1)*b(2)-a(2)*b(1)

  END FUNCTION crosprod

  FUNCTION locate_atom(atom)
    
    !USE QCFIO, ONLY:gen_io
    !USE classical, ONLY:species,nspec
    
    IMPLICIT NONE
    
    CHARACTER*1 atom   ! character type
    INTEGER i,locate_atom
    
    do i=1,nspec
       if(atom.eq.species(i)) then
          locate_atom=i
          return
       endif
    enddo
    
    !write(gen_io,*) "*** atom ",atom," is not one of the coded atoms ***"
    !stop 
    locate_atom=0
    return
    
  END FUNCTION locate_atom
    
END MODULE CLASSICAL
!-----------------------------------------------------------------------!
MODULE QUANTUM
!<id number>
! ind_nopi(i) : 
! inv_pi(i)   : 
!<charge>
! z_atom(i)   : (formally facion) The nuclear charge on each atom
!<alpha,beta,gamma>
! ibeta
! junc_flg(k)     :  if k is atom forming junction part, this variable is .true.
! inv_dih(i1,i2)  : for a couple of atoms give index of
!                   junction bond if it exist, or zero

!
! st_exp2mu(i,n)  : exp(-2*myu*(dr-dr0))   (n=1)
!                   cos(phi) cosine of torsional angle (n=1)
!                   (i = junction index (pair number of pi atoms, e.x. (i1,i2)) )   
! st_betalin(i,n) : betalin or beta at planar (n=1)
!                   cos(phi) cosine of torsional angle (n=1)
!                   (i = junction index (pair number of pi atoms, e.x. (i1,i2)) )   
! st_fbetalin(i,j) : -d(betalin)/d(r_i1)_j*cos(phi)  (j=1,2,3(=x,y,z))
!                    (betaline is for i1,i2 pair)

  IMPLICIT NONE
  SAVE
  
  !APW units, Hartree for ppp subroutines
  !DOUBLE PRECISION, PARAMETER :: alpha_mu=-0.410113 
  !DOUBLE PRECISION, PARAMETER :: umubcc=1.076 !old
  DOUBLE PRECISION, PARAMETER :: umubcc=0.73868 !new
  DOUBLE PRECISION, PARAMETER :: umubcs=0.81325
  DOUBLE PRECISION, PARAMETER :: umubcn=0.80659
  DOUBLE PRECISION, PARAMETER :: umubco=0.89987
  DOUBLE PRECISION, PARAMETER :: linb2cc=-0.019198 !old
  !DOUBLE PRECISION, PARAMETER :: linb2cc=0.027022 !new
  DOUBLE PRECISION, PARAMETER :: linb2cs=0.025506 
  DOUBLE PRECISION, PARAMETER :: linb2cn=0.029008
  DOUBLE PRECISION, PARAMETER :: linb2co=0.035262
  DOUBLE PRECISION, PARAMETER :: a0cc=-0.089586 !old
  !DOUBLE PRECISION, PARAMETER :: a0cc=-0.085765 !new
  DOUBLE PRECISION, PARAMETER :: a0cs=-0.078013
  DOUBLE PRECISION, PARAMETER :: a0cn=-0.074106
  DOUBLE PRECISION, PARAMETER :: a0co=-0.095542
  DOUBLE PRECISION, PARAMETER :: rcc = 2.639948d0 !units of bohr radii
  DOUBLE PRECISION, PARAMETER :: rcs = 2.66941d0  !units of bohr radii
  DOUBLE PRECISION, PARAMETER :: rcn = 2.639948d0 !units of bohr radii
  DOUBLE PRECISION, PARAMETER :: rco = 2.32436d0  !units of bohr radii
  !APW old values
  !paramPPP.h:      parameter (A0CC=-0.089586,A1CC=0.0)
  !paramPPP.h:c$$$      parameter (A0CC=-0.0,A1CC=0.0)
  !paramPPP.h:      parameter (A0CN=-0.074106,A1CN=0.0) 
  !paramPPP.h:      parameter (A0CO=-0.095542,A1CO=0.0)
  !AY  values are from the paper of original PI/QCFF by Warshel & Karplus (1972)
  !    and from the paper by Warshel & Lappicirella (1981)
  !    also see "Semiempirical methods of electronic structure calculation.Part A (1977)"
  DOUBLE PRECISION, PARAMETER :: Hartree2eV = 27.2116d0              ! [Hartree]->[eV]
  ! Beta'
  DOUBLE PRECISION, PARAMETER :: beta_dash_c  = 0.235d0 / Hartree2eV ! [Hartree]
!  DOUBLE PRECISION, PARAMETER :: beta_dash_c  = 0.20d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: beta_dash_s  = 0.10d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: beta_dash_n  = 0.10d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: beta_dash_np = 0.20d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: beta_dash_o  = 0.20d0  / Hartree2eV ! [Hartree]
  ! I-A   (Now, not use)
  DOUBLE PRECISION, PARAMETER :: i_a_cc  =   9.81d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_cs  =  10.84d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_cn  =  12.70d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_cnp =  12.70d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_co  =  13.50d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_nn  =  16.00d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: i_a_oo  =  17.64d0  / Hartree2eV ! [Hartree]
  ! G0   (Now, not use)
  DOUBLE PRECISION, PARAMETER :: g0_cc   =   5.14d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_cs   =   2.30d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_cn   =  11.00d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_cnp  =  11.68d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_co   =   8.00d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_nn   =  11.00d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: g0_oo   =  10.00d0  / Hartree2eV ! [Hartree]
  ! Gs
  DOUBLE PRECISION, PARAMETER :: gs_cc   =  0.69d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_cs   =  0.30d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_cn   =  0.30d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_cnp  =  0.60d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_co   =  0.60d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_nn   =  0.60d0  / Hartree2eV ! [Hartree]
  DOUBLE PRECISION, PARAMETER :: gs_oo   =  0.00d0  / Hartree2eV ! [Hartree]
  ! mu_gamma  (Now, not use)
  DOUBLE PRECISION, PARAMETER :: mu_g_cc =  0.232d0  ![1/A]
  DOUBLE PRECISION, PARAMETER :: mu_g_cs =  0.450d0  ![1/A]
  DOUBLE PRECISION, PARAMETER :: mu_g_cn =  0.450d0
  DOUBLE PRECISION, PARAMETER :: mu_g_cnp=  0.070d0
  DOUBLE PRECISION, PARAMETER :: mu_g_co =  0.430d0
  DOUBLE PRECISION, PARAMETER :: mu_g_nn =  0.240d0
  DOUBLE PRECISION, PARAMETER :: mu_g_oo =  0.350d0


  LOGICAL :: ifinish
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: junc_flg
  INTEGER :: ibeta_pairs,nci,neig,ndih_nopi,nnopi,iexstate,iexstate_l,ntully
  !parameter (ntully=7)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inv_pi,inv_dih_nopi,ind_dih_nopi,ind_nopi,&
       idxzvecq,idxzvecp
  !APW changed
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibeta,inv_dih,inv_ibeta
  !APW bond order matrix
  
  DOUBLE PRECISION :: deltae
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: z_atom,sigma,ciener,ciener_l,osc,gamma_screen,alpha_mu
  !APW gamma
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: gamma_screen_i, alpha_mu_i
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: st_betalin,st_fbetalin !APW array to store betalin each step
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bomat,mocoef,mocoef_old,osc_vec,deltaab,fdv,fdv_l,dar,aeainv
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: fdv_m
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cieig_old,cieig
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ciind,ciind_old
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: fexa,fexa_l,qr,utmp
  !AY
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: st_exp2mu
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dcosdr
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: qr_nopi
  
  COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: ctully
  !ctully(0:ntully)
  
CONTAINS

  !AY added
  FUNCTION exp2mu(dr,spec1,spec2)
    !*  +---------------------------------------------------------------+
    !*  |  Term of exp(-2u_CX (R_uv-R_CX^eq))                           |
    !*  |                                                               |  
    !*  |       dr    = R_uv                                            |
    !*  |       spec1 = atom type of u                                  |
    !*  |       spec2 = atom type of v                                  |
    !*  +---------------------------------------------------------------+
    IMPLICIT NONE
    DOUBLE PRECISION :: exp2mu
    !input variables
    CHARACTER*1, INTENT(IN) :: spec1,spec2
    DOUBLE PRECISION, INTENT(IN) :: dr
    !local variables
    DOUBLE PRECISION :: mu,r0

    call  get_mu_b_r0(spec1,spec2,mu,r0)
    exp2mu = dexp( -2.0d0 * mu * (dr-r0) )
    
  END FUNCTION exp2mu

  
  !FUNCTION betalin(dr,umub,linb2,a0c,r)
  FUNCTION betalin(dr,spec1,spec2)
    !*  +---------------------------------------------------------------+
    !*  |  Beta_uv = hbar^2/m_e * (1/R_uv) dS_uv/dR_uv                  |  
    !*  |                                                               |  
    !*  |  numerical fitting expression:                                |  
    !*  |                                                               |  
    !*  |  Beta_uv = exp [ -u_CX ( R_uv - R_CX^eq ) ]                   |  
    !*  |                [Beta_1^CX + Beta_2^CX ( R_uv - R_CX^eq ) ]    |  
    !*  +---------------------------------------------------------------+
    IMPLICIT NONE
    
    DOUBLE PRECISION :: betalin
    !input variables
    CHARACTER*1, INTENT(IN) :: spec1,spec2
    DOUBLE PRECISION, INTENT(IN) :: dr
    !local variables
    DOUBLE PRECISION :: mu,a,b,r
    LOGICAL :: hit

    hit=.true.
    
    if(spec1.eq.'A' .and. spec2.eq.'A') then
       mu=umubcc
       a=a0cc
       b=linb2cc
       r=rcc
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'S') .or. (spec1.eq.'S'.and.spec2.eq.'A')) then
       mu=umubcs
       a=a0cs
       b=linb2cs
       r=rcs
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'N') .or. (spec1.eq.'N'.and.spec2.eq.'A')) then
       mu=umubcn
       a=a0cn
       b=linb2cn
       r=rcn
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'M') .or. (spec1.eq.'M'.and.spec2.eq.'A')) then
       mu=umubcn
       a=a0cn
       b=linb2cn
       r=rcn
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'O') .or. (spec1.eq.'O'.and.spec2.eq.'A')) then
       mu=umubco
       a=a0co
       b=linb2co
       r=rco
       hit=.false.
    end if
    
    if(hit) then
       print *,'***WARNING: combo not defined in linderberg approximation ',spec1,spec2,' ***'
       stop
    end if
    
    betalin=dexp(-mu*(dr-r))*(a+b*(dr-r))
       !betalin=dexp(-umub*(dr-r))*(a0c+linb2*(dr-r))

  END FUNCTION betalin
  
  !FUNCTION fbetalin(vr,umub,linb2,r,betalin,angle)
  FUNCTION fbetalin(vr,spec1,spec2,betalin,angle)
    !* =====================================================================
    !* []  DECK BETALINF --  RETURN FORCE CONTRIBUTIONS FROM THE 
    !*                       ONE-ELECTRON CORE MATRIX ELEMENTS          
    !* =====================================================================
    !c----------------------------------------------------------------------
    !c For MEH-PPV we add variable angle in the Call. 
    !c If I and J are junction atoms, then angle is the torsion
    !c otherwise is set to 1 before the call
    !c This because the twist phenomena seems to be interesting in MEH-PPV
    !c----------------------------------------------------------------------
    !AY  This function give the gradient of beta, but not complete. 
    !    This does not include a part of gradient of cosine(torsion angle)

    IMPLICIT NONE
    
    DOUBLE PRECISION, DIMENSION(3) :: fbetalin
    !input variables
    CHARACTER*1, INTENT(IN) :: spec1,spec2
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: vr
    DOUBLE PRECISION, INTENT(IN) :: angle,betalin
    !DOUBLE PRECISION, INTENT(IN) :: umub,linb2,r,betalin,angle
    
    !local variables
    DOUBLE PRECISION :: ftemp,dr,mu,b,r
    LOGICAL :: hit
    
    hit=.true.
    
    if(spec1.eq.'A' .and. spec2.eq.'A') then
       mu=umubcc
       b=linb2cc
       r=rcc
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'S') .or. (spec1.eq.'S'.and.spec2.eq.'A')) then
       mu=umubcs
       b=linb2cs
       r=rcs
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'N') .or. (spec1.eq.'N'.and.spec2.eq.'A')) then
       mu=umubcn
       b=linb2cn
       r=rcn
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'M') .or. (spec1.eq.'M'.and.spec2.eq.'A')) then
       mu=umubcn
       b=linb2cn
       r=rcn
       hit=.false.
    end if
    if((spec1.eq.'A'.and.spec2.eq.'O') .or. (spec1.eq.'O'.and.spec2.eq.'A')) then
       mu=umubco
       b=linb2co
       r=rco
       hit=.false.
    end if
    
    if(hit) then
       print *,'***WARNING: combo not defined in lindenberg approximation',spec1,spec2
       stop
    end if    

    dr=dsqrt(vr(1)**2+vr(2)**2+vr(3)**2)
    !ftemp=linb2*dexp(-umub*(dr-r)) - umub*betalin
    ftemp=b*dexp(-mu*(dr-r)) - mu*betalin
    !c---------------------------------------------------------------------
    !c For MEH-PPV we test if I and J are involved in junction torsion
    !c If yes we set angle = torsion. Otherwise is 1
    !c---------------------------------------------------------------------
    ftemp=ftemp*angle
    
    fbetalin=-vr/dr*ftemp
  END FUNCTION fbetalin

END MODULE QUANTUM
!-----------------------------------------------------------------------!

MODULE PPV
  
  IMPLICIT NONE
  SAVE

  !APW used
  INTEGER :: ntor_junc
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: junction
  DOUBLE PRECISION :: u_conf
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: pibond,torsion
  
  
END MODULE PPV

!----------------------------------------------------------------!

MODULE QCFIO
  
  IMPLICIT NONE
  SAVE
  
  INTEGER, PARAMETER :: gen_io=99
  INTEGER, PARAMETER :: mov_io=60
  INTEGER, PARAMETER :: mo_io=61
  INTEGER, PARAMETER :: stat_io=78
  INTEGER, PARAMETER :: tor_io=79
  INTEGER, PARAMETER :: bond_io=80
  INTEGER, PARAMETER :: gap_io=81
  INTEGER, PARAMETER :: osc_io=82
  INTEGER, PARAMETER :: temp_io=85
  INTEGER, PARAMETER :: tul_io=30
  INTEGER, PARAMETER :: fdotv_io=31
  INTEGER, PARAMETER :: fssh_io=32
  INTEGER, PARAMETER :: cicoef_io=62
  INTEGER, PARAMETER :: vel_io=83
  INTEGER, PARAMETER :: parm_io=7
  INTEGER, PARAMETER :: mol_io=19
  INTEGER, PARAMETER :: orb_io=84
  INTEGER, PARAMETER :: bo_io=86
  INTEGER, PARAMETER :: hop_io=87
  INTEGER, PARAMETER :: swap_io=88

END MODULE QCFIO
   
!----------------------------------------!

!MODULE class_atom
!  TYPE ATOM
!     INTEGER :: atype
!     DOUBLE PRECISION :: mass
!     DOUBLE PRECISION, DIMENSION(3) :: pos,vel,force,deriv
!     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dderiv
!  END type ATOM
!END MODULE class_atom
     
