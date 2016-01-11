!!$c----------------------------------------------------------------------
!!$c     
!!$c     this modeule computes the phi angle cotributions to e , d and dd
!!$c     
!!$c----------------------------------------------------------------------
SUBROUTINE het_phip(eng)
  !subroutine het_phip(d_,x_,e,ph,wp_ppv,nc_ppv,iac2_ppv)
  
  USE commondata, ONLY:pi,istep,log_dderiv,idsx
  USE classical
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE 
  
  DOUBLE PRECISION :: eng
  
  LOGICAL :: iarc
  INTEGER :: i,j,k,iphi,btype,hit1,hit2,ipc,btype2,ii
  INTEGER, DIMENSION(4) :: itemp
  INTEGER, DIMENSION(4) :: ind
  INTEGER, DIMENSION(2) :: iind
  
  DOUBLE PRECISION :: ra,rb,rasq,rbsq,cthet,acthet,rij,rjk,rkl,ctijk,ctjkl,&
       stijk,stjkl,f,df,ddf
  DOUBLE PRECISION :: ct1,act1,ct2,act2,tz,dsign,stemp,add_f,add_df
  DOUBLE PRECISION, DIMENSION(4) :: ptemp,etemp,ftemp,ctemp
  DOUBLE PRECISION, DIMENSION(3) :: dxa,dxb,lvec,dftemp
  DOUBLE PRECISION, DIMENSION(3,3) :: tvec
  DOUBLE PRECISION, DIMENSION(3,4) :: dp
  DOUBLE PRECISION, DIMENSION(3,2) :: dvcos
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ddvcos
  DOUBLE PRECISION, DIMENSION(nphis) :: vphi
  DOUBLE PRECISION, DIMENSION(4,4) :: cro
  DOUBLE PRECISION, DIMENSION(3,6) :: torvec
  DOUBLE PRECISION, DIMENSION(3,4,4) :: dtp

!  AY
!  DOUBLE PRECISION :: costhe   !function
  
  
  ipc=0
  do iphi=1,nphis
     etemp=0.d0
     !d4  = 0.0
     !dd4 = 0.0
     
     cro=0.d0
     dtp=0.d0
     
     ind=phis(:,iphi)
     
     btype=pair_code(atype(ind(2)),atype(ind(3)))
     
     !APW bytype=11 corresponds to an A-H bond
     itemp(1)=locate_atom('A')
     itemp(2)=locate_atom('H')
     if(itemp(1).eq.0 .or. itemp(2).eq.0) then
        write(gen_io,*) '****WARNING: A or H are not defined in input, check het_phip.f90***'
        !stop
     end if
     if(iphi.gt.nphis_in) btype=pair_code(itemp(1),itemp(2))     

     !c:::  generate effective elements of phi vectors
     do i=1,3
        tvec(i,1)=pos(i,ind(1))-pos(i,ind(2))
        tvec(i,2)=pos(i,ind(3))-pos(i,ind(2))
        tvec(i,3)=pos(i,ind(3))-pos(i,ind(4))
     end do
    
     do i=1,3
        lvec(i)=dsqrt(length(tvec(:,i)))
     end do
     
     dxa(1)=tvec(2,1)*tvec(3,2)-tvec(3,1)*tvec(2,2)
     dxa(2)=tvec(3,1)*tvec(1,2)-tvec(1,1)*tvec(3,2)
     dxa(3)=tvec(1,1)*tvec(2,2)-tvec(2,1)*tvec(1,2)
     dxb(1)=tvec(3,2)*tvec(2,3)-tvec(2,2)*tvec(3,3)
     dxb(2)=tvec(1,2)*tvec(3,3)-tvec(3,2)*tvec(1,3)
     dxb(3)=tvec(2,2)*tvec(1,3)-tvec(1,2)*tvec(2,3)
     
     !c:::  calculate the length and angle between effective vectors
     rasq=length(dxa)
     rbsq=length(dxb)
     ra=dsqrt(rasq)
     rb=dsqrt(rbsq)
     
     cthet=costhe(dxa,dxb,ra,rb)
     
     if(cthet-1.0e-19 .ge. 1.0) cthet = 1.0-1.0e-19
     if(cthet+1.0e-19 .le. -1.0) cthet = -1.0 + 1.0e-19
     !if(cthet>1.d0) cthet=1.d0
     !if(cthet<-1.d0) cthet=-1.d0
     acthet=dacos(cthet)
     vphi(iphi)=acthet
     
     !c:::  adjust dsign of phi angle
     stemp=0.d0
     
     stemp=stemp+tvec(1,2)*(dxa(3)*dxb(2)-dxa(2)*dxb(3))
     stemp=stemp+tvec(2,2)*(dxa(1)*dxb(3)-dxa(3)*dxb(1))
     stemp=stemp+tvec(3,2)*(dxa(2)*dxb(1)-dxa(1)*dxb(2))
     
     if(stemp.lt.0.d0) acthet=-acthet
     acthet=pi-acthet
     etemp(4)=dcos(acthet)
     ftemp(1)=-ftemp(1)
     
     if(idsx.eq.0) then
        !call dpxs(bi,bk,dp,a,b_,c,dx,dx1)
        !c-----------------------------------------------------------
        !c
        !c     compute dphi/dx for small angles with wilson' s formula
        !c
        !c----------------------------------------------------------------------
        !subroutine dpxs(bi,bk,dp,a,b,c,dx,dx1)
        rij=dsqrt(lvec(1))
        rjk=dsqrt(lvec(2))
        rkl=dsqrt(lvec(3))
        
        ctijk=DOT_PRODUCT(tvec(:,1),tvec(:,2))/(rij*rjk)
        ctjkl=DOT_PRODUCT(tvec(:,2),tvec(:,3))/(rjk*rkl)
        stijk=dsqrt(1.d0-ctijk**2)
        stjkl=dsqrt(1.d0-ctjkl**2)
        do i=1,3
           dp(i,1)=dxa(i)/(ra*rij*stijk)
           dp(i,2)=-((rjk-rij*ctijk)/(rij*rjk*stijk))*dxa(i)/ra&
                -ctjkl*dxb(i)/(rb*rjk*stjkl)
           dp(i,3)=-((rjk-rkl*ctjkl)/(rkl*rjk*stjkl))*dxb(i)/rb&
                -ctijk*dxa(i)/(ra*rjk*stijk)
           dp(i,4)=dxb(i)/(rb*rkl*stjkl)
        end do
        
        !APW this hasn't been checked, previous version of code didn't allocate enough space for evaluation of normal modes for phi
        do i=1,4
           do j=1,3
              deriv_phi(i,j,iphi)=dp(j,i)
           end do
        end do
     else
        !calculate the energy and its derivatives
        !attention please. de/dx=de/d(cos(phi))*d(cos(phi)/dx ,not
        !de/dx=de/d(phi)*d(phi)/dx
!!$            e_4=dcos(ap)
!!$            d4=1.
!!$            dd4=0.
        f=0.d0
        df=0.d0
        ddf=0.d0
        !cc--> differenza sta qui
        !call het_ththph(i1,j1,k1,l1,ap,ic,d_,x_,nc_ppv,iac2_ppv,e)
        !c----------------------------------------------------------------------
        !c
        !c     this subroutine evaluates the energy and derivatives of
        !c     function of the form f(phi)*g(th1)*g(th1)
        !c
        !c----------------------------------------------------------------------
        !subroutine het_ththph(i1,j1,k1,l1,ap,ic,d_,x_,nc_,iac_,e)

        !APW evidently exclude A-H bonds from this subroutine
        if(locate_atom('A').ne.0 .and. locate_atom('H').ne.0) then
           itemp(1)=pair_code(locate_atom('A'),locate_atom('H'))
           if(btype.ne.itemp(1)) then
              iarc=.true.
              ftemp(1)=0.0
              
              ct1=costhe(tvec(:,1),tvec(:,2),lvec(1),lvec(2))
              if(ct1.gt.1.d0) ct1=1.d0
              if(ct1.lt.-1.d0) ct1=-1.d0
              act1=dacos(ct1)
              
              ct2=costhe(tvec(:,2),tvec(:,3),lvec(2),lvec(3))
              if(ct2.gt.1.d0) ct2=1.d0
              if(ct2.lt.-1.d0) ct2=-1.d0
              act2=dacos(ct2)
              
              !c-----------------------------------------------------------------------
              !c
              !c     the theta*theta*phi energy fc (f=ctt*(th1-th0)*(th2-th0)*dcos(phi)
              !c
              !c-----------------------------------------------------------------------
              hit1=0
              !hit = npp                 !note that hit = npp instead of = 0
              i=0
              do while(i.lt.nphi_type .and. hit1.eq.0)
                 i=i+1
                 if(btype.eq.phi_code(i)) hit1=i
              end do
              if(hit1.eq.0) then
                 write(gen_io,*) '***Warning: phi parameters not defined for',&
                      atype(ind(1)),atype(ind(2)),atype(ind(3)),atype(ind(4))
                 stop
              end if
              
              etemp(1)=phi_pms(2,hit1)
              if(dabs(etemp(1)).gt.1.e-4) then
                 
                 !APW for CC (btype=10) or CB (btype=19)
                 if(locate_atom('C').ne.0) then
                    itemp(1)=pair_code(locate_atom('C'),locate_atom('C'))
                    itemp(3)=atype(ind(1))
                    itemp(4)=atype(ind(4))
                    if(btype.eq.itemp(1)) then 
                       if(species(itemp(3)).eq.'H'.and.species(itemp(4)).eq.'H') then
                          !if(atype(ind(1)).eq.1 .and.atype(ind(4)).eq.1) then
                          etemp(1)=phi_pms(3,hit1)
                       else if(species(itemp(3)).eq.'H' .or. species(itemp(4)).eq.'H') then
                          etemp(1)=-dsqrt(phi_pms(2,hit1)*phi_pms(3,hit1))
                       end if
                    end if
                 end if

                 if(locate_atom('C').ne.0 .and. locate_atom('B').ne.0) then
                    itemp(1)=pair_code(locate_atom('C'),locate_atom('B'))
                    itemp(3)=atype(ind(1))
                    itemp(4)=atype(ind(4))
                    if(btype.eq.itemp(1)) then 
                       if(species(itemp(3)).eq.'H'.and.species(itemp(4)).eq.'H') then
                          !if(atype(ind(1)).eq.1 .and.atype(ind(4)).eq.1) then
                          etemp(1)=phi_pms(3,hit1)
                       else if(species(itemp(3)).eq.'H' .or. species(itemp(4)).eq.'H') then
                          etemp(1)=-dsqrt(phi_pms(2,hit1)*phi_pms(3,hit1))
                       end if
                    end if
                 end if
              
                 tz=1.911 !APW ???
                 !APW NA, NN, or AA respectively
                 if(locate_atom('N').ne.0) then
                    itemp(1)=pair_code(locate_atom('N'),locate_atom('N'))
                    if(btype.eq.itemp(1)) tz=2.094
                 end if
                 if(locate_atom('A').ne.0) then
                    itemp(1)=pair_code(locate_atom('A'),locate_atom('A'))
                    if(btype.eq.itemp(1)) tz=2.094
                 end if
                 if(locate_atom('N').ne.0 .and. locate_atom('A').ne.0) then
                    itemp(1)=pair_code(locate_atom('N'),locate_atom('N'))
                    if(btype.eq.itemp(1)) tz=2.094
                 end if
                 !if(btype.eq.13 .or. btype.eq.6 .or. btype.eq.15) tz=2.094
                 etemp(2)=act1-tz
                 etemp(3)=act2-tz
                 f=etemp(1)*(act1-tz)*(act2-tz)*etemp(4)
                 dftemp(1)=etemp(1)*(act2-tz)*etemp(4)
                 dftemp(2)=etemp(1)*(act1-tz)*etemp(4)
                 
                 eng=eng+f
                 dsign=1.0
                 
                 call difang(tvec(:,1),tvec(:,2),lvec(1),lvec(2),ct1,dvcos,ddvcos,dsign,iarc)
                 df=dftemp(1)
                 itemp(1)=ind(1)
                 itemp(2)=ind(2)
                 itemp(3)=ind(3)
                 
                 call dtd(itemp,df,ddf,dvcos,ddvcos)
                 
                 !end if
                 !APW not sure what this is
                 if(log_dderiv) then  
                    do i=1,3
                       dtp(i,1,2)=dvcos(i,1)
                       dtp(i,2,2)=-(dvcos(i,1)+dvcos(i,2))
                       dtp(i,3,2)=dvcos(i,2)
                    end do
                 end if
                 
                 call difang(tvec(:,2),tvec(:,3),lvec(2),lvec(3),ct2,dvcos,ddvcos,dsign,iarc)
                 
                 dvcos=-dvcos
                 df=dftemp(2)
                 itemp(1)=ind(2)
                 itemp(2)=ind(3)
                 itemp(3)=ind(4)
                 
                 call dtd(itemp,df,ddf,dvcos,ddvcos)
                 if(log_dderiv) then
                    do i=1,3
                       dtp(i,2,3)=dvcos(i,1)
                       dtp(i,3,3)=-(dvcos(i,1)+dvcos(i,2))
                       dtp(i,4,3)=dvcos(i,2)
                    end do
                 end if
                 !:::  cro are the cross term derivatives :
                 !1 for b 
                 !2 for th1 
                 !3 for th2 
                 !4 for phi
                 
                 cro(2,3)=etemp(1)*etemp(4)
                 cro(2,4)=etemp(1)*etemp(3)
                 cro(3,4)=etemp(1)*etemp(2)
              end if
           end if
        end if !APW if btype AH exists
        !c be carefull qua entra in gioco la matrice H. che non hai cccccc
        !c ecco perche viene diverso -------------------------------------
        !c ma non capisco
        !APW roughly translated (google translate)
        !Be carefull here that touches on the matrix H. you have not here is because divers--------but I do not understand
        
        !call ephi(ap,ic,i1,j1,k1,l1,d_,x_,nc_ppv,iac2_ppv)
        !c---------------------------------------------------------------------
        !
        !     calculates the torsional energy beetwen the atoms i1,j1,k1,l1
        !     and the corresponding pi-term f(b)*g(phi). the respective
        !     derivatives are also calculated.
        !
        !----------------------------------------------------------------------
        !subroutine ephi(ap,ic,i1,j1,k1,l1,d_,x_,nc_,iac_)
        
        !APW for C-C bonds
        itemp(1)=1000
        if(locate_atom('C').ne.0) itemp(1)=pair_code(locate_atom('C'),locate_atom('C'))
        if(btype.eq.itemp(1)) then
           itemp(1)=ipc
           ipc=ind(2)*100+ind(3)
           itemp(2)=ind(3)*100+ind(2)
           if(ipc.eq.itemp(1) .or. itemp(2).eq.itemp(1)) then
              f=0.d0
              df=0.d0
              ddf=0.d0
           end if
        else
           
           !APW below all seem aimed at het_pfunc, I'm going to bypass all the if
           !only one cccc torsion angle around a cc bond enters the phi func.,
           !but all 9 are used in theta*theta*phi functions
           !hit1=0
           !hit2=0
           !do i=1,npi
           !if(ind(2).eq.pi_iden(i)) hit1=i
           !if(ind(3).eq.pi_iden(i)) hit2=i
           !end do
           
!!$      if(ic.eq.10) then
!!$         ipc1=ipc
!!$         ipc=j1*100+k1
!!$         ipc2=k1*100+j1
!!$         if(ipc.eq.ipc1.or.ipc2.eq.ipc1) then
!!$            f   = 0.
!!$            df  = 0.
!!$            ddf = 0.
!!$            return
!!$         endif
!!$      endif
!!$
!!$      ih=0
!!$      jh=0
!!$      do i=1,nap
!!$         if(j1.eq.nc_(i)) ih=i
!!$         if(k1.eq.nc_(i)) jh=i
!!$      enddo
!!$
!!$      if(ic.eq.11) then
!!$         ic1=11
!!$         call het_pfunc(ic1,j1,k1,ap,lp)
!!$         return
!!$      endif
!!$
!!$      if(ih.eq.0.or.jh.eq.0) then
!!$         call het_pfunc(ic,j1,k1,ap,lp)
!!$         return
!!$      endif
!!$
!!$      call het_pfunc(ic,j1,k1,ap,lp)
           !-------------------------------------------------------------------
           !c     
           !c     phi potential function of the form f=cp*dcos(n*ph)
           !c------------------------------------------------------------------
           !subroutine het_pfunc(ic,j1,k1,ap,lp)
           if(p_sw.eq.1) then
              btype2=btype
              if(iphi.gt.nphis_in) then
                 if(locate_atom('H').eq.0 .or. locate_atom('A').eq.0) then
                    write(gen_io,*)'***WARNING: H or A not defined in het_phip.f90***'
                    !stop
                 end if
                 !btype2=pair_code(locate_atom('H'),locate_atom('A')) !APW H-A
              end if
              hit1=0
              i=0
              do while(i.lt.nphi_type .and. hit1.eq.0)
                 i=i+1
                 if(btype2.eq.phi_code(i)) hit1=i
              end do
                 
              if(hit1.eq.0) then
                 print *,btype2,phi_code(3),i,pair_code(locate_atom('A'),locate_atom('A'))
                 print *,iphi,nphis_in
                 write(gen_io,*) '***Warning: phi parameters not defined for ',&
                      species(atype(ind(2))),species(atype(ind(3))),'***'
                 stop
              end if
              
              itemp(1)=ncos(hit1)
              ftemp(3)=dcos(acthet)
              itemp(2)=1
              ftemp(1)=phi_pms(1,hit1)
              
              add_f=0.0d0
              add_df=0.0d0
              do ii=1,nphib
                 iind(1)=phib_pairs(1,ii)
                 iind(2)=phib_pairs(2,ii)
                 if((ind(2).eq.iind(1).and.ind(3).eq.iind(2)).or.(ind(3).eq.iind(1).and.ind(2).eq.iind(2))) then
                    !ctemp(1)=1.00*(0.3587)
                    !ctemp(2)=0.25*(3.2149)
                    !ctemp(3)=1.00*(0.0489)
                    !ctemp(4)=0.25*(0.2884)
                    !ctemp(2)=0.25*(3.2249)
                    
                    ctemp(2)=0.25*phib_pms(1,ii)
                    ctemp(4)=0.25*phib_pms(2,ii)
                    
                    add_f =add_f + ctemp(2)*dcos(acthet*2)
                    add_df=add_df+ 4.0d0*ctemp(2)*ftemp(3)
                    
                    add_f =add_f + ctemp(4)*dcos(acthet*4)
                    add_df=add_df+ ctemp(4)*16.0d0*ftemp(3)*(2*ftemp(3)**2-1)
                    
                    !if((ind(1).eq.1.and.ind(4).eq.9).or.(ind(1).eq.9.and.ind(4).eq.1)) then
                    !print *,acthet*57.29,cthet
                    !add_f =add_f + ctemp(1)*dcos(acthet)
                    !add_df=add_df+ ctemp(1)
                    
                    !add_f =add_f + ctemp(3)*dcos(acthet*3)
                    !add_df=add_df+ ctemp(3)*(12.0d0*ftemp(3)**2-3.0d0)
                    !end if
                    
                    ftemp(1)=0.0d0
                    itemp(1)=3
                 end if
              end do
              ftemp(2)=acthet*itemp(1)
              if(itemp(1).eq.2) itemp(2)=-1
              
              if(locate_atom('C').ne.0) then
                 if(btype.ne.pair_code(locate_atom('C'),locate_atom('C'))) ftemp(1)=ftemp(1)/iatg(hit1)
              end if
              !although iatg for torsion angle around cc
              !is not 1, only one angle is used in the
              !phi-functions,but all 9 enter the
              !the*the*phi funct
              ftemp(4)=ftemp(1)*dcos(ftemp(2))*itemp(2)
              f=ftemp(1)+ftemp(4)+add_f

              df=add_df
              if(itemp(1).eq.1) then
                 df=df+ftemp(1)*itemp(2)
                 ddf=0.d0
              else if(itemp(1).eq.2) then
                 df=df+ftemp(1)*itemp(2)*4.d0*ftemp(3)
                 ddf=ftemp(1)*itemp(2)*4.d0
              else if(itemp(1).eq.3) then
                 df=df+ftemp(1)*itemp(2)*(12.d0*ftemp(3)**2-3.d0)
                 ddf=24.d0*ftemp(1)*itemp(2)*ftemp(3)
              end if
              if(itemp(1).ne.3 .and. iphi.le.nphis_in) then
                 f=f+phi_pms(3,hit1)*ftemp(3)
                 df=df+phi_pms(3,hit1)
              end if
           end if
           !APW translated now enter into this comment but the electronics
           !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           !ccccc ora lo commento ma questo entra in parte elettronica ----cccc
           !cccc
           !c      call het_vphi(ap,j1,k1,ih,jh,lp,ic)                  
           !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end if
        
        eng=eng+f
        df=df+etemp(1)*etemp(2)*etemp(3)*1.0
        ddf=ddf+etemp(1)*etemp(2)*etemp(3)*0.0
        !e = e + f
        !df=df+e_1*e_2*e_3*d4
        !ddf=ddf+e_1*e_2*e_3*dd4
        
        !c     
        !c     calculate the effective derivatives
        !c
        iarc=.false.
        dsign=stemp
        call difang(dxa,dxb,ra,rb,cthet,dvcos,ddvcos,dsign,iarc)
        dvcos=-dvcos
        ddvcos=-ddvcos
        
        !call dxvec(i1,j1,k1,l1,x_,v)
        !c---------------------------------------------------------------------
        !c
        !c     build the vectors invoved in the torsion angle calculation.
        !c
        !c---------------------------------------------------------------------
        !subroutine dxvec(i1,j1,k1,l1,x,v)
        
        torvec(:,1)=pos(:,ind(3))-pos(:,ind(2))
        torvec(:,2)=pos(:,ind(1))-pos(:,ind(3))
        torvec(:,3)=pos(:,ind(3))-pos(:,ind(4))
        torvec(:,4)=pos(:,ind(2))-pos(:,ind(1))
        torvec(:,5)=pos(:,ind(4))-pos(:,ind(2))
        torvec(:,6)=-torvec(:,1)
        
        !call dpx(dp,v,dc)
        !c-------------------------------------------
        !c
        !c     calculate 'dp',i.e.d(phi)/dx, from 'dc',i.e. the derivatives of
        !c     cos(phi) w.r.t. the components of the two normal vectors that
        !c     define the cosinus.
        !c-------------------------------------------
        !subroutine dpx(dp,v,dc)
        
        !if(iphi.eq.20) then
        call dpd(ind,df,ddf,dvcos,ddvcos,torvec,cro,dtp)
        !end if
!!$            call dpd(ii,d_,dp,df,newton)
        
     end if
  end do
end SUBROUTINE het_phip

