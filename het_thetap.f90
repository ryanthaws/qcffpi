
!subroutine het_thetap (d,x,wt_ppv,th,iac2_ppv,e,dsx)
SUBROUTINE het_thetap(eng)
  
  USE commondata, ONLY:idsx
  USE classical
  USE qcfio, ONLY:gen_io
  
  IMPLICIT NONE
  
  !input variables
  DOUBLE PRECISION, INTENT(out) :: eng

  !function
!  AY
!  DOUBLE PRECISION :: costhe
  
  !local variables
  LOGICAL :: iarc
  INTEGER :: i,j,k,btype,jtype,ttype,hit,itheta
  INTEGER, DIMENSION(3) :: ind
  INTEGER, DIMENSION(2) :: itemp
  DOUBLE PRECISION :: r,rsq,tcub,f,df,ddf,dr1,ra,rasq,rb,rbsq,cthet,dsign
  
  DOUBLE PRECISION, DIMENSION(3) :: ftemp,dx,dy
  DOUBLE PRECISION, DIMENSION(3,2) :: dr
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ddr,ddvcos
  DOUBLE PRECISION, DIMENSION(nthetas) :: vtheta
  DOUBLE PRECISION, DIMENSION(3,2) :: dvcos
  
  do itheta=1,nthetas
     ind(1)=thetas(1,itheta)
     ind(2)=thetas(2,itheta)
     ind(3)=thetas(3,itheta)
     
     btype=pair_code(atype(ind(1)),atype(ind(3)))
     jtype=atype(ind(2))
     
     if(idsx.eq.1 .and. nb_sw.eq.1) then
        do k=1,3
           dx(k)=pos(k,ind(1))-pos(k,ind(3))
        end do
        
        rsq=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
        r=dsqrt(rsq)
        
        call difbond(dx,dr,ddr,r)
        
        !c----------------------------------------------------------------------
        !c     
        !c     the urey bradley potential function ( f=cub*(q-q0)**2 )
        !c     
        !c-----------------------------------------------------------------------
        !subroutine eurey(i1,j1,k1,bi,lt,iac_2)
        ttype=pair_code(theta_equiv(atype(ind(1))),theta_equiv(atype(ind(3))))*100+theta_equiv(atype(ind(2)))
        hit=0
        i=0
        do while(i.lt.ntheta_type .and. hit.eq.0)
           i=i+1
           if(theta_code(i).eq.ttype) hit=i
        end do
        if(hit.eq.0) then
           write(gen_io,*) '***Warning: theta parameters not defined for',&
                theta_equiv(atype(ind(1))),theta_equiv(atype(ind(2))),&
                theta_equiv(atype(ind(3))),'***'
           stop
        end if
        
        !APW 1004 corresponds to CCC and 1904 corresponds to CCB, why?
        if(ttype.eq.1004 .or. ttype.eq.1904) then
           tcub=0.0
        else
           tcub=theta_cub(2,hit)
        end if
        dr1=r-theta_cub(3,hit)
        f  =theta_cub(1,hit)*dr1**2+tcub*dr1
        df =2.d0*theta_cub(1,hit)*dr1+tcub
        ddf=2.d0*theta_cub(1,hit)
        
        !write(gen_io,*) ttype,f,df,ddf
        eng=eng+f
     end if
     
     if(t_sw.eq.1) then
        
        itemp(1)=ind(1)
        itemp(2)=ind(3)
        call drd(itemp,df,ddf,dr,ddr)
        
        do j=1,3
           dx(j)=pos(j,ind(1))-pos(j,ind(2))
           dy(j)=pos(j,ind(3))-pos(j,ind(2))
        end do
        rasq=length(dx)
        rbsq=length(dy)
        ra=dsqrt(rasq)
        rb=dsqrt(rbsq)
        cthet=costhe(dx,dy,ra,rb)
        if(cthet>1.d0) cthet=1.d0
        if(cthet<-1.d0) cthet=-1.d0
        ftemp(2)=dacos(cthet)
        vtheta(itheta)=ftemp(2)
        
        !c------------------------------------------------------------------
        !c     
        !c     the quadratic theta energy function ( f=ct*(th-th0)**2 )
        !c     
        !c------------------------------------------------------------------
        !subroutine het_tfunc (i1,j1,k1,at1,lt,iac_2)
    
        ttype=pair_code(theta_equiv(atype(ind(1))),theta_equiv(atype(ind(3))))*100+theta_equiv(atype(ind(2)))
        hit=0
        i=0
        do while(i.lt.ntheta_type .and. hit.eq.0)
           i=i+1
           if(theta_code(i).eq.ttype) hit=i
        end do
        
        if(hit.eq.0) then
           write(gen_io,*) '***Warning: theta parameters not defined for',&
                theta_equiv(atype(ind(1))),theta_equiv(atype(ind(2))),&
                theta_equiv(atype(ind(3))),'***'
           stop
        end if
        
        f  =theta_pms(1,hit)*(ftemp(2)-theta_pms(2,hit))**2
        df =2.d0*theta_pms(1,hit)*(ftemp(2)-theta_pms(2,hit))
        ddf=2d0*theta_pms(1,hit)
        
        !APW 1004 corresponds to CCC and 1904 corresponds to CCB
        !the theta term contains an extra linear term
        if(locate_atom('C').ne.0) then
           itemp(1)=pair_code(theta_equiv(locate_atom('C')),theta_equiv(locate_atom('C')))*100+theta_equiv(locate_atom('C'))
           if(ttype.eq.itemp(1)) then
              !if(ttype.eq.1004 .or. ttype.eq.1904) then
              f =f+theta_cub(2,hit)*(ftemp(2)-theta_pms(2,hit))
              df=df+theta_cub(2,hit)
           end if
           if(locate_atom('B').ne.0) then
              itemp(1)=pair_code(theta_equiv(locate_atom('C')),theta_equiv(locate_atom('B')))*100+theta_equiv(locate_atom('C'))
              if(ttype.eq.itemp(1)) then
                 !if(ttype.eq.1004 .or. ttype.eq.1904) then
                 f =f+theta_cub(2,hit)*(ftemp(2)-theta_pms(2,hit))
                 df=df+theta_cub(2,hit)
              end if
           end if
        end if
        eng=eng+f
        
        dsign=1.0
        
        iarc=.true.
        
        call difang(dx,dy,ra,rb,cthet,dvcos,ddvcos,dsign,iarc)
        !call difang (dx,bi,bk,ct,dc,ddc,dsign,newton,iarc)
        
        call dtd(ind,df,ddf,dvcos,ddvcos)
        
     end if
  end do

  if(idsx.eq.0) then
     deriv_theta(1,:,itheta) = dvcos(:,1)
     deriv_theta(2,:,itheta) = -(dvcos(:,1)+dvcos(:,2))
     deriv_theta(3,:,itheta) = dvcos(:,2)
  end if
  
END SUBROUTINE het_thetap

