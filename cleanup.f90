!APW cleanup memory

SUBROUTINE cleanup
  
  USE commondata
  USE qcfio
  USE PPV
  USE classical
  USE quantum
  
  IMPLICIT NONE


  !close all the output files
  close (mov_io)
  close(mo_io)
  close(stat_io)
  close(tor_io)
  close(bond_io)
  close(gap_io)
  close(osc_io)
  close(temp_io)
  close(orb_io)
  close(bo_io)
  
  if (log_exdyn) then
     close(tul_io)
     close(fdotv_io)
     close(fssh_io)
     close(hop_io)
  end if
  
  if (log_cicalc) close(cicoef_io)
  
  if(log_dumpvel) close(vel_io)

  
  

!deallocate arrays
  !ALLOCATED IN setup_vars.f90
  DEALLOCATE(deriv,dderiv,atype,amass,pos,vel,vel_l,pos_m,pos_l,fra)
  DEALLOCATE(bomat,mocoef,mocoef_old)
  if(log_cicalc) then
     DEALLOCATE(ciener,osc,ciener_l,osc_vec,cieig_old,cieig,ciind,ciind_old)
  end if
  if(log_exdyn) then
     DEALLOCATE(qr,utmp,aeainv)
  end if
  if(log_ppv) then
     DEALLOCATE(junction,inv_dih,junc_flg,inv_dih_nopi,torsion,ind_dih_nopi)
  end if

  DEALLOCATE(z_atom,deltaab,idxzvecq,idxzvecp,ibeta,pibond)
  if(log_exdyn) then
     DEALLOCATE(ctully,fexa,fexa_l,fdv,fdv_l,fdv_m,dar)
  end if

  
  !ALLOCATED IN read3d.f90
  DEALLOCATE(bonds)

  !ALLOCATED IN npar.f90
  !DEALLOCATE(theta_pms)
  DEALLOCATE(bond_pms,bond_code,theta_cub)
  DEALLOCATE(theta_code,phi_pms,phi_code,iatg,ncos,beta,beta_code)
  DEALLOCATE(phi_in,nonb_pms,nonb_code,pispec,alpha,gamma)
  DEALLOCATE(gamma_screen,alpha_mu)

  
  !ALLOCATED IN allpak.f90
  DEALLOCATE(thetas,phis,nonbs)
  !DEALLOCATE(deriv_bond,deriv_theta,deriv_phi)
  
  !ALLOCATED IN molein.f90
  DEALLOCATE(pi_iden,inv_pi,ind_nopi)
  
  write(gen_io,*) "exiting gracefully:  bye bye"
  close(gen_io)

END SUBROUTINE cleanup
