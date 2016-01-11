!c======================================================================
!c CONNECT subroutine that creates arrays of the chemical bonds
!c======================================================================
SUBROUTINE connect()
  
  USE commondata, ONLY: nnuc
  USE classical, ONLY: bonds,nbonds,atype,species,locate_atom
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  INTEGER :: i,j,k,numhydro,itemp
  INTEGER ::  iatom,itype,numbond
  
  
  do iatom=1,nnuc
     itype=atype(iatom)
     numbond=0
     do k=1,nbonds
        if(bonds(1,k) .eq. iatom) then
           numbond = numbond+1
           itemp=atype(bonds(2,k))
           if(species(itemp).eq.'X') numbond=numbond-1
           if(species(itemp).eq.'Z') numbond=numbond-1
           !APW generalize atype
           !if(atype(bonds(2,k)).eq.12) numbond=numbond-1
           !if(atype(bonds(2,k)).eq.14) numbond=numbond-1
        end if
        if(bonds(2,k) .eq. iatom) then
           numbond = numbond+1
           itemp=atype(bonds(1,k))
           if(species(itemp).eq.'X') numbond=numbond-1
           if(species(itemp).eq.'Z') numbond=numbond-1
           !APW generalize atype
           !if(atype(bonds(1,k)).eq.12) numbond=numbond-1
           !if(atype(bonds(1,k)).eq.14) numbond=numbond-1
        end if
     end do
     
     !c<<<<<<<<<<<<<<<<<<<<<   atom is n     >>>>>>>>>>>>>>>>>>>>>c
     if(species(itype).eq.'N'.or.species(itype).eq.'M') then
        if(species(itype).eq.'N') then
           !APW generalize atype
           !if(itype.eq.3 .or. itype.eq.8) then
           !if(itype .eq. 3) then
           if(numbond .eq. 3) then
              write(gen_io,*) iatom,' nitrogen is N'
           else
              write(gen_io,*) iatom,' **Unknown nitrogen type**'
              !stop
           end if
        end if
        !c<<<<<<<<<<<<<<<<<<<<<   atom is o or q    >>>>>>>>>>>>>>>>>>>>>c
     else if(species(itype).eq.'O'.or.species(itype).eq.'Q') then
        !APW generalize atype
        !else if(itype.eq.2 .or. itype.eq.9) then
        if(numbond .eq. 1) then
           write(gen_io,*) iatom,' oxygen is O'
           atype(iatom)=locate_atom('O')
           !atype(iatom)=2
        else if(numbond .eq. 2) then
           write(gen_io,*) iatom,' oxygen is Q'
           atype(iatom)=locate_atom('Q')
           !atype(iatom)=9
        else
           write(gen_io,*) iatom,' **Unknown oxygen type**'
           stop
        end if
        !c<<<<<<<<<<<<<<<<<<<<<   atom is a,b or c  >>>>>>>>>>>>>>>>>>>
        
     else if(species(itype).eq.'C'.or.species(itype).eq.'A'.or.species(itype).eq.'B') then
        
        !else if(itype.eq.4 .or. itype.eq.5 .or. itype.eq.6) then
        if(numbond .eq. 3) then
           atype(iatom) = locate_atom('A')
           write(gen_io,*) iatom, ' carbon is A but not changed'
        else if(numbond.eq.4) then
           numhydro=0
           do k=1,nbonds
              if(bonds(1,k).eq.iatom)then
                 itemp=atype(bonds(2,k))
                 if(species(itemp).eq.'H') numhydro=numhydro+1
                 !if(atype(bonds(2,k)).eq.1) numhydro = numhydro+1
              end if
              if(bonds(2,k).eq.iatom)then
                 itemp=atype(bonds(1,k))
                 if(species(itemp).eq.'H') numhydro=numhydro+1
                 !if(atype(bonds(1,k)).eq.1) numhydro = numhydro+1
              end if
           end do
           if(numhydro.eq.2) then
              write(gen_io,*) iatom, 'carbon is C (-CH2-)'
              atype(iatom) = locate_atom('C')
              !atype(iatom) = 4
           else if(numhydro.eq.3) then
              write(gen_io,*) iatom, 'carbon is B (-CH3)'
              atype(iatom) = locate_atom('B')
              !atype(iatom) = 6
           else if(numhydro.eq.1) then !AY 2011.01, allow C to one H-bond on sp3
              write(gen_io,*) iatom, 'carbon is C (-CH-)'
              atype(iatom) = locate_atom('C')
           else
              write(gen_io,*) iatom,' 1***unknown type of carbon***'
              stop
           end if
        else
           write(gen_io,*) iatom, ' 2***unknown type of carbon***'
        end if
        !c<<<<<<<<<<<<<<<<<<<< S atom >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     else if(species(itype).eq.'S' .or. species(itype).eq.'T') then
        write(gen_io,*) iatom, 'Sulfer is S'
        !c<<<<<<<<<<<<<<<<<<<< H atom >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     else if(species(itype).eq.'H' .or. species(itype).eq.'D') then
        !else if(itype.eq.1 .or. itype.eq.7) then
        !c     Do nothing
        write(gen_io,*) iatom, 'Hydrogen or D'
        !c<<<<<<<<<<<<<<<<<<<< girifalco atoms >>>>>>>>>>>>>>>>>>>>>>>>
        !c took over classification "m" (nitrogen w/ 2 bonds) to create
        !c girifalco molecule type
        !APW changed girifalco molecule to type G
     else if(species(itype).eq.'G') then
        !else if (itype .eq. 8) then
        write(gen_io,*) iatom,' girifalco molecule is G'
     else if(species(itype).eq.'E') then
        !else if (itype .eq. 8) then
        write(gen_io,*) iatom,' is coarse-grained aromatic carbond E'
     else
        write(gen_io,*) iatom, '***unknown atom type, itype=',itype,'***'
        numbond = 0
        do k=1,nbonds
           if(bonds(1,k).eq.iatom .or. bonds(2,k).eq.iatom) then
              numbond = numbond+1
           end if
        end do
!!$c:::  The following settings need to be improved.
        if(numbond.eq.1) then
           atype(iatom) = locate_atom('H')
           !atype(iatom) = 1
           write(gen_io,*) iatom, 'set to H'
        else if(numbond.eq.2) then
           atype(iatom) = locate_atom('A')
           !atype(iatom) = 4
           write(gen_io,*) iatom, 'set to carbon A'
           stop
        else if(numbond.eq.3) then
           atype(iatom) = locate_atom('A')
           !atype(iatom) = 4
           write(gen_io,*) iatom, 'set to carbon A'
           stop
        else if(numbond.eq.4) then
           atype(iatom) = locate_atom('A')
           !atype(iatom) = 4
           write(gen_io,*) iatom, 'set to carbon A'
           stop
        else
           atype(iatom) = locate_atom('A')
           !atype(iatom) = 4
           write(gen_io,*) iatom, 'set to carbon A'
           stop
        end if
     end if
  end do

END SUBROUTINE connect
