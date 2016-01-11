MODULE nrtype
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(SP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(DP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp
END MODULE nrtype


MODULE ran_state
  USE nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
  INTEGER(K4B), SAVE :: lenran=0, seq=0
  INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
  INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
  INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
       nran,mran,ranv
  REAL(SP), SAVE :: amm
  INTERFACE ran_hash
     MODULE PROCEDURE ran_hash_s, ran_hash_v
  END INTERFACE
CONTAINS
  !BL
  SUBROUTINE ran_init(length)
    USE nrtype; USE nrutil, ONLY : arth,nrerror,reallocate
    IMPLICIT NONE
    INTEGER(K4B), INTENT(IN) :: length
    INTEGER(K4B) :: new,j,hgt
    if (length < lenran) RETURN
    hgt=hg
    if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
    if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
    if (hgt+1-hgng /= 0)    call nrerror('ran_init: arith assump 3 fails')
    if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
    if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
    if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
    if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
    if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
    if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
    if (lenran > 0) then
       ranseeds=>reallocate(ranseeds,length,5)
       ranv=>reallocate(ranv,length-1)
       new=lenran+1
    else
       allocate(ranseeds(length,5))
       allocate(ranv(length-1))
       new=1
       amm=nearest(1.0_sp,-1.0_sp)/hgng
       if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            call nrerror('ran_init: arth assump 10 fails')
    end if
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
    do j=1,4
       call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    end do
    where (ranseeds(new:,1:3) < 0) &
         ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
    where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
    if (new == 1) then
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
       rans=nran0
    end if
    if (length > 1) then
       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran
    end if
    lenran=length
  END SUBROUTINE ran_init
  !BL
  SUBROUTINE ran_deallocate
    if (lenran > 0) then
       deallocate(ranseeds,ranv)
       nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
       lenran = 0
    end if
  END SUBROUTINE ran_deallocate
  !BL
  SUBROUTINE ran_seed(sequence,size,put,get)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: sequence
    INTEGER, OPTIONAL, INTENT(OUT) :: size
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
    if (present(size)) then
       size=5*lenran
    else if (present(put)) then
       if (lenran == 0) RETURN
       ranseeds=reshape(put,shape(ranseeds))
       where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
       where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
    else if (present(get)) then
       if (lenran == 0) RETURN
       ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
       get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
       call ran_deallocate
       seq=sequence
    end if
  END SUBROUTINE ran_seed
  !BL
  SUBROUTINE ran_hash_s(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), INTENT(INOUT) :: il,ir
    INTEGER(K4B) :: is,j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_s
  !BL
  SUBROUTINE ran_hash_v(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
    INTEGER(K4B), DIMENSION(size(il)) :: is
    INTEGER(K4B) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_v
END MODULE ran_state



MODULE nr
  INTERFACE ran2
     SUBROUTINE ran2_s(harvest)
       USE nrtype
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran2_s
     !BL
     SUBROUTINE ran2_v(harvest)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran2_v
  END INTERFACE
  
END MODULE nr

  SUBROUTINE ran2_s(harvest)
    USE nrtype
    USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran0,jran0,kran0,nran0,mran0,rans
    IMPLICIT NONE
    REAL(SP), INTENT(OUT) :: harvest
    !Lagged Fibonacci generator combined with a Marsaglia shift sequence and a linear congruential generator. Returns as harvest a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). This generator has the same calling and initialization conventions as Fortran 90’s random number routine. Use ran seed to initialize or reinitialize to a particular sequence. The period of this generator is about 8.5 × 1037 , and it fully vectorizes. Validity of the integer model assumed by this generator is tested at initialization.
    if (lenran < 1) call ran_init(1)
    !Initialization routine in ran state.
    rans=iran0-kran0
    !Update Fibonacci generator, which
    if (rans < 0) rans=rans+2147483579_k4b
    !has period p2 + p + 1, p = 2^31 − 69
    iran0=jran0
    jran0=kran0
    kran0=rans
    nran0=ieor(nran0,ishft(nran0,13))
    !Update Marsaglia shift sequence with period 2^32 − 1.
    nran0=ieor(nran0,ishft(nran0,-17))
    nran0=ieor(nran0,ishft(nran0,5))
    rans=iand(mran0,65535)
    !Update the sequence m ← 69069m + 820265819 mod 2^32 using shifts instead of multiplies. Wrap-around addition (tested at initialization) is used.
    mran0=ishft(3533*ishft(mran0,-16)+rans,16)+ &
         3533*rans+820265819_k4b
    rans=ieor(nran0,kran0)+mran0
    !Combine the generators.
    harvest=amm*merge(rans,not(rans), rans<0 )
    !Make the result positive definite (note that amm is negative).
  END SUBROUTINE ran2_s
  
  SUBROUTINE ran2_v(harvest)
    USE nrtype
    USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran,jran,kran,nran,mran,ranv
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
    INTEGER(K4B) :: n
    n=size(harvest)
    if (lenran < n+1) call ran_init(n+1)
    ranv(1:n)=iran(1:n)-kran(1:n)
    where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
    iran(1:n)=jran(1:n)
    jran(1:n)=kran(1:n)
    kran(1:n)=ranv(1:n)
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
    ranv(1:n)=iand(mran(1:n),65535)
    mran(1:n)=ishft(3533*ishft(mran(1:n),-16)+ranv(1:n),16)+ &
         3533*ranv(1:n)+820265819_k4b
    ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)
    harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
  END SUBROUTINE ran2_v
  

SUBROUTINE gasdev_s(harvest)
  USE nrtype
  USE nr, ONLY : ran2

  IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: harvest
  !Returns in harvest a normally distributed deviate with zero mean and unit variance, using ran2 as the source of uniform deviates.
  REAL(SP) :: rsq,v1,v2
  REAL(SP), SAVE :: g
  LOGICAL, SAVE :: gaus_stored=.false.
  if (gaus_stored) then
     !We have an extra deviate handy,
     harvest=g
     !so return it,
     gaus_stored=.false.
     !and unset the flag.
  else
     !We don’t have an extra deviate handy, so
     do
        call ran2(v1)
        !pick two uniform numbers in the square ex-
        call ran2(v2)
        !tending from -1 to +1 in each direction,
        v1=2.0_sp*v1-1.0_sp
        v2=2.0_sp*v2-1.0_sp
        rsq=v1**2+v2**2
        !see if they are in the unit circle,
        if (rsq > 0.0 .and. rsq < 1.0) exit
     end do
     !otherwise try again.
     rsq=sqrt(-2.0_sp*log(rsq)/rsq)
     !Now make the Box-Muller transformation to
     harvest=v1*rsq
     !get two normal deviates. Return one and
     g=v2*rsq
     !save the other for next time.
     gaus_stored=.true.
     !Set flag.
  end if
END SUBROUTINE gasdev_s

SUBROUTINE gasdev_v(harvest)
  USE nrtype; USE nrutil, ONLY : array_copy
  USE nr, ONLY: ran2
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
  REAL(SP), DIMENSION(size(harvest)) :: rsq,v1,v2
  REAL(SP), ALLOCATABLE, DIMENSION(:), SAVE :: g
  INTEGER(I4B) :: n,ng,nn,m
  INTEGER(I4B), SAVE :: last_allocated=0
  LOGICAL, SAVE :: gaus_stored=.false.
  LOGICAL, DIMENSION(size(harvest)) :: mask
  n=size(harvest)
  if (n /= last_allocated) then
     if (last_allocated /= 0) deallocate(g)
     allocate(g(n))
     last_allocated=n
     gaus_stored=.false.
  end if
  if (gaus_stored) then
     harvest=g
     gaus_stored=.false.
  else
     ng=1
     do
        if (ng > n) exit
        call ran2(v1(ng:n))
        call ran2(v2(ng:n))
        v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
        v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
        rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
        mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
        call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
        v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
        rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
        ng=ng+nn
     end do
     rsq=sqrt(-2.0_sp*log(rsq)/rsq)
     harvest=v1*rsq
     g=v2*rsq
     gaus_stored=.true.
  end if
END SUBROUTINE gasdev_v

