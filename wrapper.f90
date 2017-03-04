
module wrapper

use utils
use chebyshev
use odesolve
use kummer
use, intrinsic :: iso_c_binding


integer                       :: nints, nintsq,k
double precision, allocatable :: abq(:,:),valsq(:,:)

double precision, allocatable :: xscheb(:),whtscheb(:),chebintl(:,:),chebintr(:,:), &
   ucheb(:,:),vcheb(:,:)

double precision, allocatable :: ab(:,:),alpha(:,:),alphap(:,:),alphapp(:,:)

contains

subroutine qfun(t,val)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: t
double precision, intent(out) :: val

call chebpw_eval(nintsq,abq,k,xscheb,valsq,t,val)

end subroutine



subroutine  phase_function(kin,nintsin,abin,qvals)
implicit double precision (a-h,o-z)

integer, intent(in)          :: kin, nintsin
double precision, intent(in) :: abin(2,nintsin),qvals(kin,nintsin)



k = kin
call chebexps(k,xscheb,whtscheb,ucheb,vcheb,chebintl,chebintr)


allocate(abq(2,nintsin),valsq(k,nintsin))
nintsq = nintsin
abq    = abin
valsq  = qvals



eps = 1.0d-13
a   = abin(1,1)
b   = abin(2,nintsin)


call kummer_adap(eps,a,b,qfun,k,xscheb,chebintl,chebintr,ucheb, &
   nints,ab,alphap,alphapp)

 call prin2("after kummer_adap, ab = ",ab)
 call prin2("after kummer_adap, alphap = ",alphap)

ifleft = 1

call kummer_phase(ifleft,k,xscheb,chebintl,chebintr,ucheb, &
   nints,ab,alpha,alphap,alphapp)

call prin2("after kummer_phase, alpha = ",alpha)

end subroutine



function phase_intervals() result(nintsout)
implicit double precision (a-h,o-z)
integer :: nintsout
nintsout = nints
end function


subroutine phase_data(nintsout,kout,about,alphaout,alphapout,alphappout)
implicit double precision (a-h,o-z)

integer          :: nintsout,kout
double precision :: about(2,nintsout),alphaout(kout,nintsout),alphapout(kout,nintsout)
double precision :: alphappout(kout,nintsout)



do int=1,nintsout
about(1,int) = ab(1,int)
about(2,int) = ab(2,int)

do i=1,kout
alphappout(i,int) = alphapp(i,int)
alphapout(i,int) = alphap(i,int)
alphaout(i,int)  = alpha(i,int)

end do

end do



end subroutine

end module wrapper
