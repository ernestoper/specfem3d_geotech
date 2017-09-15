! collection of solvers
! AUTHOR
!   Hom Nath Gharti
! REVISION
!   HNG, Jul 12,2011; HNG, Apr 09,2010
module solver
use set_precision
use global, only : cg_maxiter,cg_tol,g_num,nedof
use math_constants, only : zero

contains

! diagonally preconditioned conjuate-gradient solver
subroutine pcg_solver(neq,nelmt,k,u,f,dprecon,gdof_elmt,cg_iter,errcode,errtag)
use,intrinsic :: iso_c_binding
implicit none

interface
  subroutine ebe_matvec_c(nelmt,nedof,neq,gdof_elmt,k,vin,vout) bind(c,NAME='ebe_matvec_c')
  use,intrinsic :: iso_c_binding
  integer(c_int),value :: nelmt,nedof,neq
  integer(c_int),dimension(:,:) :: gdof_elmt
  real(c_double),dimension(:,:,:) :: k
  real(c_double),dimension(:) :: vin,vout
  end subroutine ebe_matvec_c
end interface

integer,intent(in) :: neq,nelmt
! nelmt (for intact) may not be same as global nelmt
real(kind=kreal),dimension(nedof,nedof,nelmt),intent(in) :: k
! only for intact elements
real(kind=kreal),dimension(0:neq),intent(inout) :: u
real(kind=kreal),dimension(0:neq),intent(in) :: f,dprecon
!integer,dimension(nedof,nelmt),intent(in) :: gdof_elmt
integer,dimension(:,:),intent(in) :: gdof_elmt
! only for intact elements
integer,intent(out) :: cg_iter
integer,intent(out) :: errcode
character(len=250),intent(out) :: errtag

integer :: i_elmt
integer,dimension(nedof) :: egdof
real(kind=kreal) :: alpha,beta,rz
real(kind=kreal),dimension(0:neq) :: kp,p,r,z
real(kind=kreal),dimension(nedof,nedof) :: km

errtag="ERROR: unknown!"
errcode=-1

kp=zero
if(maxval(abs(u)).gt.zero)then
  !do i_elmt=1,nelmt
  !  egdof=gdof_elmt(:,i_elmt)
  !  km=k(:,:,i_elmt)
  !  kp(egdof)=kp(egdof)+matmul(km,u(egdof))
  !enddo
  !kp(0)=zero
  call ebe_matvec_c(nelmt,nedof,neq,gdof_elmt,k,u,kp)
endif
r=f-kp
z=dprecon*r

p=z
! pcg iteration
pcg: do cg_iter=1,cg_maxiter
  !kp=zero
  !do i_elmt=1,nelmt
  !  egdof=gdof_elmt(:,i_elmt)
  !  km=k(:,:,i_elmt)
  !  kp(egdof)=kp(egdof)+matmul(km,p(egdof))
  !enddo
  !kp(0)=zero
  call ebe_matvec_c(nelmt,nedof,neq,gdof_elmt,k,p,kp)

  rz=dot_product(r,z)
  alpha=rz/dot_product(p,kp)
  u=u+alpha*p

  if(abs(alpha)*maxval(abs(p))/maxval(abs(u)).le.cg_tol)then
    errcode=0
    return
  endif

  r=r-alpha*kp
  z=dprecon*r
  beta=dot_product(r,z)/rz
  p=z+beta*p
  !write(*,'(i3,f25.18,f25.18,f25.18)')cg_iter,alpha,beta,rz

enddo pcg
write(errtag,'(a)')'ERROR: PCG solver doesn''t converge!'
return
end subroutine pcg_solver
!===============================================================================

end module solver
!===============================================================================
