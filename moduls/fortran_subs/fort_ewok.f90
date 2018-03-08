!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      potential solver routine            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine explicit_scheme_loop(potent, binsy,binsx,k_y_loc,rhs, dt,c1)

implicit none

integer, intent(in) :: binsx, binsy
real (kind=8), intent(in) :: dt, c1, k_y_loc(0:binsx,0:binsy)
complex (kind=8), dimension(0:binsx,0:binsy), intent(in) :: rhs
complex (kind=8), dimension(0:binsx,0:binsy), intent(inout) :: potent
real (kind=8) :: pi=3.1415926535897931
complex (kind=8) :: ii=(0.0,1.0), c2

!f2py intent(in,out) :: potent
!f2py intent(in) :: k_y_loc, rhs, dt
!f2py intent(hide) :: binsx, binsy

c2 = (pi*dt/ii)

potent(1:binsx,:) = (potent(0:binsx-1,:)*c1 + potent(1:binsx,:)*(1.-c1)+ &
                     rhs(0:binsx-1,:)*c2)/( 1.0 - c2*k_y_loc(1:binsx,:)**2)

end subroutine



subroutine push_them_relativ(potent_estat, potent_em, particles, left, bottom, &
      x_ij, y_ij, dx, dy, dt, binsx, binsy,length)

!!!! calculates the interpolated X,Y forces at particle; !!!
!!!! calculates the relativistic factor of each particle !!!
!!!! with account for oscillations in e.m. field;        !!!
!!!! evaluates the particles motion and returns modified !!!
!!!! coordinates, velocities and gamma factor            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
integer, intent(in) :: binsx, binsy, length
real (kind=8), dimension(0:binsx,0:binsy), intent(in) :: potent_estat, potent_em, x_ij, y_ij
real (kind=8) , intent(in) :: left, bottom, dx, dy, dt
real (kind=8), dimension(length,5), intent(inout) :: particles

integer :: ii, jj,i
real (kind=8) :: f00, f10, f01, f11,f0m1,fm10,f1m1,fm11,f02,f20,fm1m1,f12,f21,fm12,f2m1,f22, potent_proj,x,y
real (kind=8) :: a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,a11,a12,a13,a14,a15,a16
real (kind=8) :: ax01,ax02,ax03,ax04,ax05,ax06,ax07,ax08,ax09,ax10,ax11,ax12,ax13,ax14,ax15,ax16
real (kind=8) :: ay01,ay02,ay03,ay04,ay05,ay06,ay07,ay08,ay09,ay10,ay11,ay12,ay13,ay14,ay15,ay16
real (kind=8), dimension(2) :: force_em, force_estat

!f2py intent(in,out) :: particles
!f2py intent(in) :: potent_estat, potent_em, x_ij, y_ij, dx, dy, dt,left, bottom
!f2py intent(hide) :: binsx, binsy, length
do i=1,length
   ii = INT(FLOOR( (particles(i,1)-left)/dx   ))
   jj = INT(FLOOR( (particles(i,2)-bottom)/dy ))

   x = (particles(i,1) - x_ij(ii,jj) )/dx
   y = (particles(i,2) - y_ij(ii,jj))/dy

   if (ii<0 .OR. jj<0 .OR. ii> binsx .OR. jj> binsy) then 
!     write(*,*) 'oups..',ii,jj, 'bins', binsx, binsy
     write(*,*) 'oups..',particles(i,1), particles(i,2),left,bottom
     CYCLE
   endif

   if(ii<=1 .OR. jj<=1 .OR. ii>= binsx-1 .OR. jj>= binsy-1) then

   f00 = potent_em(ii,jj)
   f10 = potent_em(ii+1,jj)
   f01 = potent_em(ii,jj+1)
   f11 = potent_em(ii+1,jj+1)

   potent_proj = (1.-x)*(1.-y)*f00 + x*(1.-y)*f10 + &
                          (1.-x)*y*f01 + x*y*f11

   force_em(1) = -(1.-y)*f00 + (1.-y)*f10 - y*f01 + y*f11
   force_em(2) = -(1.-x)*f00 - x*f10 + (1.-x)*f01 + x*f11

   f00 = potent_estat(ii,jj)
   f10 = potent_estat(ii+1,jj)
   f01 = potent_estat(ii,jj+1)
   f11 = potent_estat(ii+1,jj+1)

   force_estat(1) = -(1.-y)*f00 + (1.-y)*f10 - y*f01 + y*f11
   force_estat(2) = -(1.-x)*f00 - x*f10 + (1.-x)*f01 + x*f11

   else

   a01 = 0.25*(x-1.)*(x-2.)*(x+1.)*(y-1.)*(y-2.)*(y+1.)
   a02 = -0.25*x*(x+1.)*(x-2.)*(y-1.)*(y-2.)*(y+1.)
   a03 = -0.25*(x-1.)*(x-2.)*(x+1.)*y*(y+1.)*(y-2.)
   a04 = 0.25*x*(x+1.)*(x-2.)*y*(y+1.)*(y-2.)
   a05 = -(1./12.)*x*(x-1.)*(x-2.)*(y-1.)*(y-2.)*(y+1.)
   a06 = -(1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y-1.)*(y-2.)
   a07 = (1./12.)*x*(x-1.)*(x-2.)*y*(y+1.)*(y-2.)
   a08 = (1./12.)*x*(x+1.)*(x-2.)*y*(y-1.)*(y-2.)
   a09 = (1./12.)*x*(x-1.)*(x+1.)*(y-1.)*(y-2.)*(y+1.)
   a10 = (1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y-1.)*(y+1.)
   a11 = (1./36.)*x*(x-1.)*(x-2.)*y*(y-1.)*(y-2.)
   a12 = -(1./12.)*x*(x-1.)*(x+1.)*y*(y+1.)*(y-2.)
   a13 = -(1./12.)*x*(x+1.)*(x-2.)*y*(y-1.)*(y+1.)
   a14 = -(1./36.)*x*(x-1.)*(x+1.)*y*(y-1.)*(y-2.)
   a15 = -(1./36.)*x*(x-1.)*(x-2.)*y*(y-1.)*(y+1.)
   a16 = (1./36.)*x*(x-1.)*(x+1.)*y*(y-1.)*(y+1.)

   ax01 = 0.25*(x-2.)*(x+1.)*(y-1.)*(y-2.)*(y+1.)+0.25*(x-1.)*(x+1.)*(y-1.)*(y-2.)*(y+1.)+0.25*(x-1.)*(x-2.)*(y-1.)*(y-2.)*(y+1.)
   ax02 = -0.25*(x+1.)*(x-2.)*(y-1.)*(y-2.)*(y+1.)-0.25*x*(x-2.)*(y-1.)*(y-2.)*(y+1.)-0.25*x*(x+1.)*(y-1.)*(y-2.)*(y+1.)
   ax03 = -0.25*(x-2.)*(x+1.)*y*(y+1.)*(y-2.)-0.25*(x-1.)*(x+1.)*y*(y+1.)*(y-2.)-0.25*(x-1.)*(x-2.)*y*(y+1.)*(y-2.)
   ax04 = 0.25*(x+1.)*(x-2.)*y*(y+1.)*(y-2.)+0.25*x*(x-2.)*y*(y+1.)*(y-2.)+0.25*x*(x+1.)*y*(y+1.)*(y-2.)
   ax05 = -(1./12.)*(x-1.)*(x-2.)*(y-1.)*(y-2.)*(y+1.)-(1./12.)*x*(x-2.)*(y-1.)*(y-2.)*(y+1.)-(1./12.)*x*(x-1.)*(y-1.)*(y-2.)*(y+1.)
   ax06 = -(1./12.)*(x-2.)*(x+1.)*y*(y-1.)*(y-2.)-(1./12.)*(x-1.)*(x+1.)*y*(y-1.)*(y-2.)-(1./12.)*(x-1.)*(x-2.)*y*(y-1.)*(y-2.)
   ax07 = (1./12.)*(x-1.)*(x-2.)*y*(y+1.)*(y-2.)+(1./12.)*x*(x-2.)*y*(y+1.)*(y-2.)+(1./12.)*x*(x-1.)*y*(y+1.)*(y-2.)
   ax08 = (1./12.)*(x+1.)*(x-2.)*y*(y-1.)*(y-2.)+(1./12.)*x*(x-2.)*y*(y-1.)*(y-2.)+(1./12.)*x*(x+1.)*y*(y-1.)*(y-2.)
   ax09 = (1./12.)*(x-1.)*(x+1.)*(y-1.)*(y-2.)*(y+1.)+(1./12.)*x*(x+1.)*(y-1.)*(y-2.)*(y+1.)+(1./12.)*x*(x-1.)*(y-1.)*(y-2.)*(y+1.)
   ax10 = (1./12.)*(x-2.)*(x+1.)*y*(y-1.)*(y+1.)+(1./12.)*(x-1.)*(x+1.)*y*(y-1.)*(y+1.)+(1./12.)*(x-1.)*(x-2.)*y*(y-1.)*(y+1.)
   ax11 = (1./36.)*(x-1.)*(x-2.)*y*(y-1.)*(y-2.)+(1./36.)*x*(x-2.)*y*(y-1.)*(y-2.)+(1./36.)*x*(x-1.)*y*(y-1.)*(y-2.)
   ax12 = -(1./12.)*(x-1.)*(x+1.)*y*(y+1.)*(y-2.)-(1./12.)*x*(x+1.)*y*(y+1.)*(y-2.)-(1./12.)*x*(x-1.)*y*(y+1.)*(y-2.)
   ax13 = -(1./12.)*(x+1.)*(x-2.)*y*(y-1.)*(y+1.)-(1./12.)*x*(x-2.)*y*(y-1.)*(y+1.)-(1./12.)*x*(x+1.)*y*(y-1.)*(y+1.)
   ax14 = -(1./36.)*(x-1.)*(x+1.)*y*(y-1.)*(y-2.)-(1./36.)*x*(x+1.)*y*(y-1.)*(y-2.)-(1./36.)*x*(x-1.)*y*(y-1.)*(y-2.)
   ax15 = -(1./36.)*(x-1.)*(x-2.)*y*(y-1.)*(y+1.)-(1./36.)*x*(x-2.)*y*(y-1.)*(y+1.)-(1./36.)*x*(x-1.)*y*(y-1.)*(y+1.)
   ax16 = (1./36.)*(x-1.)*(x+1.)*y*(y-1.)*(y+1.)+(1./36.)*x*(x+1.)*y*(y-1.)*(y+1.)+(1./36.)*x*(x-1.)*y*(y-1.)*(y+1.)

   ay01 = 0.25*(x-1.)*(x-2.)*(x+1.)*(y-2.)*(y+1.)+0.25*(x-1.)*(x-2.)*(x+1.)*(y-1.)*(y+1.)+0.25*(x-1.)*(x-2.)*(x+1.)*(y-1.)*(y-2.)
   ay02 = -0.25*x*(x+1.)*(x-2.)*(y-2.)*(y+1.)-0.25*x*(x+1.)*(x-2.)*(y-1.)*(y+1.)-0.25*x*(x+1.)*(x-2.)*(y-1.)*(y-2.)
   ay03 = -0.25*(x-1.)*(x-2.)*(x+1.)*(y+1.)*(y-2.)-0.25*(x-1.)*(x-2.)*(x+1.)*y*(y-2.)-0.25*(x-1.)*(x-2.)*(x+1.)*y*(y+1.)
   ay04 = 0.25*x*(x+1.)*(x-2.)*(y+1.)*(y-2.)+0.25*x*(x+1.)*(x-2.)*y*(y-2.)+0.25*x*(x+1.)*(x-2.)*y*(y+1.)
   ay05 = -(1./12.)*x*(x-1.)*(x-2.)*(y-2.)*(y+1.)-(1./12.)*x*(x-1.)*(x-2.)*(y-1.)*(y+1.)-(1./12.)*x*(x-1.)*(x-2.)*(y-1.)*(y-2.)
   ay06 = -(1./12.)*(x-1.)*(x-2.)*(x+1.)*(y-1.)*(y-2.)-(1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y-2.)-(1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y-1.)
   ay07 = (1./12.)*x*(x-1.)*(x-2.)*(y+1.)*(y-2.)+(1./12.)*x*(x-1.)*(x-2.)*y*(y-2.)+(1./12.)*x*(x-1.)*(x-2.)*y*(y+1.)
   ay08 = (1./12.)*x*(x+1.)*(x-2.)*(y-1.)*(y-2.)+(1./12.)*x*(x+1.)*(x-2.)*y*(y-2.)+(1./12.)*x*(x+1.)*(x-2.)*y*(y-1.)
   ay09 = (1./12.)*x*(x-1.)*(x+1.)*(y-2.)*(y+1.)+(1./12.)*x*(x-1.)*(x+1.)*(y-1.)*(y+1.)+(1./12.)*x*(x-1.)*(x+1.)*(y-1.)*(y-2.)
   ay10 = (1./12.)*(x-1.)*(x-2.)*(x+1.)*(y-1.)*(y+1.)+(1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y+1.)+(1./12.)*(x-1.)*(x-2.)*(x+1.)*y*(y-1.)
   ay11 = (1./36.)*x*(x-1.)*(x-2.)*(y-1.)*(y-2.)+(1./36.)*x*(x-1.)*(x-2.)*y*(y-2.)+(1./36.)*x*(x-1.)*(x-2.)*y*(y-1.)
   ay12 = -(1./12.)*x*(x-1.)*(x+1.)*(y+1.)*(y-2.)-(1./12.)*x*(x-1.)*(x+1.)*y*(y-2.)-(1./12.)*x*(x-1.)*(x+1.)*y*(y+1.)
   ay13 = -(1./12.)*x*(x+1.)*(x-2.)*(y-1.)*(y+1.)-(1./12.)*x*(x+1.)*(x-2.)*y*(y+1.)-(1./12.)*x*(x+1.)*(x-2.)*y*(y-1.)
   ay14 = -(1./36.)*x*(x-1.)*(x+1.)*(y-1.)*(y-2.)-(1./36.)*x*(x-1.)*(x+1.)*y*(y-2.)-(1./36.)*x*(x-1.)*(x+1.)*y*(y-1.)
   ay15 = -(1./36.)*x*(x-1.)*(x-2.)*(y-1.)*(y+1.) -(1./36.)*x*(x-1.)*(x-2.)*y*(y+1.) -(1./36.)*x*(x-1.)*(x-2.)*y*(y-1.)
   ay16 = (1./36.)*x*(x-1.)*(x+1.)*(y-1.)*(y+1.)+(1./36.)*x*(x-1.)*(x+1.)*y*(y+1.)+(1./36.)*x*(x-1.)*(x+1.)*y*(y-1.)

   f00 = potent_em(ii,jj)
   f01 = potent_em(ii+1,jj)
   f10 = potent_em(ii,jj+1)
   f11 = potent_em(ii+1,jj+1)
   f0m1 = potent_em(ii-1,jj)
   fm10 = potent_em(ii,jj-1)
   f1m1 = potent_em(ii-1,jj+1)
   fm11 = potent_em(ii+1,jj-1)
   f02 = potent_em(ii+2,jj)
   f20 = potent_em(ii,jj+2)
   fm1m1 = potent_em(ii-1,jj-1)
   f12 = potent_em(ii+2,jj+1)
   f21 = potent_em(ii+1,jj+2)
   fm12 = potent_em(ii+2,jj-1)
   f2m1 = potent_em(ii-1,jj+2)
   f22 = potent_em(ii+2,jj+2)

   potent_proj = a01*f00 + a02*f01 + a03*f10 + a04*f11 + a05*f0m1 + a06*fm10 + a07*f1m1 + a08*fm11 + &
                 a09*f02 + a10*f20 + a11*fm1m1 + a12*f12 + a13*f21 + a14*fm12 + a15*f2m1 + a16*f22
   force_em(1) = ax01*f00 + ax02*f01 + ax03*f10 + ax04*f11 + ax05*f0m1 + ax06*fm10 + ax07*f1m1 + ax08*fm11 + &
                 ax09*f02 + ax10*f20 + ax11*fm1m1 + ax12*f12 + ax13*f21 + ax14*fm12 + ax15*f2m1 + ax16*f22
   force_em(2) = ay01*f00 + ay02*f01 + ay03*f10 + ay04*f11 + ay05*f0m1 + ay06*fm10 + ay07*f1m1 + ay08*fm11 + &
                 ay09*f02 + ay10*f20 + ay11*fm1m1 + ay12*f12 + ay13*f21 + ay14*fm12 + ay15*f2m1 + ay16*f22

   f00 = potent_estat(ii,jj)
   f01 = potent_estat(ii+1,jj)
   f10 = potent_estat(ii,jj+1)
   f11 = potent_estat(ii+1,jj+1)
   f0m1 = potent_estat(ii-1,jj)
   fm10 = potent_estat(ii,jj-1)
   f1m1 = potent_estat(ii-1,jj+1)
   fm11 = potent_estat(ii+1,jj-1)
   f02 = potent_estat(ii+2,jj)
   f20 = potent_estat(ii,jj+2)
   fm1m1 = potent_estat(ii-1,jj-1)
   f12 = potent_estat(ii+2,jj+1)
   f21 = potent_estat(ii+1,jj+2)
   fm12 = potent_estat(ii+2,jj-1)
   f2m1 = potent_estat(ii-1,jj+2)
   f22 = potent_estat(ii+2,jj+2)

   force_estat(1) = ax01*f00 + ax02*f01 + ax03*f10 + ax04*f11 + ax05*f0m1 + ax06*fm10 + ax07*f1m1 + ax08*fm11 + &
                 ax09*f02 + ax10*f20 + ax11*fm1m1 + ax12*f12 + ax13*f21 + ax14*fm12 + ax15*f2m1 + ax16*f22
   force_estat(2) = ay01*f00 + ay02*f01 + ay03*f10 + ay04*f11 + ay05*f0m1 + ay06*fm10 + ay07*f1m1 + ay08*fm11 + &
                 ay09*f02 + ay10*f20 + ay11*fm1m1 + ay12*f12 + ay13*f21 + ay14*fm12 + ay15*f2m1 + ay16*f22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   endif
   particles(i,5) = dsqrt(1.0 + potent_proj/2. + particles(i,3)**2 + particles(i,4)**2)
   if (particles(i,5)<1.) then
      particles(i,5) = 1.
   end if

   force_em(1) = (-0.25/dx/particles(i,5))*force_em(1)
   force_em(2) = (-0.25/dx/particles(i,5))*force_em(2)

   force_estat(1) = (-1.0/dx)*force_estat(1)
   force_estat(2) = (-1.0/dx)*force_estat(2)

   particles(i,3:4) = particles(i,3:4) + dt*(force_em(:) + force_estat(:))
   particles(i,1:2) = particles(i,1:2) +dt*particles(i,3:4)/particles(i,5)

end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! finds and sorts migrating and staying particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bound_and_migration(particles, left_p, right_p, bottom_p, & 
    top_p, leng ,  part_2left, part_2right, part_2stay, Npart_2left, & 
    Npart_2right, Npart_2stay)

implicit none

real (kind=8), dimension(leng,7), intent(inout) :: particles
real (kind=8) , intent(in) :: left_p,right_p,bottom_p, top_p
integer, intent(in) :: leng

real (kind=8), dimension(leng,7), intent(out) :: part_2left, part_2right, part_2stay
integer, intent(out) :: Npart_2left, Npart_2right, Npart_2stay

!f2py intent(in) :: particles, left_p,right_p,bottom_p, top_p
!f2py intent(out) :: part_2left, part_2right, part_2stay, Npart_2left, Npart_2right, Npart_2stay
!f2py intent(hide) :: leng

integer :: i

Npart_2left = 0
Npart_2right = 0 
Npart_2stay = 0
part_2left=0
part_2right=0
part_2stay = 0

!write(*,*) top_p,  bottom_p, right_p, left_p

do i=1,leng
   if(particles(i,2) >= top_p) then
      particles(i,2) = particles(i,2) - (top_p-bottom_p)
   else if(particles(i,2) < bottom_p) then 
      particles(i,2) = particles(i,2) + (top_p-bottom_p)
   endif

   if(particles(i,1) > right_p) then
      Npart_2right = Npart_2right+1
      part_2right(Npart_2right,:) = particles(i,:)
   else if(particles(i,1) <= left_p) then
      Npart_2left = Npart_2left+1
      part_2left(Npart_2left,:) = particles(i,:)
   else
      Npart_2stay = Npart_2stay+1
      part_2stay(Npart_2stay,:) = particles(i,:)
   endif

end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! finds and sorts migrating and staying particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bound_migration_absorbtion(particles, left_p, right_p, bottom_p, &
        top_p, leng ,  part_2left, part_2right, part_2stay, Npart_2left, &
        Npart_2right, Npart_2stay,Npart_absorbed)

implicit none

real (kind=8), dimension(leng,7), intent(inout) :: particles
real (kind=8) , intent(in) :: left_p,right_p,bottom_p, top_p
integer, intent(in) :: leng

real (kind=8), dimension(leng,7), intent(out) :: part_2left, part_2right, part_2stay
integer, intent(out) :: Npart_2left, Npart_2right, Npart_2stay, Npart_absorbed

!f2py intent(in) :: particles, left_p,right_p,bottom_p, top_p
!f2py intent(out) :: part_2left, part_2right, part_2stay, Npart_2left, Npart_2right, Npart_2stay, Npart_absorbed
!f2py intent(hide) :: leng

integer :: i

Npart_2left = 0
Npart_2right = 0
Npart_2stay = 0
Npart_absorbed = 0

part_2left=0
part_2right=0
part_2stay = 0

do i=1,leng
   if(particles(i,1) > right_p .AND. particles(i,2) <= top_p .AND. particles(i,2) >= bottom_p) then
      Npart_2right = Npart_2right+1
      part_2right(Npart_2right,:) = particles(i,:)
   else if(particles(i,1) <= left_p .AND. particles(i,2) <= top_p .AND. particles(i,2) >= bottom_p) then
      Npart_2left = Npart_2left+1
      part_2left(Npart_2left,:) = particles(i,:)
   else if( particles(i,2) <= top_p .AND. particles(i,2) >= bottom_p) then
      Npart_2stay = Npart_2stay+1
      part_2stay(Npart_2stay,:) = particles(i,:)
   else
      Npart_absorbed = Npart_absorbed+1
   endif

end do
end subroutine


!======================================================================================
SUBROUTINE DENSITY_2x(coord,wght,origx, origy, dlt_xg, dlt_yg,bins_x,bins_y,dens, N_part)
!-------------------------------------------------------------------------------------
!     Calculates the 2D density ditribution (histograming subroutine)
!--------------------------------------------------------------------------------------

IMPLICIT NONE
REAL(kind=8), INTENT(IN), dimension(N_part,2) :: coord
REAL(kind=8), INTENT(IN), dimension(N_part) :: wght
REAL(kind=8), INTENT(IN)   :: origx, origy, dlt_xg, dlt_yg
INTEGER, INTENT(IN)        :: bins_x, bins_y, N_part
REAL(kind=8), INTENT(OUT)  :: dens(-2:bins_x+2,-2:bins_y+2)
REAL(kind=8)               :: dx, dy, dlt_xg_inv, dlt_yg_inv, xg(-2:bins_x+2), yg(-2:bins_y+2)
INTEGER                    :: jp, kx, ky, j
REAL(kind=8)               :: S0(-3:3,1:2),x_max,y_max

!f2py intent(in) :: coord,wght, origx, origy, dlt_xg, dlt_yg, bins_x, bins_y
!f2py intent(out) :: dens
!f2py intent(hide) :: N_part

dens = 0.d0
dlt_xg_inv = 1.0d0/dlt_xg

dlt_yg_inv= 1.0d0/dlt_yg


x_max=origx+dlt_xg*bins_x
y_max=origy+dlt_yg*bins_y

DO j  =  -2,bins_x+2
   xg(j)  =  dlt_xg*REAL(j) + origx
ENDDO
DO j  =  -2,bins_y+2
   yg(j)  =  dlt_yg*REAL(j) + origy
ENDDO


DO jp = 1, N_part
 if((coord(jp,1)>=origx).AND.(coord(jp,1)<=x_max).AND. &
      (coord(jp,2)>=origy).AND.(coord(jp,2)<=y_max)) then
   S0(:,:) = 0.d0
   kx     = FLOOR((coord(jp,1)-origx  )*dlt_xg_inv+0.5d0)
   ky     = FLOOR((coord(jp,2)-origy  )*dlt_yg_inv+0.5d0)
   dx     =       (coord(jp,1)-xg(kx ))*dlt_xg_inv
   dy     =       (coord(jp,2)-yg(ky ))*dlt_yg_inv
   S0(-1,1) = (0.25d0-0.5d0*dx+(dx**3)/3.d0)
   S0( 0,1) = (0.5d0-dabs(dx**3)/3.d0)      
   S0( 1,1) = (0.25d0+0.5d0*dx-(dx**3)/3.d0)
   if(dx >= 0.0) then
       S0(2,1)  = dabs(dx**3)/3.d0
   else
       S0(-2,1) = dabs(dx**3)/3.d0
   endif
   S0(-1,2) = (0.25d0-0.5d0*dy+(dy**3)/3.d0)
   S0( 0,2) = (0.5d0-dabs(dy**3)/3.d0)      
   S0( 1,2) = (0.25d0+0.5d0*dy-(dy**3)/3.d0)
   if(dy >= 0.0) then
      S0(2,2)  = dabs(dy**3)/3.d0
   else
      S0(-2,2) = dabs(dy**3)/3.d0
   endif
   DO j = -2, 2
      dens(kx-2:kx+2,ky+j) = dens(kx-2:kx+2,ky+j) + wght(jp)*S0(-2:2,1)*S0(j,2)
   ENDDO
 end if
ENDDO

dens = dens * dlt_xg_inv * dlt_yg_inv

END SUBROUTINE
