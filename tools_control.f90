!*******************************************************************
!
subroutine test_scalar_min_max(phi)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI


implicit none

integer :: i,j,k
real(mytype) :: phimax,phimin,cfl
real(mytype) :: phimax1,phimin1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi
integer :: code

phimax=-1609.
phimin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   if (phi(i,j,k).gt.phimax) phimax=phi(i,j,k)
   if (phi(i,j,k).lt.phimin) phimin=phi(i,j,k)
enddo
enddo
enddo


call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
if (nrank==0) then
   print *,'PHI max=',phimax1
   print *,'PHI min=',phimin1
endif


return
end subroutine test_scalar_min_max

!*******************************************************************
!
subroutine test_speed_min_max(ux,uy,uz)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI


implicit none

integer :: i,j,k
real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin,cfl
real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer :: code

uxmax=-1609.
uymax=-1609.
uzmax=-1609.
uxmin=1609.
uymin=1609.
uzmin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   if (ux(i,j,k).gt.uxmax) uxmax=ux(i,j,k)
   if (uy(i,j,k).gt.uymax) uymax=uy(i,j,k)
   if (uz(i,j,k).gt.uzmax) uzmax=uz(i,j,k)
   if (ux(i,j,k).lt.uxmin) uxmin=ux(i,j,k)
   if (uy(i,j,k).lt.uymin) uymin=uy(i,j,k)
   if (uz(i,j,k).lt.uzmin) uzmin=uz(i,j,k)
enddo
enddo
enddo


call MPI_REDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymax,uymax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymin,uymin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
if (nrank==0) then
   print *,'U,V,W max=',uxmax1,uymax1,uzmax1
   print *,'U,V,W min=',uxmin1,uymin1,uzmin1
endif


return
end subroutine test_speed_min_max

!*******************************************************************
!
subroutine restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,px1,py1,pz1,phis1,&
                   hx1,hy1,hz1,phiss1,phG,irestart)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE MPI

implicit none

TYPE(DECOMP_INFO) :: phG
integer :: i,j,k,irestart,nzmsize,fh,ierror,code
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: phi1, phis1, phiss1
real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
integer (kind=MPI_OFFSET_KIND) :: filesize, disp
real(mytype) :: xdt
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

if (iscalar==0) then
   if (nscheme.ne.4) then
      if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)   
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   else !AB3
     if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,hx1)
         call decomp_2d_write_var(fh,disp,1,hy1)
         call decomp_2d_write_var(fh,disp,1,hz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)   
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,hx1)
         call decomp_2d_read_var(fh,disp,1,hy1)
         call decomp_2d_read_var(fh,disp,1,hz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      endif 
   endif
else !SCALAR 
if (nscheme.ne.4) then
      if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call decomp_2d_write_var(fh,disp,1,phi1)
         call decomp_2d_write_var(fh,disp,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)   
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call decomp_2d_read_var(fh,disp,1,phi1)
         call decomp_2d_read_var(fh,disp,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   else !SCALAR + AB3
     if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,hx1)
         call decomp_2d_write_var(fh,disp,1,hy1)
         call decomp_2d_write_var(fh,disp,1,hz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call decomp_2d_write_var(fh,disp,1,phi1)
         call decomp_2d_write_var(fh,disp,1,phis1)
         call decomp_2d_write_var(fh,disp,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)   
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,hx1)
         call decomp_2d_read_var(fh,disp,1,hy1)
         call decomp_2d_read_var(fh,disp,1,hz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call decomp_2d_read_var(fh,disp,1,phi1)
         call decomp_2d_read_var(fh,disp,1,phis1)
         call decomp_2d_read_var(fh,disp,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)
      endif 
   endif
endif

if (irestart==0) then
! reconstruction of the dp/dx, dp/dy and dp/dz from px1,py1 and pz1
! Temporal scheme (1:AB2, 2: RK3, 3:RK4, 4:AB3)
   if (nscheme==1) xdt=adt(1)+bdt(1)
   if (nscheme==2) xdt=(3./4.)*dt + (-5./12.)*dt
   if (nscheme==3) xdt=0.041717869325*dt
   if (nscheme==4) xdt=adt(1)+bdt(1)+cdt(1)
   
   do k=1,xsize(3)
   do j=1,xsize(2)
      dpdyx1(j,k)=py1(1,j,k)/xdt
      dpdzx1(j,k)=pz1(1,j,k)/xdt
      dpdyxn(j,k)=py1(nx,j,k)/xdt
      dpdzxn(j,k)=pz1(nx,j,k)/xdt
   enddo
   enddo
   
   if (xsize(3)==1) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxz1(i,j)=px1(i,j,1)/xdt
         dpdyz1(i,j)=py1(i,j,1)/xdt
      enddo
      enddo
   endif
   if (xsize(3)==nz) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxzn(i,j)=px1(i,j,nz)/xdt
         dpdyzn(i,j)=py1(i,j,nz)/xdt
      enddo
      enddo
   endif

   ! determine the processor grid in use
   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, code)

   if (dims(1)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxy1(i,k)=px1(i,1,k)/xdt
      dpdzy1(i,k)=pz1(i,1,k)/xdt
      enddo
      enddo
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
      dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
      enddo
      enddo
   else
!find j=1 and j=ny
      if (xstart(2)==1) then
         do k=1,xsize(3)
         do i=1,xsize(1)
         dpdxy1(i,k)=px1(i,1,k)/xdt
         dpdzy1(i,k)=pz1(i,1,k)/xdt
         enddo
         enddo
      endif
!      print *,nrank,xstart(2),ny-(nym/p_row)
       if (ny-(nym/dims(1))==xstart(2)) then
         do k=1,xsize(3)
         do i=1,xsize(1)
         dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
         dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
         enddo
         enddo
      endif

   endif


   if (nrank==0) print *,'reconstruction pressure gradients done!'
endif

if (irestart==0) then
if (ivirt==2) then
   call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
        MPI_MODE_RDONLY, MPI_INFO_NULL, &
        fh, ierror)
   disp = 0_MPI_OFFSET_KIND
   call decomp_2d_read_var(fh,disp,1,ep1) 
   call MPI_FILE_CLOSE(fh,ierror)
   if (nrank==0) print *,'read epsilon file done from restart'
endif
endif

end subroutine restart
!
!*******************************************************************
subroutine stretching()
!
!*******************************************************************
!
USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
integer :: j

yinf=-yly/2.
den=2.*beta*yinf
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
alpha=abs(xnum/den)
xcx=1./beta/alpha
if (alpha.ne.0.) then
   if (istret.eq.1) yp(1)=0.
   if (istret.eq.2) yp(1)=0.
   if (istret.eq.1) yeta(1)=0.
   if (istret.eq.2) yeta(1)=-1./2.
   if (istret.eq.3) yp(1)=0.
   if (istret.eq.3) yeta(1)=-1./2.
   do j=2,ny!-1
      if (istret==1) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      if ((ncly.ne.0).and.(j==ny).and.(istret==1)) then
         xnum1=0.
      else
         xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      endif
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
         if (yeta(j).eq.0.5) yp(j)=0.-yinf
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
         if (yeta(j).eq.0.5) yp(j)=0.+yly
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
         if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
         if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
      endif
   enddo

endif
!if (nrank==0) then
!do j=1,ny
!print *,j,yp(j),yeta(j)
!enddo
!endif
!stop
if (alpha.eq.0.) then
   yp(1)=-1.e10
   do j=2,ny
      yeta(j)=(j-1.)*(1./ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
   enddo
endif
if (alpha.ne.0.) then
   do j=1,ny
      if (istret==1) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./2./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./2./(ny-1.))-0.5
      endif
      
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
         if (yetai(j).eq.0.5) ypi(j)=0.-yinf
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst+yly
         if (yetai(j).eq.0.5) ypi(j)=0.+yly
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yetai(j).lt.0.5) ypi(j)=(xnum1-cst+yly)*2.
         if (yetai(j).eq.0.5) ypi(j)=(0.+yly)*2.
         if (yetai(j).gt.0.5) ypi(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   ypi(1)=-1.e10
   do j=2,ny
      yetai(j)=(j-1.)*(1./ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
   enddo
endif

do j=1,ny
   ppy(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
        sin(pi*yeta(j)))
   pp2y(j)=ppy(j)*ppy(j)
   pp4y(j)=(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
enddo
do j=1,ny
   ppyi(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yetai(j))* &
        sin(pi*yetai(j)))
   pp2yi(j)=ppyi(j)*ppyi(j)
   pp4yi(j)=(-2./beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
enddo

! if(nrank.eq.0) then
!    open(69,file='ygrid.txt',status='unknown')
!    do j=1,ny
!       write(69,*) yp(j)
!    end do
! end if
   
end subroutine stretching

!*****************************************************************
!
subroutine inversion5_v1(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
   do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
         if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
         if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
         sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
         eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
              aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
      enddo
      enddo
      do jc=ja(i),jb(i)
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
         aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
              aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
      enddo
      enddo
      enddo
   enddo
enddo


do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
        aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
   else
      tmp1=0.
      tmp3=0.
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
   else
      tmp2=0.
      tmp4=0.
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)
   
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
        aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
enddo
enddo

do i=ny/2-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,i,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v1 

!*****************************************************************
!
subroutine inversion5_v2(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,nym-2
do i=1,2
   mi=m+i
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1) 
      if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
      if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
      sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
           aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
   enddo
   enddo
   do jc=ja(i),jb(i)
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1) 
      aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
           aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
   enddo
   enddo
   enddo
enddo
enddo
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,nym,k,2), kind=mytype)/real(aaa(j,nym-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,nym,k,2))/aimag(aaa(j,nym-1,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,nym,k,3), kind=mytype)-tmp1*real(aaa(j,nym-1,k,4), kind=mytype),&
        aimag(aaa(j,nym,k,3))-tmp2*aimag(aaa(j,nym-1,k,4)), kind=mytype)
   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,nym,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,nym-1,k), kind=mytype)
   else
      tmp1=0.
      tmp3=0.
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,nym,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,nym-1,k))
   else
      tmp2=0.
      tmp4=0.
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,nym,k)=cmplx(tmp3,tmp4, kind=mytype)
   
   if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,nym-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,nym-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,nym-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,nym-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,nym-1,k)=cmplx(real(eee(j,nym-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,nym,k), kind=mytype),&
        aimag(eee(j,nym-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,nym,k)), kind=mytype)
enddo
enddo

do i=nym-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,i,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v2

!********************************************************************
!
subroutine channel (ux)
!
!********************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux

integer :: j,i,k,code
real(mytype) :: can,ut3,ut,ut4

ut3=0.
do k=1,ysize(3)
do i=1,ysize(1)
   ut=0.
   do j=1,ny-1
      if (istret.ne.0) ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
      if (istret.eq.0) ut=ut+(yly/(ny-1))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
   enddo
   ut=ut/yly
   ut3=ut3+ut
enddo
enddo
ut3=ut3/ysize(1)/ysize(3)

call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
ut4=ut4/nproc

!can=-2.*xnu*gdt(itr) ! Poisseuille    
can=-(2./3.-ut4) ! constant flow rate

if (nrank==0) print *,nrank,'UT',ut4,can

do k=1,ysize(3)
do i=1,ysize(1)
do j=2,ny-1
   ux(i,j,k)=-can+ux(i,j,k)
enddo
enddo
enddo

return
end subroutine channel



!*****************************************************************
!
subroutine collect_data()
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

integer :: i,j,imin,ii,code
real(mytype), dimension(200) :: ta
integer,dimension(200) :: tai,tbi

!TOTAL TIME FOR COLLECTION T=50000*0.01=500
!I want 100 files
!100 numbers from 0 to 500

if (nrank==0) then
call random_number(ta)
do i=1,200
   tai(i)=int(ta(i)*2000./dt)
enddo
do j=1,200
imin=999998
ii=0
do i=1,200
   if (tai(i).le.imin) then
      imin=tai(i)
      ii=i
   endif
enddo
tbi(j)=imin
tai(ii)=999999
enddo
idata=tbi
do i=1,200
print *,i,idata(i),'VALUE COLLECTION DATA'
enddo
endif

call MPI_BCAST(idata,200,MPI_INTEGER,0,MPI_COMM_WORLD,code)

end subroutine collect_data

subroutine load_disks()

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE var
USE MPI

implicit none

integer :: i,j,k,fhx,fhz,ierrorx,ierrorz
integer :: fhx1,fhz1,fhx2,fhz2,fhx3,fhz3,fhx4,fhz4,fhx5,fhz5,fhx6,fhz6
integer :: fhx7,fhz7,fhx8,fhz8,fhx9,fhz9,fhx10,fhz10,fhx11,fhz11,fhx12,fhz12
integer :: ierrorx1,ierrorz1,ierrorx2,ierrorz2,ierrorx3,ierrorz3
integer :: ierrorx4,ierrorz4,ierrorx5,ierrorz5,ierrorx6,ierrorz6
integer :: ierrorx7,ierrorz7,ierrorx8,ierrorz8,ierrorx9,ierrorz9
integer :: ierrorx10,ierrorz10,ierrorx11,ierrorz11,ierrorx12,ierrorz12
integer (kind=MPI_OFFSET_KIND) :: dispx,dispz
integer (kind=MPI_OFFSET_KIND) :: dispx1,dispz1,dispx2,dispz2
integer (kind=MPI_OFFSET_KIND) :: dispx3,dispz3,dispx4,dispz4
integer (kind=MPI_OFFSET_KIND) :: dispx5,dispz5,dispx6,dispz6
integer (kind=MPI_OFFSET_KIND) :: dispx7,dispz7,dispx8,dispz8
integer (kind=MPI_OFFSET_KIND) :: dispx9,dispz9,dispx10,dispz10
integer (kind=MPI_OFFSET_KIND) :: dispx11,dispz11,dispx12,dispz12

dispx = 0_MPI_OFFSET_KIND  ;dispz = 0_MPI_OFFSET_KIND
dispx1 = 0_MPI_OFFSET_KIND ;dispz1 = 0_MPI_OFFSET_KIND
dispx2 = 0_MPI_OFFSET_KIND ;dispz2 = 0_MPI_OFFSET_KIND
dispx3 = 0_MPI_OFFSET_KIND ;dispz3 = 0_MPI_OFFSET_KIND
dispx4 = 0_MPI_OFFSET_KIND ;dispz4 = 0_MPI_OFFSET_KIND
dispx5 = 0_MPI_OFFSET_KIND ;dispz5 = 0_MPI_OFFSET_KIND
dispx6 = 0_MPI_OFFSET_KIND ;dispz6 = 0_MPI_OFFSET_KIND
dispx7 = 0_MPI_OFFSET_KIND ;dispz7 = 0_MPI_OFFSET_KIND
dispx8 = 0_MPI_OFFSET_KIND ;dispz8 = 0_MPI_OFFSET_KIND
dispx9 = 0_MPI_OFFSET_KIND ;dispz9 = 0_MPI_OFFSET_KIND
dispx10 = 0_MPI_OFFSET_KIND;dispz10 = 0_MPI_OFFSET_KIND
dispx11 = 0_MPI_OFFSET_KIND;dispz11 = 0_MPI_OFFSET_KIND
dispx12 = 0_MPI_OFFSET_KIND;dispz12 = 0_MPI_OFFSET_KIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! all disks for control
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'disks_opac_x.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx, ierrorx)
call decomp_2d_read_var(fhx,dispx,1,tabx)
call MPI_FILE_CLOSE(fhx,ierrorx)

call MPI_FILE_OPEN(MPI_COMM_WORLD, 'disks_opac_z.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz, ierrorz)
call decomp_2d_read_var(fhz,dispz,1,tabz)
call MPI_FILE_CLOSE(fhz,ierrorz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 1
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_1/disks_opac_x1.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx1, ierrorx1)
call decomp_2d_read_var(fhx1,dispx1,1,tabx1)
call MPI_FILE_CLOSE(fhx1,ierrorx1)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_1/disks_opac_z1.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz1, ierrorz1)
call decomp_2d_read_var(fhz1,dispz1,1,tabz1)
call MPI_FILE_CLOSE(fhz1,ierrorz1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 2
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_2/disks_opac_x2.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx2, ierrorx2)
call decomp_2d_read_var(fhx2,dispx2,1,tabx2)
call MPI_FILE_CLOSE(fhx2,ierrorx2)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_2/disks_opac_z2.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz2, ierrorz2)
call decomp_2d_read_var(fhz2,dispz2,1,tabz2)
call MPI_FILE_CLOSE(fhz2,ierrorz2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 3
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_3/disks_opac_x3.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx3, ierrorx3)
call decomp_2d_read_var(fhx3,dispx3,1,tabx3)
call MPI_FILE_CLOSE(fhx3,ierrorx3)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_3/disks_opac_z3.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz3, ierrorz3)
call decomp_2d_read_var(fhz3,dispz3,1,tabz3)
call MPI_FILE_CLOSE(fhz3,ierrorz3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 4
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_4/disks_opac_x4.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx4, ierrorx4)
call decomp_2d_read_var(fhx4,dispx4,1,tabx4)
call MPI_FILE_CLOSE(fhx4,ierrorx4)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_4/disks_opac_z4.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz4, ierrorz4)
call decomp_2d_read_var(fhz4,dispz4,1,tabz4)
call MPI_FILE_CLOSE(fhz4,ierrorz4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 5
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_5/disks_opac_x5.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx5, ierrorx5)
call decomp_2d_read_var(fhx5,dispx5,1,tabx5)
call MPI_FILE_CLOSE(fhx5,ierrorx5)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_5/disks_opac_z5.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz5, ierrorz5)
call decomp_2d_read_var(fhz5,dispz5,1,tabz5)
call MPI_FILE_CLOSE(fhz5,ierrorz5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 6
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_6/disks_opac_x6.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx6, ierrorx6)
call decomp_2d_read_var(fhx6,dispx6,1,tabx6)
call MPI_FILE_CLOSE(fhx6,ierrorx6)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_6/disks_opac_z6.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz6, ierrorz6)
call decomp_2d_read_var(fhz6,dispz6,1,tabz6)
call MPI_FILE_CLOSE(fhz6,ierrorz6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 7
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_7/disks_opac_x7.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx7, ierrorx7)
call decomp_2d_read_var(fhx7,dispx7,1,tabx7)
call MPI_FILE_CLOSE(fhx7,ierrorx7)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_7/disks_opac_z7.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz7, ierrorz7)
call decomp_2d_read_var(fhz7,dispz7,1,tabz7)
call MPI_FILE_CLOSE(fhz7,ierrorz7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 8
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_8/disks_opac_x8.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx8, ierrorx8)
call decomp_2d_read_var(fhx8,dispx8,1,tabx8)
call MPI_FILE_CLOSE(fhx8,ierrorx8)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_8/disks_opac_z8.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz8, ierrorz8)
call decomp_2d_read_var(fhz8,dispz8,1,tabz8)
call MPI_FILE_CLOSE(fhz8,ierrorz8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 9
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_9/disks_opac_x9.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx9, ierrorx9)
call decomp_2d_read_var(fhx9,dispx9,1,tabx9)
call MPI_FILE_CLOSE(fhx9,ierrorx9)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_9/disks_opac_z9.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz9, ierrorz9)
call decomp_2d_read_var(fhz9,dispz9,1,tabz9)
call MPI_FILE_CLOSE(fhz9,ierrorz9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 10
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_10/disks_opac_x10.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx10, ierrorx10)
call decomp_2d_read_var(fhx10,dispx10,1,tabx10)
call MPI_FILE_CLOSE(fhx10,ierrorx10)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_10/disks_opac_z10.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz10, ierrorz10)
call decomp_2d_read_var(fhz10,dispz10,1,tabz10)
call MPI_FILE_CLOSE(fhz10,ierrorz10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 11
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_11/disks_opac_x11.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx11, ierrorx11)
call decomp_2d_read_var(fhx11,dispx11,1,tabx11)
call MPI_FILE_CLOSE(fhx11,ierrorx11)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_11/disks_opac_z11.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz11, ierrorz11)
call decomp_2d_read_var(fhz11,dispz11,1,tabz11)
call MPI_FILE_CLOSE(fhz11,ierrorz11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 12
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_12/disks_opac_x12.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx12, ierrorx12)
call decomp_2d_read_var(fhx12,dispx12,1,tabx12)
call MPI_FILE_CLOSE(fhx12,ierrorx12)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_12/disks_opac_z12.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz12, ierrorz12)
call decomp_2d_read_var(fhz12,dispz12,1,tabz12)
call MPI_FILE_CLOSE(fhz12,ierrorz12)

if (nrank.eq.0) then
   write(6,*) "###############################################"
   write(6,*) "       read disk_opac_x.dat file done          "
   write(6,*) "       read disk_opac_z.dat file done          "
   write(6,*) "###############################################"
   write(6,*) "       read disk_opac_x*.dat file done     "
   write(6,*) "###############################################"
end if

do k=1,xsize(3)
   do i=1,xsize(1)
      uxd(i,k)=tabx(i,1,k)
   end do
end do

do k=1,xsize(3)
   do i=1,xsize(1)
      uzd(i,k)=tabz(i,1,k)
   end do
end do

end subroutine load_disks


subroutine load_disks_ratio()

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE var
USE MPI

implicit none

integer :: i,j,k
integer :: fhx_rt1,fhz_rt1,fhx_rt2,fhz_rt2,fhx_rt3,fhz_rt3,fhx_rt4,fhz_rt4,fhx_rt5,fhz_rt5,fhx_rt6,fhz_rt6
integer :: fhx_rt7,fhz_rt7,fhx_rt8,fhz_rt8,fhx_rt9,fhz_rt9,fhx_rt10,fhz_rt10,fhx_rt11,fhz_rt11,fhx_rt12,fhz_rt12
integer :: ierrorx_rt1,ierrorz_rt1,ierrorx_rt2,ierrorz_rt2,ierrorx_rt3,ierrorz_rt3
integer :: ierrorx_rt4,ierrorz_rt4,ierrorx_rt5,ierrorz_rt5,ierrorx_rt6,ierrorz_rt6
integer :: ierrorx_rt7,ierrorz_rt7,ierrorx_rt8,ierrorz_rt8,ierrorx_rt9,ierrorz_rt9
integer :: ierrorx_rt10,ierrorz_rt10,ierrorx_rt11,ierrorz_rt11,ierrorx_rt12,ierrorz_rt12
integer (kind=MPI_OFFSET_KIND) :: dispx_rt1,dispz_rt1,dispx_rt2,dispz_rt2
integer (kind=MPI_OFFSET_KIND) :: dispx_rt3,dispz_rt3,dispx_rt4,dispz_rt4
integer (kind=MPI_OFFSET_KIND) :: dispx_rt5,dispz_rt5,dispx_rt6,dispz_rt6
integer (kind=MPI_OFFSET_KIND) :: dispx_rt7,dispz_rt7,dispx_rt8,dispz_rt8
integer (kind=MPI_OFFSET_KIND) :: dispx_rt9,dispz_rt9,dispx_rt10,dispz_rt10
integer (kind=MPI_OFFSET_KIND) :: dispx_rt11,dispz_rt11,dispx_rt12,dispz_rt12

dispx_rt1 = 0_MPI_OFFSET_KIND ;dispz_rt1 = 0_MPI_OFFSET_KIND
dispx_rt2 = 0_MPI_OFFSET_KIND ;dispz_rt2 = 0_MPI_OFFSET_KIND
dispx_rt3 = 0_MPI_OFFSET_KIND ;dispz_rt3 = 0_MPI_OFFSET_KIND
dispx_rt4 = 0_MPI_OFFSET_KIND ;dispz_rt4 = 0_MPI_OFFSET_KIND
dispx_rt5 = 0_MPI_OFFSET_KIND ;dispz_rt5 = 0_MPI_OFFSET_KIND
dispx_rt6 = 0_MPI_OFFSET_KIND ;dispz_rt6 = 0_MPI_OFFSET_KIND
dispx_rt7 = 0_MPI_OFFSET_KIND ;dispz_rt7 = 0_MPI_OFFSET_KIND
dispx_rt8 = 0_MPI_OFFSET_KIND ;dispz_rt8 = 0_MPI_OFFSET_KIND
dispx_rt9 = 0_MPI_OFFSET_KIND ;dispz_rt9 = 0_MPI_OFFSET_KIND
dispx_rt10 = 0_MPI_OFFSET_KIND;dispz_rt10 = 0_MPI_OFFSET_KIND
dispx_rt11 = 0_MPI_OFFSET_KIND;dispz_rt11 = 0_MPI_OFFSET_KIND
dispx_rt12 = 0_MPI_OFFSET_KIND;dispz_rt12 = 0_MPI_OFFSET_KIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 1
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_1/disks_opac_x1_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt1, ierrorx_rt1)
call decomp_2d_read_var(fhx_rt1,dispx_rt1,1,tabx_rt1)
call MPI_FILE_CLOSE(fhx_rt1,ierrorx_rt1)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_1/disks_opac_z1_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt1, ierrorz_rt1)
call decomp_2d_read_var(fhz_rt1,dispz_rt1,1,tabz_rt1)
call MPI_FILE_CLOSE(fhz_rt1,ierrorz_rt1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 2
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_2/disks_opac_x2_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt2, ierrorx_rt2)
call decomp_2d_read_var(fhx_rt2,dispx_rt2,1,tabx_rt2)
call MPI_FILE_CLOSE(fhx_rt2,ierrorx_rt2)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_2/disks_opac_z2_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt2, ierrorz_rt2)
call decomp_2d_read_var(fhz_rt2,dispz_rt2,1,tabz_rt2)
call MPI_FILE_CLOSE(fhz_rt2,ierrorz_rt2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 3
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_3/disks_opac_x3_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt3, ierrorx_rt3)
call decomp_2d_read_var(fhx_rt3,dispx_rt3,1,tabx_rt3)
call MPI_FILE_CLOSE(fhx_rt3,ierrorx_rt3)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_3/disks_opac_z3_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt3, ierrorz_rt3)
call decomp_2d_read_var(fhz_rt3,dispz_rt3,1,tabz_rt3)
call MPI_FILE_CLOSE(fhz_rt3,ierrorz_rt3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 4
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_4/disks_opac_x4_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt4, ierrorx_rt4)
call decomp_2d_read_var(fhx_rt4,dispx_rt4,1,tabx_rt4)
call MPI_FILE_CLOSE(fhx_rt4,ierrorx_rt4)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_4/disks_opac_z4_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt4, ierrorz_rt4)
call decomp_2d_read_var(fhz_rt4,dispz_rt4,1,tabz_rt4)
call MPI_FILE_CLOSE(fhz_rt4,ierrorz_rt4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 5
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_5/disks_opac_x5_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt5, ierrorx_rt5)
call decomp_2d_read_var(fhx_rt5,dispx_rt5,1,tabx_rt5)
call MPI_FILE_CLOSE(fhx_rt5,ierrorx_rt5)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_5/disks_opac_z5_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt5, ierrorz_rt5)
call decomp_2d_read_var(fhz_rt5,dispz_rt5,1,tabz_rt5)
call MPI_FILE_CLOSE(fhz_rt5,ierrorz_rt5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 6
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_6/disks_opac_x6_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt6, ierrorx_rt6)
call decomp_2d_read_var(fhx_rt6,dispx_rt6,1,tabx_rt6)
call MPI_FILE_CLOSE(fhx_rt6,ierrorx_rt6)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_6/disks_opac_z6_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt6, ierrorz_rt6)
call decomp_2d_read_var(fhz_rt6,dispz_rt6,1,tabz_rt6)
call MPI_FILE_CLOSE(fhz_rt6,ierrorz_rt6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 7
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_7/disks_opac_x7_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt7, ierrorx_rt7)
call decomp_2d_read_var(fhx_rt7,dispx_rt7,1,tabx_rt7)
call MPI_FILE_CLOSE(fhx_rt7,ierrorx_rt7)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_7/disks_opac_z7_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt7, ierrorz_rt7)
call decomp_2d_read_var(fhz_rt7,dispz_rt7,1,tabz_rt7)
call MPI_FILE_CLOSE(fhz_rt7,ierrorz_rt7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 8
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_8/disks_opac_x8_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt8, ierrorx_rt8)
call decomp_2d_read_var(fhx_rt8,dispx_rt8,1,tabx_rt8)
call MPI_FILE_CLOSE(fhx_rt8,ierrorx_rt8)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_8/disks_opac_z8_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt8, ierrorz_rt8)
call decomp_2d_read_var(fhz_rt8,dispz_rt8,1,tabz_rt8)
call MPI_FILE_CLOSE(fhz_rt8,ierrorz_rt8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 9
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_9/disks_opac_x9_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt9, ierrorx_rt9)
call decomp_2d_read_var(fhx_rt9,dispx_rt9,1,tabx_rt9)
call MPI_FILE_CLOSE(fhx_rt9,ierrorx_rt9)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_9/disks_opac_z9_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt9, ierrorz_rt9)
call decomp_2d_read_var(fhz_rt9,dispz_rt9,1,tabz_rt9)
call MPI_FILE_CLOSE(fhz_rt9,ierrorz_rt9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 10
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_10/disks_opac_x10_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt10, ierrorx_rt10)
call decomp_2d_read_var(fhx_rt10,dispx_rt10,1,tabx_rt10)
call MPI_FILE_CLOSE(fhx_rt10,ierrorx_rt10)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_10/disks_opac_z10_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt10, ierrorz_rt10)
call decomp_2d_read_var(fhz_rt10,dispz_rt10,1,tabz_rt10)
call MPI_FILE_CLOSE(fhz_rt10,ierrorz_rt10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 11
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_11/disks_opac_x11_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt11, ierrorx_rt11)
call decomp_2d_read_var(fhx_rt11,dispx_rt11,1,tabx_rt11)
call MPI_FILE_CLOSE(fhx_rt11,ierrorx_rt11)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_11/disks_opac_z11_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt11, ierrorz_rt11)
call decomp_2d_read_var(fhz_rt11,dispz_rt11,1,tabz_rt11)
call MPI_FILE_CLOSE(fhz_rt11,ierrorz_rt11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! disk 12
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_12/disks_opac_x12_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhx_rt12, ierrorx_rt12)
call decomp_2d_read_var(fhx_rt12,dispx_rt12,1,tabx_rt12)
call MPI_FILE_CLOSE(fhx_rt12,ierrorx_rt12)
call MPI_FILE_OPEN(MPI_COMM_WORLD, 'twelve_disk_12/disks_opac_z12_rt.dat', &
     MPI_MODE_RDONLY, MPI_INFO_NULL, &
     fhz_rt12, ierrorz_rt12)
call decomp_2d_read_var(fhz_rt12,dispz_rt12,1,tabz_rt12)
call MPI_FILE_CLOSE(fhz_rt12,ierrorz_rt12)

if (nrank.eq.0) then
   write(6,*) "###############################################"
   write(6,*) "       read disk_opac_x*_rt.dat file done     "
   write(6,*) "###############################################"
end if

end subroutine load_disks_ratio

subroutine pw_output(prd_in,nd)

USE param
USE decomp_2d, only : nrank

implicit none

integer :: nd
real(mytype) :: prd_in,coeff_w,coeff_psp,u_bulk,re_tau0,re_p

re_p=4200.0
re_tau0=179.5
u_bulk=2.0/3.0
coeff_w=(-1.0/re_p)
coeff_psp=coeff_w/(2.0*u_bulk*(re_tau0/re_p)*(re_tau0/re_p))

! 991 format('vort',i6.6,'.dat')
! write(filename, 991) itime/imodulo

if(nrank.eq.0) then
   write(6,*) "##############################################################################################"
!    write(6,*) "W_",nd," = ", coeff_w*prd_in
   write(6,"(a,i2.2,x,a,x,f20.14)") "Psp_",nd," = ", 100.0*coeff_psp*prd_in
!    write(6,*) "##############################################################################################"
end if

end subroutine pw_output
