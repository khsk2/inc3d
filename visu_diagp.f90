!############################################################################
!
subroutine CALC_VORT_GRADVEL(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,&
     th1,ti1,di1,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,vrtx,vrty,vrtz,&
     ddudx1,ddudy1,ddudz1,ddvdx1,ddvdy1,ddvdz1,ddwdx1,ddwdy1,ddwdz1)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI

implicit none

TYPE(DECOMP_INFO) :: phG
integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: vrtx,vrty,vrtz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: ddudx1,ddudy1,ddudz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: ddvdx1,ddvdy1,ddvdz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: ddwdx1,ddwdy1,ddwdz1

call perfon('calc_vort_gradvel')
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

!############################################################################
!VORTICITY - x-y-z component
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         vrtx(i,j,k)=tf1(i,j,k)-th1(i,j,k)
         vrty(i,j,k)=tg1(i,j,k)-tc1(i,j,k)
         vrtz(i,j,k)=tb1(i,j,k)-td1(i,j,k)
      end do
   end do
end do

!############################################################################
! VELOCITY GRADIENT TENSOR elements
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         ddudx1(i,j,k)=ta1(i,j,k);ddudy1(i,j,k)=td1(i,j,k);ddudz1(i,j,k)=tg1(i,j,k)
         ddvdx1(i,j,k)=tb1(i,j,k);ddvdy1(i,j,k)=te1(i,j,k);ddvdz1(i,j,k)=th1(i,j,k)
         ddwdx1(i,j,k)=tc1(i,j,k);ddwdy1(i,j,k)=tf1(i,j,k);ddwdz1(i,j,k)=ti1(i,j,k)
      end do
   end do
end do

!############################################################################

call perfoff
return
end subroutine CALC_VORT_GRADVEL

!############################################################################
!
subroutine DISSIPATION(ux1,ddudx1,ddudy1,ddudz1,ddvdx1,&
     ddvdy1,ddvdz1,ddwdx1,ddwdy1,ddwdz1)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: dissp1_tavg,dubardy_tavg

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,ddudx1,ddudy1,ddudz1,diss1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddvdx1,ddvdy1,ddvdz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddwdx1,ddwdy1,ddwdz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tdiss1,tux1,tmp_tux1_spavg,dtmp_tux1_spavg,di2
real(mytype),dimension(ysize(2)) :: tux1_spavg,dissp1_spavg

call perfon('dissipation')

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         diss1(i,j,k)=2*ddudx1(i,j,k)*ddudx1(i,j,k)+2*ddvdy1(i,j,k)*ddvdy1(i,j,k)+&
              2*ddwdz1(i,j,k)*ddwdz1(i,j,k)+(ddudy1(i,j,k)+ddvdx1(i,j,k))*(ddudy1(i,j,k)&
              +ddvdx1(i,j,k))+(ddudz1(i,j,k)+ddwdx1(i,j,k))*(ddudz1(i,j,k)+ddwdx1(i,j,k))+&
              (ddvdz1(i,j,k)+ddwdy1(i,j,k))*(ddvdz1(i,j,k)+ddwdy1(i,j,k))
      end do
   end do
end do

call transpose_x_to_y(diss1,tdiss1)
call space_avg(tdiss1,dissp1_spavg)

do j=1,ny
   dissp1_tavg(j)=dissp1_tavg(j)+dissp1_spavg(j)
end do

call transpose_x_to_y(ux1,tux1)
call space_avg(tux1,tux1_spavg)

do k=1,ysize(3)
   do j=1,ysize(2)
      do i=1,ysize(1)
         tmp_tux1_spavg(i,j,k)=tux1_spavg(j)
      end do
   end do
end do

call dery (dtmp_tux1_spavg,tmp_tux1_spavg,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),0)

do j=1,ny
   dubardy_tavg(j)=dubardy_tavg(j)+dtmp_tux1_spavg(1,j,1)
end do

call perfoff
end subroutine DISSIPATION

!############################################################################
!
subroutine TURB_CNV(ux1,uy1,uz1)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: turb_conv

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,uuv,uv,vvv,wwv
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu,tuuv,tuv,tvvv,twwv,di2,tmp2,tmp3
real(mytype),dimension(ysize(2)) :: tmp1,tu_spavg,tuuv_spavg,tuv_spavg,tvvv_spavg,twwv_spavg

call perfon('turb_conv')

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         uuv(i,j,k)=ux1(i,j,k)*ux1(i,j,k)*uy1(i,j,k)
         uv(i,j,k)=ux1(i,j,k)*uy1(i,j,k)
         vvv(i,j,k)=uy1(i,j,k)*uy1(i,j,k)*uy1(i,j,k)
         wwv(i,j,k)=uz1(i,j,k)*uz1(i,j,k)*uy1(i,j,k)
      end do
   end do
end do

call transpose_x_to_y(ux1,tu)
call transpose_x_to_y(uuv,tuuv);call transpose_x_to_y(uv,tuv)
call transpose_x_to_y(vvv,tvvv);call transpose_x_to_y(wwv,twwv)

call space_avg(tu,tu_spavg)
call space_avg(tuuv,tuuv_spavg);call space_avg(tuv,tuv_spavg)
call space_avg(tvvv,tvvv_spavg);call space_avg(twwv,twwv_spavg)

do j=1,ny
   tmp1(j)=-0.5*(tuuv_spavg(j)-2.0*tu_spavg(j)*tuv_spavg(j)+&
        tvvv_spavg(j)+twwv_spavg(j))
end do

do k=1,ysize(3)
   do j=1,ysize(2)
      do i=1,ysize(1)
         tmp2(i,j,k)=tmp1(j)
      end do
   end do
end do

call dery (tmp3,tmp2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),0)

do j=1,ny
   turb_conv(j)=turb_conv(j)+tmp3(10,j,15)
end do

call perfoff
end subroutine TURB_CNV

!############################################################################
!
subroutine VISC_DFF(ux1,uy1,uz1)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: visc_diff

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,uu,vv,ww
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu,tuu,tvv,tww,tmp2,vscdf,di2
real(mytype),dimension(ysize(2)) :: tu_spavg,tuu_spavg,tvv_spavg,tww_spavg,tmp1

call perfon('visc_diff')

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         uu(i,j,k)=ux1(i,j,k)*ux1(i,j,k)
         vv(i,j,k)=uy1(i,j,k)*uy1(i,j,k)
         ww(i,j,k)=uz1(i,j,k)*uz1(i,j,k)
      end do
   end do
end do

call transpose_x_to_y(ux1,tu)
call transpose_x_to_y(uu,tuu)
call transpose_x_to_y(vv,tvv)
call transpose_x_to_y(ww,tww)

call space_avg(tu,tu_spavg)
call space_avg(tuu,tuu_spavg)
call space_avg(tvv,tvv_spavg)
call space_avg(tww,tww_spavg)

do j=1,ny
   tmp1(j)=0.5*(tuu_spavg(j)-tu_spavg(j)*tu_spavg(j)+&
        tvv_spavg(j)+tww_spavg(j))
end do

do k=1,ysize(3)
   do j=1,ysize(2)
      do i=1,ysize(1)
         tmp2(i,j,k)=tmp1(j)
      end do
   end do
end do

call deryy (vscdf,tmp2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)

do j=1,ny
   visc_diff(j)=visc_diff(j)+vscdf(10,j,15)
end do

call perfoff
end subroutine VISC_DFF

!############################################################################
!
subroutine STATISTIC(ux1,uy1,uz1,vortx1,vorty1,vortz1,&
     ddudx1,ddudy1,ddudz1,ddvdx1,ddvdy1,ddvdz1,&
     ddwdx1,ddwdy1,ddwdz1,dirout)
!
!############################################################################
USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only: S_dudy_lw_tot,S_dudy_uw_tot,S_dudy_tot,c_f_tot,&
     S_ubar,S_vbar,S_wbar,S_urms,S_vrms,S_wrms,S_uv,S_uw,S_vw,&
     S_vortxbar,S_vortybar,S_vortzbar,S_vortxrms,S_vortyrms,S_vortzrms,&
     ud,vd,wd,udrms,vdrms,wdrms,&
     vortxd,vortyd,vortzd,vortxdrms,vortydrms,vortzdrms,udvdrey,&
     uturb,vturb,wturb,vortxturb,vortyturb,vortzturb,&
     u_tavg,v_tavg,w_tavg,vortx_tavg,vorty_tavg,vortz_tavg,&
     u2_tavg,v2_tavg,w2_tavg,vortx2_tavg,vorty2_tavg,vortz2_tavg,&
     uv_tavg

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddudx1,ddudy1,ddudz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddvdx1,ddvdy1,ddvdz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ddwdx1,ddwdy1,ddwdz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tux1,tuy1,tuz1,tddudy1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortx1,vorty1,vortz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tvortx1,tvorty1,tvortz1
real(mytype),dimension(ysize(2)) :: ubar_red,urms_red
real(mytype),dimension(ysize(2)) :: vbar_red,vrms_red
real(mytype),dimension(ysize(2)) :: wbar_red,wrms_red
real(mytype),dimension(ysize(2)) :: vortxbar,vortxbar_red,vortxrms_red
real(mytype),dimension(ysize(2)) :: vortybar,vortybar_red,vortyrms_red
real(mytype),dimension(ysize(2)) :: vortzbar,vortzbar_red,vortzrms_red
real(mytype),dimension(ysize(2)) :: uv_red,uw,uw_red,vw,vw_red
real(mytype) :: dudy_mean_lw,dudy_mean_uw,dudy_mean,c_f
real(mytype) :: dudy_mean_lw_red,dudy_mean_uw_red
real(mytype),dimension(5,ilast-ifirst+1) :: drg_inst
integer :: i,j,k,ss,code
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods
character(len=128) :: dirout

call perfon('statistic')
if(itime.eq.ifirst.and.nrank.eq.0) then
   open(500,file=trim(dirout)//'/drg.txt',access='append')
end if

call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
     dims, dummy_periods, dummy_coords, code)

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         u_tavg(i,j,k)=u_tavg(i,j,k)+ux1(i,j,k)
         v_tavg(i,j,k)=v_tavg(i,j,k)+uy1(i,j,k)
         w_tavg(i,j,k)=w_tavg(i,j,k)+uz1(i,j,k)
         u2_tavg(i,j,k)=u2_tavg(i,j,k)+ux1(i,j,k)*ux1(i,j,k)
         v2_tavg(i,j,k)=v2_tavg(i,j,k)+uy1(i,j,k)*uy1(i,j,k)
         w2_tavg(i,j,k)=w2_tavg(i,j,k)+uz1(i,j,k)*uz1(i,j,k)

         vortx_tavg(i,j,k)=vortx_tavg(i,j,k)+vortx1(i,j,k)
         vorty_tavg(i,j,k)=vorty_tavg(i,j,k)+vorty1(i,j,k)
         vortz_tavg(i,j,k)=vortz_tavg(i,j,k)+vortz1(i,j,k)
         vortx2_tavg(i,j,k)=vortx2_tavg(i,j,k)+vortx1(i,j,k)*vortx1(i,j,k)
         vorty2_tavg(i,j,k)=vorty2_tavg(i,j,k)+vorty1(i,j,k)*vorty1(i,j,k)
         vortz2_tavg(i,j,k)=vortz2_tavg(i,j,k)+vortz1(i,j,k)*vortz1(i,j,k)

         uv_tavg(i,j,k)=uv_tavg(i,j,k)+ux1(i,j,k)*uy1(i,j,k)
      end do
   end do
end do

!!!!!!!!!!!!!!!!!!!
! mean quantities !
!!!!!!!!!!!!!!!!!!!

call transpose_x_to_y(ux1,tux1)
call transpose_x_to_y(uy1,tuy1)
call transpose_x_to_y(uz1,tuz1)

call transpose_x_to_y(vortx1,tvortx1)
call transpose_x_to_y(vorty1,tvorty1)
call transpose_x_to_y(vortz1,tvortz1)

call space_avg(tux1,ubar_red)
call space_avg(tuy1,vbar_red)
call space_avg(tuz1,wbar_red)

call space_avg(tvortx1,vortxbar_red)
call space_avg(tvorty1,vortybar_red)
call space_avg(tvortz1,vortzbar_red)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rms quantities and Reynolds stresses !            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         ud(i,j,k)=ud(i,j,k)+(tux1(i,j,k)-ubar_red(j))
         vd(i,j,k)=vd(i,j,k)+(tuy1(i,j,k)-vbar_red(j))
         wd(i,j,k)=wd(i,j,k)+(tuz1(i,j,k)-wbar_red(j))

         vortxd(i,j,k)=vortxd(i,j,k)+(tvortx1(i,j,k)-vortxbar_red(j))
         vortyd(i,j,k)=vortyd(i,j,k)+(tvorty1(i,j,k)-vortybar_red(j))
         vortzd(i,j,k)=vortzd(i,j,k)+(tvortz1(i,j,k)-vortzbar_red(j))

         uturb(i,j,k)=uturb(i,j,k)+(tux1(i,j,k)-ud(i,j,k)-ubar_red(j))
         vturb(i,j,k)=vturb(i,j,k)+(tuy1(i,j,k)-vd(i,j,k)-vbar_red(j))
         wturb(i,j,k)=wturb(i,j,k)+(tuz1(i,j,k)-wd(i,j,k)-wbar_red(j))

         vortxturb(i,j,k)=vortxturb(i,j,k)+(tvortx1(i,j,k)-vortxd(i,j,k)-vortxbar_red(j))
         vortyturb(i,j,k)=vortyturb(i,j,k)+(tvorty1(i,j,k)-vortyd(i,j,k)-vortybar_red(j))
         vortzturb(i,j,k)=vortzturb(i,j,k)+(tvortz1(i,j,k)-vortzd(i,j,k)-vortzbar_red(j))
      end do
   end do
end do

call rms_comp(tux1,ubar_red,urms_red)
call rms_comp(tuy1,vbar_red,vrms_red)
call rms_comp(tuz1,wbar_red,wrms_red)
call rms_comp(tvortx1,vortxbar_red,vortxrms_red)
call rms_comp(tvorty1,vortybar_red,vortyrms_red)
call rms_comp(tvortz1,vortzbar_red,vortzrms_red)

call rey_comp(tux1,tuy1,ubar_red,vbar_red,uv_red)
call rey_comp(tux1,tuz1,ubar_red,wbar_red,uw_red)
call rey_comp(tuy1,tuz1,vbar_red,wbar_red,vw_red)

do j=1,ny
   S_ubar(j)=S_ubar(j)+ubar_red(j)
   S_vbar(j)=S_vbar(j)+vbar_red(j)
   S_wbar(j)=S_wbar(j)+wbar_red(j)
   S_vortxbar(j)=S_vortxbar(j)+vortxbar_red(j)
   S_vortybar(j)=S_vortybar(j)+vortybar_red(j)
   S_vortzbar(j)=S_vortzbar(j)+vortzbar_red(j)
   S_urms(j)=S_urms(j)+urms_red(j)
   S_vrms(j)=S_vrms(j)+vrms_red(j)
   S_wrms(j)=S_wrms(j)+wrms_red(j)
   S_uv(j)=S_uv(j)+uv_red(j)
   S_uw(j)=S_uw(j)+uw_red(j)
   S_vw(j)=S_vw(j)+vw_red(j)
   S_vortxrms(j)=S_vortxrms(j)+vortxrms_red(j)
   S_vortyrms(j)=S_vortyrms(j)+vortyrms_red(j)
   S_vortzrms(j)=S_vortzrms(j)+vortzrms_red(j)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! drag-related quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

call transpose_x_to_y(ddudy1,tddudy1)
dudy_mean_lw=0.0
dudy_mean_uw=0.0

do k=1,ysize(3)
   do i=1,ysize(1)
      dudy_mean_lw=dudy_mean_lw+tddudy1(i,1,k)
      dudy_mean_uw=dudy_mean_uw+tddudy1(i,ny,k)
   end do
end do
dudy_mean_lw=dudy_mean_lw/ysize(1)/ysize(3)
dudy_mean_uw=dudy_mean_uw/ysize(1)/ysize(3)

call mpi_allreduce(dudy_mean_lw,dudy_mean_lw_red,1,real_type,&
     MPI_SUM,MPI_COMM_WORLD,code)
call mpi_allreduce(dudy_mean_uw,dudy_mean_uw_red,1,real_type,&
     MPI_SUM,MPI_COMM_WORLD,code)

dudy_mean_lw_red=dudy_mean_lw_red/nproc
dudy_mean_uw_red=dudy_mean_uw_red/nproc

dudy_mean=(dudy_mean_lw_red-dudy_mean_uw_red)/(2.0)
c_f=4.5*dudy_mean/re

if(nrank.eq.0) then
   ss=itime+(ilast-ifirst+1)-ilast
   drg_inst(1,ss)=dble(itime*imodulo*dt)
   drg_inst(2,ss)=dudy_mean_lw_red
   drg_inst(3,ss)=dudy_mean_uw_red
   drg_inst(4,ss)=dudy_mean
   drg_inst(5,ss)=c_f
   write(500,*) drg_inst(1,ss),drg_inst(2,ss),drg_inst(3,ss),&
        drg_inst(4,ss),drg_inst(5,ss)
end if

S_dudy_lw_tot=S_dudy_lw_tot+dudy_mean_lw_red
S_dudy_uw_tot=S_dudy_uw_tot+dudy_mean_uw_red
S_dudy_tot=S_dudy_tot+dudy_mean
c_f_tot=c_f_tot+c_f

call perfoff
end subroutine STATISTIC

!############################################################################
!
subroutine STAT_DISK(u_to_avg,v_to_avg,w_to_avg,vortx_to_avg,vorty_to_avg,vortz_to_avg,&
     utavg,vtavg,wtavg,vortxtavg,vortytavg,vortztavg,&
     u2tavg,v2tavg,w2tavg,vortx2tavg,vorty2tavg,vortz2tavg,&
     uvtavg,u_b,v_b,w_b,vortx_b,vorty_b,vortz_b)
!
!############################################################################
USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE FR_perf_mod
USE MPI
USE var, only : udrms,vdrms,wdrms,vortxdrms,vortydrms,vortzdrms,&
     uturbrms,vturbrms,wturbrms,vortxturbrms,vortyturbrms,vortzturbrms,&
     udvdrey,utvtrey

implicit none

integer :: i,j,k,code
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tu_tavg,dtu_tavg
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tv_tavg,dtv_tavg
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tw_tavg,dtw_tavg

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tvortx_tavg,dtvortx_tavg
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tvorty_tavg,dtvorty_tavg
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tvortz_tavg,dtvortz_tavg

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: utavg,utavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vtavg,vtavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wtavg,wtavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u2tavg,v2tavg,w2tavg

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortxtavg,vortxtavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortytavg,vortytavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortztavg,vortztavgsqrd
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortx2tavg,vorty2tavg,vortz2tavg

real(mytype),dimension(ysize(2)) :: u_tavgsqrd_spavg,u_tavgsqrd_spavg_red
real(mytype),dimension(ysize(2)) :: v_tavgsqrd_spavg,v_tavgsqrd_spavg_red
real(mytype),dimension(ysize(2)) :: w_tavgsqrd_spavg,w_tavgsqrd_spavg_red

real(mytype),dimension(ysize(2)) :: vortx_tavgsqrd_spavg,vortx_tavgsqrd_spavg_red
real(mytype),dimension(ysize(2)) :: vorty_tavgsqrd_spavg,vorty_tavgsqrd_spavg_red
real(mytype),dimension(ysize(2)) :: vortz_tavgsqrd_spavg,vortz_tavgsqrd_spavg_red

real(mytype),dimension(ysize(2)) :: u2_tavg_spavg,u2_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: v2_tavg_spavg,v2_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: w2_tavg_spavg,w2_tavg_spavg_red

real(mytype),dimension(ysize(2)) :: vortx2_tavg_spavg,vortx2_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: vorty2_tavg_spavg,vorty2_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: vortz2_tavg_spavg,vortz2_tavg_spavg_red
                                    
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u_to_avg,v_to_avg,w_to_avg
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: vortx_to_avg,vorty_to_avg,vortz_to_avg

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu_tavgsqrd,tv_tavgsqrd,tw_tavgsqrd
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tvortx_tavgsqrd,tvorty_tavgsqrd,tvortz_tavgsqrd

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu2_tavg,tv2_tavg,tw2_tavg
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tvortx2_tavg,tvorty2_tavg,tvortz2_tavg

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u_vtavg,uvtavg
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu_v_tavg,tuv_tavg
real(mytype),dimension(ysize(2)) :: u_v_tavg_spavg,u_v_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: uv_tavg_spavg,uv_tavg_spavg_red
real(mytype),dimension(ysize(2)) :: u_b,v_b,w_b,vortx_b,vorty_b,vortz_b

real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
integer :: ijk,nvect1
character(len=128) :: filename

call perfon('stat_disk')

nvect1=xsize(1)*xsize(2)*xsize(3)
call transpose_y_to_x(u_to_avg,tu_tavg)
call transpose_y_to_x(v_to_avg,tv_tavg)
call transpose_y_to_x(w_to_avg,tw_tavg)

call transpose_y_to_x(vortx_to_avg,tvortx_tavg)
call transpose_y_to_x(vorty_to_avg,tvorty_tavg)
call transpose_y_to_x(vortz_to_avg,tvortz_tavg)

dtu_tavg=0.0;dtv_tavg=0.0;dtw_tavg=0.0
dtvortx_tavg=0.0;dtvorty_tavg=0.0;dtvortz_tavg=0.0

do ijk=1,nvect1
   dtu_tavg(ijk,1,1)=tu_tavg(ijk,1,1)/(dble(ilast-ifirst+1))
   dtv_tavg(ijk,1,1)=tv_tavg(ijk,1,1)/(dble(ilast-ifirst+1))
   dtw_tavg(ijk,1,1)=tw_tavg(ijk,1,1)/(dble(ilast-ifirst+1))

   dtvortx_tavg(ijk,1,1)=tvortx_tavg(ijk,1,1)/(dble(ilast-ifirst+1))
   dtvorty_tavg(ijk,1,1)=tvorty_tavg(ijk,1,1)/(dble(ilast-ifirst+1))
   dtvortz_tavg(ijk,1,1)=tvortz_tavg(ijk,1,1)/(dble(ilast-ifirst+1))
enddo

!!!!!!!!!!!!
! velocity !
!!!!!!!!!!!!
uvisu=0.0
call fine_to_coarseV(1,dtu_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/velx.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

uvisu=0.0
call fine_to_coarseV(1,dtv_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/vely.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

uvisu=0.0
call fine_to_coarseV(1,dtw_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/velz.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

!!!!!!!!!!!!!
! vorticity !
!!!!!!!!!!!!!
uvisu=0.0
call fine_to_coarseV(1,dtvortx_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/vortx.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

uvisu=0.0
call fine_to_coarseV(1,dtvorty_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/vorty.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

uvisu=0.0
call fine_to_coarseV(1,dtvortz_tavg,uvisu)
write(filename, "(a)") trim(outdir)//"/vortz.dat"
call decomp_2d_write_one(1,uvisu,filename,2)

!###########################################################################################
utavgsqrd=(utavg/(dble(ilast-ifirst+1)))**2.0
vtavgsqrd=(vtavg/(dble(ilast-ifirst+1)))**2.0
wtavgsqrd=(wtavg/(dble(ilast-ifirst+1)))**2.0

vortxtavgsqrd=(vortxtavg/(dble(ilast-ifirst+1)))**2.0
vortytavgsqrd=(vortytavg/(dble(ilast-ifirst+1)))**2.0
vortztavgsqrd=(vortztavg/(dble(ilast-ifirst+1)))**2.0

u2tavg=u2tavg/(dble(ilast-ifirst+1))
v2tavg=v2tavg/(dble(ilast-ifirst+1))
w2tavg=w2tavg/(dble(ilast-ifirst+1))

vortx2tavg=vortx2tavg/(dble(ilast-ifirst+1))
vorty2tavg=vorty2tavg/(dble(ilast-ifirst+1))
vortz2tavg=vortz2tavg/(dble(ilast-ifirst+1))

u_vtavg=(utavg/(dble(ilast-ifirst+1)))*(vtavg/(dble(ilast-ifirst+1)))
uvtavg=uvtavg/(dble(ilast-ifirst+1))

call transpose_x_to_y(utavgsqrd,tu_tavgsqrd)
call transpose_x_to_y(vtavgsqrd,tv_tavgsqrd)
call transpose_x_to_y(wtavgsqrd,tw_tavgsqrd)

call transpose_x_to_y(vortxtavgsqrd,tvortx_tavgsqrd)
call transpose_x_to_y(vortytavgsqrd,tvorty_tavgsqrd)
call transpose_x_to_y(vortztavgsqrd,tvortz_tavgsqrd)

call transpose_x_to_y(u2tavg,tu2_tavg)
call transpose_x_to_y(v2tavg,tv2_tavg)
call transpose_x_to_y(w2tavg,tw2_tavg)

call transpose_x_to_y(vortx2tavg,tvortx2_tavg)
call transpose_x_to_y(vorty2tavg,tvorty2_tavg)
call transpose_x_to_y(vortz2tavg,tvortz2_tavg)

call transpose_x_to_y(u_vtavg,tu_v_tavg)
call transpose_x_to_y(uvtavg,tuv_tavg)

do j=1,ny
   u_tavgsqrd_spavg(j)=0.0
   v_tavgsqrd_spavg(j)=0.0
   w_tavgsqrd_spavg(j)=0.0

   vortx_tavgsqrd_spavg(j)=0.0
   vorty_tavgsqrd_spavg(j)=0.0
   vortz_tavgsqrd_spavg(j)=0.0

   u2_tavg_spavg(j)=0.0
   v2_tavg_spavg(j)=0.0
   w2_tavg_spavg(j)=0.0

   vortx2_tavg_spavg(j)=0.0
   vorty2_tavg_spavg(j)=0.0
   vortz2_tavg_spavg(j)=0.0

   u_v_tavg_spavg(j)=0.0
   uv_tavg_spavg(j)=0.0
end do

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         u_tavgsqrd_spavg(j)=u_tavgsqrd_spavg(j)+tu_tavgsqrd(i,j,k)
         v_tavgsqrd_spavg(j)=v_tavgsqrd_spavg(j)+tv_tavgsqrd(i,j,k)
         w_tavgsqrd_spavg(j)=w_tavgsqrd_spavg(j)+tw_tavgsqrd(i,j,k)

         vortx_tavgsqrd_spavg(j)=vortx_tavgsqrd_spavg(j)+tvortx_tavgsqrd(i,j,k)
         vorty_tavgsqrd_spavg(j)=vorty_tavgsqrd_spavg(j)+tvorty_tavgsqrd(i,j,k)
         vortz_tavgsqrd_spavg(j)=vortz_tavgsqrd_spavg(j)+tvortz_tavgsqrd(i,j,k)

         u2_tavg_spavg(j)=u2_tavg_spavg(j)+tu2_tavg(i,j,k)
         v2_tavg_spavg(j)=v2_tavg_spavg(j)+tv2_tavg(i,j,k)
         w2_tavg_spavg(j)=w2_tavg_spavg(j)+tw2_tavg(i,j,k)

         vortx2_tavg_spavg(j)=vortx2_tavg_spavg(j)+tvortx2_tavg(i,j,k)
         vorty2_tavg_spavg(j)=vorty2_tavg_spavg(j)+tvorty2_tavg(i,j,k)
         vortz2_tavg_spavg(j)=vortz2_tavg_spavg(j)+tvortz2_tavg(i,j,k)

         u_v_tavg_spavg(j)=u_v_tavg_spavg(j)+tu_v_tavg(i,j,k)
         uv_tavg_spavg(j)=uv_tavg_spavg(j)+tuv_tavg(i,j,k)
      end do
   end do
end do

do j=1,ny
   u_tavgsqrd_spavg(j)=u_tavgsqrd_spavg(j)/ysize(1)/ysize(3)
   v_tavgsqrd_spavg(j)=v_tavgsqrd_spavg(j)/ysize(1)/ysize(3)
   w_tavgsqrd_spavg(j)=w_tavgsqrd_spavg(j)/ysize(1)/ysize(3)

   vortx_tavgsqrd_spavg(j)=vortx_tavgsqrd_spavg(j)/ysize(1)/ysize(3)
   vorty_tavgsqrd_spavg(j)=vorty_tavgsqrd_spavg(j)/ysize(1)/ysize(3)
   vortz_tavgsqrd_spavg(j)=vortz_tavgsqrd_spavg(j)/ysize(1)/ysize(3)

   u2_tavg_spavg(j)=u2_tavg_spavg(j)/ysize(1)/ysize(3)
   v2_tavg_spavg(j)=v2_tavg_spavg(j)/ysize(1)/ysize(3)
   w2_tavg_spavg(j)=w2_tavg_spavg(j)/ysize(1)/ysize(3)

   vortx2_tavg_spavg(j)=vortx2_tavg_spavg(j)/ysize(1)/ysize(3)
   vorty2_tavg_spavg(j)=vorty2_tavg_spavg(j)/ysize(1)/ysize(3)
   vortz2_tavg_spavg(j)=vortz2_tavg_spavg(j)/ysize(1)/ysize(3)

   u_v_tavg_spavg(j)=u_v_tavg_spavg(j)/ysize(1)/ysize(3)
   uv_tavg_spavg(j)=uv_tavg_spavg(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(u_tavgsqrd_spavg,u_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(v_tavgsqrd_spavg,v_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(w_tavgsqrd_spavg,w_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

call MPI_ALLREDUCE(vortx_tavgsqrd_spavg,vortx_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vorty_tavgsqrd_spavg,vorty_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortz_tavgsqrd_spavg,vortz_tavgsqrd_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

call MPI_ALLREDUCE(u2_tavg_spavg,u2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(v2_tavg_spavg,v2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(w2_tavg_spavg,w2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

call MPI_ALLREDUCE(vortx2_tavg_spavg,vortx2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vorty2_tavg_spavg,vorty2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(vortz2_tavg_spavg,vortz2_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

call MPI_ALLREDUCE(u_v_tavg_spavg,u_v_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(uv_tavg_spavg,uv_tavg_spavg_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   u_tavgsqrd_spavg_red(j)=u_tavgsqrd_spavg_red(j)/nproc
   v_tavgsqrd_spavg_red(j)=v_tavgsqrd_spavg_red(j)/nproc
   w_tavgsqrd_spavg_red(j)=w_tavgsqrd_spavg_red(j)/nproc

   vortx_tavgsqrd_spavg_red(j)=vortx_tavgsqrd_spavg_red(j)/nproc
   vorty_tavgsqrd_spavg_red(j)=vorty_tavgsqrd_spavg_red(j)/nproc
   vortz_tavgsqrd_spavg_red(j)=vortz_tavgsqrd_spavg_red(j)/nproc

   u2_tavg_spavg_red(j)=u2_tavg_spavg_red(j)/nproc
   v2_tavg_spavg_red(j)=v2_tavg_spavg_red(j)/nproc
   w2_tavg_spavg_red(j)=w2_tavg_spavg_red(j)/nproc

   vortx2_tavg_spavg_red(j)=vortx2_tavg_spavg_red(j)/nproc
   vorty2_tavg_spavg_red(j)=vorty2_tavg_spavg_red(j)/nproc
   vortz2_tavg_spavg_red(j)=vortz2_tavg_spavg_red(j)/nproc

   u_v_tavg_spavg_red(j)=u_v_tavg_spavg_red(j)/nproc
   uv_tavg_spavg_red(j)=uv_tavg_spavg_red(j)/nproc
end do

do j=1,ny
   udrms(j)=sqrt(u_tavgsqrd_spavg_red(j)-u_b(j)*u_b(j))
   vdrms(j)=sqrt(v_tavgsqrd_spavg_red(j)-v_b(j)*v_b(j))
   wdrms(j)=sqrt(w_tavgsqrd_spavg_red(j)-w_b(j)*w_b(j))

   vortxdrms(j)=sqrt(vortx_tavgsqrd_spavg_red(j)-vortx_b(j)*vortx_b(j))
   vortydrms(j)=sqrt(vorty_tavgsqrd_spavg_red(j)-vorty_b(j)*vorty_b(j))
   vortzdrms(j)=sqrt(vortz_tavgsqrd_spavg_red(j)-vortz_b(j)*vortz_b(j))

   uturbrms(j)=sqrt(u2_tavg_spavg_red(j)-u_tavgsqrd_spavg_red(j))
   vturbrms(j)=sqrt(v2_tavg_spavg_red(j)-v_tavgsqrd_spavg_red(j))
   wturbrms(j)=sqrt(w2_tavg_spavg_red(j)-w_tavgsqrd_spavg_red(j))

   vortxturbrms(j)=sqrt(vortx2_tavg_spavg_red(j)-vortx_tavgsqrd_spavg_red(j))
   vortyturbrms(j)=sqrt(vorty2_tavg_spavg_red(j)-vorty_tavgsqrd_spavg_red(j))
   vortzturbrms(j)=sqrt(vortz2_tavg_spavg_red(j)-vortz_tavgsqrd_spavg_red(j))

   udvdrey(j)=u_v_tavg_spavg_red(j)-u_b(j)*v_b(j)
   utvtrey(j)=uv_tavg_spavg_red(j)-u_v_tavg_spavg_red(j)
end do

call perfoff
return

end subroutine STAT_DISK

!############################################################################
!
  subroutine WRITE_STATISTICS(Subar,Svbar,Swbar,Svortxbar,Svortybar,Svortzbar,&
     Suv,Suw,Svw,Surms,Svrms,Swrms,Svortxrms,Svortyrms,Svortzrms,&
     Sdudy_lw_tot,Sdudy_uw_tot,Sdudy_tot,cf_tot,dissp1,dubardyt,turbconv,&
     viscdiff,dirout)
!
!############################################################################

USE param
USE variables
USE MPI
USE FR_perf_mod
USE var, only : udrms,vdrms,wdrms,vortxdrms,vortydrms,vortzdrms,&
     uturbrms,vturbrms,wturbrms,vortxturbrms,vortyturbrms,vortzturbrms,&
     udvdrey,utvtrey

implicit none

real(mytype), dimension(ny) :: Subar,Svbar,Swbar,Svortxbar,Svortybar,Svortzbar
real(mytype), dimension(ny) :: Surms,Svrms,Swrms,Svortxrms,Svortyrms,Svortzrms
real(mytype), dimension(ny) :: Suv,Suw,Svw
real(mytype), dimension(ny) :: dissp1,dubardyt,turbconv,viscdiff
real(mytype), dimension(ny) :: convtrspx123,convtrspy123,convtrspz123
real(mytype) :: Sdudy_lw_tot,Sdudy_uw_tot,Sdudy_tot,cf_tot
integer :: j
character(len=128) :: dirout

call perfon('write_statistics')

open(120,file=trim(dirout)//'/velbar.txt',status='unknown')
open(130,file=trim(dirout)//'/vortbar.txt',status='unknown')
open(140,file=trim(dirout)//'/velrms.txt',status='unknown')
open(150,file=trim(dirout)//'/vortrms.txt',status='unknown')
open(160,file=trim(dirout)//'/rey.txt',status='unknown')
open(170,file=trim(dirout)//'/veldrms.txt',status='unknown')
open(180,file=trim(dirout)//'/veltrms.txt',status='unknown')
open(190,file=trim(dirout)//'/rey_dt.txt',status='unknown')
open(200,file=trim(dirout)//'/vortdrms.txt',status='unknown')
open(210,file=trim(dirout)//'/vorttrms.txt',status='unknown')
open(220,file=trim(dirout)//'/balance.txt',status='unknown')

open(300,file=trim(dirout)//'/drg_avg.txt',status='unknown')

do j=1,ny 
   write(120,*) yp(j),Subar(j),Svbar(j),Swbar(j)
   write(130,*) yp(j),Svortxbar(j),Svortybar(j),Svortzbar(j)
   write(140,*) yp(j),Surms(j),Svrms(j),Swrms(j)
   write(150,*) yp(j),Svortxrms(j),Svortyrms(j),Svortzrms(j)
   write(160,*) yp(j),Suv(j),Suw(j),Svw(j)
   write(170,*) yp(j),udrms(j),vdrms(j),wdrms(j)
   write(180,*) yp(j),uturbrms(j),vturbrms(j),wturbrms(j)
   write(190,*) yp(j),udvdrey(j),utvtrey(j)
   write(200,*) yp(j),vortxdrms(j),vortydrms(j),vortzdrms(j)
   write(210,*) yp(j),vortxturbrms(j),vortyturbrms(j),vortzturbrms(j)
   ! dissipation,production,turb_conv,visc_diff
   write(220,*) yp(j),-(dissp1(j)-dubardyt(j)*dubardyt(j))*(re**2/re_tau**4),-Suv(j)*dubardyt(j)*(re**3/re_tau**4),&
        turbconv(j)*(re**3/re_tau**4),viscdiff(j)*(re/re_tau**2)
end do

write(300,*) Sdudy_lw_tot,Sdudy_uw_tot,&
     Sdudy_tot,cf_tot,sqrt(re*Sdudy_tot)

close(120);close(130);close(140);close(150);close(160);close(170)
close(180);close(190);close(200);close(210);close(220);close(300)

call perfoff
end subroutine WRITE_STATISTICS

subroutine time_div()

USE param
USE var
USE FR_perf_mod

implicit none

integer :: j

call perfon('time_div')

do j=1,ny
   S_ubar(j)=S_ubar(j)/(dble(ilast-ifirst+1))
   S_vbar(j)=S_vbar(j)/(dble(ilast-ifirst+1))
   S_wbar(j)=S_wbar(j)/(dble(ilast-ifirst+1))
   S_vortxbar(j)=S_vortxbar(j)/(dble(ilast-ifirst+1))
   S_vortybar(j)=S_vortybar(j)/(dble(ilast-ifirst+1))
   S_vortzbar(j)=S_vortzbar(j)/(dble(ilast-ifirst+1))
   S_urms(j)=S_urms(j)/(dble(ilast-ifirst+1))
   S_urms(j)=sqrt(S_urms(j))
   S_vrms(j)=S_vrms(j)/(dble(ilast-ifirst+1))
   S_vrms(j)=sqrt(S_vrms(j))
   S_wrms(j)=S_wrms(j)/(dble(ilast-ifirst+1))
   S_wrms(j)=sqrt(S_wrms(j))
   S_uv(j)=S_uv(j)/(dble(ilast-ifirst+1))
   S_uw(j)=S_uw(j)/(dble(ilast-ifirst+1))
   S_vw(j)=S_vw(j)/(dble(ilast-ifirst+1))
   S_vortxrms(j)=S_vortxrms(j)/(dble(ilast-ifirst+1))
   S_vortxrms(j)=sqrt(S_vortxrms(j))
   S_vortyrms(j)=S_vortyrms(j)/(dble(ilast-ifirst+1))
   S_vortyrms(j)=sqrt(S_vortyrms(j))
   S_vortzrms(j)=S_vortzrms(j)/(dble(ilast-ifirst+1))
   S_vortzrms(j)=sqrt(S_vortzrms(j))

   dissp1_tavg(j)=dissp1_tavg(j)/(dble(ilast-ifirst+1))
   dubardy_tavg(j)=dubardy_tavg(j)/(dble(ilast-ifirst+1))
   turb_conv(j)=turb_conv(j)/(dble(ilast-ifirst+1))
   visc_diff(j)=visc_diff(j)/(dble(ilast-ifirst+1))
end do

S_dudy_lw_tot=S_dudy_lw_tot/(dble(ilast-ifirst+1))
S_dudy_uw_tot=S_dudy_uw_tot/(dble(ilast-ifirst+1))
S_dudy_tot=S_dudy_tot/(dble(ilast-ifirst+1))
c_f_tot=c_f_tot/(dble(ilast-ifirst+1))

call perfoff
end subroutine time_div

subroutine space_avg(tu_1,u_1_xz_red)

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI
USE FR_perf_mod

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu_1
real(mytype),dimension(ysize(2)) :: u_1_xz
real(mytype),dimension(ysize(2)), intent(out) :: u_1_xz_red
integer :: i,j,k,code

call perfon('space_avg')

do j=1,ny
   u_1_xz(j)=0.0
end do

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         u_1_xz(j)=u_1_xz(j)+tu_1(i,j,k)
      end do
   end do
end do

do j=1,ny
   u_1_xz(j)=u_1_xz(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(u_1_xz,u_1_xz_red,ny,&
     real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   u_1_xz_red(j)=u_1_xz_red(j)/nproc
end do

call perfoff
end subroutine space_avg

subroutine rms_comp(tu_1,u_1_bar,u_1_rms_red)

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI
USE FR_perf_mod

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu_1
real(mytype),dimension(ysize(2)) :: u_1_bar,u_1_rms
real(mytype),dimension(ysize(2)), intent(out) :: u_1_rms_red
integer :: i,j,k,code

call perfon('rms_comp')

do j=1,ny
   u_1_rms(j)=0.0
end do

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         u_1_rms(j)=u_1_rms(j)+(tu_1(i,j,k)-u_1_bar(j))**2
      end do
   end do
end do

do j=1,ny
   u_1_rms(j)=u_1_rms(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(u_1_rms,u_1_rms_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   u_1_rms_red(j)=u_1_rms_red(j)/nproc
end do

call perfoff
end subroutine rms_comp

subroutine rey_comp(tu_1,tu_2,u_1_bar,u_2_bar,u_12_rey_red)

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI
USE FR_perf_mod

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tu_1,tu_2
real(mytype),dimension(ysize(2)) :: u_1_bar,u_2_bar,u_12_rey
real(mytype),dimension(ysize(2)), intent(out) :: u_12_rey_red
integer :: i,j,k,code

call perfon('rey_comp')

do j=1,ny
   u_12_rey(j)=0.0
end do

do j=1,ny
   do k=1,ysize(3)
      do i=1,ysize(1)
         u_12_rey(j)=u_12_rey(j)+(tu_1(i,j,k)-u_1_bar(j))*(tu_2(i,j,k)-u_2_bar(j))
      end do
   end do
end do

do j=1,ny
   u_12_rey(j)=u_12_rey(j)/ysize(1)/ysize(3)
end do

call MPI_ALLREDUCE(u_12_rey,u_12_rey_red,ny,real_type,MPI_SUM,MPI_COMM_WORLD,code)

do j=1,ny
   u_12_rey_red(j)=u_12_rey_red(j)/nproc
end do

call perfoff
end subroutine rey_comp

subroutine read_file(ff,fh,dtd,nstr,it,disp,tab)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI
  USE FR_perf_mod
  
  implicit none

  integer :: nstr,it,fh,ierror
  integer (kind=MPI_OFFSET_KIND) :: disp
  character(len=200) :: ff,dtd
  character (len=2), dimension(3) :: arr_str 
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)), &
       intent(out) :: tab(xsize(1),xsize(2),xsize(3))

  call perfon('read_file')

  arr_str=(/'ux','uy','uz'/)

  write(ff, "(a,i6.6,a)") trim(dtd)//"/"//trim(arr_str(nstr)), it ,".dat"
  call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(ff), &
       MPI_MODE_RDONLY, MPI_INFO_NULL, &
       fh, ierror)
  call decomp_2d_read_var(fh,disp,1,tab)
  call MPI_FILE_CLOSE(fh,ierror)

  call perfoff
end subroutine read_file
