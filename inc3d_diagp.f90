PROGRAM incompact3d

  USE decomp_2d
  use decomp_2d_io
  USE variables
  USE param
  USE var
  USE FR_perf_mod
  USE MPI
  USE derivX
  USE derivZ

  implicit none

  integer :: code,nlock,j,bcx,bcy,bcz
  integer :: fhx,fhy,fhz
  double precision :: t1,t2
  character(len=200) :: filename
  integer (kind=MPI_OFFSET_KIND) :: dispx,dispy,dispz,status
  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  CALL MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)
  call parameter()

  call perfinit
  call perfon('Incompact3d')
  call init_variables

  call schemes()

  if (nclx==0) then
     bcx=0
  else
     bcx=1
  endif
  if (ncly==0) then
     bcy=0
  else
     bcy=1
  endif
  if (nclz==0) then
     bcz=0
  else
     bcz=1
  endif

  call decomp_info_init(nxm,nym,nzm,phG)

  if(nrank.eq.0) then
     S_dudy_lw_tot=0.0;S_dudy_uw_tot=0.0;S_dudy_tot=0.0;c_f_tot=0.0
     do j=1,ny
        S_ubar(j)=0.0;S_vbar(j)=0.0;S_wbar(j)=0.0
        S_uv(j)=0.0;S_uw(j)=0.0;S_vw(j)=0.0
        S_urms(j)=0.0;S_vrms(j)=0.0;S_wrms(j)=0.0
        S_vortxbar(j)=0.0;S_vortybar(j)=0.0;S_vortzbar(j)=0.0
        S_vortxrms(j)=0.0;S_vortyrms(j)=0.0;S_vortzrms(j)=0.0

        dissp1_tavg(j)=0.0;dubardy_tavg(j)=0.0
        turb_conv(j)=0.0;visc_diff(j)=0.0
     end do
  end if
  ud=0.0;vd=0.0;wd=0.0
  vortxd=0.0;vortyd=0.0;vortzd=0.0

  u_tavg=0.0;v_tavg=0.0;w_tavg=0.0
  vortx_tavg=0.0;vorty_tavg=0.0;vortz_tavg=0.0

  u2_tavg=0.0;v2_tavg=0.0;w2_tavg=0.0
  vortx2_tavg=0.0;vorty2_tavg=0.0;vortz2_tavg=0.0

  t1 = MPI_WTIME()

  do itime=ifirst,ilast

     dispx = 0_MPI_OFFSET_KIND
     dispy = 0_MPI_OFFSET_KIND
     dispz = 0_MPI_OFFSET_KIND

     call read_file(filename,fhx,datdir,1,itime,dispx,tabx)
     call read_file(filename,fhy,datdir,2,itime,dispy,taby)
     call read_file(filename,fhz,datdir,3,itime,dispz,tabz)
  
     call CALC_VORT_GRADVEL(tabx,taby,tabz,ta1,tb1,tc1,td1,te1,tf1,tg1,&
          th1,ti1,di1,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
          ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,&
          tabvox_w_vort,tabvoy_w_vort,tabvoz_w_vort,&
          tabdudx,tabdudy,tabdudz,tabdvdx,tabdvdy,tabdvdz,&
          tabdwdx,tabdwdy,tabdwdz)
     call STATISTIC(tabx,taby,tabz,tabvox_w_vort,tabvoy_w_vort,tabvoz_w_vort,&
          tabdudx,tabdudy,tabdudz,tabdvdx,tabdvdy,&
          tabdvdz,tabdwdx,tabdwdy,tabdwdz,outdir)

     call DISSIPATION(tabx,tabdudx,tabdudy,tabdudz,tabdvdx,&
          tabdvdy,tabdvdz,tabdwdx,tabdwdy,tabdwdz)

     call TURB_CNV(tabx,taby,tabz)
     call VISC_DFF(tabx,taby,tabz)
  enddo

  if(nrank.eq.0) then
     call time_div()
  end if

  call STAT_DISK(ud,vd,wd,vortxd,vortyd,vortzd,u_tavg,v_tavg,w_tavg,&
       vortx_tavg,vorty_tavg,vortz_tavg,&
       u2_tavg,v2_tavg,w2_tavg,vortx2_tavg,vorty2_tavg,vortz2_tavg,&
       uv_tavg,S_ubar,S_vbar,S_wbar,S_vortxbar,S_vortybar,S_vortzbar)

  if(nrank.eq.0) then
     call WRITE_STATISTICS(S_ubar,S_vbar,S_wbar,S_vortxbar,S_vortybar,S_vortzbar,&
          S_uv,S_uw,S_vw,S_urms,S_vrms,S_wrms,S_vortxrms,S_vortyrms,S_vortzrms,&
          S_dudy_lw_tot,S_dudy_uw_tot,S_dudy_tot,c_f_tot,dissp1_tavg,dubardy_tavg,&
          turb_conv,visc_diff,outdir)
  end if

  call perfoff
  if(nrank.eq.0) then
     call perfout('Incompact3d')
  end if

  t2=MPI_WTIME()-t1
  call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
       MPI_COMM_WORLD,code)
  if (nrank.eq.0) print *,'time per time_step: ', &
       t1/float(nproc)/(ilast-ifirst+1),' seconds'
  if (nrank.eq.0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
  if (nrank.eq.0) print *,'Mapping p_row*p_col=',p_row,p_col
  if (nrank.eq.0) print *,'Total runtime=',t1/float(nproc)

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
