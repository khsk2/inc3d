!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module var

use decomp_2d
USE variables
USE param

! define all major arrays here

real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3,po3,dv3,pp3
real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
real(mytype), save, allocatable, dimension(:,:,:) :: phi1, phi2, phi3
real(mytype), save, allocatable, dimension(:,:,:) :: gx1, gy1, gz1, hx1, hy1, hz1, phis1,phiss1
real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
real(mytype), save, allocatable, dimension(:,:,:) :: ep1

!arrays for statistic collection
real(mytype), save, allocatable, dimension(:,:,:) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype), save, allocatable, dimension(:,:,:) :: phimean, phiphimean

!arrays for visualization
real(mytype), save, allocatable, dimension(:,:,:) :: uvisu

! define all work arrays here
real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
     te1,tf1,tg1,th1,ti1,di1
real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
     te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
     te3,tf3,tg3,th3,ti3,di3
! for reading files
real(mytype), save, allocatable, dimension(:,:,:) :: tabx,taby,tabz
real(mytype), save, allocatable, dimension(:,:,:) :: tabvox,tabvoy,tabvoz
real(mytype), save, allocatable, dimension(:,:,:) :: tabvox_w_vort,tabvoy_w_vort,tabvoz_w_vort

real(mytype), save, allocatable, dimension(:,:,:) :: ud,vd,wd
real(mytype), save, allocatable, dimension(:,:,:) :: vortxd,vortyd,vortzd

real(mytype), save, allocatable, dimension(:,:,:) :: uturb,vturb,wturb
real(mytype), save, allocatable, dimension(:,:,:) :: vortxturb,vortyturb,vortzturb

real(mytype), save, allocatable, dimension(:,:,:) :: u_tavg,v_tavg,w_tavg
real(mytype), save, allocatable, dimension(:,:,:) :: vortx_tavg,vorty_tavg,vortz_tavg

real(mytype), save, allocatable, dimension(:,:,:) :: u2_tavg,v2_tavg,w2_tavg
real(mytype), save, allocatable, dimension(:,:,:) :: vortx2_tavg,vorty2_tavg,vortz2_tavg

real(mytype), save, allocatable, dimension(:,:,:) :: uv_tavg
real(mytype), save, allocatable, dimension(:,:,:) :: tabdudx,tabdudy,tabdudz
real(mytype), save, allocatable, dimension(:,:,:) :: tabdvdx,tabdvdy,tabdvdz
real(mytype), save, allocatable, dimension(:,:,:) :: tabdwdx,tabdwdy,tabdwdz

real(mytype), save, allocatable, dimension(:) :: dissp1_tavg,dubardy_tavg,&
     turb_conv,visc_diff

integer, save :: nxmsize, nymsize, nzmsize 

real(mytype) :: S_dudy_lw_tot,S_dudy_uw_tot,S_dudy_tot,c_f_tot
real(mytype), allocatable, dimension(:) :: S_ubar,S_vbar,S_wbar,&
     S_vortxbar,S_vortybar,S_vortzbar
real(mytype), allocatable, dimension(:) :: S_urms,S_vrms,S_wrms,&
     S_vortxrms,S_vortyrms,S_vortzrms
real(mytype), allocatable, dimension(:) :: S_uv,S_uw,S_vw
real(mytype), allocatable, dimension(:) :: udrms,vdrms,wdrms
real(mytype), allocatable, dimension(:) :: vortxdrms,vortydrms,vortzdrms
real(mytype), allocatable, dimension(:) :: uturbrms,vturbrms,wturbrms
real(mytype), allocatable, dimension(:) :: vortxturbrms,vortyturbrms,vortzturbrms
real(mytype), allocatable, dimension(:) :: udvdrey,utvtrey
real(mytype), allocatable, dimension(:) :: dubardy

contains

  
  subroutine init_variables

    TYPE(DECOMP_INFO), save :: ph  ! decomposition object

    if (nclx==0) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly==0) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz==0) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)

!X PENCILS
    call alloc_x(ux1, opt_global=.true.)
    call alloc_x(uy1, opt_global=.true.)
#ifndef TWOD
    call alloc_x(uz1, opt_global=.true.)
    call alloc_x(pz1, opt_global=.true.)
#else
    allocate (uz1(1,1,1))
    allocate (pz1(1,1,1))
#endif
    call alloc_x(px1, opt_global=.true.)
    call alloc_x(py1, opt_global=.true.)
    call alloc_x(phi1, opt_global=.true.)
    call alloc_x(gx1);call alloc_x(gy1);call alloc_x(gz1);call alloc_x(phis1) 
    call alloc_x(hx1);call alloc_x(hy1);call alloc_x(hz1);call alloc_x(phiss1)
    call alloc_x(ta1);call alloc_x(tb1);call alloc_x(tc1)
    call alloc_x(td1);call alloc_x(te1);call alloc_x(tf1)
    call alloc_x(tg1);call alloc_x(th1);call alloc_x(ti1)
    call alloc_x(di1);call alloc_x(ep1)
    allocate(sx(xsize(2),xsize(3)),vx(xsize(2),xsize(3)))
    !inflow/ouflow 2d arrays
    allocate(bxx1(xsize(2),xsize(3)),bxy1(xsize(2),xsize(3)))
    allocate(bxz1(xsize(2),xsize(3)),bxxn(xsize(2),xsize(3)))
    allocate(bxyn(xsize(2),xsize(3)),bxzn(xsize(2),xsize(3)))
    allocate(bxo(xsize(2),xsize(3)),byo(xsize(2),xsize(3)))
    allocate(bzo(xsize(2),xsize(3)))
    allocate(byx1(xsize(1),xsize(3)),byy1(xsize(1),xsize(3)))
    allocate(byz1(xsize(1),xsize(3)),byxn(xsize(1),xsize(3)))
    allocate(byyn(xsize(1),xsize(3)),byzn(xsize(1),xsize(3)))   
    allocate(bzx1(xsize(1),xsize(2)),bzy1(xsize(1),xsize(2)))
    allocate(bzz1(xsize(1),xsize(2)),bzxn(xsize(1),xsize(2)))
    allocate(bzyn(xsize(1),xsize(2)),bzzn(xsize(1),xsize(2)))
    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))

    allocate(S_ubar(ny),S_vbar(ny),S_wbar(ny),S_vortxbar(ny),S_vortybar(ny),S_vortzbar(ny))
    allocate(S_urms(ny),S_vrms(ny),S_wrms(ny),S_vortxrms(ny),S_vortyrms(ny),S_vortzrms(ny))
    allocate(S_uv(ny),S_uw(ny),S_vw(ny))

    allocate(udrms(ny),vdrms(ny),wdrms(ny))
    allocate(vortxdrms(ny),vortydrms(ny),vortzdrms(ny))

    allocate(uturbrms(ny),vturbrms(ny),wturbrms(ny))
    allocate(vortxturbrms(ny),vortyturbrms(ny),vortzturbrms(ny))

    allocate(udvdrey(ny),utvtrey(ny))
    allocate(dissp1_tavg(ny),dubardy_tavg(ny))
    allocate(turb_conv(ny),visc_diff(ny))

!arrays for statistic collection!pay attention to the size!
    allocate (umean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uumean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (wwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (uwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (vwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    allocate (tmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))    
    if (iscalar==1) then
       allocate (phimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate (phiphimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    else
       allocate (phimean(1,1,1))
       allocate (phiphimean(1,1,1))
    endif

    allocate(tabx(xsize(1),xsize(2),xsize(3)))
    allocate(taby(xsize(1),xsize(2),xsize(3)))
    allocate(tabz(xsize(1),xsize(2),xsize(3)))

    allocate(tabvox(xsize(1),xsize(2),xsize(3)))
    allocate(tabvoy(xsize(1),xsize(2),xsize(3)))
    allocate(tabvoz(xsize(1),xsize(2),xsize(3)))

    allocate(tabvox_w_vort(xsize(1),xsize(2),xsize(3)))
    allocate(tabvoy_w_vort(xsize(1),xsize(2),xsize(3)))
    allocate(tabvoz_w_vort(xsize(1),xsize(2),xsize(3)))

    allocate(tabdudx(xsize(1),xsize(2),xsize(3)))
    allocate(tabdudy(xsize(1),xsize(2),xsize(3)))
    allocate(tabdudz(xsize(1),xsize(2),xsize(3)))

    allocate(tabdvdx(xsize(1),xsize(2),xsize(3)))
    allocate(tabdvdy(xsize(1),xsize(2),xsize(3)))
    allocate(tabdvdz(xsize(1),xsize(2),xsize(3)))

    allocate(tabdwdx(xsize(1),xsize(2),xsize(3)))
    allocate(tabdwdy(xsize(1),xsize(2),xsize(3)))
    allocate(tabdwdz(xsize(1),xsize(2),xsize(3)))

    allocate(ud(ysize(1),ysize(2),ysize(3)))
    allocate(vd(ysize(1),ysize(2),ysize(3)))
    allocate(wd(ysize(1),ysize(2),ysize(3)))

    allocate(vortxd(ysize(1),ysize(2),ysize(3)))
    allocate(vortyd(ysize(1),ysize(2),ysize(3)))
    allocate(vortzd(ysize(1),ysize(2),ysize(3)))

    allocate(uturb(ysize(1),ysize(2),ysize(3)))
    allocate(vturb(ysize(1),ysize(2),ysize(3)))
    allocate(wturb(ysize(1),ysize(2),ysize(3)))

    allocate(vortxturb(ysize(1),ysize(2),ysize(3)))
    allocate(vortyturb(ysize(1),ysize(2),ysize(3)))
    allocate(vortzturb(ysize(1),ysize(2),ysize(3)))

    allocate(u_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(v_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(w_tavg(xsize(1),xsize(2),xsize(3)))

    allocate(vortx_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(vorty_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(vortz_tavg(xsize(1),xsize(2),xsize(3)))

    allocate(u2_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(v2_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(w2_tavg(xsize(1),xsize(2),xsize(3)))

    allocate(vortx2_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(vorty2_tavg(xsize(1),xsize(2),xsize(3)))
    allocate(vortz2_tavg(xsize(1),xsize(2),xsize(3)))

    allocate(uv_tavg(xsize(1),xsize(2),xsize(3)))

!arrays for visualization!pay attention to the size!
    allocate (uvisu(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))

!Y PENCILS
    call alloc_y(ux2);call alloc_y(uy2);call alloc_y(uz2)
    call alloc_y(ta2);call alloc_y(tb2);call alloc_y(tc2)
    call alloc_y(td2);call alloc_y(te2);call alloc_y(tf2)
    call alloc_y(tg2);call alloc_y(th2);call alloc_y(ti2)
    call alloc_y(tj2)
    call alloc_y(di2);call alloc_y(phi2)
    allocate(sy(ysize(1),ysize(3)),vy(ysize(1),ysize(3)))
!Z PENCILS
    call alloc_z(ux3);call alloc_z(uy3);call alloc_z(uz3)
    call alloc_z(ta3);call alloc_z(tb3);call alloc_z(tc3)
    call alloc_z(td3);call alloc_z(te3);call alloc_z(tf3)
    call alloc_z(tg3);call alloc_z(th3);call alloc_z(ti3)
    call alloc_z(di3);call alloc_z(phi3)
    allocate(sz(zsize(1),zsize(2)),vz(zsize(1),zsize(2)))

 ! if all periodic
 !   allocate (pp3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
 !   allocate (dv3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
 !   allocate (po3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    call alloc_z(pp3,ph,.true.)
    call alloc_z(dv3,ph,.true.)
    call alloc_z(po3,ph,.true.)

    return
  end subroutine init_variables

end module var

