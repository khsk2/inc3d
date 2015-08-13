#ifdef WITHFUTILS
       if(write_h5) then
          call creatf(trim(diagdir)//'/nrg'//trim(file_extension)//'.h5', &
               fidnrg_h5, "phase space averages", 'd')
          call creatg(fidnrg_h5, '/nrg')
          allocate(hbufnrg_h5(0:n_spec-1))
          do n=0, n_spec-1
             call creatg(fidnrg_h5, '/nrg'//trim(spec(n)%name))
             call htable_init(hbufnrg_h5(n), bufsize)
             call set_htable_fileid(hbufnrg_h5(n), fidnrg_h5, '/nrg'//trim(spec(n)%name))
          end do
       end if
#endif
    endif
#ifdef WITHFUTILS
    integer:: ng, isnap_nrg
#endif
       
#ifdef WITHFUTILS
       if ((write_h5).and.(mype).eq.0) then
          isnap_nrg = itime/istep_nrg + 1
          call attach(fidnrg_h5, "/nrg", "n_steps", isnap_nrg)    
          do n=ln1,ln2
             ng = pes*ln0 + n-ln1
             call add_record(hbufnrg_h5(ng), "time", "simulation time", time)
             do o=1,8
                call add_record(hbufnrg_h5(ng), trim(nrg_label1(o)), &
                     &trim(nrg_label2(o)),rat(o,n))
             enddo
             call htable_endstep(hbufnrg_h5(ng))
             call htable_hdf5_flush(hbufnrg_h5(ng))
          end do
          call flushh5(fidnrg_h5)
          
       end if
#endif
#ifdef WITHFUTILS
         if(write_h5) then
            call closef(fidnrg_h5)
            deallocate(hbufnrg_h5)
         end if
#endif
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank, o
#endif

#ifdef WITHFUTILS
    if (write_h5.and.((my_pev+my_pew+my_pespec).eq.0)) then
       isnap_field = 0
       call creatf(trim(diagdir)//'/field'//trim(file_extension)//'.h5', &
            fidfield_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidfield_h5, '/field')
       rank = 0
       call creatd(fidfield_h5, rank, dims, "/field/time", "time")
       
       do o=1,n_fields
          call creatg(fidfield_h5, '/field/'//trim(field_label(o)))
       enddo
    end if
#endif

#ifdef WITHFUTILS
    character(len=128) :: dset_name_field
#endif

#ifdef WITHFUTILS
      if(write_h5) then
         call append(fidfield_h5, "/field/time", time)
         call attach(fidfield_h5, "/field/time", "n_steps", isnap_field+1)
         
         do o=1,n_fields
            write(dset_name_field, "(A, '/', i10.10)") "/field/"//trim(field_label(o)), isnap_field
            !to avoid temporary array in call to putarrnd:
            tmp_field = emfields(li1:li2,lj1:lj2,lk1:lk2,o)
            IF ((o .EQ. 2) .AND. (Apar0_antenna .NE. 0.0)) tmp_field = Apar_pre_antenna
            call putarrnd(fidfield_h5, dset_name_field, tmp_field, (/3, 2, 1/))
            call attach(fidfield_h5, dset_name_field, "time", time)      
         enddo

         isnap_field = isnap_field+1
         call flushh5(fidfield_h5)
      end if
#endif


#ifdef WITHFUTILS
    if(write_h5) then
       if ((my_pev+my_pew+my_pespec).eq.0) then

          call closef(fidfield_h5)
       end if
    end if
#endif
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank, o, trap_level
#endif

#ifdef WITHFUTILS
    if ((write_h5).and.(my_pev+my_pew+my_pespec).eq.0) then
       isnap_mom = 0
       allocate(fidmom_h5(0:n_spec-1))
       allocate(momtrp_label(0:diag_trap_levels+1))
       if (diag_trap_levels.gt.0) THEN
          momtrp_label(:) = '_trap'
          momtrp_label(diag_trap_levels) = '_pass'
          momtrp_label(diag_trap_levels+1) = '_FLR'
       else
          momtrp_label(:) = ''
       endif

       do n=0,n_spec-1
          call creatf(trim(diagdir)//'/mom_'//trim(spec(n)%name)//trim(file_extension)//'.h5', &
               fidmom_h5(n), "moments diagnostics", 'd', mpi_comm_xyz)
          call creatg(fidmom_h5(n), '/mom_'//trim(spec(n)%name))
          rank = 0
          call creatd(fidmom_h5(n), rank, dims, "/mom_"//trim(spec(n)%name)//"/time", "time")

          do o = 1, 6 + n_moms_Bpar
             do trap_level=0,diag_trap_levels
                call creatg(fidmom_h5(n), '/mom_'//trim(spec(n)%name)//'/'//&
                     &trim(mom_label(o))//trim(momtrp_label(trap_level)))
             enddo
             if (diag_trap_levels.gt.0) call creatg(fidmom_h5(n), '/mom_'//&
                  &trim(spec(n)%name)//'/'//trim(mom_label(o))//&
                  &trim(momtrp_label(diag_trap_levels+1)))
          enddo
        end do
    endif
#endif       

#ifdef WITHFUTILS
    complex, dimension(:,:,:,:,:),allocatable:: tmp_buf,tmpcorr_buf
    character(len=128) :: dset_name_mom
    integer:: ng, pes, ierr
#else

#ifdef WITHFUTILS
    if(write_h5) then
       allocate(tmp_buf(li1:li2,lj1:lj2,lk1:lk2,6+n_moms_Bpar,0:n_procs_s-1))
       allocate(tmpcorr_buf(li1:li2,lj1:lj2,lk1:lk2,n_corrs,0:n_procs_s-1))
    end if
#endif
#ifdef WITHFUTILS
       if (write_h5.and.((my_pev+my_pew) .eq. 0))  then
          if(my_pespec .eq. 0) then
             do pes = 0, n_procs_s-1
                ng = pes*ln0 + n-ln1
                call append(fidmom_h5(ng), "/mom_"//trim(spec(ng)%name)//"/time", time)
                call attach(fidmom_h5(ng), "/mom_"//trim(spec(ng)%name)//"/time", "n_steps", isnap_mom+1)
             end do
          endif

          do trap_level=0,diag_trap_levels
             call mpi_gather(vmom(:,:,:,:,trap_level), (6+n_moms_Bpar)*lijk0, MPI_COMPLEX_TYPE, tmp_buf,&
                  &(6+n_moms_Bpar)*lijk0, MPI_COMPLEX_TYPE, 0, mpi_comm_spec, ierr)

             if(my_pespec .eq. 0) then
                do pes = 0, n_procs_s-1
                   ng = pes*ln0 + n-ln1
                   do o = 1, 6 + n_moms_Bpar
                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(trap_level)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom),&
                           & tmp_buf(:,:,:,o,pes), (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)

                   enddo
                   call flushh5(fidmom_h5(ng))
                 end do
             end if
          enddo

          if (diag_trap_levels .gt. 0) then
             call mpi_gather(vmom_corr(:,:,:,:), n_corrs*lijk0, MPI_COMPLEX_TYPE, tmpcorr_buf,&
                  &n_corrs*lijk0, MPI_COMPLEX_TYPE, 0, mpi_comm_spec, ierr)
             
             if(my_pespec .eq. 0) then
                do pes = 0, n_procs_s-1
                   ng = pes*ln0 + n-ln1

                   do o = 1, n_corrs

                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(diag_trap_levels+1)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom), tmpcorr_buf(:,:,:,o,pes), (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)
                   enddo

                   do o=n_corrs+1,6
                      ! no flr corrections, writing zeros for compability with IDL diagnostic
                      write(dset_name_mom, "(A, '/', i10.10)") "/mom_"&
                           &//trim(spec(ng)%name)//"/"//trim(mom_label(o))//&
                           &trim(momtrp_label(diag_trap_levels+1)), isnap_mom
                      call putarrnd(fidmom_h5(ng), trim(dset_name_mom), dummy_corr, (/3, 2, 1/))
                      call attach(fidmom_h5(ng), trim(dset_name_mom), "time", time)
                   enddo

                   call flushh5(fidmom_h5(ng))
                end do
             end if
          end if
       endif
#endif

#ifdef WITHFUTILS
    if(write_h5) then
       isnap_mom = isnap_mom+1
       deallocate(tmp_buf,tmpcorr_buf)
    end if
#endif

#ifdef WITHFUTILS
    if (write_h5.and.(my_pev+my_pew+my_pespec).eq.0) then
       do n=0, n_spec-1
          call closef(fidmom_h5(n))
       end do
       deallocate(fidmom_h5,momtrp_label)
    end if
#endif
  End Subroutine finalize_diag_mom


!!!******************************************************************!!!
!!!******************************************************************!!!
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: o, rank
#endif

    if (mype==0) then
       if(write_std) then    
          call get_unit_nr(VSPFILE)
          OPEN(VSPFILE, file=trim(diagdir)//&
               &'/vsp'//trim(file_extension), form='unformatted', &
               status=filestat, position=filepos)
       END IF

#ifdef WITHFUTILS
       if (write_h5) then
          isnap_vsp = 0
          call creatf(trim(diagdir)//'/vsp'//trim(file_extension)//'.h5', &
               fidvsp_h5, "Velocity space diagnostics", 'd')
          call creatg(fidvsp_h5, '/vsp')
          do o = 1,5
             call creatg(fidvsp_h5, '/vsp/'//trim(vsp_label(o)))
          enddo
          rank = 0
          call creatd(fidvsp_h5, rank, dims, "/vsp/time", "time")
       end if
#endif
    endif
#ifdef WITHFUTILS
    character(len=128):: dset_name_vsp
#endif

#ifdef WITHFUTILS
    if(write_h5) then
      if(mype.eq.0) then          
          call append(fidvsp_h5, "/vsp/time", time)
          call attach(fidvsp_h5, "/vsp/time", "n_steps", isnap_vsp+1)
          do o = 1,5
             write(dset_name_vsp, "(A, '/', i10.10)") "/vsp/"//&
                  &trim(vsp_label(o)), isnap_vsp
             call putarr(fidvsp_h5, dset_name_vsp, fullarr(:,:,:,:,o)*fnorm)
             call attach(fidvsp_h5, dset_name_vsp, "time", time)             
          enddo

          call flushh5(fidvsp_h5)
       end if
      isnap_vsp = isnap_vsp+1
   end if
#endif

#ifdef WITHFUTILS
       if (write_h5) call closef(fidvsp_h5)
#endif
