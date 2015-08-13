###########################################################################
### The source files	                                                ###
###########################################################################
F90PSRC = \
	decomp_2d.f90\
	glassman.f90\
	fft_$(FFTLIB).f90\
	module_param.f90\
	io.f90\
	variables.f90\
	schemes.f90\
	FR_perf.f90\
	inc3d.f90\
	derive.f90\
	parameters.f90\
	tools.f90\
	visu.f90

F90FULLSRC = $(addprefix $(SRCDIR)/,$(F90PSRC))
F90PRENAMES = $(addprefix F,$(subst .f90,.f,$(F90PSRC)))
F90POBJ = $(addprefix $(OBJDIR)/,$(subst .f90,.o,$(F90PSRC)))

OBJLIST += $(F90POBJ)

RPL_STR =-WF,
CPREPROC += $(subst $(RPL_STR),,$(PREPROC))

#For libraries
LIBELEMENTS = $(OBJLIST) $(OBJDIR)/$(LIBN).o
LIBN2=$(subst lib,,$(LIBN))

# to make some conditionals in the rules simpler, we map the case
# with USE_PERFLIB=none on the case when USE_PERFLIB is undefined or empty
USE_PERFLIB := $(if $(findstring none,$(USE_PERFLIB)),,$(strip $(USE_PERFLIB)))

###########################################################################
### The dependencies  (might be checked/updated using ../tools/makedep) ###
###########################################################################

$(OBJLIST):		$(MKFILE) $(BASEDIR)/makefile

$(OBJDIR)/derive.o:	$(OBJDIR)/module_param.o

$(OBJDIR)/fft_generic.o:	$(OBJDIR)/decomp2d.o\
				$(OBJDIR)/glassman.o

$(OBJDIR)/glassman.o:	$(OBJDIR)/decomp2d.o

$(OBJDIR)/inc3d.o:	$(OBJDIR)/module_param.o\
			$(OBJDIR)/decomp2d.o\
			$(OBJDIR)/variables.o\
			$(OBJDIR)/derive.o\
			$(OBJDIR)/FR_perf.o

$(OBJDIR)/io.o:	$(OBJDIR)/decomp2d.o\
		$(OBJDIR)/variables.o\
		$(OBJDIR)/module_param.o

$(OBJDIR)/module_param.o:	$(OBJDIR)/decomp2d.o

$(OBJDIR)/parameters.o:	$(OBJDIR)/decomp2d.o\
			$(OBJDIR)/variables.o\
			$(OBJDIR)/module_param.o

$(OBJDIR)/schemes.o:	$(OBJDIR)/module_param.o\
			$(OBJDIR)/derive.o\
			$(OBJDIR)/variables.o

$(OBJDIR)/tools.o:	$(OBJDIR)/decomp2d.o\
			$(OBJDIR)/variables.o\
			$(OBJDIR)/module_param.o

$(OBJDIR)/variables.o:	$(OBJDIR)/decomp2d.o\
			$(OBJDIR)/module_param.o

$(OBJDIR)/visu.o:	$(OBJDIR)/decomp2d.o\
			$(OBJDIR)/variables.o\
			$(OBJDIR)/module_param.o

$(OBJDIR)/alloc.o:	
$(OBJDIR)/decomp2d.o:	
$(OBJDIR)/factor.o:	
$(OBJDIR)/fft_common3d.o:	
$(OBJDIR)/fft_common.o:		
$(OBJDIR)/fft_fftw3.o:		
$(OBJDIR)/FR_perf.o:	
$(OBJDIR)/halo.o:	
$(OBJDIR)/halo_common.o:	
$(OBJDIR)/io_new.o:	
$(OBJDIR)/io_old.o:	
$(OBJDIR)/io_read_one.o:	
$(OBJDIR)/io_read_var.o:	
$(OBJDIR)/io_write_every.o:	
$(OBJDIR)/io_write_one.o:	
$(OBJDIR)/io_write_plane.o:	
$(OBJDIR)/io_write_var.o:	
$(OBJDIR)/transpose_x_to_y.o:		
$(OBJDIR)/transpose_y_to_x.o:		
$(OBJDIR)/transpose_y_to_z.o:		
$(OBJDIR)/transpose_z_to_y.o:		

################################
##### FUTILS dep. parts ########
################################
ifeq ($(FUTILS),yes)

$(OBJDIR)/inc3d.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/visu.o:			$(FUTILSDIR)/libfutils.a

ifeq ($(FUTILSDIR),$(EXTDIR)/futils/src)
export FUTILSDIR FC FFLAGS LDFLAGS CC CFLAGS HDF5PATH OPT ARCHIVE HDF5_LIBS HDF5VAR
$(FUTILSDIR)/libfutils.a: 
	@make -C $(EXTDIR) futils
endif

endif

####################################
##### PERFormance libraries ########
####################################
ifneq ($(strip $(USE_PERFLIB)),)
 ifeq ($(USE_PERFLIB),FR)
  $(OBJDIR)/libmyperf.a: $(OBJDIR)/FR_perf.o
 endif
 $(EXECDIR)/$(EXEC): $(OBJDIR)/libmyperf.a
endif

show_srcnames:
	@echo $(F90FULLSRC)
	for x in $(F90FULLSRC);\
	do\
		echo $$x >> flist.txt;\
	done

# --- create executable
$(EXECDIR)/$(EXEC): $(OBJLIST) $(OBJDIR)/inc3d.o
	@echo "Linking "
	$(LD) -o $(OBJDIR)/inc3d.o -o $@ $(OBJLIST) $(LDFLAGS) $(LIBS)

# --- create library
$(EXECDIR)/$(LIBN).a:	$(LIBELEMENTS) $(OBJDIR)/test_$(LIBN).o
	$(ARCHIVE) $@ $(LIBELEMENTS)
	ranlib $@
	$(FC) $(FFLAGS) $(INCPATHS) $(PREPROC) -c -o $(OBJDIR)/test_$(LIBN).o $(SRCDIR)/test_$(LIBN).f90
	$(LD) -o $(EXECDIR)/test_$(LIBN) $(LIBELEMENTS) $(OBJDIR)/test_$(LIBN).o -L$(EXECDIR) -l$(LIBN2) $(LDFLAGS) $(LIBS)
	@echo ""
	@echo "--- Link $(LIBN) to your program with ---"
	@echo "$(LD) $(FFLAGS) $(INCPATHS) -o YOUR_PROGRAM YOUR_PROGRAM.f90 -L$(EXECDIR) -l$(LIBN2) $(LIBS)"
	@echo ""


###########################################################################
### Compiling rules (FLAGS should be set in machine dep. makefiles)     ###
###########################################################################
#note: adding the @ sign in front of $(CC) and $(FC) will suppress the 
#      full compiler call which increases readability

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo "Compiling " $<
	$(CC) $(CFLAGS) $(INCPATHS) $(CPREPROC) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@echo "Compiling " $<
	$(FC) $(FFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

#use the following line for profiling of compilation
# @time $(FC) $(FFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

###########################################################################
