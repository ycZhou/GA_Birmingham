SHELL=/bin/sh
SOURCE = ./source
OBJDIR = ./build
LIBDIR = ./lib
BINDIR = ./bin
DOCDIR = ./doc
AR = ar crl 
RANLIB = ranlib
RM = rm -rf
MKDIR = mkdir -p
LIB = $(LIBDIR)/libga.a
EXEC = $(BINDIR)/ga

include make.flags

ga :  dir $(OBJS) lib
	$(F77) $(OBJS) -o $(EXEC) $(LDFLAGS)


OBJS = $(OBJDIR)/ga.o

LOBJS =	$(OBJDIR)/cluster.o \
	$(OBJDIR)/fitness.o \
	$(OBJDIR)/gen_pop.o \
	$(OBJDIR)/gen_pop_ionic.o \
	$(OBJDIR)/gen_pop_silica.o \
	$(OBJDIR)/ionic.o \
	$(OBJDIR)/silica.o \
	$(OBJDIR)/mate.o \
	$(OBJDIR)/minimiser.o \
	$(OBJDIR)/morse.o \
	$(OBJDIR)/murrmott.o \
	$(OBJDIR)/gupta.o \
	$(OBJDIR)/mutate.o \
	$(OBJDIR)/pop.o \
	$(OBJDIR)/print_sum.o \
	$(OBJDIR)/ran3.o \
	$(OBJDIR)/read_iparams.o \
	$(OBJDIR)/read_mmparams.o \
	$(OBJDIR)/read_morse.o \
	$(OBJDIR)/read_gupta.o \
	$(OBJDIR)/read_cutoff.o \
	$(OBJDIR)/routines.o \
	$(OBJDIR)/select.o \
	$(OBJDIR)/utils.o \
	$(OBJDIR)/xyz.o \
	$(OBJDIR)/z_heapsort.o \
	$(OBJDIR)/token.o \
	$(OBJDIR)/inter.o \
	$(OBJDIR)/DFT_module.o \
	$(OBJDIR)/write_DFT.o \
	$(OBJDIR)/read_DFT.o \
        $(OBJDIR)/output_DFT.o \
        $(OBJDIR)/restart.o \
        $(OBJDIR)/moi.o \
        $(OBJDIR)/ga_commons.o \
	$(OBJDIR)/ga.o


$(OBJDIR)/cluster.o:	$(SOURCE)/cluster.f $(SOURCE)/ran3.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/cluster.f

$(OBJDIR)/fitness.o:	$(SOURCE)/fitness.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/fitness.f

$(OBJDIR)/gen_pop.o:	$(SOURCE)/gen_pop.f $(SOURCE)/ran3.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/gen_pop.f

$(OBJDIR)/gen_pop_ionic.o:	$(SOURCE)/gen_pop_ionic.f $(SOURCE)/ran3.f 
				$(F77) $(FFLAGS) -o $@ \
				$(SOURCE)/gen_pop_ionic.f

$(OBJDIR)/gen_pop_silica.o:	$(SOURCE)/gen_pop_silica.f $(SOURCE)/ran3.f 
				$(F77) $(FFLAGS) -o $@ \
				$(SOURCE)/gen_pop_silica.f

$(OBJDIR)/ionic.o:	$(SOURCE)/ionic.f $(OBJDIR)/ga_commons.o
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/ionic.f

$(OBJDIR)/silica.o:	$(SOURCE)/silica.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/silica.f

$(OBJDIR)/mate.o:	$(SOURCE)/mate.f $(SOURCE)/ran3.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/mate.f

$(OBJDIR)/minimiser.o:	$(SOURCE)/minimiser.f $(SOURCE)/routines.f $(SOURCE)/write_DFT.f $(OBJDIR)/read_DFT.o
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/minimiser.f

$(OBJDIR)/morse.o:	$(SOURCE)/morse.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/morse.f

$(OBJDIR)/murrmott.o:	$(SOURCE)/murrmott.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/murrmott.f

$(OBJDIR)/gupta.o:	$(SOURCE)/gupta.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/gupta.f

$(OBJDIR)/mutate.o:	$(SOURCE)/mutate.f $(SOURCE)/z_heapsort.f \
			$(SOURCE)/ran3.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/mutate.f

$(OBJDIR)/pop.o:	$(SOURCE)/pop.f $(SOURCE)/minimiser.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/pop.f

$(OBJDIR)/print_sum.o:	$(SOURCE)/print_sum.f $(SOURCE)/utils.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/print_sum.f

$(OBJDIR)/ran3.o:	$(SOURCE)/ran3.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/ran3.f

$(OBJDIR)/read_iparams.o:	$(SOURCE)/read_iparams.f
				$(F77) $(FFLAGS) -o $@ $(SOURCE)/read_iparams.f

$(OBJDIR)/read_mmparams.o:	$(SOURCE)/read_mmparams.f
				$(F77) $(FFLAGS) -o $@ \
				$(SOURCE)/read_mmparams.f

$(OBJDIR)/read_morse.o:	$(SOURCE)/read_morse.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/read_morse.f

$(OBJDIR)/read_gupta.o:	$(SOURCE)/read_gupta.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/read_gupta.f

$(OBJDIR)/read_cutoff.o:	$(SOURCE)/read_cutoff.f
				$(F77) $(FFLAGS) -o $@ $(SOURCE)/read_cutoff.f

$(OBJDIR)/routines.o: 	$(SOURCE)/routines.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/routines.f

$(OBJDIR)/select.o:	$(SOURCE)/select.f $(SOURCE)/ran3.f 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/select.f

$(OBJDIR)/utils.o:	$(SOURCE)/utils.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/utils.f

$(OBJDIR)/xyz.o:	$(SOURCE)/xyz.f $(SOURCE)/utils.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/xyz.f

$(OBJDIR)/z_heapsort.o:	$(SOURCE)/z_heapsort.f
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/z_heapsort.f

$(OBJDIR)/token.o:	$(SOURCE)/token.F 
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/token.F

$(OBJDIR)/inter.o:	$(SOURCE)/inter.f $(OBJDIR)/ga_commons.o
			$(F77) $(FFLAGS) -o $@ $(SOURCE)/inter.f

$(OBJDIR)/DFT_module.o:  $(SOURCE)/DFT_module.f90
	                $(F77) $(FFLAGS) -o $@ $(SOURCE)/DFT_module.f90

$(OBJDIR)/write_DFT.o:  $(SOURCE)/write_DFT.f $(OBJDIR)/read_DFT.o $(OBJDIR)/DFT_module.o
	                $(F77) $(FFLAGS) -o $@ $(SOURCE)/write_DFT.f

$(OBJDIR)/read_DFT.o:   $(SOURCE)/read_DFT.f $(OBJDIR)/DFT_module.o $(OBJDIR)/ga_commons.o
	                $(F77) $(FFLAGS) -o $@ $(SOURCE)/read_DFT.f

$(OBJDIR)/output_DFT.o: $(SOURCE)/output_DFT.f $(OBJDIR)/DFT_module.o
	                $(F77) $(FFLAGS) -o $@ $(SOURCE)/output_DFT.f

$(OBJDIR)/restart.o:$(SOURCE)/restart.f
	                $(F77) $(FFLAGS) -o $@ $(SOURCE)/restart.f

$(OBJDIR)/moi.o: $(SOURCE)/moi.f
			 $(F77) $(FFLAGS) -o $@ $(SOURCE)/moi.f

$(OBJDIR)/ga_commons.o: $(SOURCE)/ga_commons.f90
			 $(F77) $(FFLAGS) -o $@ $(SOURCE)/ga_commons.f90

$(OBJDIR)/ga.o: $(SOURCE)/ga.f
			 $(F77) $(FFLAGS) -o $@ $(SOURCE)/ga.f

dir :
	if [ ! -d $(OBJDIR) ] ; then $(MKDIR) $(OBJDIR) ;fi
	if [ ! -d $(BINDIR) ] ; then $(MKDIR) $(BINDIR) ;fi
	if [ ! -d $(LIBDIR) ] ; then $(MKDIR) $(LIBDIR) ;fi

clean :
	rm -rv $(LIBDIR)/* $(OBJDIR)/* $(BINDIR)/* *.mod

lib : $(LOBJS)
	$(AR) $(LIB) $(LOBJS)
	$(RANLIB) $(LIB)
	
.PHONY : doc
doc : 
	cd $(DOCDIR) ;\
	latex manual ;\
	latex manual ;\
	dvips -o manual.ps manual.dvi ;\
	ps2pdf manual.ps manual.pdf
