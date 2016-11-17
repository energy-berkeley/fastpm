-include ../Makefile.local
-include ../Makefile.rules

.PHONY: all clean
BIGFILE_LIBS = libbigfile.a libbigfile-mpi.a
MPSORT_LIBS = libradixsort.a libmpsort-mpi.a
PFFT_LIBS = libfftw3.a \
			libfftw3_mpi.a \
			libfftw3_omp.a \
			libfftw3f.a \
			libfftw3f_mpi.a \
			libfftw3f_omp.a \
			libpfft_omp.a \
			libpfftf_omp.a \

BIGFILE_H = libbigfile.h libbigfile-mpi.h
MPSORT_H = libradixsort.h libmpsort-mpi.h
PFFT_H = libfftw3.h \
			libfftw3_mpi.h \
			libfftw3_omp.h \
			libfftw3f.h \
			libfftw3f_mpi.h \
			libfftw3f_omp.h \
			libpfft_omp.h \
			libpfftf_omp.h \

.PHONY: all

all : libfastpm-dep.a .header

libfastpm-dep.a: $(BIGFILE_LIBS:%.a=install/lib/%.a) \
	  $(MPSORT_LIBS:%.a=install/lib/%.a) \
	  $(PFFT_LIBS:%.a=install/lib/%.a)
	../ar-join $@ $^ #create a lib file

.header: $(BIGFILE_H:%.h=install/include/%.h) \
	  $(MPSORT_H:%.h=install/include/%.h) \
	  $(PFFT_H:%.h=install/include/%.h)
	#header no need a library file

.INTERMEDIATE: .bigfile .mpsort .pfft

$(BIGFILE_LIBS:%.a=install/lib/%.a) $(BIGFILE_H:%.h=install/include/%.h): .bigfile

.bigfile:
	(cd bigfile; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC)" "MPICC=$(CC)")

$(MPSORT_LIBS:%.a=install/lib/%.a) $(MPSORT_H:%.h=install/include/%.h): .mpsort

.mpsort:
	(cd mpsort; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC)" "MPICC=$(CC)")

$(PFFT_LIBS:%.a=install/lib/%.a) $(PFFT_H:%.h=install/include/%.h): .pfft

.pfft:
	($(MAKE) "CC=$(CC)" -f Makefile.pfft)






