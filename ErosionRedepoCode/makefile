sources := geompack3.f90 sputter_mod.f90 iomesh_mod.f90 Redepo_v3.f90 
objects := $(patsubst %.f90,%.o,$(sources))
# LDIR = ${MKLROOT}/lib/intel64 /opt/intel/compilers_and_libraries_2020.4.301/mac/compiler/lib
# LIBS = ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl

Redepo_v3 : $(objects)
	ifort -o Redepo_v3 $(objects) # $(LIBS) -I${MKLROOT}/include -fpp  -parallel

$(objects) : $(sources)
	ifort -c $(sources) #-I${MKLROOT}/include -fpp  -parallel



.PHONY: x=clean
clean:
	rm -f $(objects)
