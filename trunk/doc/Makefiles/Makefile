
#Use GNU C compiler
#CC = gcc

#CFLAGS = -O3 
#CFLAGS = -O3 -static 
#CFLAGS = -O0 -Wall
#CFLAGS = -g -O0 


#Use Intel C compiler
CC = icc

#wd: disable: remark #981: operands are evaluated in unspecified order
CFLAGS = -O3 -traceback
#CFLAGS = -O3 -static -traceback
#CFLAGS = -O3 -static -traceback -Wcheck -Wall -wd981 
CFLAGS = -O0 -traceback -Wcheck -Wall -wd981 
#CFLAGS = -g -O0  -traceback  -Wcheck -Wall -wd981 





DEF_MCMCOBJ = SPINspiral_main.o SPINspiral_mcmc.o SPINspiral_templates.o SPINspiral_routines.o SPINspiral_data.o SPINspiral_parameters.o SPINspiral_signal.o SPINspiral_lal.o SPINspiral_3rdparty.o 
DEF_CFLAGS = -I/usr/include -I./include -I ~/bin/include   -I ~/local/programs/lscsoft/lal/include/ -I ~/local/programs/lscsoft/libmetaio/include/
DEF_LFLAGS = -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
LAL_LFLAGS = -L ~/bin/lib -llal -llalmetaio 






%.o: %.c Makefile include/SPINspiral.h
	$(CC) $(CFLAGS) $(DEF_CFLAGS) -c $<

all: SPINspiral

exit:


SPINspiral_3rdparty.o: include/remez.h


SPINspiral: $(DEF_MCMCOBJ)
	$(CC) -o SPINspiral $(DEF_MCMCOBJ) $(DEF_LFLAGS) $(LAL_LFLAGS)




clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f SPINspiral

cleanoutput:
	rm -f mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

