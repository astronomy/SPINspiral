
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
#CFLAGS = -O0 -traceback -Wcheck -Wall -wd981 
#CFLAGS = -g -O0  -traceback 



DEF_MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o 
DEF_CFLAGS = -I/usr/include -I./include -I ~/bin/include   -I ~/local/programs/lscsoft/lal/include/
DEF_LFLAGS = -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
LAL_LFLAGS = -L ~/bin/lib -llal 





%.o: %.c Makefile include/mcmc.h
	$(CC) $(CFLAGS) $(DEF_CFLAGS) -c $<

exit:

all: mcmc_nolal mcmc_lal


mcmc_lal.o: include/mcmc_lal.h
mcmc_3rdparty.o: include/remez.h



mcmc_nolal: $(DEF_MCMCOBJ) mcmc_nolal.o
	$(CC) -o mcmc_nolal $(DEF_MCMCOBJ) mcmc_nolal.o $(DEF_LFLAGS)

mcmc_lal: $(DEF_MCMCOBJ) mcmc_lal.o
	$(CC) -o mcmc_lal $(DEF_MCMCOBJ) mcmc_lal.o $(DEF_LFLAGS) $(LAL_LFLAGS)




clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f mcmc 

cleanoutput:
	rm -f mcmc.input.?????? mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

