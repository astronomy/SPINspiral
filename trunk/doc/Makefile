
#Use GNU C compiler
#CC = gcc
#CFLAGS = -O3 
#CFLAGS = -O3 -static 
#CFLAGS = -O0 -Wall
#CFLAGS = -g -O0 


#Use Intel C compiler
CC = icc

#wd: disable: remark #981: operands are evaluated in unspecified order
#CFLAGS = -O3 -traceback
#CFLAGS = -O3 -static -traceback
#CFLAGS = -O3 -static -traceback -Wcheck -Wall -wd981 
CFLAGS = -O0 -traceback -Wcheck -Wall -wd981 
#CFLAGS = -g -O0  #-traceback -Wcheck -Wall -wd981 



#Use LAL
DEFCFLAGS = -I/usr/include -I./include -I ~/bin/include -I ~/local/programs/lscsoft/lal/include/
MCMCLFLAGS = -L ~/bin/lib -llal -lFrame -lfftw3 -lgsl -lgslcblas -lm
MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o mcmc_lal.o 

#No LAL
##DEFCFLAGS = -ipo -parallel -I/usr/include -I./include -I ~/bin/include
DEFCFLAGS = -I/usr/include -I./include -I ~/bin/include 
##MCMCLFLAGS = -ipo -parallel -static -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
MCMCLFLAGS = -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
##MCMCLFLAGS = -static -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o mcmc_nolal.o 


WAVELFLAGS = -ipo -static -L ~/bin/lib -llal -lfftw3 -lgsl
WAVEOBJ = waveform_main.o mcmc_templates.o mcmc_routines.o waveform_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o mcmc_nolal.o



%.o: %.c include/mcmc.h
	$(CC) $(CFLAGS) $(DEFCFLAGS) -c $<
mcmc_lal.o: include/mcmc_lal.h
mcmc_3rdparty.o: include/remez.h



mcmc: $(MCMCOBJ) 
	$(CC) -o mcmc $(MCMCOBJ) $(MCMCLFLAGS)

waveform: $(WAVEOBJ) 
	$(CC) -o waveform $(WAVEOBJ) $(WAVELFLAGS)





clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f mcmc waveform

cleanoutput:
	rm -f mcmc.input.?????? mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

