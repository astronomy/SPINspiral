
#Use GNU C compiler
CC = gcc
CFLAGS = -O3 -static -I/usr/include -I ./include
#CFLAGS = -g -O0 -I/usr/include -I ./include

#Use Intel C compiler
#CC = icc
#CFLAGS = -O3 -static -I/usr/include -I ./include -traceback -Wcheck -Wall -wd981 #wd: disable: remark #981: operands are evaluated in unspecified order
#CFLAGS = -g -O0 -I/usr/include -I ./include -traceback -Wcheck -Wall -wd981 #wd: disable: remark #981: operands are evaluated in unspecified order


#MCMCLFLAGS = -L ~/bin/lib -lFrame -lm -lgsl -lgslcblas -lfftw3
MCMCLFLAGS = -static -L ~/bin/lib -lFrame -lm -lgsl -lgslcblas -lfftw3
MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o

#WAVELFLAGS = -lfftw3
WAVELFLAGS = -static -lfftw3
WAVEOBJ = waveform_main.o mcmc_templates.o mcmc_routines.o waveform_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o



.c.o:
	$(CC) $(CFLAGS) -c $?

exit:


mcmc: $(MCMCOBJ) 
	$(CC) -o mcmc $(MCMCOBJ) $(MCMCLFLAGS)

waveform: $(WAVEOBJ) 
	$(CC) -o waveform $(WAVEOBJ) $(WAVELFLAGS)

clean:
	rm -f *.o 

cleanoutput:
	rm -f mcmc.input.* mcmc.output.*



