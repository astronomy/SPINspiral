COMP = ifort
#CFLAGS = -vec-report0 -free -O0
#CFLAGS = -g -vec-report0 -free -O0
#CFLAGS = -vec-report0 -free -O0  -warn alignments -warn declarations -warn uncalled -warn usage #-warn unused
CFLAGS = -vec-report0 -free -O0  -warn alignments -warn declarations -warn uncalled -warn usage -check bounds -check format -check output_conversion -traceback #-warn unused
#LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -L/opt/local/lib -laquaterm -Wl,-framework -Wl,Foundation -lifcore -limf -lirc -lsvml -lifport -Xlinker -multiply_defined suppress -bind_at_load 
#LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -L/opt/local/lib -Wl,-framework -Wl,Foundation -lifcore -limf -lirc -lsvml -lifport -Xlinker -multiply_defined suppress -bind_at_load 
LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -lifcore -limf -lirc -lsvml -lifport -Xlinker 

#gfortran needs ~60% more time (mostly for reading data in plotspins, it seems)
#COMP = gfortran
#CFLAGS = -O0 -ffree-form -ffree-line-length-none
#LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -L/opt/intel/fc/10.0.026/lib -lifcore -limf -lirc -lsvml -lifport -Xlinker 


#Comp time:
# O0: 1.2s
# O1: 7.4s
# O2: 5.8s
# O3: 8.0s

STDOBJ =  


.f.o: 
	$(COMP) $(CFLAGS) -c $?

exit:


#convert_mcmcoutput: $(STDOBJ) convert_mcmcoutput.o
#	$(COMP) $(LFLAGS) -o convert_mcmcoutput $(STDOBJ) convert_mcmcoutput.o

convergence: $(STDOBJ) convergence.o
	$(COMP) $(LFLAGS) -o convergence $(STDOBJ) convergence.o

plotdata: $(STDOBJ) plotdata.o
	$(COMP) $(LFLAGS) -o plotdata $(STDOBJ) plotdata.o

plotspectrum: $(STDOBJ) plotspectrum.o
	$(COMP) $(LFLAGS) -o plotspectrum $(STDOBJ) plotspectrum.o

plotsignal: $(STDOBJ) plotsignal.o
	$(COMP) $(LFLAGS) -o plotsignal $(STDOBJ) plotsignal.o

plotspins: $(STDOBJ) plotspins.o
	$(COMP) -o plotspins $(LFLAGS) $(STDOBJ) plotspins.o

plotspins_old: $(STDOBJ) plotspins_old.o
	$(COMP) $(LFLAGS) -o plotspins_old $(STDOBJ) plotspins_old.o

plotstats: $(STDOBJ) plotstats.o
	$(COMP) $(LFLAGS) -o plotstats $(STDOBJ) plotstats.o

clean:
	rm *.o
