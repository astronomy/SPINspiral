!Read and plot the data output from the spinning MCMC code.

program plotspins
  implicit none
  integer, parameter :: narr1=2.01e5+2,npar0=13,npar1=15,nchs=10,nbin=100,nbin2dx=60,nbin2dy=40,nival1=5,nr1=5,nstat1=10,ndets=3
  integer :: n(nchs),ntot(nchs),n0,n1,n2,nd,i,j,j1,j2,k,nburn(nchs),nburn0(nchs),iargc,io,pgopen,system,thin,narr,maxdots,reverseread
  integer :: niter(nchs),seed(nchs),ndet(nchs)
  integer :: index(npar1,nchs*narr1),index1(nchs*narr1)
  real :: is(nchs,narr1),isburn(nchs),jumps(nchs,npar1,narr1),startval(nchs,npar1,2)
  real :: sig(npar1,nchs,narr1),acc(npar1,nchs,narr1),avgtotthin
  real :: sn,dlp,xbin(nchs,nbin+1),ybin(nchs,nbin+1),xbin1(nbin+1),ybin1(nbin+1),ybin2(nbin+1),x(nchs,nchs*narr1),y(nchs,narr1),ybintot,xx(nchs*narr1),yy(nchs*narr1),zz(nchs*narr1)
  real :: ysum(nbin+1),yconv(nbin+1),ycum(nbin+1),a,b,r2d,r2h,rat,plx,ply
  real*8 :: t,t0,pi,tpi,dvar,dvar1,dvar2,ra,nullh
  real, allocatable :: dat(:,:,:),alldat(:,:,:),pldat(:,:,:)
  character :: js*4,varnames(npar1)*5,pgunits(npar1)*99,pgvarns(npar1)*99,pgvarnss(npar1)*99,pgorigvarns(npar1)*99,infile*100,infiles(nchs)*100,str*99,str1*99,str2*99,fmt*99,bla*10
  
  integer :: nn,nn1,lowvar(npar1),nlowvar,highvar(npar1),nhighvar,ntotrelvar,nlogl1,nlogl2,ksn1,ksn2
  real*8 :: chmean(nchs,npar1),totmean(npar1),chvar(npar1),chvar1(nchs,npar1),totvar(npar1),rhat(npar1),totrelvar,ksdat1(narr1),ksdat2(narr1),ksd,ksprob
  character :: ch
  
  integer :: samplerate(nchs,ndets),samplesize(nchs,ndets),FTsize(nchs,ndets),detnr(nchs,ndets)
  real :: snr(nchs,ndets),flow(nchs,ndets),fhigh(nchs,ndets),t_before(nchs,ndets),t_after(nchs,ndets),deltaFT(nchs,ndets)
  real*8 :: FTstart(nchs,ndets)
  character :: detname(nchs,ndets)*16,string*99
  
  integer :: nfx,nfy,fx,fy,file,update,plvars(npar1),nplvar,smooth,fillpdf,quality,combinechainplots,chainsymbol
  real :: x0,x1,x2,y1,y2,dx,dy,xmin,xmax,ymin,ymax,xmin1,xmax1,xpeak,ymin1,ymax1,ymaxs(nchs+2),z(nbin2dx+1,nbin2dy+1),zs(nchs,nbin2dx+1,nbin2dy+1)
  real :: coefs(100),coefs1(100),cont(11),tr(6)
  
  integer :: nchains,nchains0,ic,ip,i0,i1,i2,lw,lw1,lw2,rdsigacc
  integer :: prprogress,prruninfo,prinitial,prstat,prcorr,prival,prconv,prvalues
  integer :: placorr,plot,pllogl,plchain,pljump,plsigacc,plpdf1d,savepdf,plpdf2d,plotsky,plmovie,chainpli,logjump,logsig,pltrue,plstart,plmedian,plrange
  integer :: moviescheme,savestats
  real :: sch
  character :: header*1000,outputname*99,psclr*4,colournames(15)*20
  
  integer :: o,o1,o2,do12,p,p1,p2,p11,p12,p21,p22,par1,par2,nr,c,c0,nstat,wrap(nchs,npar1),nival,wrapdata,changevar,mergechains,npar,colour,ncolours,colours(10),defcolour,plotthis,tempintarray(99)
  real :: range,minrange,range1,range2,drange,maxgap,ranges(nchs,nival1+1,npar1,nr1),ival,ival0,ivals(nival1+1),centre,shift(nchs,npar1),plshift
  real :: median,medians(npar1),mean(npar1),stdev1(npar1),stdev2(npar1),var1(npar1),var2(npar1),absvar1(npar1),absvar2(npar1)
  real :: stats(nchs,npar1,nstat1),corrs(npar1,npar1),acorrs(nchs,0:npar1,0:narr1)
  real :: norm
  
  integer :: xpanels,ypanels
  real :: scrsz,scrrat,bmpsz,bmprat,rev360,rev24,rev2pi
  character :: bmpxpix*99
  
  integer :: nframes,iframe,nplt
  character :: framename*99
  
  pi = 4*datan(1.d0)
  tpi = 2*pi
  r2d = real(180.d0/pi)
  r2h = real(12.d0/pi)
  
  thin = 10       !If >1, 'thin' the output; read every thin-th line 
  nburn = 4e6       !If >=0: override length of the burn-in phase, for all chains! This is now the ITERATION number, but it becomes the line number later on in the code.  Nburn > Nchain sets Nburn = 0.1*Nchain
  file = 0          !Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf
  colour = 1        !Use colours: 0-no (grey scales), 1-yes
  quality = 0       !'Quality' of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster
  reverseread = 1   !Read files reversely (anti-alphabetically), to plot coolest chain last so that it becomes better visible: 0-no, 1-yes, 2-use colours in reverse order too
  update = 0        !Update screen plot every 10 seconds: 0-no, 1-yes
  mergechains = 1   !Merge the data from different files into one chain: 0-no (treat separately), 1-yes
  wrapdata = 1      !Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is around 0)
  changevar = 1     !Change variables (e.g. logd->d, kappa->theta_SL, rad->deg)
  prprogress = 1    !Print general messages about the progress of the program: 0-no, 1-yes
  prruninfo = 0     !Print run info at read (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-yes.
  prinitial = 0     !Print true values, starting values and their difference
  prstat = 0        !Print statistics: 0-no, 1-yes
  prcorr = 0        !Print correlations: 0-no, 1-yes
  prival = 0        !Print interval info: 0-no, 1-yes
  prconv = 0        !Print convergence information for multiple chains to screen and chains plot: 0-no, 1-yes: 1 summary line, 2-yes: medians, stdevs etc too.
  savestats = 0     !Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS
  savepdf = 0       !Save the binned data for 1d and/or 2d pdfs (depending on plpdf1d and plpdf2d).  This causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying
  
  plot = 1          !0: plot nothing at all, 1: plot the items selected below
  combinechainplots = 0  !Combine logL, chain, sigma and acc plots into one multipage file
  pllogl = 1        !Plot log L chains: 0-no, 1-yes
  plchain = 0       !Plot parameter chains: 0-no, 1-yes
  pljump = 0        !Plot actual jump sizes
  logjump = 1       !Plot the log of the jump size: 0-no, 1-yes
  rdsigacc = 0      !Read sigma and acceptance rate: 0-no, 1-yes   (0-Don't read these data, save 40% read-in time).  0 can give problems with large scale, or high-temperature chains
  plsigacc = 0      !Plot sigma and acceptance rate: 0-no, 1-yes   (Sets rdsigacc to 1)
  logsig = 1        !Plot the log of sigma: 0-no, 1-yes
  plpdf1d = 0       !Plot 1d posterior distributions: 0-no, 1-yes. If plot=0 and savepdf=1, this determines whether to write the pdfs to file or not.
  plpdf2d = 0       !Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours only. If plot=0 and savepdf=1, this determines whether to write the pdfs to file (>0) or not (=0).
  placorr = 0e4     !Plot autocorrelations: 0-no, >0-yes: plot placorr steps
  plotsky = 0       !Plot 2d pdf with stars, implies plpdf2d=1
  plmovie = 0       !Plot movie frames
  
  chainsymbol = 1   !Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default), 2: plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle
  chainpli = 0      !Plot every chainpli-th point in chains, logL, jump plots:  chainpli=0: autodetermine, chainpli>0: use this chainpli.  All states in between *are* used for statistics, pdf generation, etc.
  pltrue = 1        !Plot true values in the chains and pdfs
  plstart = 1       !Plot starting values in the chains and pdfs
  plmedian = 1      !Plot median values in the pdfs
  plrange = 1       !Plot the probability range in the pdfs
  prvalues = 1      !Print values (true, median, range) in pdfs
  smooth = 3        !Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).   This is 1D only for now, and can introduce artefacts on narrow peaks!
  fillpdf = 1       !Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched
  ival0 = 0.90      !Standard probability interval, e.g. 0.90, 0.95
  nframes = 1       !Number of frames for the movie
  
  par1 = 1          !First parameter to treat (stats, plot): 0-all
  par2 = 15         !Last parameter to treat (0: use npar)
  
  scrsz  = 10.8     !Screen size for X11 windows:  MacOS: 16.4, Gentoo: 10.8
  scrrat = 0.57     !Screen ratio for X11 windows, MacBook: 0.57
  bmpsz  = 12.      !Size for bitmap:  10.6
  bmprat = 0.75     !Ratio for bitmap: 0.75
  !bmprat = 1.25     !Ratio for bitmap: 0.75
  bmpxpix = '850'     !Final size of converted bitmap output in pixels (string)
  
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  !Choose plot variables.  Number of variables to plot: 1,2,3,4,5,6,8,9,10,12,15,16
  !nplvar = 1;  plvars(1:nplvar) = (/2/)
  !nplvar = 1;  plvars(1:nplvar) = (/14/)
  !nplvar = 2;  plvars(1:nplvar) = (/2,3/)
  !nplvar = 4;  plvars(1:nplvar) = (/2,3,14,15/) !All masses
  !nplvar = 4;  plvars(1:nplvar) = (/2,3,6,7/) !Masses, spins
  !nplvar = 4;  plvars(1:nplvar) = (/4,5,8,9/) !tc, d, position
  !nplvar = 8;  plvars(1:nplvar) = (/4,5,8,9,10,11,12,13/) !All except masses, spins
  !nplvar = 9;  plvars(1:nplvar) = (/2,3,4,6,7,5,10,8,9/)
  nplvar = 12;  plvars(1:nplvar) = (/2,3,4,5, 6,7,8,9, 10,11,12,13/) !Default: all parameters
  !nplvar = 12;  plvars(1:nplvar) = (/2,3,4, 6,7,5, 8,9,10, 11,12,13/) !For poster (3x4 rather than 4x3)
  !nplvar = 14;  plvars(1:nplvar) = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
  !nplvar = 15;  plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  
  xpanels = 0 ! 0 - use default values, >0 have xpanels panels in the horizontal direction   \
  ypanels = 0 ! 0 - use default values, >0 have ypanels panels in the vertical direction     / Default values are use if either xpanels or ypanels = 0
  
  psclr = '/cps'
  if(colour.eq.0) psclr = '/ps '
  
  ncolours = 5; colours(1:ncolours)=(/4,2,3,6,5/) !Paper
  if(quality.eq.2) then !Beamer
     ncolours = 5
     colours(1:ncolours)=(/4,2,5,11,15/)
  end if
  if(colour.ne.1) then
     ncolours=2
     colours(1:ncolours)=(/14,15/)
  end if
  if(quality.eq.0.and.nchs.gt.5) then
     ncolours = 10
     colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  end if
  !Overrule
  ncolours = 10
  colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  !ncolours = 1
  !colours(1:ncolours)=(/6/)
  !defcolour = 2 !Red e.g. in case of 1 chain
  defcolour = colours(1)
  
  if(reverseread.ge.2) then !Reverse colours too
     do i=1,ncolours
        tempintarray(i) = colours(i)
     end do
     do i=1,ncolours
        colours(i) = tempintarray(ncolours-i+1) !Reverse colours too
     end do
  end if
  
  !Sort out implicit options:
  if(plot.eq.0) then
     pllogl = 0
     plchain = 0
     pljump = 0
     plsigacc = 0
     if(savepdf.eq.0) then
        plpdf1d = 0
        plpdf2d = 0
     end if
     plmovie = 0
  end if
  if(savepdf.eq.1) then
     nplvar = 14; plvars(1:nplvar) = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     wrapdata = 0
  end if
  !if(par1.lt.2) par1 = 2
  if(par1.lt.1) par1 = 1 !Include log(L)
  if(plsigacc.eq.1.or.plmovie.ge.1) rdsigacc = 1
  if(file.eq.1) combinechainplots = 0
  if(file.ge.1) update = 0
  if(plmovie.eq.1) update = 0
  if(plotsky.ge.1) plpdf2d = 1
  
  colournames(1:15) = (/'white','red','dark green','dark blue','very light blue','light purple','yellow','orange','light green','light blue-green','light blue','dark purple','red-purple','dark grey','light grey'/)
  if(file.ge.2) colournames(1) = 'black'

  
  
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  varnames(1:15) = (/'logL','Mc','eta','tc','log dl','spin','kappa','RA','sin dec','phase','sin thJo','phJo','alpha','M1','M2'/)
  pgvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
       'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','R.A. (rad)            ', &
       'sin dec.              ','\(2147)\dc\u (rad)    ','sin \(2134)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
       '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
  pgvarnss(1:15)  = (/'log L    ','M\dc\u ','\(2133)','t\dc\u','log d\dL\u','a\dspin\u','\(2136)','R.A.','sin dec.','\(2147)\dc\u', &
       'sin \(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  pgorigvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
       'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','R.A. (rad)            ', &
       'sin dec.              ','\(2147)\dc\u (rad)    ','sin \(2134)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
       '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
  !pgorigvarns(1:15)  = (/'log L    ','M\dc\u ','\(2133)','t\dc\u','log d\dL\u','a\dspin\u','\(2136)','R.A.','sin dec.','\(2147)\dc\u', &
  !     'sin \(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','rad','rad','','rad','','rad','rad','M\d\(2281)\u','M\d\(2281)\u'/)
  
  
  !nival = 1; ivals(1:nival) = (/0.9/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 2; ivals(1:nival) = (/0.9,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 3; ivals(1:nival) = (/0.683,0.9,0.997/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  nival = 3; ivals(1:nival) = (/0.90,0.95,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 4; ivals(1:nival) = (/0.6827,0.9,0.9545,0.9973/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 4; ivals(1:nival) = (/0.8,0.9,0.95,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  ivals(nival+1) = 1.
  
  
  
  
  
  j=0
  do i=1,nival
     if(abs(ivals(i)-ival0).lt.1.e-6) j=1
  end do
  if(j.eq.0) write(*,'(A44,F5.3,A35)')' !!! Error:  standard probability interval (',ival0,') is not present in array ivals !!!'

  npar = 13
  nchains0 = iargc()
  if(nchains0.lt.1) then
     write(*,'(A)')'  Syntax: plotspins <file1> [file2] ...'
     goto 9999
  end if
  if(prprogress.ge.1) then
     if(nchains0.gt.nchs) then
        write(*,'(A,I3,A)')' Too many chains, please increase nchs. Only',nchs,' files will be read.'
     else
        write(*,'(I5,A)')nchains0,' input files will be read.'
     end if
  end if
  nchains0 = min(nchains0,nchs)
  nchains = nchains0
  
  
  
  
  
  
  
  !*******************************************************************************************************************************
  !***   READ INPUT FILE(S)   ****************************************************************************************************
  !*******************************************************************************************************************************

101 continue
  narr = narr1
  allocate(dat(npar1,nchains,narr1))
  
  do ic = 1,nchains0
     if(reverseread.eq.0) then
        call getarg(ic,infile) !Read file name from the command-line arguments
     else
        call getarg(nchains0-ic+1,infile) !Read file name from the command-line arguments in reverse order
     end if
     infiles(ic) = infile
     
     open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io)
     if(io.ne.0) then
        write(*,'(A)')'File not found: '//trim(infile)//'. Quitting the programme.'
        goto 9999
     end if
     rewind(10)
     
     !if(update.ne.1) 
     if(prprogress.ge.1) write(*,'(A19,I3,A,I3,A)')'Reading input file',ic,':  '//trim(infile)//'...    Using colour',colours(mod(ic-1,ncolours)+1),': '//trim(colournames(colours(mod(ic-1,ncolours)+1)))
     
     !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
     !Read the headers
     read(10,*,end=199,err=199)bla
     if(prruninfo.eq.1.and.update.eq.0) write(*,'(6x,A10,A12,A8,A22,A8)')'niter','nburn','seed','null likelihood','ndet'
     read(10,'(I10,I12,I8,F22.10,I8)')niter(ic),nburn0(ic),seed(ic),nullh,ndet(ic)
     if(prruninfo.eq.1.and.update.eq.0) write(*,'(6x,I10,I12,I8,F22.10,I8)')niter(ic),nburn0(ic),seed(ic),nullh,ndet(ic)
     read(10,*,end=199,err=199)bla
     !if(prprogress.ge.1.and.update.eq.0) write(*,*)''
     if(prruninfo.eq.1.and.update.eq.0) write(*,'(A16,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
     do i=1,ndet(ic)
        read(10,'(A16,F18.8,4F12.2,F22.8,F17.7,3I14)') detname(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        if(detname(ic,i).eq.'         Hanford') detnr(ic,i) = 1
        if(detname(ic,i).eq.'      Livingston') detnr(ic,i) = 2
        if(detname(ic,i).eq.'            Pisa') detnr(ic,i) = 3
        if(prruninfo.eq.1.and.update.eq.0) write(*,'(A16,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(ic,i),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
     end do
     !if(prprogress.ge.1.and.update.eq.0) write(*,*)''
     read(10,*,end=199,err=199)bla
     read(10,*,end=199,err=199)bla
     
     !Read the data
     if(rdsigacc.eq.1) then
        !do i=1,narr1
        i=1
        do while(i.le.narr1)
           !read(10,*,end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=2,npar0)
           read(10,*,end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=2,3),  t,sig(4,ic,i),acc(4,ic,i),  (dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=5,npar0)
           is(ic,i) = real(i1)
           if(ic.eq.1.and.i.eq.1) t0 = dble(floor(t/10.d0)*10)
           dat(4,ic,i) = real(t - t0)
           if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
              do j=1,thin-1
                 read(10,*,end=199,err=198)bla
              end do
           end if
           i = i+1
        end do !i
     else !Read 40% quicker
        !do i=1,narr1
        i=1
        do while(i.le.narr1)
           !read(10,'(I12,F21.10,2(F17.10,17x),F22.10,17x,9(F17.10,17x))',end=199,err=198)i1,dat(1:npar0,ic,i)
           read(10,'(I12,F21.10,2(F17.10,17x),F22.10,17x,9(F17.10,17x))',end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),j=2,3),  t,  (dat(j,ic,i),j=5,npar0)
           is(ic,i) = real(i1)
           if(ic.eq.1.and.i.eq.1) t0 = dble(floor(t/10.d0)*10)
           dat(4,ic,i) = real(t - t0)
           !if(abs(dat(7,ic,i)).gt.0.995) i=i-1
           if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
              do j=1,thin-1
                 read(10,*,end=199,err=198)bla
              end do
           end if
           i = i+1
        end do !i
     end if
     write(*,'(A,$)')'   *** WARNING ***   Not all lines in this file were read    '
     goto 199
198  write(*,'(A,I7)')' Read error in line',i
     i = i-1
199  close(10)
     ntot(ic) = i-1
     n(ic) = ntot(ic) !n can be changed in rearranging chains, ntot can't
     !write(*,'(A,I,A,ES7.1,A,I,A,ES7.1,A,I4)')' Total number of lines: ',n(ic),'  (',real(n(ic)),'),   last iteration: ',nint(is(ic,n(ic))),'  (',is(ic,n(ic)),'),   thinning in file: ',nint(is(ic,n(ic))/real(n(ic)*max(thin,1)))
     if(prruninfo.eq.1.and.update.eq.0) write(*,*)''
  end do !do ic = 1,nchains0
  
  narr = maxval(n(1:nchains0))
  
  
  !*** It seems nburn is now the iteration number, but becomes the line number in the next bit, after isburn has taken that role
  do ic=1,nchains0
     if(nburn(ic).le.0) nburn(ic) = nburn0(ic)
     if(nburn(ic).ge.nint(is(ic,n(ic)))) then
        !print*,nburn(ic),nint(is(ic,n(ic)))
        if(nburn0(ic).ge.nint(is(ic,n(ic)))) then
           write(*,'(A,I3)')'   *** WARNING ***  Nburn larger than Nchain, setting nburn to 10% for chain',ic
           nburn(ic) = nint(is(ic,n(ic))*0.1)
        else
           nburn(ic) = nburn0(ic)
        end if
     end if
  end do
  
  do ic=1,nchains0
     isburn(ic) = real(nburn(ic))
     do i=1,ntot(ic)
        if(is(ic,i).le.isburn(ic)) nburn(ic) = i   !isburn is the true iteration number at which the burnin ends
     end do
     if(prruninfo.eq.1.and.update.ne.1) write(*,'(I3,A,I9,A,ES7.1,A,I9,A,ES7.1,A,I4,A,ES9.1,A,ES9.1,A,I4)')ic,':  Lines:',n(ic),'  (',real(n(ic)),'),  Iterations:',nint(is(ic,n(ic))),'  (',is(ic,n(ic)),'),  thinning in file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),'     Burn-in:  line:',real(nburn(ic)),', iteration:',isburn(ic),', total thinning:',nint(isburn(ic)/real(nburn(ic)))
  end do
  avgtotthin = sum(isburn(1:nchains0))/real(sum(nburn(1:nchains0))) !Total thinning, averaged over all chains
  
  
  
  if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Changing some variables...   '
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  !Calculate the masses from Mch and eta:
  do ic=1,nchains0
     do i=1,ntot(ic)
        dvar = dsqrt(0.25d0-dat(3,ic,i))
        dvar1 = (0.5d0-dvar)/(0.5d0+dvar)
        dvar2 = dvar1**0.6d0
        dat(14,ic,i) = dat(2,ic,i) * ((1.d0+dvar1)**0.2d0 / dvar2)
        dat(15,ic,i) = dat(2,ic,i) * ((1.d0+1.d0/dvar1)**0.2d0 * dvar2)
     end do
  end do
  npar = 15
  if(par2.eq.0) par2 = npar
  
  acc = acc*0.25  !Transfom back to the actual acceptance rate
    
  !*** Put plot data in pldat, startval and jumps
  allocate(pldat(nchains,npar,narr))
  jumps = 0.
  if(prinitial.ne.0) write(*,'(8x,A20,12A10)')varnames(1:13)
  do ic=1,nchains
     pldat(ic,1:npar,1:n(ic)) = real(dat(1:npar,ic,1:n(ic))) !Note the change of order of indices!!!  Pldat has the same structure as alldat.  Pldat contains all data that's in dat (also before the burnin), but in single precision.  Use it to plot log(L), jumps, chains, etc., but not for PDF creation.
     startval(ic,1:npar,1:2) = real(dat(1:npar,ic,1:2)) !True value and starting value
     jumps(ic,1:npar,2:n(ic)) = real(dat(1:npar,ic,2:n(ic)) -  dat(1:npar,ic,1:n(ic)-1))
     if(prinitial.ne.0) then 
        write(*,'(A8,F20.5,12F10.5)')' True:  ',startval(ic,1:13,1)
        if(abs((sum(startval(ic,1:13,1))-sum(startval(ic,1:13,2)))/sum(startval(ic,1:13,1))).gt.1.e-10) then
           write(*,'(A8,F20.5,12F10.5)')' Start: ',startval(ic,1:13,2)
           write(*,'(A8,F20.5,12F10.5)')' Diff:  ',abs(startval(ic,1:13,1)-startval(ic,1:13,2))!/abs(startval(ic,1:13,1))
           write(*,*)''
        end if
     end if
  end do
  
  
  !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  Done.'
  if(prprogress.ge.1.and.update.eq.0) write(*,'(A,I12)')' t0:',nint(t0)
  
  
  !Construct output file name
  ic = 1
  outputname = ' '
  !Number of detectors
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,I1,A1)')trim(outputname),ndet(ic),'d'
     if(ic.gt.1 .and. (ndet(ic).ne.ndet(1) .or. sum(detnr(ic,1:ndet(ic))).ne.sum(detnr(1,1:ndet(1))))) write(outputname,'(A,I1,A1)')trim(outputname)//'-',ndet(ic),'d'
     if(ic.eq.1 .or. (ic.gt.1.and. (ndet(ic).ne.ndet(1) .or. sum(detnr(ic,1:ndet(ic))).ne.sum(detnr(1,1:ndet(1)))) )) then 
        do i=1,ndet(ic)
           write(outputname,'(A,I1)')trim(outputname),detnr(ic,i)
        end do
     end if
  end do
  !write(outputname,'(A,F4.2,A1,I3.3)')trim(outputname)//'_',startval(ic,6,1),'_',nint(acos(startval(ic,7,1))*r2d)
  !Spin magnitudes
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,F4.2)')trim(outputname)//'_',startval(ic,6,1)
     if(ic.gt.1.and.startval(ic,6,1).ne.startval(1,6,1)) write(outputname,'(A,F4.2)')trim(outputname)//'-',startval(ic,6,1)
  end do
  !Theta_SLs
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,I3.3)')trim(outputname)//'_',nint(acos(startval(ic,7,1))*r2d)
     if(ic.gt.1.and.startval(ic,7,1).ne.startval(1,7,1)) write(outputname,'(A,I3.3)')trim(outputname)//'-',nint(acos(startval(ic,7,1))*r2d)
  end do
  !print*,outputname
  
  
  
  !*** Put data in alldat
  if(mergechains.eq.1) then  !Merge chains, leave out burnin (then nchains = 1)
     allocate(alldat(1,npar1,nchains*narr))
     j = 1
     do ic=1,nchains
        do i=nburn(ic)+1,n(ic)
           alldat(1,1:npar,j) = real(dat(1:npar,ic,i))  !Note the change of order of indices!!!  Alldat has the same structure as pldat, but contains only info AFTER the burnin.
           j = j+1
        end do
     end do
     nchains = 1
     n(1) = j-1
     if(prprogress.ge.1) write(*,'(A,I8,A,ES7.1,A)')' Data points in combined chains: ',n(1),'  (',real(n(1)),')'
  else
     allocate(alldat(nchains,npar1,narr))
     do ic=1,nchains
        alldat(ic,1:npar,1:n(ic)-nburn(ic)) = real(dat(1:npar,ic,nburn(ic)+1:n(ic)))  !Note the change of order of indices!!!  Alldat has the same structure as pldat, but contains only info AFTER the burnin.
        n(ic) = n(ic)-nburn(ic)
     end do
     if(prprogress.ge.1) write(*,'(A,I8)')' Datapoints in combined chains: ',sum(n(1:nchains))
  end if
  
  maxdots = 10000 !~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much
  !if(file.ge.2) maxdots = 10000 !Smaller max for eps,pdf to reduce file size (?)
  !if(maxval(n(1:nchains)).gt.maxdots) chainpli = nint(real(maxval(n(1:nchains)))/real(maxdots).)  !Change the number of points plotted in chains
  !if(maxval(n(1:nchains)).gt.maxdots .and. (file.ge.2.and.chainpli.eq.0)) then  !Change the number of points plotted in chains, for eps,pdf, to reduce file size
  if(sum(n(1:nchains)).gt.maxdots.and.chainpli.eq.0) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
     chainpli = nint(real(sum(n(1:nchains)))/real(maxdots))
     write(*,'(A,I5,A,I5,A,I6,A)')' Plotting every',chainpli,'-th state in likelihood, chains, jumps, etc. plots.  Average total thinning is',nint(avgtotthin),', for these plots it is',nint(avgtotthin*chainpli),'.'
  else
     chainpli = 1
     write(*,'(A,I5,A)')' Plotting *every* state in likelihood, chains, jumps, etc. plots.  Average total thinning remais',nint(avgtotthin),' for these plots.'
  end if
  
  
  
  
  
  ! **********************************************************************************************************************************
  ! ***  DO STATISTICS   *************************************************************************************************************
  ! **********************************************************************************************************************************
  
  
  !Sort all data and find the 90% interval limits for the wrapable parameters
  if(prprogress.ge.1) write(*,*)''
  shift = 0.
  wrap = 0
  do ic=1,nchains
     ival = ival0
     index = 0
     if(prprogress.ge.1.and.wrapdata.ge.1) write(*,'(A)')' Wrapping data...'
     do p=par1,par2
        if(wrapdata.eq.0 .or. (p.ne.8.and.p.ne.10.and.p.ne.12.and.p.ne.13)) then
           call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
           index(p,1:n(ic)) = index1(1:n(ic))
           cycle
        end if
        !if(p.ne.8.and.p.ne.10.and.p.ne.12.and.p.ne.13) cycle
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
        !Make sure data are between 0 and 2pi to start with:
        do i=1,n(ic)
           alldat(ic,p,i) = rev2pi(alldat(ic,p,i))
        end do
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
        index(p,1:n(ic)) = index1(1:n(ic))
        minrange = 1.e30
        
        
        
        do i=1,n(ic)
           x1 = alldat(ic,p,index(p,i))
           x2 = alldat(ic,p,index(p,mod(i+nint(n(ic)*ival)-1,n(ic))+1))
           range = mod(x2 - x1 + real(20*pi),real(tpi))
           if(range.lt.minrange) then
              minrange = range
              y1 = x1
              y2 = x2
              !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*ival),n(ic)),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end if
           !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*ival),n(ic)),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        if(y1.gt.y2) then
           wrap(ic,p) = 1
           centre = mod(centre + pi, tpi) !Then distribution peaks close to 0/2pi, shift centre by pi
        end if
        !if(p.eq.8) write(*,'(3I6,7F10.5)')ic,p,wrap(ic,p),y1,y2,minrange,centre

        !See whether there's a gap in the data  WHY is necessary, should it work like this???
        if(wrap(ic,p).eq.0 .and. 1.eq.2) then
           !ymin = minval(alldat(ic,p,1:n(ic)))
           !ymax = maxval(alldat(ic,p,1:n(ic)))
           i0 = -1
           maxgap = -1.e30
           !write(*,'(2I3,I8)')ic,p,n(ic)
           do i=1,n(ic)-1
              x1 = alldat(ic,p,index(p,i))
              x2 = alldat(ic,p,index(p,i+1))
              !write(*,'(2I3,2I8,4F10.5)')ic,p,i,i0,x1,x2,x2-x2,maxgap
              if(x2-x1.gt.maxgap) then
                 maxgap = x2-x1
                 i0 = i
              end if
           end do !i
           x1 = alldat(ic,p,index(p,i0))
           x2 = alldat(ic,p,index(p,i0+1))
           !if(maxgap.gt.2*tpi/sqrt(real(n(ic)))) then
           if(maxgap.gt.0.1) then 
              x0 = (x1+x2)/2.
              !write(*,'(10F10.5)')x1,x2,(x1+x2)/2.,maxgap,ymin,ymax,centre,minrange,y1,y2
              !if(y1.lt.y2.and.(x0.lt.y1.or.x0.gt.y2)) wrap(ic,p) = 1  !If centre of max gap is outside 90% range  WHY???
              !if(y1.gt.y2.and.(x0.gt.y2.and.x0.lt.y1)) wrap(ic,p) = 1
           end if
        end if
        !if(p.eq.8) write(*,'(3I6,9F10.5)')ic,p,wrap(ic,p),y1,y2,x1/pi*12,x2/pi*12,x0/pi*12

        !Now, wrap around anticentre
        shift(ic,p) = 0.
        if(wrap(ic,p).eq.1) shift(ic,p) = tpi - mod(centre + pi, tpi)
        alldat(ic,p,1:n(ic)) = mod(alldat(ic,p,1:n(ic))+shift(ic,p),tpi)-shift(ic,p)
        pldat(ic,p,1:ntot(ic)) = mod(pldat(ic,p,1:ntot(ic))+shift(ic,p),tpi)-shift(ic,p) !Original data
        y1 = mod(y1+shift(ic,p),tpi)-shift(ic,p)
        y2 = mod(y2+shift(ic,p),tpi)-shift(ic,p)
        centre = mod(centre+shift(ic,p),tpi)-shift(ic,p)
        minrange = y2-y1
        !call rindexx(n(ic),alldat(ic,p,1:n(ic)),index(p,1:n(ic)))  !Re-sort
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))  !Re-sort
        index(p,1:n(ic)) = index1(1:n(ic))
        !if(p.eq.8) write(*,'(I3,A8,4x,6F10.5,I4)')ic,varnames(p),y1,y2,minrange,centre,minval(alldat(ic,p,1:n(ic))),maxval(alldat(ic,p,1:n(ic))),wrap(ic,p)
        if(abs(abs(minval(alldat(ic,p,1:n(ic)))-maxval(alldat(ic,p,1:n(ic))))-2*pi).lt.1.e-3) wrap(ic,p)=1 !If centre is around pi, still needs to be flagged 'wrap' to plot PDF
     end do !p







     !Do statistics
     if(prprogress.ge.1) write(*,'(A)')' Calculating statistics...'
     do p=par1,par2
        !Determine the median
        if(mod(n(ic),2).eq.0) medians(p) = 0.5*(alldat(ic,p,index(p,n(ic)/2)) + alldat(ic,p,index(p,n(ic)/2+1)))
        if(mod(n(ic),2).eq.1) medians(p) = alldat(ic,p,index(p,(n(ic)+1)/2))
        
        !Mean:
        mean(p) = sum(alldat(ic,p,1:n(ic)))/real(n(ic))
        
        !Variances, etc:
        var1(p)=0.; var2(p)=0.; absvar1(p)=0.; absvar2(p)=0.; stdev1(p)=0.; stdev2(p)=0.
        do i=1,n(ic)
           var1(p) = var1(p) + (alldat(ic,p,i) - medians(p))**2
           var2(p) = var2(p) + (alldat(ic,p,i) - mean(p))**2
           absvar1(p) = absvar1(p) + abs(alldat(ic,p,i) - medians(p))
           absvar2(p) = absvar2(p) + abs(alldat(ic,p,i) - mean(p))
           stdev1(p) = stdev1(p) + (alldat(ic,p,i) - medians(p))*(alldat(ic,p,i) - medians(p))
           stdev2(p) = stdev2(p) + (alldat(ic,p,i) - mean(p))*(alldat(ic,p,i) - mean(p))
        end do
        
        absvar1(p) = absvar1(p)/real(n(ic))
        absvar2(p) = absvar2(p)/real(n(ic))
        stdev1(p)  = sqrt(stdev1(p)/real(n(ic)-1))
        stdev2(p)  = sqrt(stdev2(p)/real(n(ic)-1))
        
        !Save statistics:
        nstat = 6
        stats(ic,p,1) = medians(p)
        stats(ic,p,2) = mean(p)
        stats(ic,p,3) = absvar1(p)
        stats(ic,p,4) = absvar2(p)
        stats(ic,p,5) = stdev1(p)
        stats(ic,p,6) = stdev2(p)
     end do
     
     
     !Correlations:
     if(prcorr.gt.0.or.savestats.gt.0) then
        write(*,'(A)')' Calculating correlations...'
        do p1=par1,par2
           !do p2=par1,par2
           do p2=p1,par2
              corrs(p1,p2) = 0.
              do i=1,n(ic)
                 !corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - medians(p1))*(alldat(ic,p2,i) - medians(p2))
                 corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - mean(p1))*(alldat(ic,p2,i) - mean(p2)) !Hardly differs from median method
              end do
              !corrs(p1,p2) = corrs(p1,p2) / (stdev1(p1)*stdev1(p2)*(n(ic)-1))
              corrs(p1,p2) = corrs(p1,p2) / (stdev2(p1)*stdev2(p2)*(n(ic)-1))
           end do !p2
        end do !p1
     end if
     
     
     !Autocorrelations:
     if(placorr.gt.0) then
        write(*,'(A)')' Calculating autocorrelations...'
        j1 = placorr/100 !Step size to get 100 autocorrelations per var
        do p=par1,par2
           acorrs(ic,p,:) = 0.
           !do j=1,ntot(ic)-1
           !do j=1,min(placorr,ntot(ic)-1)
           do j=0,min(100,ntot(ic)-1)
              do i=1,ntot(ic)-j*j1
                 acorrs(ic,p,j) = acorrs(ic,p,j) + (pldat(ic,p,i) - medians(p))*(pldat(ic,p,i+j*j1) - medians(p))
                 !acorrs(p,j) = acorrs(ic,p,j) + (pldat(ic,p,i) - mean(p))*(pldat(ic,p,i+j*j1) - mean(p))
              end do
              !if(j.eq.0) write(*,'(3I6,A,4F9.3)')j1,j,j*j1,'  '//varnames(p),acorrs(ic,0,j),acorrs(ic,p,j),(stdev1(p)*stdev1(p)*(ntot(ic)-j*j1)),acorrs(ic,p,0)
              acorrs(ic,0,j) = real(j*j1)
              acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev1(p)*stdev1(p)*(ntot(ic)-j*j1))
              !acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev2(p)*stdev2(p)*(ntot(ic)-j*j1))
              !if(j.eq.0) write(*,'(3I6,A,4F9.3)')j1,j,j*j1,'  '//varnames(p),acorrs(ic,0,j),acorrs(ic,p,j),(stdev1(p)*stdev1(p)*(ntot(ic)-j*j1)),acorrs(ic,p,0)
           end do !j
           !write(*,*)''
        end do !p
     end if
     
     
     !Determine interval ranges
     if(prprogress.ge.1) write(*,'(A29,$)')' Determining interval levels: '
     c0 = 0
     do c=1,nival
        ival = ivals(c)
        if(abs(ival-ival0).lt.0.001) c0 = c
        if(c.ne.c0.and.prival.eq.0.and.savestats.eq.0) cycle
        
        if(prprogress.ge.1) write(*,'(F8.4,$)')ival
        do p=par1,par2
           minrange = 1.e30
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
           do i=1,floor(n(ic)*(1.-ival))
              x1 = alldat(ic,p,index(p,i))
              x2 = alldat(ic,p,index(p,i+floor(n(ic)*ival)))
              range = abs(x2 - x1)
              !range = x2 - x1
              if(range.lt.minrange) then
                 minrange = range
                 y1 = x1
                 y2 = x2
              end if
              !write(*,'(I6,7F10.5)')i,x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)

           !Save ranges:
           nr = 4
           ranges(ic,c,p,1) = y1
           ranges(ic,c,p,2) = y2
           ranges(ic,c,p,3) = (y1+y2)/2.
           ranges(ic,c,p,4) = y2-y1
           ranges(ic,c,p,5) = ranges(ic,c,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,c,p,5) = ranges(ic,c,p,4)/ranges(ic,c,p,3)
        end do !p
     end do !c
     if(prprogress.ge.1) write(*,'(A34,F8.4)')'.  Standard probability interval: ',ival0
     
     
     
     
     
     !Change variables
     !Columns in alldat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl,   6:longi, 7:sinlati:, 8:phase, 9:spin,   10:kappa,     11:sinthJ0, 12:phiJ0, 13:alpha
     if(changevar.eq.1) then
        if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' Changing some variables...   '
        do p=par1,par2
           if(p.eq.5) then
              alldat(ic,p,1:n(ic)) = exp(alldat(ic,p,1:n(ic)))     !logD -> Distance
              !startval(ic,p,1:2) = exp(startval(ic,p,1:2))
              startval(1:nchains0,p,1:2) = exp(startval(1:nchains0,p,1:2))
              stats(ic,p,1:nstat) = exp(stats(ic,p,1:nstat))
              ranges(ic,1:nival,p,1:nr) = exp(ranges(ic,1:nival,p,1:nr))
              !print*,ic,p
           end if
           if(p.eq.9.or.p.eq.11) then
              alldat(ic,p,1:n(ic)) = asin(alldat(ic,p,1:n(ic)))*r2d
              !startval(ic,p,1:2) = asin(startval(ic,p,1:2))*r2d
              startval(1:nchains0,p,1:2) = asin(startval(1:nchains0,p,1:2))*r2d
              stats(ic,p,1:nstat) = asin(stats(ic,p,1:nstat))*r2d
              ranges(ic,1:nival,p,1:nr) = asin(ranges(ic,1:nival,p,1:nr))*r2d
           end if
           if(p.eq.7) then
              alldat(ic,p,1:n(ic)) = acos(alldat(ic,p,1:n(ic)))*r2d
              !startval(ic,p,1:2) = acos(startval(ic,p,1:2))*r2d
              startval(1:nchains0,p,1:2) = acos(startval(1:nchains0,p,1:2))*r2d
              stats(ic,p,1:nstat) = acos(stats(ic,p,1:nstat))*r2d
              ranges(ic,1:nival,p,1:nr) = acos(ranges(ic,1:nival,p,1:nr))*r2d
              do c=1,nival
                 y1 = ranges(ic,c,p,2)
                 ranges(ic,c,p,2) = ranges(ic,c,p,1)  !acos is monotonously decreasing
                 ranges(ic,c,p,1) = y1
              end do
           end if
           if(p.eq.8) then
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*r2h
              !startval(ic,p,1:2) = startval(ic,p,1:2)*r2h
              startval(1:nchains0,p,1:2) = startval(1:nchains0,p,1:2)*r2h
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*r2h
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*r2h
           end if
           if(p.eq.10.or.p.eq.12.or.p.eq.13) then
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*r2d
              !startval(ic,p,1:2) = startval(ic,p,1:2)*r2d
              startval(1:nchains0,p,1:2) = startval(1:nchains0,p,1:2)*r2d
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*r2d
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*r2d
           end if
           ranges(ic,1:nival,p,3) = 0.5*(ranges(ic,1:nival,p,1) + ranges(ic,1:nival,p,2))
           ranges(ic,1:nival,p,4) = ranges(ic,1:nival,p,2) - ranges(ic,1:nival,p,1)
           ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,5)/ranges(ic,1:nival,p,3)
        end do !p
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:thJ0, 12:phiJ0, 13:alpha
        varnames(1:15) = (/'logL','Mc','eta','tc','dl','spin','th_SL','RA','Dec','phase','thJo','phJo','alpha','M1','M2'/)
        pgvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',  &
             'd\dL\u (Mpc)          ','a\dspin\u             ','\(2134)\dSL\u(\(2218))','R.A. (h)              ','Dec. (\(2218))        ', &
             '\(2147)\dc\u (\(2218))','\(2134)\dJ\u (\(2218))','\(2147)\dJ\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u) ','M\d2\u(M\d\(2281)\u)  '/)
        pgvarnss(1:15)  = (/'log L','M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u',  &
             '\(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
        !Include units
        pgvarnss(1:15)  = (/'log L','M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
             '\(2134)\dJ0\u (\(2218))','\(2147)\dJ0\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
        !Units only
        pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','\(2218)','\uh\d','\(2218)','\(2218)','\(2218)','\(2218)','\(2218)','M\d\(2281)\u','M\d\(2281)\u'/)

        if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  Done.'
     end if !if(changevar.eq.1)
     
     
     !Find 100% range
     do p=par1,par2
        if(p.eq.1) cycle
        ranges(ic,nival+1,p,1) = minval(alldat(ic,p,1:n(ic)))
        ranges(ic,nival+1,p,2) = maxval(alldat(ic,p,1:n(ic)))
        ranges(ic,nival+1,p,3) = 0.5*(ranges(ic,nival+1,p,1) + ranges(ic,nival+1,p,2))
        ranges(ic,nival+1,p,4) = ranges(ic,nival+1,p,2) - ranges(ic,nival+1,p,1)
        ranges(ic,nival+1,p,5) = ranges(ic,nival+1,p,4)
        if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,nival+1,p,5) = ranges(ic,nival+1,p,5)/ranges(ic,nival+1,p,3)
        !write(*,'(I3,2F10.4)')p,ranges(ic,nival+1,p,1),ranges(ic,nival+1,p,2)
     end do
     
     
     
     
     
     !**********************************************************************************************
     !******   PRINT STATISTICS   ******************************************************************
     !**********************************************************************************************
     
     !Print statistics to screen
     o=6
     if(prstat.gt.0) then
        write(o,*)''
        c = c0
        write(o,'(A8,13A10,A4,A12,I3,A2)')'param.','model','median','mean','stdev1','stdev2','abvar1','abvar2',  &
             'rng_c','rng1','rng2','drng','d/drng','delta','ok?','result (',nint(ivals(c0)*100),'%)'
        do p=par1,par2
           if(stdev1(p).lt.1.d-20) cycle !Parameter was probably not fitted
           write(o,'(A8,13F10.4,$)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),stdev1(p),stdev2(p),absvar1(p),  &
                absvar2(p),ranges(ic,c,p,3),ranges(ic,c,p,1),ranges(ic,c,p,2),ranges(ic,c,p,4),  &
                !abs(startval(ic,p,1)-stats(ic,p,1))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt median
                2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt centre of range
           if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
              write(o,'(A4,$)')'y '
           else
              write(o,'(A4,$)')'*N*'
           end if
           write(o,'(F10.4,A3,F9.4)')ranges(ic,c,p,3),'+-',0.5*ranges(ic,c,p,4)
        end do
        write(o,*)''
     end if
     
     !Print correlations:
     if(prcorr.gt.0) then
        write(o,'(A8,$)')''
        do p=par1,par2
           if(stdev1(p).gt.1.d-20) write(o,'(A7,$)')trim(varnames(p))
        end do
        write(o,*)''
        do p1=par1,par2
           if(stdev1(p1).lt.1.d-20) cycle
           write(o,'(A8,$)')trim(varnames(p1))
           do p2=par1,par2
              if(stdev1(p2).lt.1.d-20) cycle
              if(abs(corrs(p1,p2)).gt.0.5) then 
              !if(abs(corrs(p1,p2)).gt.-0.5) then 
                 write(o,'(F7.2,$)')corrs(p1,p2)
              else
                 write(o,'(A7,$)')''
              end if
           end do
           write(o,'(A)')'   '//varnames(p1)
        end do
        write(o,*)''
     end if
     
     !Print intervals:
     if(prival.gt.0) then
        write(o,'(A20,A8,$)')'Interval:',''
        do c=1,nival+1
           write(o,'(F20.4,A9,$)')ivals(c),''
        end do
        write(o,*)''
        
        write(o,'(A8,2x,2A9,$)')'param.','model','median'
        do c=1,nival+1
           !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
           write(o,'(2x,3A9,$)')'centre','delta','in rnge'
        end do
        write(o,*)''
        do p=par1,par2
           if(stdev1(p).lt.1.d-20) cycle
           write(o,'(A8,2x,2F9.4,$)')varnames(p),startval(ic,p,1),stats(ic,p,1)
           do c=1,nival+1
              !write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-stats(ic,p,1))/ranges(ic,c,p,4) !Defined with median
              !write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
                 write(o,'(A3,$)')'y '
              else
                 write(o,'(A3,$)')'N*'
              end if
           end do
           write(o,*)''
        end do
     end if
     
     
  end do !ic
  
  
     
     
  
  !Write statistics to file
  if(savestats.ge.1) write(*,*)''
  if(savestats.ge.1.and.nchains.gt.1) write(*,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
  if(savestats.ge.1.and.nchains.eq.1) then
     ic = 1 !Use chain 1
     o = 20 !Output port
     open(unit=o, form='formatted', status='replace',file='bulk/'//trim(outputname)//'__statistics.dat')
     write(o,'(A)')trim(outputname)
     write(o,*)''
     !write(o,'(3(A,I))')'Npoints: ',n(ic),'  Nburn: ',nburn(ic),'  Ndet: ',ndet(ic)
     !write(o,*)''
     write(o,'(6x,A10,A12,A8,A22,A8)')'niter','nburn','seed','null likelihood','ndet'
     write(o,'(6x,I10,I12,I8,F22.10,I8)')n(ic),nint(isburn(ic)),seed(ic),nullh,ndet(ic)
     write(o,*)''
     write(o,'(A16,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
     do i=1,ndet(ic)
        write(o,'(A16,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(ic,i),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
     end do
     write(o,*)''
     
     write(o,'(A,I11)')' t0:',nint(t0)
     write(o,*)''
     
     !Print statistics
     write(o,'(A,2I3)')'Npar,ncol:',par2-par1+1,7
     write(o,'(A8,7A12)')'param.','model','median','mean','stdev1','stdev2','abvar1','abvar2'
     
     do p=par1,par2
        write(o,'(A8,7F12.6)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),stdev1(p),stdev2(p),absvar1(p),absvar2(p)
     end do
     write(o,*)''
     
     
     !Print correlations:
     write(o,'(A,I3)')'Npar:',par2-par1+1
     write(o,'(A9,$)')''
     do p=par1,par2
        write(o,'(A10,$)')trim(varnames(p))
     end do
     write(o,*)''
     do p1=par1,par2
        write(o,'(A9,$)')trim(varnames(p1))
        do p2=par1,par2
           write(o,'(F10.5,$)')corrs(p1,p2)
        end do
        write(o,'(A)')'   '//varnames(p1)
     end do
     write(o,*)''
     
     
     !Print intervals:
     write(o,'(A,I3)')'Nival:',nival
     write(o,'(A22,$)')'Interval:'
     do c=1,nival+1
        write(o,'(F20.5,A14,$)')ivals(c),''
     end do
     write(o,*)''
     
     write(o,'(A8,2x,$)')'param.'
     do c=1,nival+1
        !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
        write(o,'(2x,2A12,A8,$)')'centre','delta','in rnge'
     end do
     write(o,*)''
     do p=par1,par2
        write(o,'(A8,2x,$)')varnames(p)
        do c=1,nival+1
           !write(o,'(2x,2F11.6,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
           write(o,'(2x,2F12.6,F6.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
           if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
              write(o,'(A2,$)')'y'
           else
              write(o,'(A2,$)')'N'
           end if
        end do
        write(o,*)''
     end do
     
     close(o) !Statistics output file
     if(savestats.eq.2) i = system('a2ps -1rf7 bulk/'//trim(outputname)//'__statistics.dat -o bulk/'//trim(outputname)//'__statistics.ps')
     write(*,*)''
     if(savestats.eq.1) write(*,'(A)')' Statistics saved in '//trim(outputname)//'__statistics.dat'
     if(savestats.eq.2) write(*,'(A)')' Statistics saved in '//trim(outputname)//'__statistics.dat,ps'
  end if !if(savestats.ge.1.and.nchains.eq.1) then
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !Check convergence for multiple chains. this works only for fixed chain length, so take the min N
  if(nchains0.gt.1 .and. (prconv.ge.1.or.savestats.ge.1)) then
     chmean = 1.d-30
     totmean = 1.d-30
     nn = minval(ntot(1:nchains0))/2
     
     do p=par1,par2
        do ic=1,nchains0
           do i=nn+1,2*nn
              chmean(ic,p) = chmean(ic,p) + dat(p,ic,i) !Can't use pldat, because it may get wrapped above
              totmean(p) = totmean(p) + dat(p,ic,i)
           end do
        end do
     end do
     chmean = chmean/dble(nn)
     totmean = totmean/dble(nn*nchains0)
     
     chvar = 1.d-30
     chvar1 = 1.d-30
     totvar = 1.d-30
     do p=par1,par2
        do ic=1,nchains0
           do i=nn+1,2*nn
              dx = (dat(p,ic,i) - chmean(ic,p))**2 !Can't use pldat, because it may get wrapped above
              chvar(p) = chvar(p) + dx
              chvar1(ic,p) = chvar1(ic,p) + dx !Keep track of the variance per chain
           end do
           totvar(p) = totvar(p) + (chmean(ic,p) - totmean(p))**2
           chvar1(ic,p) = chvar1(ic,p)/dble(nn-1)
        end do
        chvar(p) = chvar(p)/dble(nchains0*(nn-1))
        totvar(p) = totvar(p)/dble(nchains0-1)
        
        !rhat(p) = ( dble(nn-1)/dble(nn) * chvar(p)  +  totvar(p) * (1.d0 + 1.d0/dble(nchains0)) ) / chvar(p)
        rhat(p) = min( dble(nn-1)/dble(nn)  +  totvar(p)/chvar(p) * (1.d0 + 1.d0/dble(nchains0)), 99.d0)
     end do
     
     if(prconv.ge.1) then
        write(*,*)''
        if(prconv.ge.2) write(*,'(A,I7,A)')'  Convergence criterion for',nn,' parameters:'
        write(*,'(18x,20A12)')varnames(par1:min(par2,13))
        if(prconv.ge.2) then
           write(*,'(A)')'  Means:'
           do ic=1,nchains0
              write(*,'(I16,A2,20F12.6)')ic,': ',chmean(ic,par1:min(par2,13))
           end do
           write(*,'(A18,20F12.6)')'           Total: ',totmean(par1:min(par2,13))
           
           write(*,*)''
           write(*,'(A)')'  Variances:'
        end if !if(prconv.ge.2)
     end if !if(prconv.ge.1)
     do ic=1,nchains0
        !write(*,'(I16,A2,20F12.6)')ic,': ',chvar1(ic,par1:min(par2,13))
        !if(chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,3).lt.0.5*chvar(3).and.chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,2).lt.0.5*chvar(2)) then
        lowvar = 0
        highvar = 0
        totrelvar = 1.d0
        ntotrelvar = 0
        do p=par1,min(par2,13)
           if(abs(chvar1(ic,p)).gt.1.e-30) then  !The parameters that were not fitted for have a variance of 0
              if(chvar1(ic,p).lt.0.5*chvar(p)) lowvar(p) = 1  !Too (?) low variance, mark it
              if(chvar1(ic,p).gt.2*chvar(p))  highvar(p) = 1  !Too (?) high variance, mark it
              totrelvar = totrelvar * chvar1(ic,p)/chvar(p) !Take geometric mean
              ntotrelvar = ntotrelvar + 1
           end if
        end do
        nlowvar = lowvar(2)+lowvar(3)+lowvar(6)+lowvar(7)  !Sum of 2 masses and 2 spin parameters
        nhighvar = highvar(2)+highvar(3)+highvar(6)+highvar(7)  !Sum of 2 masses and 2 spin parameters
        !totrelvar = totrelvar**(1.d0/dble(abs(min(par2,13)-par1+1))) !Take geometric mean of (the variance of each chain, relative to the total variance)
        totrelvar = totrelvar**(1.d0/dble(ntotrelvar)) !Take geometric mean of (the variance of each chain, relative to the total variance)
        if(prconv.ge.2) then
           ch = ' '
           if(nlowvar.eq.4) ch = '*'
           if(nhighvar.eq.4) ch = '#'
           write(*,'(I7,A3,$)')ic,': '//ch
           ch = ' '
           if(totrelvar.lt.0.5) ch = '*'
           if(totrelvar.gt.2.0) ch = '#'
           write(*,'(F8.3,A1,$)')totrelvar,ch
           do p=par1,min(par2,13)
              ch = ' '
              if(lowvar(p).eq.1) ch = '*'
              if(highvar(p).eq.1) ch = '#'
              write(*,'(F11.7,A1,$)')chvar1(ic,p),ch
           end do
           write(*,*)''
        end if !if(prconv.ge.2)
     end do
     if(prconv.ge.2) then
        write(*,'(A9,9x,20F12.7)')'  Total: ',chvar(par1:min(par2,13))
        
        write(*,*)''
        write(*,'(A)')'  Variances:'
        write(*,'(A18,20ES12.4)')'   Within chains: ',chvar(par1:min(par2,13))
        write(*,'(A18,20ES12.4)')'  Between chains: ',totvar(par1:min(par2,13))
     end if
     
     if(prconv.ge.1) write(*,'(A18,20F12.6)')'     Convergence: ',rhat(par1:min(par2,13)),sum(rhat(par1:min(par2,13)))/dble(min(par2,13)-par1+1)
     !write(*,*)''
  end if
  
  
  
  
  !Test: get mean and stdev for log(L)
  if(1.eq.1) then
     write(*,*)''
     
     nn = minval(ntot(1:nchains0))/2
     nlogl1 = nn+1
     nlogl2 = 2*nn
     nn = abs(nlogl2-nlogl1)
     
     write(*,'(A,I7,A)')'  Convergence criterion for',nn,' parameters:'
     write(*,'(16x,16x,20A12)')'Mean','Stddev','M-S','M+S','KS d','KS prob'
     
     do ic=1,nchains0
        nn1 = ntot(ic)/20
        ksn2 = 0
        ksd = 1.
        ksprob = 0.
        do nlogl1 = 1,ntot(ic),nn1
           nlogl2 = min(nlogl1+nn1,ntot(ic))
           nn = abs(nlogl2-nlogl1)+1
           if(nn.lt.nn1) cycle
           
           chmean = 1.d-30
           totmean = 1.d-30
           chvar = 1.d-30
           chvar1 = 1.d-30
           totvar = 1.d-30
           
           p=1
           do i=nlogl1,nlogl2
              chmean(ic,p) = chmean(ic,p) + dat(p,ic,i) !Can't use pldat, because it may get wrapped above
              totmean(p) = totmean(p) + dat(p,ic,i)
           end do
           chmean = chmean/dble(nn)
           totmean = totmean/dble(nn*nchains0)
           
           do i=nlogl1,nlogl2
              dx = (dat(p,ic,i) - chmean(ic,p))**2 !Can't use pldat, because it may get wrapped above
              chvar(p) = chvar(p) + dx
              chvar1(ic,p) = chvar1(ic,p) + dx !Keep track of the variance per chain
           end do
           totvar(p) = totvar(p) + (chmean(ic,p) - totmean(p))**2
           chvar1(ic,p) = chvar1(ic,p)/dble(nn-1)
           chvar(p) = chvar(p)/dble(nchains0*(nn-1))
           totvar(p) = totvar(p)/dble(nchains0-1)
           
           !write(*,'(I16,2I8,20F12.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob
           !!print*,ic,p,nlogl1,nlogl2,nn
           ksdat1(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn1   = nn
           !!call kstwo(data1,n1,data2,n2,d,prob)
           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
           !
           !ksdat2 = ksdat1
           ksdat2(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn2 = ksn1
           
           write(*,'(I16,2I8,20F12.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob,dlog10(ksprob+1.d-100)
        end do
        write(*,*)''
     end do
     write(*,*)''
     
     !KS test
     call pgbegin(1,'21/xs',1,1)
     call pgpap(scrsz,scrrat)
     call pgsch(1.5)
     call pgsvp(0.07,0.99,0.10,0.96)
     call pgswin(0.,real(maxval(ntot(1:nchains0))),-100.,0.)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     call pgmtxt('B',2.2,0.5,0.5,'i')
     call pgmtxt('L',1.8,0.5,0.5,'log(d\dKS\u)')
     
     do ic=1,nchains0
        call pgsci(colours(mod(ic-1,ncolours)+1))
        nn = ntot(ic)/10
        ksd = 1.
        ksprob = 0.
        do nlogl1 = 1,ntot(ic),nn
           nlogl2 = nlogl1+nn-1
           if(nlogl2.gt.ntot(ic)) cycle
           
           ksn1   = nn
           ksdat1(1:ksn1) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn2   = ntot(ic)-nlogl1+1
           ksdat2(1:ksn2) = dble(dat(p,ic,nlogl1:ntot(ic)))
           
           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
           !write(*,'(I16,2I8,20F12.6)')ic,nlogl1,nlogl2,ksd,ksprob,dlog10(ksprob+1.d-100)
           call pgpoint(1,real(nlogl1+nlogl2)/2.,real(dlog10(ksprob+1.d-100)),2)
        end do
        write(*,*)''
     end do
     
     call pgend
     write(*,*)''
     
  end if
  
  
  
  
  
  
  
  
  
  !Change the original chain data
  if(changevar.eq.1) then
     do ic=1,nchains0
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:thJ0, 12:phiJ0, 13:alpha
        !if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')'Changing some variables...   '
        do p=par1,par2
           if(p.eq.5) pldat(ic,p,1:ntot(ic)) = exp(pldat(ic,p,1:ntot(ic)))
           if(p.eq.9.or.p.eq.11) pldat(ic,p,1:ntot(ic)) = asin(pldat(ic,p,1:ntot(ic)))*r2d
           if(p.eq.7) pldat(ic,p,1:ntot(ic)) = acos(pldat(ic,p,1:ntot(ic)))*r2d
           if(p.eq.8) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2h
           if(p.eq.10.or.p.eq.12.or.p.eq.13) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2d
        end do !p
     end do
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  Done.'
  end if !if(changevar.eq.1)




  deallocate(dat)



  ! **********************************************************************************************************************************
  ! ***  CREATE PLOTS   **************************************************************************************************************
  ! **********************************************************************************************************************************
  
  
  if(prprogress.ge.1) write(*,*)''
  if(combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.eq.1)) then
     io = pgopen('chaininfo.eps'//trim(psclr))
     call pginitl(colour,file)
  end if
    
  
  
  
  !***********************************************************************************************************************************      
  !Plot likelihood chain
  if(pllogl.eq.1) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting chain likelihood...'
     if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('logL.ppm/ppm')
        if(file.ge.2) io = pgopen('logL.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file)
     !call pgsubp(1,2)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.94) !To make room for title
     
     ic = 1
     p=1
     !call pgpage
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nchains0
        xmin = 0.
        xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
        ymin = min(ymin,real(minval(pldat(ic,p,10:ntot(ic))))) !Use pldat iso alldat, to also take into account burn-in, ic should loop to nchains0 iso nchains
        ymax = max(ymax,real(maxval(pldat(ic,p,10:ntot(ic)))))
     end do
     ic = 1
     p = 1
     if(ymax.gt.0.) then !This is log(L)-log(Lo) (which we started saving later on), so that nullh=Lo=0
        ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),0.)
        ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),0.)
     else !This is log(L), so that Lo = nullh
        ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),nullh)
        ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),nullh)
     end if
     dx = abs(xmax-xmin)*0.01
     dy = abs(ymax-ymin)*0.05
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(abs(startval(1,1,1)-startval(1,1,2))/abs(startval(1,1,1)).gt.1.e-10) then
        call pgsls(4)
        call pgbox('',0.0,0,'G',0.0,0)
        call pgsls(1)
     end if
     
     do ic=1,nchains0
        call pgsci(defcolour)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !if(thin.le.1) then
        !   do i=1,ntot(ic),chainpli
        !      call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        !   end do
        !else
        !   call pgpoint(ntot(ic),is(ic,1:ntot(ic)),pldat(ic,p,1:ntot(ic)),1)
        !end if
        !do i=1,ntot(ic),chainpli
        do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
     end do
     
     do ic=1,nchains0
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))
        call pgsci(6)
        !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
        call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(4)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))
        call pgsci(6)
        call pgline(2,(/-1.e20,1.e20/),real((/nullh,nullh/)))
     end do
     call pgsci(1)
     call pgsls(1)
     call pgmtxt('T',0.5,0.1,0.1,trim(pgvarns(p)))
     
     if(quality.eq.0) then
        !call pgsubp(1,1)
        !call pgsvp(0.,1.,0.,1.)
        !call pgswin(-1.,1.,-1.,1.)
     
        !call pgsch(sch*0.8)
        call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  !Print title
        !call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf logL.eps  -o bulk/'//trim(outputname)//'__logL.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f logL.eps bulk/'//trim(outputname)//'__logL.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8  -resize '//trim(bmpxpix)//' logL.ppm  bulk/'//trim(outputname)//'__logL.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f logL.ppm')
        !if(i.ne.0) write(*,'(A)')'  Error removing file',i
     end if
  end if !if(pllogl.eq.1) then
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot chain for each parameter
  if(plchain.eq.1) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting parameter chains...'
     if(file.eq.0) then
        io = pgopen('13/xs')
        sch = 1.5
        lw = 1
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('chains.ppm/ppm')
        if(file.ge.2) io = pgopen('chains.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.le.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.gt.12) then
              sch = sch*1.5
              lw = 1
           end if
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           end if
           !if(nplvar.le.6) then
           !   !sch = sch*4
           !   !lw = 4
           !end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
              sch = sch*2.7
              lw = 3
           else
              !sch = sch*1.25
              sch = sch*1.5
              lw = 2
           end if
        end if
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgscf(2)
     !if(file.eq.1) call pgsch(1.5)
     
     call pgsch(sch)
     call pgslw(lw)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file)
        xmin = 0.
        !xmax = real(maxval(ntot(1:nchains0)))
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           ymin = min(ymin,minval(pldat(ic,p,1:ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,p,1:ntot(ic))))
        end do
        
        if(p.eq.8) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,24.)
           !ymin = max(rev24(ymin),0.)
           !ymax = min(rev24(ymax),24.)
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        end if
        if(p.eq.10.or.p.eq.12.or.p.eq.13) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,360.)
           !write(*,'(I5,2F10.5,$)')p,ymin,ymax
           !ymin = max(rev360(ymin),0.)
           !ymax = min(rev360(ymax),360.)
           !write(*,'(2F10.5)')ymin,ymax
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
        end if
        dx = abs(xmax-xmin)*0.01
        dy = abs(ymax-ymin)*0.05
        if(dx.eq.0) then
           xmin = 0.5*xmin
           xmax = 2*xmax
        end if
        if(dy.eq.0) then
           ymin = 0.5*ymin
           ymax = 2*ymax
        end if
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        call pgsch(1.)
        call pgslw(1)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !print*,ic,colours(mod(ic-1,ncolours)+1)
           !if(thin.le.1) then
           !   do i=1,ntot(ic),chainpli
           !      call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
           !   end do
           !else
           !   call pgpoint(ntot(ic),is(ic,1:ntot(ic)),pldat(ic,p,1:ntot(ic)),1)
           !end if
           !if(thin.le.1) then
           if(chainsymbol.eq.0) then !Plot lines rather than symbols
              call pgline(ntot(ic),is(ic,1:ntot(ic)),pldat(ic,p,1:ntot(ic)))
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 ply = pldat(ic,p,i)
                 if(p.eq.8) ply = rev24(ply)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
                 !call pgpoint(1,is(ic,i),ply,1) !Plot small dots
                 call pgpoint(1,is(ic,i),ply,chainsymbol) !Plot symbols
              end do
           end if
        end do
        call pgsch(sch)
        call pgslw(lw)

        do ic=1,nchains0
           call pgsls(2)
           call pgsci(6)
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))   !Burnin phase
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           call pgsci(1)
           
           !Plot true values in chains
           if(pltrue.eq.1) then
              if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !True value
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                 call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
                 if(p.eq.8) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
                 end if
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
                    !print*,p,plx,plx-360.,plx+360.
                 end if
              end if
           end if
           
           !Plot starting values in chains
           if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)) .gt. 1.e-10) then
              call pgsls(4)
              if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              plx = startval(ic,p,2) !Initial value
              if(p.eq.8) plx = rev24(plx)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
              call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
              if(p.eq.8) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
              end if
              if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
              end if
           end if
        end do
        
        call pgsci(1)
        call pgsls(1)
        write(string,'(F6.3)')rhat(p)
        !call pgmtxt('T',1.,0.,0.,'Chain: '//trim(pgvarns(p)))
        call pgmtxt('T',-1.,0.,0.,' '//trim(pgvarns(p)))
        if(nchains0.gt.1.and.prconv.ge.1) call pgmtxt('T',1.,1.,1.,'Conv: '//trim(string))
     end do !do j=1,nplvar
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf chains.eps  -o bulk/'//trim(outputname)//'__chains.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f chains.eps bulk/'//trim(outputname)//'__chains.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8 -resize '//trim(bmpxpix)//' chains.ppm  bulk/'//trim(outputname)//'__chains.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f chains.ppm')
     end if
  end if !if(plchain.eq.1)








  !***********************************************************************************************************************************            
  !Plot jump sizes
  if(pljump.eq.1) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting jump sizes...'
     if(file.eq.0) then
        io = pgopen('18/xs')
        sch=1.5
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('jumps.ppm/ppm')
        if(file.ge.2) io = pgopen('jumps.eps'//trim(psclr))
        sch=1.2
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) sch=1.5
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) sch=1.5
     if(quality.eq.0) then !Draft
        !sch = sch*1.75
        sch = sch*1.4
        lw = 2
        call pgsvp(0.1,0.95,0.06,0.87) !To make room for title
     end if
     
     call pgsch(sch)
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file)
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(logjump.eq.0) then
              ymin = min(ymin,minval(jumps(ic,p,2:ntot(ic))))
              ymax = max(ymax,maxval(jumps(ic,p,2:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(logjump.eq.1.and.jumps(ic,p,i).gt.1.e-20) then
                 ymin = min(ymin,log10(abs(jumps(ic,p,i))))
                 ymax = max(ymax,log10(abs(jumps(ic,p,i))))
              end if
              ymin = -6.
              ymax = 1.
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(logjump.eq.0) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(logjump.eq.1) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log
        
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(logjump.eq.0) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(logjump.eq.0) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),jumps(ic,p,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(abs(jumps(ic,p,1:ntot(ic)))+1.e-30),1)
           !   end if
           !end if
           if(logjump.eq.0) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',1.,0.5,0.5,'Jumps: '//trim(pgorigvarns(p)))
        call pgmtxt('T',-1.2,0.05,0.0,trim(pgorigvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf jumps.eps  -o bulk/'//trim(outputname)//'__jumps.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f jumps.eps bulk/'//trim(outputname)//'__jumps.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8 -resize '//trim(bmpxpix)//' jumps.ppm  bulk/'//trim(outputname)//'__jumps.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f jumps.ppm')
     end if
  end if !if(pljump.eq.1)









  !***********************************************************************************************************************************            
  !Plot sigma values ('jump proposal width')
  if(plsigacc.eq.1) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting sigma...'
     if(file.eq.0) then
        io = pgopen('16/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('sigs.ppm/ppm')
        if(file.ge.2) io = pgopen('sigs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file)

     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(logsig.eq.0) then
              ymin = min(ymin,minval(sig(p,ic,10:ntot(ic))))
              ymax = max(ymax,maxval(sig(p,ic,10:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(logsig.eq.1.and.sig(p,ic,i).gt.1.e-20) then
                 ymin = min(ymin,log10(sig(p,ic,i)))
                 ymax = max(ymax,log10(sig(p,ic,i)))
              end if
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(logsig.eq.0) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(logsig.eq.1) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(logsig.eq.0) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),sig(p,ic,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(logsig.eq.0) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),sig(p,ic,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(sig(p,ic,1:ntot(ic))+1.e-30),1)
           !   end if
           !end if
           if(logsig.eq.0) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),sig(p,ic,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Sigma: '//trim(pgvarns(p)))
        !call pgmtxt('T',1.,0.5,0.5,'log Sigma: '//trim(pgvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf sigs.eps  -o bulk/'//trim(outputname)//'__sigs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f sigs.eps bulk/'//trim(outputname)//'__sigs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8 -resize '//trim(bmpxpix)//' sigs.ppm  bulk/'//trim(outputname)//'__sigs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f sigs.ppm')
     end if
  end if !if(plsigacc.eq.1)









  !***********************************************************************************************************************************      
  !Plot acceptance rates
  if(plsigacc.eq.1) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting acceptance rates...'
     if(file.eq.0) then
        io = pgopen('17/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('accs.ppm/ppm')
        if(file.ge.2) io = pgopen('accs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file)

     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           !xmax = max(xmax,real(ntot(ic)))
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           do i=1,ntot(ic)
              if(acc(p,ic,i).gt.1.e-10 .and. acc(p,ic,i).lt.1.-1.e-10) then
                 n0 = i
                 exit
              end if
           end do
           n0 = n0+10
           ymin = min(ymin,minval(acc(p,ic,n0:ntot(ic))))
           ymax = max(ymax,maxval(acc(p,ic,n0:ntot(ic))))
           dy = abs(ymax-ymin)*0.05
        end do

        call pgsci(1)
        call pgsls(1)
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)


        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Acceptance: '//trim(pgvarns(p)))

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   do i=1,ntot(ic),chainpli
           !      call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           !   end do
           !else
           !   call pgpoint(ntot(ic),is(ic,1:ntot(ic)),acc(p,ic,1:ntot(ic)),1)
           !end if
           !do i=1,ntot(ic),chainpli
           do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           end do
        end do
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     !if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf accs.eps  -o bulk/'//trim(outputname)//'__accs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f accs.eps bulk/'//trim(outputname)//'__accs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8 -resize '//trim(bmpxpix)//' accs.ppm  bulk/'//trim(outputname)//'__accs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f accs.ppm')
     end if
  end if !if(plsigacc.eq.1)
  
  
  
  if(file.ge.1.and.combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.eq.1)) then
     call pgend
     if(file.eq.3) i = system('eps2pdf chaininfo.eps -o bulk/'//trim(outputname)//'__chaininfo.pdf  >& /dev/null')
     i = system('mv -f chaininfo.eps bulk/'//trim(outputname)//'__chaininfo.eps')
  end if
  !***********************************************************************************************************************************        
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot autocorrelations for each parameter
  if(placorr.gt.0) then
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting autocorrelations...'
     if(file.eq.0) then
        io = pgopen('19/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('acorrs.ppm/ppm')
        if(file.ge.2) io = pgopen('acorrs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file)
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file)
        xmin = 0.
        xmin = minval(acorrs(1,0,0:100))
        xmax = maxval(acorrs(1,0,0:100))
        dx = abs(xmax-xmin)*0.01
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           ymin = min(ymin,minval(acorrs(ic,p,0:100)))
           ymax = max(ymax,maxval(acorrs(ic,p,0:100)))
        end do
        dy = abs(ymax-ymin)*0.05
        !write(*,'(I3,5F10.2)')p,xmin,xmax,ymin,ymax,acorrs(1,0,100)
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        do ic=1,nchains0
           call pgsci(defcolour)
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           do i=1,100!placorr!ntot(ic),chainpli
              call pgpoint(1,acorrs(ic,0,i),acorrs(ic,p,i),1)
           end do
        end do
        
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.,0./))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Autocorrelation: '//trim(pgvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf acorrs.eps  -o bulk/'//trim(outputname)//'__acorrs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f acorrs.eps bulk/'//trim(outputname)//'__acorrs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -depth 8 -resize '//trim(bmpxpix)//' acorrs.ppm  bulk/'//trim(outputname)//'__acorrs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f acorrs.ppm')
     end if
  end if !if(placorrs.gt.0)
  !***********************************************************************************************************************************      







  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot pdfs (1d)
  if(plpdf1d.eq.1) then
     if(plot.eq.0.and.savepdf.eq.1) write(*,'(A,$)')' Saving 1D pdfs...   '
     if(plot.eq.1) then
        if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' Plotting 1D pdfs...   '
        if(file.eq.0) then
           io = pgopen('14/xs')
           sch = 1.5
           lw = 1
        end if
        if(file.ge.1) then
           if(file.eq.1) io = pgopen('pdfs.ppm/ppm')
           if(file.ge.2) io = pgopen('pdfs.eps'//trim(psclr))
           lw = 3
           if(nplvar.ge.10) lw = 2
           if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
           sch = 2.5
           if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
           if(quality.eq.0) then !Draft
              sch = sch*1.75
              lw = 2
           end if
           if(quality.eq.1) then !Paper
              if(nplvar.eq.12) then
                 sch = sch*1.75
                 lw = 2
              else
                 sch = sch*1.25
                 lw = 1
              end if
           end if
           if(quality.eq.2) then !Talk
              if(nplvar.le.12) then
                 sch = sch*2
                 lw = 2
              else
                 sch = sch*1.5
                 lw = 1
              end if
           end if
           if(quality.eq.3) then !Poster
              if(nplvar.eq.12.and.file.ge.2) then
                 sch = sch*2.7
                 lw = 3
              else
                 !sch = sch*1.25
                 !lw = 1
                 sch = sch*1.5
                 lw = 2
              end if
           end if
        end if
        if(io.le.0) then
           write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
           goto 9999
        end if
        if(file.eq.0) call pgpap(scrsz,scrrat)
        if(file.eq.1) call pgpap(bmpsz,bmprat)
        if(file.ge.2.and.quality.eq.3.and.nplvar.eq.12) call pgpap(10.6,0.925)
        if(file.ge.2) call pgscf(2)
        !call pgscr(3,0.,0.5,0.)
        !call pginitl(colour,file)
        call pgslw(lw)
        call pgsch(sch)
        call pgsfs(fillpdf)
        
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        
        if(xpanels*ypanels.lt.1) then
           if(nplvar.eq.2) call pgsubp(2,1)
           if(nplvar.eq.3) call pgsubp(3,1)
           if(nplvar.eq.4) call pgsubp(2,2)
           if(nplvar.eq.6) call pgsubp(3,2)
           if(nplvar.eq.8) call pgsubp(4,2)
           if(nplvar.eq.9) call pgsubp(3,3)
           if(nplvar.eq.10) call pgsubp(5,2)
           if(nplvar.eq.12) call pgsubp(4,3)
           if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
           if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
           if(nplvar.eq.16) call pgsubp(4,4)
        else
           call pgsubp(xpanels,ypanels)
        end if
     end if !if(plot.eq.1)
     
     !Save 1D PDF data
     if(savepdf.eq.1) then
        open(unit=30,action='write',form='formatted',status='replace',file='bulk/'//trim(outputname)//'__pdf1d.dat')
        write(30,'(3I6,T100,A)')nplvar,nchains,nbin,'Total number of plot variables, total number of chains, number of bins'
     end if
     
     !do p=par1,par2
     do j=1,nplvar
        p = plvars(j)
        if(plot.eq.1) then
           call pgpage
           if(j.eq.1) call pginitl(colour,file)
        end if

        !Set x-ranges, bin the data and get y-ranges
        xmin = 1.e30
        xmax = -1.e30
        do ic=1,nchains
           xmin = min(xmin,minval(alldat(ic,p,1:n(ic))))
           xmax = max(xmax,maxval(alldat(ic,p,1:n(ic))))
        end do
        dx = xmax - xmin
        !dx = max(xmax - xmin,1.e-30)
        
        do ic=1,nchains
           x(ic,1:n(ic)) = alldat(ic,p,1:n(ic))
           xmin1 = minval(alldat(ic,p,1:n(ic)))
           xmax1 = maxval(alldat(ic,p,1:n(ic)))
           
           call bindata(n(ic),x(ic,1:n(ic)),1,nbin,xmin1,xmax1,xbin1,ybin1)
           
           !Weigh with likelihood.  I should probably do something like this at the start, to get updated ranges etc.
           !y(ic,1:n(ic)) = alldat(ic,1,1:n(ic))
           !call bindata1(n(ic),x(ic,1:n(ic)),y(ic,1:n(ic)),1,nbin,xmin1,xmax1,xbin1,ybin1) !Measure the amount of likelihood in each bin
           
           !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:d_l, 6:spin, 7:theta_SL, 8: RA, 9:dec,10:phase, 11:thetaJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
           if(p.eq.5.or.p.eq.7.or.p.eq.9.or.p.eq.11) then  !Do something about chains 'sticking to the wall'
              if(ybin1(1).gt.ybin1(2)) ybin1(1)=0.
              if(ybin1(nbin).gt.ybin1(nbin-1)) ybin1(nbin)=0.
           end if

           
           !Normalise the SURFACE, not the height (because of different bin size)
           norm = 0.
           do i=1,nbin+1
              norm = norm + ybin1(i)
           end do
           norm = norm*(xmax1-xmin1)
           ybin1 = ybin1/norm
           
           !Smoothen
           ybin2 = ybin1
           if(smooth.gt.1) then
              !i0 = nbin/10
              !print*,nbin/10,nint(min(max(real(nbin)/real(smooth),1.0),real(nbin)/2.))
              !i0 = nint(min(max(real(nbin)/real(smooth),1.0),real(nbin)/2.))
              i0 = min(max(smooth,1),floor(real(nbin)/2.))
              do i=1+i0,nbin+1-i0
                 coefs1(1:2*i0+1) = ybin1(i-i0:i+i0)
                 call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
                 do i1=1,i0+1
                    coefs(i0-i1+2) = coefs1(i1)
                 end do
                 do i1 = i0+2,2*i0+1
                    coefs(3*i0+3-i1) = coefs1(i1)
                 end do
                 ybin2(i) = 0.
                 do i1=1,2*i0+1
                    ybin2(i) = ybin2(i) + coefs(i1) * ybin1(i+i1-i0-1)
                 end do
              end do
              ybin1 = ybin2
           end if !if(smooth.gt.1)
           xbin(ic,1:nbin+1) = xbin1(1:nbin+1)
           ybin(ic,1:nbin+1) = ybin1(1:nbin+1)
           
           !Save binned data
           if(savepdf.eq.1) then
              write(30,'(3I6,T100,A)')ic,p,wrap(ic,p),'Chain number, variable number, and wrap'
              write(30,'(2ES15.7,T100,A)')startval(ic,p,1:2),'True and starting value'
              write(30,'(6ES15.7,T100,A)')stats(ic,p,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p,1:5),'Ranges: lower,upper limit, centre, width, relative width'
              write(30,'(2ES15.7,T100,A)')xmin1,xmax1,'Xmin and Xmax of PDF'
              do i=1,nbin+1
                 write(30,'(2ES15.7)')xbin1(i),ybin1(i)
              end do
           end if
        end do !ic
        
        if(plot.eq.1) then
           xmin = xmin - 0.1*dx
           xmax = xmax + 0.1*dx
           ymin = 0.
           ymax = -1.e30
           do ic=1,nchains
              !ymax = max(ymax,maxval(ybin(ic,1:nbin+1)))
              do i=1,nbin+1
                 if(ybin(ic,i).gt.ymax) then
                    ymax = ybin(ic,i)
                    xpeak = xbin(ic,i)
                 end if
              end do
           end do
           ymax = ymax*1.1
           if(dx.eq.0) then
              xmin = 0.5*xmin
              xmax = 2*xmax
           end if
           
           call pgsch(sch)
           call pgswin(xmin,xmax,ymin,ymax)
           if(abs(dx).lt.1.e-30) then !So that the programme doesn't hang if a parameter is kept constant
              xbin = 0.
              ybin = 0.
           end if
           
           call pgsci(1)
           if(file.ge.2) call pgslw(lw)
           !Plot 1D PDF
           do ic=1,nchains
              if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              xbin1(1:nbin+1) = xbin(ic,1:nbin+1)
              ybin1(1:nbin+1) = ybin(ic,1:nbin+1)
              if(wrap(ic,p).eq.0) then
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin+2,(/xbin1(1),xbin1(1:nbin+1)/),(/0.,ybin1(1:nbin+1)/))
                 !Plot pdf contour
                 !if(nchains.eq.1) call pgsci(1)
                 !call pgsci(1)
                 if(fillpdf.eq.1) call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin+1,xbin1(1:nbin+1),ybin1(1:nbin+1)) !:nbin) ?
                 
                 !Fix the loose ends
                 call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/))
                 call pgline(2,(/xbin1(nbin+1),xbin1(nbin+1)/),(/ybin1(nbin+1),0./))
              else !If parameter is wrapped
                 plshift = real(2*pi)
                 if(changevar.eq.1) plshift = 360.
                 if(changevar.eq.1.and.p.eq.8) plshift = 24.  !RA in hours
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin+3,(/xbin1(1),xbin1(1:nbin),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin),ybin1(1),0./))
                 
                 !Plot pdf contour
                 !call pgsci(1)
                 !if(fillpdf.ne.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 if(fillpdf.eq.1) call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin,xbin1(1:nbin),ybin1(1:nbin))
                 
                 !Plot dotted lines outside the pdf for wrapped periodic variables
                 call pgsls(4)
                 !if(file.ge.2) call pgslw(2)
                 call pgline(nbin+1,(/xbin1(1:nbin)-plshift,xbin1(1)/),(/ybin1(1:nbin),ybin1(1)/))
                 call pgline(nbin,xbin1+plshift,ybin1)
                 
                 !Fix the loose end
                 call pgsls(1)
                 if(file.ge.2) call pgslw(lw)
                 call pgline(2,(/xbin1(nbin),xbin1(1)+plshift/),(/ybin1(nbin),ybin1(1)/))
              end if
           end do !ic
           !Plot lines again over surface of overlapping distributions
           if(nchains.gt.1.and.fillpdf.eq.1) then
              call pgsls(4)
              do ic=1,nchains
                 call pgsci(1)
                 !call pgsci(colours(mod(ic-1,ncolours)+1))
                 xbin1(1:nbin+1) = xbin(ic,1:nbin+1)
                 ybin1(1:nbin+1) = ybin(ic,1:nbin+1)
                 if(wrap(ic,p).eq.0) then
                    call pgline(nbin+1,xbin1(1:nbin+1),ybin1(1:nbin+1))
                 else
                    call pgline(nbin,xbin1(1:nbin),ybin1(1:nbin))
                 end if
              end do
              call pgsls(4)
           end if
           
           
           
           
           !Plot median and model value
           call pgsch(sch)
           
           do ic=1,nchains
              !Draw white lines
              if(nchains.gt.1) then
                 call pgslw(lw)
                 call pgsls(1); call pgsci(0)
                 call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
                 call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
                 call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
                 if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
              end if
              
              call pgslw(lw+1)
              !Draw coloured lines over the white ones
              !Median
              if(plmedian.eq.1) then
                 call pgsls(2); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))
              end if
              
              !Plot ranges in PDF
              if(plrange.eq.1) then
                 call pgsls(4); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
                 !if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
              end if
              
              !Plot true value in PDF
              if(pltrue.eq.1) then !Plot true values
                 if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                    call pgsls(2); call pgsci(1)
                    plx = startval(ic,p,1)
                    if(p.eq.8) plx = rev24(plx)
                    if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                    call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !True value
                    if(p.eq.8) then
                       call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !True value
                    end if
                    if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                       call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !True value
                    end if
                 end if
              end if
              
              !Plot starting value in PDF
              !if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
              !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
              !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
              !end if
              
              call pgsls(1)
              call pgsci(1)
           end do !ic
           
           
           !Print median, model value and range widths
           call pgslw(lw)
           call pgsci(1)
           ic = 1
           !if(nplvar.lt.7.or.nplvar.eq.9) then  !Three or less columns
           if(quality.ne.2.and.quality.ne.3) then  !Not a talk/poster
              if(nplvar.lt.5) then
                 if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                    write(str,'(A,F7.3,A5,F7.3,A9,F6.2,A1)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))/startval(ic,p,1)*100,'%'
                         ' \(2030):',ranges(ic,c0,p,5)*100,'%'
                 else
                    write(str,'(A,F7.3,A5,F7.3,A9,F7.3)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))  
                         ' \(2030):',ranges(ic,c0,p,5)
                 end if
              else  !if nplvar>=5
                 if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                    !write(str,'(A,F7.3,A9,F6.2,A1)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),
                    write(str,'(A4,F7.3,A9,F6.2,A1)')'mdl:',startval(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))/startval(ic,p,1)*100,'%'  
                         ' \(2030):',ranges(ic,c0,p,5)*100,'%'
                 else
                    !write(str,'(A,F7.3,A9,F7.3)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),
                    write(str,'(A4,F7.3,A9,F7.3)')'mdl:',startval(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))  
                         ' \(2030):',ranges(ic,c0,p,5)
                 end if
                 call pgsch(sch*1.2)
                 call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              end if
           end if
           
           
           if(quality.eq.2.or.quality.eq.3) then  !Talk/poster
              !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
              !   write(str,'(A9,F6.2,A1)')' \(2030):',ranges(ic,c0,p,5)*100,'%'
              !else
              !   write(str,'(A9,F7.3)')' \(2030):',ranges(ic,c0,p,5)
              !end if
              x0 = ranges(ic,c0,p,5)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) x0 = x0*100
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F6.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F5.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F4.2)')x0
              if(x0.ge.1.and.x0.lt.9.95) write(str,'(F3.1)')x0
              if(x0.ge.9.95.and.x0.lt.99.5) write(str,'(I2)')nint(x0)
              if(x0.ge.99.5) write(str,'(I3)')nint(x0)
              write(str,'(A)')'\(2030): '//trim(str)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                 write(str,'(A)')trim(str)//'%'
              else
                 write(str,'(A)')trim(str)//trim(pgunits(p))
              end if
              call pgsch(sch*1.2)
              !call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              if(abs(xmin-xpeak).lt.abs(xmax-xpeak)) then !peak is at left, put varname at right
                 call pgptxt(xmax-0.05*dx,ymax*0.9,0.,1.,trim(pgvarnss(p)))
              else
                 call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              end if
           end if
           
           
           !if(nchains.gt.1) write(str,'(A,F7.3,A)')trim(pgvarns(p))//'  mdl: ',startval(ic,p,1),''
           !if(nchains.gt.1) write(str,'(A4,F7.3,A)')'mdl:',startval(ic,p,1),''
           !if(nchains.eq.2) write(str,'(A,3(A6,F7.3))')trim(pgvarns(p)),' mdl:',startval(ic,p,1),
           !                       ' med1:',stats(1,p,1),' med2:',stats(2,p,1)
           !if(nchains.eq.2.and.abs((startval(1,p,1)-startval(2,p,1))/startval(1,p,1)).gt.1.e-10)  &
           !                        !         write(str,'(A,A5,F7.3,A2,F7.3)')trim(pgvarns(p)),'mdl:',startval(1,p,1),', ',startval(2,p,1)  
           !     write(str,'(A4,F6.2,A2,F6.2)')'mdl:',startval(1,p,1),', ',startval(2,p,1)
           
           !Write the deltas of the two pdfs
           if(nchains.eq.2.) then
              write(str,'(A8)')'\(2030)'
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                 write(str1,'(A8,F6.2,A1)')'\(2030):',ranges(1,c0,p,5)*100.,'%'
                 write(str2,'(A8,F6.2,A1)')'\(2030):',ranges(2,c0,p,5)*100.,'%'
              else
                 write(str1,'(A8,F7.3)')'\(2030):',ranges(1,c0,p,5)
                 write(str2,'(A8,F7.3)')'\(2030):',ranges(2,c0,p,5)
              end if
           end if
           
           !call pgsch(sch*0.9)
           call pgsch(sch*1.1)
           if(prvalues.eq.1) then
              if(nchains.eq.2) then
                 call pgsci(colours(mod(0,ncolours)+1))
                 call pgmtxt('T',0.5,0.25,0.5,trim(str1))
                 call pgsci(colours(mod(1,ncolours)+1))
                 call pgmtxt('T',0.5,0.75,0.5,trim(str2))
              else
                 if(quality.eq.2.or.quality.eq.3) call pgsci(2)
                 !call pgmtxt('T',0.2,0.5,0.5,trim(str))
                 !call pgptxt(ranges(ic,c0,p,3),ymax,0.,0.5,trim(str)) !Align with centre of 90%-probability range
                 call pgptxt((xmin+xmax)/2.,ymax,0.,0.5,trim(str)) !Centre
                 !print*,trim(str)
                 !call pgarro(ranges(ic,c0,p,3),0.95*ymax,ranges(ic,c0,p,1),0.95*ymax)
                 !call pgarro(ranges(ic,c0,p,3),0.95*ymax,ranges(ic,c0,p,2),0.95*ymax)
                 call pgsci(2)
                 call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,2)/),(/0.99*ymax,0.99*ymax/))  !Plot line at top over 90%-probability width
                 call pgsci(1)
              end if
           end if
           
           call pgsci(1)
           call pgsch(sch)
           !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
           call pgbox('BNTS',0.0,0,'',0.0,0)
           !print*,sch,lw,trim(str)
        end if !if(plot.eq.1) 
     end do !p
     
     if(savepdf.eq.1) close(30)
     
     if(plot.eq.1) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        if(quality.eq.0) then
           !Remove also the pgsvp at the beginning of the plot
           string=' '
           do ic=1,nchains
              write(string,'(A,I7,A,I6)')trim(string)//'n:',ntot(ic),', nburn:',nburn(ic)
           end do
           call pgsch(sch*0.7)
           call pgmtxt('T',-0.7,0.5,0.5,trim(outputname)//'  '//trim(string))  !Print title
           call pgsch(sch)
        end if
        
        
        !Make sure gv auto-reloads on change
        !      if(file.ge.2) then
        !        call pgpage
        !        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        !      end if
        
        call pgend
        
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf pdfs.eps -o bulk/'//trim(outputname)//'__pdfs.pdf >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f pdfs.eps bulk/'//trim(outputname)//'__pdfs.eps')
        else if(file.eq.1) then
           i = system('convert -resize '//trim(bmpxpix)//' -depth 8 pdfs.ppm bulk/'//trim(outputname)//'__pdfs.png')
           if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           i = system('rm -f pdfs.ppm')
        end if
     end if !if(plot.eq.1)
     write(*,*)''   
  end if !if(plpdf1d.eq.1)












  !***********************************************************************************************************************************      
  if(plpdf2d.ge.1) then
     ic = 1 !Can only do one chain
     if(plot.eq.0.and.savepdf.eq.1) write(*,'(A)')' Saving 2D pdfs...    '
     if(plot.eq.1) then
        if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting 2D pdfs...    '
        if(file.eq.0) then
           io = pgopen('15/xs')
           lw = 1
           lw2 = 1
           sch = 1.5
        end if
        if(file.ge.1) then
           !if(file.eq.1) io = pgopen('pdf2d.ppm/ppm')
           if(file.ge.2) io = pgopen('pdf2d.eps'//trim(psclr))
           lw = 3
           lw2 = 2 !Font lw
           sch = 1.5
           if(quality.eq.3) then !Poster
              lw = 4
              lw2 = 3 !Font lw
              sch = 2.
           end if
        end if
        if(io.le.0) then
           write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
           goto 9999
        end if
        if(file.eq.0) call pgpap(scrsz,scrrat)
        !if(file.eq.1) call pgpap(bmpsz,bmprat)
        if(file.ge.2) call pgscf(2)
        if(file.ne.1) call pginitl(colour,file)
     end if !if(plot.eq.1)
     
     !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
     j1 = 2
     j2 = npar
     !j2 = 13 !Don't do M1,M2
     !j1 = 6
     !j2 = 7
     if(plotsky.eq.1) then
        j1 = 8
        j2 = 9
     end if
     
     if(savepdf.eq.1) then
        open(unit=30,action='write',form='formatted',status='replace',file='bulk/'//trim(outputname)//'__pdf2d.dat')
        write(30,'(5I6,T100,A)')j1,j2,1,nbin2dx,nbin2dy,'Plot variable 1,2, total number of chains, number of bins x,y'
     end if
     
     
     do p1=j1,j2-1
        do p2=p1+1,j2
     !do p1=3,3
        !do p2=6,6
           
           !Skip some 2d pdfs to save time:
           !if(p1.eq.4.or.p1.eq.5.or.p1.eq.10.or.p1.eq.11.or.p1.eq.12.or.p1.eq.13) cycle
           !if(p2.eq.4.or.p2.eq.5.or.p2.eq.10.or.p2.eq.11.or.p2.eq.12.or.p2.eq.13) cycle
           
           !if(p1.eq.2.and.p2.ne.3 .or. p2.eq.2.or.p1.eq.3) cycle !Mc only with eta
           !if(p1.eq.6.and.p2.ne.7 .or. p2.eq.6.or.p1.eq.7) cycle !a_spin only with theta
           !if(p1.eq.8.and.p2.ne.9 .or. p2.eq.8.or.p1.eq.9) cycle !a_spin only with theta
           !if(p1.eq.14.and.p2.ne.15 .or. p2.eq.14.or.p1.eq.15) cycle !M1 only with M2
           
           !if(p1.ne.5.and.p1.ne.8.and.p1.ne.9) cycle  !Only position and distance
           !if(p2.ne.5.and.p2.ne.8.and.p2.ne.9) cycle  !Only position and distance
           
           plotthis = 0  !Determine to plot or save this combination of j1/j2 or p1/p2
           !if(p1.eq.2.and.p2.eq.3) plotthis = 1  !Mc-eta
           !if(p1.eq.6.and.p2.eq.7) plotthis = 1  !a-theta
           if(p1.eq.8.and.p2.eq.9) plotthis = 1  !RA-dec
           !if(p1.eq.11.and.p2.eq.12) plotthis = 1  !theta/phi_Jo
           !if(p1.eq.14.and.p2.eq.15) plotthis = 1  !M1-M2
           
           
           if(plotthis.eq.0) cycle
           
           !print*,p1,p2
           write(*,'(A)')'   PDF 2D:  '//trim(varnames(p1))//'-'//trim(varnames(p2))
           
           if(plot.eq.1) then
              if(file.eq.1) then
                 io = pgopen('pdf2d.ppm/ppm')
                 call pgpap(bmpsz,bmprat)
                 call pginitl(colour,file)
              end if
              !call pgscr(3,0.,0.5,0.)
              call pgsch(sch)
           end if
           
           xmin = minval(alldat(ic,p1,1:n(ic)))
           xmax = maxval(alldat(ic,p1,1:n(ic)))
           ymin = minval(alldat(ic,p2,1:n(ic)))
           ymax = maxval(alldat(ic,p2,1:n(ic)))
           dx = xmax - xmin
           dy = ymax - ymin
           !print*,xmin,xmax
           !print*,ymin,ymax
           !print*,ic

           xx(1:n(ic)) = alldat(ic,p1,1:n(ic)) !Parameter 1
           yy(1:n(ic)) = alldat(ic,p2,1:n(ic)) !Parameter 2
           zz(1:n(ic)) = alldat(ic,1,1:n(ic))   !Likelihood

           xmin = xmin - 0.05*dx
           xmax = xmax + 0.05*dx
           ymin = ymin - 0.05*dy
           ymax = ymax + 0.05*dy
           
           
           if(plot.eq.1.and.plotsky.eq.1) then !Plot a cute sky map
              xmax = 18.2!14.
              xmin = 14.8!20.
              ymin = 20.!0.
              ymax = 50.!70.
              rat = 0.75
              call pgpap(11.,rat)
              dx = xmax - xmin
              dy = ymax - ymin
              if(abs(dx)*15.lt.dy/rat) then !Expand x
                 dx = dy/(15*rat)
                 a = (xmin+xmax)*0.5
                 xmin = a - 0.5*dx
                 xmax = a + 0.5*dx
              end if
              if(abs(dx)*15.gt.dy/rat) then !Expand y
                 dy = abs(dx)*rat*15
                 a = (ymin+ymax)*0.5
                 ymin = a - 0.5*dy
                 ymax = a + 0.5*dy
              end if
           end if
           
           
           call bindata2d(n(ic),xx(1:n(ic)),yy(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Count number of chain elements in each bin
           
           !Weigh with likelihood
           !call bindata2da(n(ic),xx(1:n(ic)),yy(1:n(ic)),zz(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Measure amount of likelihood in each bin
           
           !z = max(0.,log10(z + 1.e-30))
           
           if(p1.eq.8.and.p2.eq.9) then !Swap RA boundaries
              a = xmin
              xmin = xmax
              xmax = a
              dx = -dx
           end if
           
           
           z = z/(maxval(z)+1.e-30)
           
           if(plot.eq.1.and.plotsky.eq.1.and.file.ge.2) z = 1. - z !Invert grey scales
           
           if(plot.eq.1) then
              call pgsch(sch)
              !call pgsvp(0.12,0.95,0.12,0.95)
              call pgsvp(0.08*sch,0.95,0.08*sch,0.95)
              call pgswin(xmin,xmax,ymin,ymax)
              if(plotsky.eq.1.and.file.ge.2) then !Need dark background
                 !call pgsvp(0.,1.,0.,1.)
                 !call pgswin(0.,1.,0.,1.)
                 call pgsci(1)
                 call pgrect(xmin,xmax,ymin,ymax)
                 !call pgsci(0)
              end if
              
              if(plpdf2d.eq.1.or.plpdf2d.eq.2) call pggray(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,1.,0.,tr)
              
              
              !Plot stars (over the grey scales, but underneath contours, lines, etc)
              if(plotsky.eq.1) then
                 call pgswin(xmin*15,xmax*15,ymin,ymax) !Map works in degrees
                 call plotthesky(xmin*15,xmax*15,ymin,ymax)
                 call pgswin(xmin,xmax,ymin,ymax)
              end if
              call pgsci(1)
           end if !if(plot.eq.1)
           
           if((plpdf2d.eq.1.or.plpdf2d.eq.3) .and. plot.eq.1) then
              do i=1,11
                 cont(i) = 0.01 + 2*real(i-1)/10.
                 if(plotsky.eq.1) cont(i) = 1.-cont(i)
              end do
           
              call pgsls(1)
              if(plotsky.eq.0) then !First in black
                 call pgslw(2*lw)
                 call pgsci(0)
                 call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont,4,tr)
              end if
              call pgslw(lw)
              call pgsci(1)
              if(plotsky.eq.1) call pgsci(7)
              call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont,4,tr)
           end if
           
           
           !Save binned data
           if(savepdf.eq.1) then
              write(30,'(3I6,T100,A)')ic,p1,p2,'Chain number and variable number 1,2'
              write(30,'(2ES15.7,T100,A)')startval(ic,p1,1:2),'True and starting value p1'
              write(30,'(2ES15.7,T100,A)')startval(ic,p2,1:2),'True and starting value p2'
              write(30,'(6ES15.7,T100,A)')stats(ic,p1,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
              write(30,'(6ES15.7,T100,A)')stats(ic,p2,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p1,1:5),'Ranges: lower,upper limit, centre, width, relative width for p1'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p2,1:5),'Ranges: lower,upper limit, centre, width, relative width for p2'
              write(30,'(4ES15.7,T100,A)')xmin,xmax,ymin,ymax,'Xmin,Xmax,Ymin,Ymax of PDF'
              write(30,'(6ES15.7,T100,A)')tr,'Tr'              
              do i=1,nbin2dx+1
                 do j=1,nbin2dy+1
                    write(30,'(ES15.7,$)')z(i,j)
                 end do
                 write(30,'(1x)')
              end do
           end if
           
           
           
           if(plot.eq.1) then
              !Plot true value, median, ranges, etc.
              call pgsci(1)
              if(plotsky.eq.1) call pgsci(0)
              call pgsls(2)
              !True value
              if(plotsky.eq.0) then
                 call pgline(2,(/startval(ic,p1,1),startval(ic,p1,1)/),(/-1.e20,1.e20/))
                 call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,1),startval(ic,p2,1)/))
              end if
              call pgsci(1)
              call pgsls(4)
              !Starting value
              !         call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/-1.e20,1.e20/))
              !         call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,2),startval(ic,p2,2)/))
              call pgsci(2)
              !Interval ranges
              !         call pgline(2,(/ranges(ic,c0,p1,1),ranges(ic,c0,p1,1)/),(/-1.e20,1.e20/))
              !         call pgline(2,(/ranges(ic,c0,p1,2),ranges(ic,c0,p1,2)/),(/-1.e20,1.e20/))
              !         call pgline(2,(/-1.e20,1.e20/),(/ranges(ic,c0,p2,1),ranges(ic,c0,p2,1)/))
              !         call pgline(2,(/-1.e20,1.e20/),(/ranges(ic,c0,p2,2),ranges(ic,c0,p2,2)/))
              call pgsls(1)
              call pgsch(sch*0.6)
              call pgsah(1,45.,0.1)
              a = 0.0166667*sch
              call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,1),ymin+dy*a)
              call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,2),ymin+dy*a)
              a = 0.0333333*sch
              call pgptxt(ranges(ic,c0,p1,3),ymin+dy*a,0.,0.5,'\(2030)\d90%\u')
              a = 0.0233333*sch
              call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,1))
              call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,2))
              a = 0.01*sch
              call pgptxt(xmin+dx*a,ranges(ic,c0,p2,3),90.,0.5,'\(2030)\d90%\u')
              call pgsls(2)
              call pgsch(sch)
              
              !Median
              call pgline(2,(/stats(ic,p1,1),stats(ic,p1,1)/),(/-1.e20,1.e20/))
              call pgline(2,(/-1.e20,1.e20/),(/stats(ic,p2,1),stats(ic,p2,1)/))
              call pgsls(1)
              
              
              !Big star at true position
              if(plotsky.eq.1) then
                 call pgsch(sch*2)
                 call pgsci(9)
                 call pgpoint(1,startval(ic,p1,1),startval(ic,p2,1),18)
                 call pgsch(sch)
                 call pgsci(1)
              end if
              
              
              
              
              
              
              call pgsls(1)
              call pgslw(lw2)
              if(plotsky.eq.1) then
                 call pgsci(0)
                 call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !Box, ticks, etc in white
                 call pgsci(1)
                 call pgbox('N',0.0,0,'N',0.0,0) !Number labels in black
              else
                 call pgsci(1)
                 call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
              end if
              call pgmtxt('B',2.2,0.5,0.5,trim(pgvarns(p1)))
              call pgmtxt('L',1.7,0.5,0.5,trim(pgvarns(p2)))
              
              if(file.eq.1) then
                 call pgend
                 i = system('convert -resize '//trim(bmpxpix)//' -depth 8 pdf2d.ppm bulk/'//trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.png')
                 if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
                 i = system('rm -f pdf2d.ppm')
              end if
              if(file.ge.2) call pgpage
           end if !if(plot.eq.1)
           
        end do !p2
     end do !p1
        
     if(savepdf.eq.1) close(30)
     
     if(plot.eq.1) then
        if(file.ne.1) call pgend
        if(file.ge.2) then
           if(abs(j2-j1).le.1) then
              if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o bulk/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.pdf  >& /dev/null')
              i = system('mv -f pdf2d.eps bulk/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.eps')
           else
              if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o bulk/'//trim(outputname)//'__pdf2d.pdf  >& /dev/null')
              i = system('mv -f pdf2d.eps bulk/'//trim(outputname)//'__pdf2d.eps')
           end if
        end if
     end if
     
  end if !if(plpdf2d.eq.1)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  
  
  if(plmovie.eq.1) then
     p = par1
     
     !moviescheme = 1  !Left column: Chain, sigma and acceptance, right column: numbers and PDF
     moviescheme = 2  !Upper panel: numbers and PDF, lower panel: chain
     
     do iframe = 0,nframes
     nplt = nint(real(iframe)/real(nframes)*maxval(ntot(1:nchains)))-1
     
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting movie frame...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,I5,A1,I5,A,I7,A)')' Plotting movie frame',iframe,'/',nframes,'  (',nplt,' points)'
     !print*,iframe,nplt
     write(framename,'(A20,I4.4,A4)')'movies/frames/frame_',iframe,'.ppm'

     ic = 1
     
     !Determine median and ranges
     if(nplt.gt.nburn(ic)) then
        
        !Determine the median
        call rindexx(nplt-nburn(ic),pldat(ic,p,1:nplt-nburn(ic)),index(p,1:nplt-nburn(ic)))  !Sort
        if(mod(nplt-nburn(ic),2).eq.0) median = 0.5*(pldat(ic,p,index(p,(nplt-nburn(ic))/2)) + pldat(ic,p,index(p,(nplt-nburn(ic))/2+1)))  !Centre = nb + (n-nb)/2 = (n+nb)/2
        if(mod(nplt-nburn(ic),2).eq.1) median = pldat(ic,p,index(p,(nplt-nburn(ic)+1)/2))
        !print*,'median:',median
        
        
        !Determine interval ranges
        c = c0
        ival = ivals(c)
        
        minrange = 1.e30
        !do i=1,floor(nplt*(1.-ival))
        do i=1,floor((nplt-nburn(ic))*(1.-ival))
           x1 = pldat(ic,p,index(p,i))
           x2 = pldat(ic,p,index(p,i+floor((nplt-nburn(ic))*ival)))
           range = abs(x2 - x1)
           if(range.lt.minrange) then
              minrange = range
              y1 = x1
              y2 = x2
           end if
           !write(*,'(2I6,7F12.8)')i,i+floor(nplt*ival),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
        
        !Save ranges:
        range1 = y1
        range2 = y2
        drange = y2-y1
        if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) drange = drange/((y1+y2)/2.)
        
        !print*,'ranges:',range1,range2,drange
     end if
     
     
     
     
     
     if(file.eq.0) io = pgopen('17/xs')
     if(file.eq.1) io = pgopen(trim(framename)//'/ppm')
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     !if(file.eq.0) call pgpap(12.,0.72)
     if(file.eq.0) call pgpap(scrsz*0.75,scrrat)
     !if(file.eq.1) call pgpap(20.,0.72)   !for 850x612, change convert below.  1:1.388 for mpeg, make the output image twice as big, rescale in the end
     if(file.eq.1) call pgpap(24.08,0.72) !for 1024x738, change convert below. 1:1.388 for mpeg, make the output image twice as big, rescale in the end
     call pgsch(1.)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file)
     lw = 2
     if(file.eq.1) lw = 5
     
     
     
     
     
  !***********************************************************************************************************************************      
     !Plot chain for this parameter
     
     
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - parameter chains'
     call pgslw(lw)
     if(moviescheme.eq.1) call pgsvp(0.05,0.35,0.65,0.95)
     if(moviescheme.eq.2) call pgsvp(0.05,0.95,0.05,0.25)
     ic = 1
     xmin = 0.
     xmax = real(maxval(ntot(1:nchains0)))
     dx = abs(xmax-xmin)*0.01
     if(moviescheme.eq.2) dx = 0.
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nchains0
        ymin = min(ymin,minval(pldat(ic,p,10:ntot(ic))))
        ymax = max(ymax,maxval(pldat(ic,p,10:ntot(ic))))
     end do
     dy = abs(ymax-ymin)*0.05
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(moviescheme.eq.1) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
     if(moviescheme.eq.2) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     
     do ic=1,nchains0
        !call pgsci(mod(ic*2,10))
        if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !do i=1,ntot(ic),chainpli
        do i=1,nplt,chainpli
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
     end do
     
     do ic=1,nchains0
        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))
        call pgsci(6)
        !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
        call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        call pgsci(3)
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))
     end do
     call pgsci(1)
     call pgsls(1)
     !call pgmtxt('T',1.,0.5,0.5,'Chain')
     call pgmtxt('T',-1.5,0.05,0.0,'Chain')
     
     
     
     
     !***********************************************************************************************************************************            
     !Plot sigma values ('jump size')
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - sigma'
     
     if(moviescheme.eq.1) then
        call pgsvp(0.05,0.35,0.35,0.65)
        
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,real(ntot(ic)))
           dx = abs(xmax-xmin)*0.01
           ymin = min(ymin,minval(sig(p,ic,10:ntot(ic))))
           ymax = max(ymax,maxval(sig(p,ic,10:ntot(ic))))
           dy = abs(ymax-ymin)*0.05
        end do
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
        
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !        !do i=1,ntot(ic),chainpli
           !do i=1,nplt,chainpli
           do i=ic,nplt,chainpli !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),sig(p,ic,i),1)
           end do
        end do
        
        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',1.,0.5,0.5,'Sigma')
        call pgmtxt('T',-1.5,0.05,0.0,'\(2144)')
        
     end if
  

     !***********************************************************************************************************************************            
     !Plot acceptance rate
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - acceptance rate'
     
     if(moviescheme.eq.1) then
        call pgsvp(0.05,0.35,0.05,0.35)
        
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,real(ntot(ic)))
           dx = abs(xmax-xmin)*0.01
           do i=1,ntot(ic)
              if(acc(p,ic,i).gt.1.e-10 .and. acc(p,ic,i).lt.1.-1.e-10) then
                 n0 = i
                 exit
              end if
           end do
           n0 = n0+10
           ymin = min(ymin,minval(acc(p,ic,n0:ntot(ic))))
           ymax = max(ymax,maxval(acc(p,ic,n0:ntot(ic))))
           dy = abs(ymax-ymin)*0.05
        end do
        
        call pgsci(1)
        call pgsls(1)
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        
        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !do i=1,ntot(ic),chainpli
           !do i=1,nplt,chainpli
           do i=ic,nplt,chainpli !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           end do
        end do
        
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',1.,0.5,0.5,'Acceptance')
        call pgmtxt('T',-1.5,0.05,0.0,'Acceptance')
        
     end if






     !***********************************************************************************************************************************            
     !Plot 1D pdf
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - pdf'
     
     if(moviescheme.eq.1) call pgsvp(0.45,0.95,0.05,0.8)
     if(moviescheme.eq.2) call pgsvp(0.25,0.95,0.32,0.999)
     
     !Set x-ranges, bin the data and get y-ranges
     xmin = 1.e30
     xmax = -1.e30
     p = par1
     do ic=1,nchains
        xmin = min(xmin,minval(pldat(ic,p,1:n(ic))))
        xmax = max(xmax,maxval(pldat(ic,p,1:n(ic))))
     end do
     dx = xmax - xmin
     xmin = xmin - 0.1*dx
     xmax = xmax + 0.1*dx
     
     !if(iframe.gt.0) then
     if(nplt.gt.nburn(ic)) then
        do ic=1,nchains
           x(ic,1:nplt-nburn(ic)) = pldat(ic,p,nburn(ic)+1:nplt)
           xmin1 = minval(x(ic,1:nplt-nburn(ic)))
           xmax1 = maxval(x(ic,1:nplt-nburn(ic)))
           call bindata(nplt-nburn(ic),x(ic,1:nplt-nburn(ic)),1,nbin,xmin1,xmax1,xbin1,ybin1) !Count the number of points in each bin
           
           !Normalise the SURFACE, not the height (because of different bin size)
           norm = 0.
           do i=1,nbin+1
              norm = norm + ybin1(i)
           end do
           norm = norm*(xmax1-xmin1)
           ybin1 = ybin1/norm
           
           !Smoothen
           ybin2 = ybin1
           if(smooth.gt.1) then
              !i0 = nbin/10
              !i0 = nint(min(max(real(nbin)/real(smooth),1.0),real(nbin)/2.))
              i0 = min(max(smooth,1),floor(real(nbin)/2.))
              do i=1+i0,nbin+1-i0
                 coefs1(1:2*i0+1) = ybin1(i-i0:i+i0)
                 call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
                 do i1=1,i0+1
                    coefs(i0-i1+2) = coefs1(i1)
                 end do
                 do i1 = i0+2,2*i0+1
                    coefs(3*i0+3-i1) = coefs1(i1)
                 end do
                 ybin2(i) = 0.
                 do i1=1,2*i0+1
                    ybin2(i) = ybin2(i) + coefs(i1) * ybin1(i+i1-i0-1)
                 end do
              end do
              ybin1 = ybin2
           end if
           xbin(ic,1:nbin+1) = xbin1(1:nbin+1)
           ybin(ic,1:nbin+1) = ybin1(1:nbin+1)
        end do
        
        ymin = 0.
        ymax = -1.e30
        do ic=1,nchains
           ymax = max(ymax,maxval(ybin(ic,1:nbin+1)))
        end do
        ymax = ymax*1.05
        
        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
        
        
        call pgsci(1)
        call pgslw(lw)
        
        
        
        !print*,'Plot PDF'
        !Plot 1D PDF
        do ic=1,nchains
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           xbin1(1:nbin+1) = xbin(ic,1:nbin+1)
           ybin1(1:nbin+1) = ybin(ic,1:nbin+1)
           if(wrap(ic,p).eq.0) then
              if(nchains.eq.1) call pgsci(15)
              call pgpoly(nbin+2,(/xbin1(1),xbin1(1:nbin+1)/),(/0.,ybin1(1:nbin+1)/))
              call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              call pgline(nbin+1,xbin1(1:nbin+1),ybin1(1:nbin+1)) !:nbin) ?
              
              call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/)) !Fix the loose ends
              call pgline(2,(/xbin1(nbin+1),xbin1(nbin+1)/),(/ybin1(nbin+1),0./))
           else
              plshift = real(2*pi)
              if(changevar.eq.1) plshift = 360.
              if(changevar.eq.1.and.p.eq.8) plshift = 24. !RA in hours
              if(nchains.eq.1) call pgsci(15)
              call pgpoly(nbin+3,(/xbin1(1),xbin1(1:nbin),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin),ybin1(1),0./))
              call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              call pgline(nbin,xbin1(1:nbin),ybin1(1:nbin))
              
              call pgsls(4)
              call pgline(nbin+1,(/xbin1(1:nbin)-plshift,xbin1(1)/),(/ybin1(1:nbin),ybin1(1)/))
              call pgline(nbin,xbin1+plshift,ybin1)
              call pgsls(1)
              call pgline(2,(/xbin1(nbin),xbin1(1)+plshift/),(/ybin1(nbin),ybin1(1)/))
           end if
        end do !ic
        
        !print*,'Plot PDF lines'
        !Plot lines again over surface of overlapping distributions
        if(nchains.gt.1) then
           call pgsls(4)
           do ic=1,nchains
              call pgsci(1)
              xbin1(1:nbin+1) = xbin(ic,1:nbin+1)
              ybin1(1:nbin+1) = ybin(ic,1:nbin+1)
              if(wrap(ic,p).eq.0) then
                 call pgline(nbin+1,xbin1(1:nbin+1),ybin1(1:nbin+1))
              else
                 call pgline(nbin,xbin1(1:nbin),ybin1(1:nbin))
              end if
           end do
           call pgsls(4)
        end if
        
        
        
        
        !print*,'Plot median'
        !Plot median and model value
        call pgsch(sch)
        
        do ic=1,nchains
           !Draw white lines
           if(nchains.gt.1) then
              call pgslw(lw)
              call pgsls(1); call pgsci(0)
              call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
              call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
              call pgline(2,(/median,median/),(/-1.e20,1.e20/))
              call pgline(2,(/range1,range1/),(/-1.e20,1.e20/)) !Left limit of 90% interval
              call pgline(2,(/range2,range2/),(/-1.e20,1.e20/)) !Right limit of 90% interval
           end if
           
           call pgslw(lw+1)
           !Draw coloured lines over the white ones
           !Median
           call pgsls(2); call pgsci(defcolour); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           call pgline(2,(/median,median/),(/-1.e20,1.e20/))
           
           !Ranges
           call pgsls(4); call pgsci(defcolour); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           call pgline(2,(/range1,range1/),(/-1.e20,1.e20/)) !Left limit of 90% interval
           call pgline(2,(/range2,range2/),(/-1.e20,1.e20/)) !Right limit of 90% interval
           
           !True value
           call pgsls(2); call pgsci(1); if(nchains.gt.1) call pgsci(1)
           call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
           
           !Starting value
           !if(abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
           !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
           !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
           !end if
           
           call pgsls(1)
           call pgsci(1)
        end do
        
        
        !print*,'Print median'
        !Print median, model value and range widths
        call pgslw(lw)
        call pgsci(1)
        ic = 1
        if(moviescheme.eq.1) then
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
              write(str,'(A,F7.3,A5,F7.3,A9,F6.2,A1)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange*100,'%'
           else
              write(str,'(A,F7.3,A5,F7.3,A9,F7.3)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange
           end if
        end if
        if(moviescheme.eq.2) then
           call pgslw(lw)
           call pgsch(1.)
           call pgbox('BNTS',0.0,0,'',0.0,0) !Box for the PDF
           
           call pgsvp(0.05,0.25,0.35,0.999)
           call pgswin(0.,1.,1.,0.)
           call pgsch(2.)
           call pgslw(lw+2)
           
           call pgptxt(0.5,0.2,0.0,0.5,trim(pgvarns(p)))
           call pgslw(lw)
           call pgsch(1.4)
           
           if(nplt.gt.nburn(ic)) then
              write(str,'(A,F7.3)')' Model:',startval(ic,p,1)
              call pgptxt(0.05,0.45,0.,0.,trim(str)) 
              write(str,'(A,F7.3)')'Median:',median
              call pgptxt(0.05,0.6,0.,0.,trim(str)) 
              write(str,'(A,F7.3)')'      \(2030):',drange
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) write(str,'(A,F6.2,A1)')'      \(2030):',drange*100,'%' 
              call pgptxt(0.05,0.75,0.,0.,trim(str)) 
           end if
           call pgsch(1.)
        end if
        
     else
        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
     !end if  !if(iframe.gt.0)
     end if  !if(nplt.gt.nburn(ic)) 
     
     if(moviescheme.eq.1) then
        call pgsch(1.5)
        call pgslw(lw+2)
        call pgmtxt('T',2.5,0.5,0.5,trim(pgvarns(p)))
        call pgslw(lw)
        call pgsch(1.)
        !if(iframe.gt.0) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
        if(nplt.gt.nburn(ic)) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
        call pgbox('BNTS',0.0,0,'',0.0,0)
     end if

!************************************************************************************************************************************

     
     
     
     
     
     
     call pgend
     
     !if(file.eq.1) i = system('convert -depth 8 -resize 850x612 '//trim(framename)//' '//trim(framename))  !Rescale the output frame
     if(file.eq.1) then
        i = system('convert -depth 8 -resize 1024x738 '//trim(framename)//' '//trim(framename))  !Rescale the output frame
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
     end if
     end do
     
  end if
  
  
  
  
  
  
  
  
  
  if(update.eq.1) then
     deallocate(pldat,alldat)
     call sleep(5)
     if(sum(ntot).gt.1.e4) call sleep(5)
     if(sum(ntot).gt.1.e5) call sleep(10)
     if(sum(ntot).gt.1.e6) call sleep(20)
     goto 101
  end if
  
  !write(*,'(A)')'  Waiting for you to finish me off...'
  !pause
  
9999 continue
  deallocate(pldat,alldat)
  if(prprogress.ge.1) write(*,*)''
end program plotspins
!************************************************************************************************************************************






!************************************************************************************************************************************
subroutine pginitl(colour,file)  !Initialise pgplot
  implicit none
  integer :: colour,file,i
  if(1.eq.2) then
     call pgscr(0,1.,1.,1.) !Background colour always white (also on screen, bitmap)
     call pgscr(1,0.,0.,0.) !Default foreground colour always black
     if(file.le.1) then !png: create white background
        call pgsvp(-100.,100.,-100.,100.)
        call pgswin(0.,1.,0.,1.)
        call pgsci(0)
        call pgrect(-1.,2.,-1.,2.)
        call pgsvp(0.08,0.95,0.06,0.87) !Default viewport size (?)
        call pgsci(1)
     end if
  end if
  if(colour.eq.0) then
     do i=0,99
        call pgscr(i,0.,0.,0.)
     end do
     call pgscr(0,1.,1.,1.)
     call pgscr(14,0.3,0.3,0.3)
     call pgscr(15,0.8,0.8,0.8)
  else
     call pgscr(3,0.,0.5,0.) !Use darker green
  end if
end subroutine pginitl
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine bindata(n,x,norm,nbin,xmin1,xmax1,xbin,ybin)  !Count the number of points in each bin
  ! x - input: data, n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!

  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,ybintot,xmin1,xmax1

  xmin = xmin1
  xmax = xmax1

  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nbin)

  do k=1,nbin+1
     !        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  ybintot=0.
  do k=1,nbin
     ybin(k) = 0.
     do i=1,n
        if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + 1.
     end do
     ybintot = ybintot + ybin(k)
  end do
  if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)

  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if

end subroutine bindata
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata1(n,x,y,norm,nbin,xmin1,xmax1,xbin,ybin)  !Measure the amount of likelihood in each bin
  ! x - input: data, n points
  ! y - input: "weight" (likelihood), n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!

  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),y(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,ybintot,xmin1,xmax1,ymin
  
  xmin = xmin1
  xmax = xmax1
  ymin = minval(y)
  !print*,n,nbin,xmin1,xmax1
  !print*,minval(y),maxval(y)

  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nbin)

  do k=1,nbin+1
     !        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  ybintot=0.
  do k=1,nbin
     ybin(k) = 0.
     do i=1,n
        !if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + 1.
        if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + exp(y(i) - ymin)
     end do
     ybintot = ybintot + ybin(k)
  end do
  if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)

  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if

end subroutine bindata1
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2d(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)  !Count the number of points in each bin
  !x - input: data, n points
  !norm - input: normalise (1) or not (0)
  !nbin - input: number of bins
  !xmin, xmax - in/output: set xmin=xmax to auto-determine
  !xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!
  
  implicit none
  integer :: i,k,n,ix,iy,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),xbin(nxbin+1),ybin(nybin+1),z(nxbin+1,nybin+1),ztot,xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1
  real :: tr(6),zmax
  
  !write(*,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(*,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(*,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(*,'(50F5.2)'),x(1:50)
  !write(*,'(50F5.2)'),y(1:50)
  !write(*,'(20F8.5)'),xbin
  !write(*,'(20F8.5)'),ybin
  
  z = 0.
  ztot=0.
  !print*,xmin,xmax
  !print*,ymin,ymax
  do bx=1,nxbin
     !print*,bx,xbin(bx),xbin(bx+1)
     do by=1,nybin
        z(bx,by) = 0.
        do i=1,n
           !do ix=1,nx,10000
           !do iy=1,ny,10000
           if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) z(bx,by) = z(bx,by) + 1.
           !  write(*,'(2I4,2I9,7F10.5)')bx,by,ix,iy,x(ix),xbin(bx),xbin(bx+1),y(iy),ybin(by),ybin(by+1),z(bx,by)
           !end do
        end do
        ztot = ztot + z(bx,by) 
        !write(*,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(z(bx,by))
     end do
     !write(*,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(z(bx,1:nybin))
     end do
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2d
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2da(n,x,y,z,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,zz,tr)  !Measure the amount of likelihood in each bin
  !x,y - input: data, n points
  !z - input: amount for each point (x,y)
  !norm - input: normalise (1) or not (0)
  !nxbin,nybin - input: number of bins in each dimension
  !xmin1,xmax1 - in/output: ranges in x dimension, set xmin=xmax as input to auto-determine
  !ymin1,ymax1 - in/output: ranges in y dimension, set ymin=ymax as input to auto-determine
  !zz - output: binned data zz(x,y).  The x,y values are the left side of the bin(?)
  !tr - output: transformation elements for pgplot (pggray, pgcont)
  
  implicit none
  integer :: i,k,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),z(n),xbin(nxbin+1),ybin(nybin+1),zz(nxbin+1,nybin+1),zztot,xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1
  real :: tr(6),zmin
  
  !write(*,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(*,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(*,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  zmin = minval(z)
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(*,'(50F5.2)'),x(1:50)
  !write(*,'(50F5.2)'),y(1:50)
  !write(*,'(20F8.5)'),xbin
  !write(*,'(20F8.5)'),ybin
  
  zz = 0.
  zztot = 0.
  !print*,xmin,xmax
  !print*,ymin,ymax
  do bx=1,nxbin
     !print*,bx,xbin(bx),xbin(bx+1)
     do by=1,nybin
        zz(bx,by) = 0.
        do i=1,n
           !if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + 1.
           if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + exp(z(i) - zmin)
           !write(*,'(2I4,8F10.5)')bx,by,x(i),xbin(bx),xbin(bx+1),y(i),ybin(by),ybin(by+1),zz(bx,by),z(i)
        end do
        zztot = zztot + zz(bx,by) 
        !write(*,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(zz(bx,by))
     end do
     !write(*,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(zz(bx,1:nybin))
     end do
  !if(norm.eq.1) z = z/(zztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2da
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine verthist(n,x,y)  !x is the left of the bin!
  implicit none
  integer :: j,n
  real :: x(n+1),y(n+1)

  !      call pgline(2,(/0.,x(1)/),(/y(1),y(1)/))
  !      do j=1,n
  !        call pgline(2,(/x(j),x(j)/),y(j:j+1))
  !        call pgline(2,x(j:j+1),(/y(j+1),y(j+1)/))
  !      end do

  x(n+1) = x(n) + (x(n)-x(n-1))
  y(n+1) = 0.
  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do

end subroutine verthist
!************************************************************************************************************************************

!************************************************************************************************************************************
subroutine horzhist(n,x,y)
  implicit none
  integer :: j,n
  real :: x(n),y(n)

  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n-2
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do
  call pgline(2,x(n-1:n),(/y(n-1),y(n-1)/))
end subroutine horzhist
!************************************************************************************************************************************



!************************************************************************************************************************************
function ra(lon, GPSsec)
  !/* Derives right ascension (in radians!) from longitude given GMST (radians). */
  !/* Declination == latitude for equatorial coordinates.                        */


  !/* Derive the `Greenwich Mean Sidereal Time' (in radians!) */
  !/* from GPS time (in seconds).                              */
  !/* (see K.R.Lang(1999), p.80sqq.)                           */
  implicit none
  real*8 :: ra,lon,gmst,seconds,days,centuries,secCurrentDay
  real*8 :: gps0,leapseconds,GPSsec,tpi
  tpi = 8*datan(1.d0)

  gps0 = 630720013.d0 !GPS time at 1/1/2000 at midnight
  leapseconds = 32.d0 !At Jan 1st 2000
  if(GPSsec.gt.(gps0 + 189388800.d0)) leapseconds = leapseconds + 1.d0 !One more leapsecond after 1/1/2006
  if(GPSsec.lt.630720013.d0) write(*,'(A)')'WARNING: GMSTs before 1.1.2000 are inaccurate!'
  !Time since 1/1/2000 midnight
  seconds       = (GPSsec - gps0) + (leapseconds - 32.d0)
  days          = floor(seconds/(24.d0*3600.d0)) - 0.5d0
  secCurrentDay = mod(seconds, 24.d0*3600.d0)
  centuries     = days/36525.d0
  gmst = 24110.54841d0 + (centuries*(8640184.812866d0 + centuries*(0.093104d0 + centuries*6.2d-6)))
  gmst = gmst + secCurrentDay * 1.002737909350795d0   !UTC day is 1.002 * MST day
  gmst = mod(gmst/(24.d0*3600.d0),1.d0)
  gmst = gmst * tpi

  ra = mod(lon + gmst + 10*tpi,tpi)
end function ra
!************************************************************************************************************************************



!************************************************************************************************************************************
SUBROUTINE dindexx(n,arr,indx)
  INTEGER :: n,indx(n),M,NSTACK
  REAL*8 :: arr(n),a
  PARAMETER (M=7,NSTACK=50)
  INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
     end if
     i=l+1
     j=ir
     indxt=indx(l+1)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
     if(jstack.gt.NSTACK) write(*,'(A)')' NSTACK too small in dindexx'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
END SUBROUTINE dindexx
!************************************************************************************************************************************


!************************************************************************************************************************************
SUBROUTINE rindexx(n,arr,indx)
  INTEGER :: n,indx(n),M,NSTACK
  REAL :: arr(n),a
  PARAMETER (M=7,NSTACK=50)
  INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
     end if
     i=l+1
     j=ir
     indxt=indx(l+1)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
     if(jstack.gt.NSTACK) write(*,'(A)')' NSTACK too small in rindexx'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
END SUBROUTINE rindexx
!************************************************************************************************************************************


!************************************************************************************************************************************
SUBROUTINE savgol(c,np,nl,nr,ld,m)
  INTEGER :: ld,m,nl,np,nr,MMAX
  REAL :: c(np)
  PARAMETER (MMAX=6)
  !     USES lubksb,ludcmp
  INTEGER :: imj,ipj,j,k,kk,mm,indx(MMAX+1)
  REAL :: d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
  !if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m) pause 'bad args in savgol'
  if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m) write(*,'(A)')' Bad args in savgol'
  do ipj=0,2*m
     sum=0.
     if(ipj.eq.0)sum=1.
     do k=1,nr
        sum=sum+float(k)**ipj
     end do
     do k=1,nl
        sum=sum+float(-k)**ipj
     end do
     mm=min(ipj,2*m-ipj)
     do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
     end do
  end do
  call ludcmp(a,m+1,MMAX+1,indx,d)
  do j=1,m+1
     b(j)=0.
  end do
  b(ld+1)=1.
  call lubksb(a,m+1,MMAX+1,indx,b)
  do kk=1,np
     c(kk)=0.
  end do
  do k=-nl,nr
     sum=b(1)
     fac=1.
     do mm=1,m
        fac=fac*k
        sum=sum+b(mm+1)*fac
     end do
     kk=mod(np-k,np)+1
     c(kk)=sum
  end do
  return
END SUBROUTINE savgol
!************************************************************************************************************************************


!************************************************************************************************************************************
SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER :: n,np,indx(n)
  REAL :: a(np,np),b(n)
  INTEGER :: i,ii,j,ll
  REAL :: sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     end if
     b(i)=sum
  end do
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do
  return
end SUBROUTINE lubksb
!************************************************************************************************************************************



!************************************************************************************************************************************
SUBROUTINE ludcmp(a,n,np,indx,d)
 INTEGER :: n,np,indx(n),NMAX
 REAL :: d,a(np,np),TINY
 PARAMETER (NMAX=500,TINY=1.0e-20)
 INTEGER :: i,imax,j,k
 REAL :: aamax,dum,sum,vv(NMAX)
 d=1.
 do i=1,n
    aamax=0.
    do j=1,n
       if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
    if(aamax.eq.0.) write(*,'(A)')' Singular matrix in ludcmp'
    vv(i)=1./aamax
 end do
 do j=1,n
    do i=1,j-1
       sum=a(i,j)
       do k=1,i-1
          sum=sum-a(i,k)*a(k,j)
       end do
       a(i,j)=sum
    end do
    aamax=0.
    do i=j,n
       sum=a(i,j)
       do k=1,j-1
          sum=sum-a(i,k)*a(k,j)
       end do
       a(i,j)=sum
       dum=vv(i)*abs(sum)
       if (dum.ge.aamax) then
          imax=i
          aamax=dum
       end if
    end do
    if (j.ne.imax)then
       do k=1,n
          dum=a(imax,k)
          a(imax,k)=a(j,k)
          a(j,k)=dum
       end do
       d=-d
       vv(imax)=vv(j)
    end if
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=TINY
    if(j.ne.n)then
       dum=1./a(j,j)
       do i=j+1,n
          a(i,j)=a(i,j)*dum
       end do
    end if
 end do
 return
end SUBROUTINE ludcmp
!************************************************************************************************************************************






!************************************************************************************************************************************
subroutine plotthesky(bx1,bx2,by1,by2)
  implicit none
  integer, parameter :: ns=9110, nsn=80
  integer :: i,j,c(100,35),nc,snr(nsn),plcst,plstar,cf,spld,n,prslbl,rv
  real*8 :: ra(ns),dec(ns),d2r,r2d,r2h,pi,dx,dx1,dx2,dy,dy1,ra1,dec1,rev,l0,b0,par
  real :: pma,pmd,vm(ns),dist(ns),x1,y1,x2,y2,constx(99),consty(99),r1,g1,b1,r4,g4,b4
  real :: schcon,sz1,schfac,schlbl,prinf,snlim,sllim,schmag,getmag,mag,bx1,bx2,by1,by2,x,y,mlim
  character :: cn(100)*3,con(100)*20,name*10,vsopdir*99,sn(ns)*10,snam(nsn)*10,sni*10,getsname*10,mult,var*9
  
  mlim = 6.
  cf = 2
  schmag = 0.07
  schlbl = 1.
  schfac = 1.
  schcon = 1.
  plstar = 4
  plcst = 1
  
  prinf = 150.**2
  
  x = 0.
  call pgqcr(1,r1,g1,b1) !Store colours
  call pgqcr(4,r4,g4,b4)
  call pgscr(1,1.,1.,1.) !'White' (for stars)
  call pgscr(4,x,x,1.) !Blue (for constellations)
  
  pi = 4*datan(1.d0)
  d2r = pi/180.d0
  r2d = 180.d0/pi
  r2h = 12.d0/pi
  r2h = r2d
  
  
  if(bx1.gt.bx2) then
     x = bx1
     bx1 = bx2
     bx2 = x
  end if
  
  !Read BSC
  vsopdir = '/home/sluys/diverse/popular/fortran/VSOP87/'           !Linux pc
  open(unit=20,form='formatted',status='old',file=trim(vsopdir)//'data/bsc.dat')
  rewind(20)
  do i=1,ns
     read(20,320)name,ra(i),dec(i),pma,pmd,rv,vm(i),par,mult,var
320  format(A10,1x,2F10.6,1x,2F7.3,I5,F6.2,F6.3,A2,A10)
     sn(i) = getsname(name)
  end do
  close(20)


  !Read Constellation figure data
  open(unit=40,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_const.dat')
  do i=1,ns
     read(40,'(I4)',end=340,advance='no')c(i,1)
     do j=1,c(i,1)
        read(40,'(I5)',advance='no')c(i,j+1)
     end do
     read(40,'(1x,A3,A20)')cn(i),con(i)
     !Get mean star position to place const. name
     dx1 = 0.d0
     dx2 = 0.d0
     dy = 0.d0
     do j=2,c(i,1)
        dx1 = dx1 + dsin(ra(c(i,j)))
        dx2 = dx2 + dcos(ra(c(i,j)))
        dy = dy + dec(c(i,j))
     end do
     dx1 = (dx1 + dsin(ra(c(i,j))))/real(c(i,1))
     dx2 = (dx2 + dcos(ra(c(i,j))))/real(c(i,1))
     ra1 = rev(datan2(dx1,dx2))
     dec1 = (dy + dec(c(i,j)))/real(c(i,1))
     !call eq2xy(ra1,dec1,l0,b0,x1,y1)
     !constx(i) = x1
     !consty(i) = y1
     constx(i) = real(ra1*r2h)
     consty(i) = real(dec1*r2d)
  end do
340 close(40)
  nc = i-1
  
  !Read Star names
  open(unit=50,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_names.dat')
  do i=1,nsn
     read(50,'(I4,2x,A10)',end=350)snr(i),snam(i)
  end do
350 close(50)
  
  
  !!Read Milky Way data
  !do f=1,5
  !   write(mwfname,'(A10,I1,A4)')'milkyway_s',f,'.dat'
  !   open(unit=60,form='formatted',status='old',file=trim(vsopdir)//'data/'//mwfname)
  !   do i=1,mwn(f)
  !      read(60,'(F7.5,F9.5)')mwa(f,i),mwd(f,i)
  !      if(maptype.eq.1) call eq2az(mwa(f,i),mwd(f,i),agst)
  !      if(maptype.eq.2) call eq2ecl(mwa(f,i),mwd(f,i),eps)
  !   end do
  !end do
  !close(60)
  
  
  !Plot constellation figures
  if(plcst.gt.0) then
     !schcon = min(max(40./sz1,0.7),3.)
     call pgsch(schfac*schcon*schlbl)
     call pgscf(cf)
     call pgsci(4)
     call pgslw(2)
     do i=1,nc
        do j=2,c(i,1)
           !call eq2xy(ra(c(i,j)),dec(c(i,j)),l0,b0,x1,y1)
           !call eq2xy(ra(c(i,j+1)),dec(c(i,j+1)),l0,b0,x2,y2)
           x1 = real(ra(c(i,j))*r2h)
           y1 = real(dec(c(i,j))*r2d)
           x2 = real(ra(c(i,j+1))*r2h)
           y2 = real(dec(c(i,j+1))*r2d)
           !if((x1*x1+y1*y1.le.prinf.or.x2*x2+y2*y2.le.prinf).and.(x2-x1)**2+(y2-y1)**2.le.90.**2) & !Not too far from centre and each other 
           if((x2-x1)**2+(y2-y1)**2.le.90.**2) & !Not too far from centre and each other 
                call pgline(2,(/x1,x2/),(/y1,y2/))
	end do
        if(constx(i).lt.bx1.or.constx(i).gt.bx2.or.consty(i).lt.by1.or.consty(i).gt.by2) cycle
        if(plcst.eq.2) call pgptext(constx(i),consty(i),0.,0.5,cn(i))
        if(plcst.eq.3) call pgptext(constx(i),consty(i),0.,0.5,con(i))
     end do
     call pgsch(schfac)
     call pgscf(cf)
  end if !if(plcst.gt.0) then
  
  !Plot stars: BSC
  spld = 0
  if(plstar.gt.0) then
     n = 0
     do i=1,ns
        if(vm(i).lt.mlim.and.vm(i).ne.0.) then
           !call eq2xy(ra(i),dec(i),l0,b0,x,y)
           x = real(ra(i)*r2h)
           y = real(dec(i)*r2d)
           if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
           call pgsci(1)
           mag = getmag(vm(i),mlim)*schmag
           call pgcirc(x,y,mag)
           !write(*,'(3F10.3)')x,y,mag
           call pgsch(schfac*schlbl)
           sni = sn(i)
           !if(sni(1:1).eq.'\') call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font
           if(sni(1:1).eq.char(92)) call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font.  Char(92) is a \, but this way it doesn't mess up emacs' parentheses count
	   call pgsci(14)
           if(vm(i).lt.sllim) then
              if((plstar.eq.2.or.plstar.eq.5)) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              if(plstar.eq.4) then !Check if the name will be printed
                 prslbl = 1
                 if(vm(i).lt.snlim) then
                    do j=1,nsn
                       if(snr(j).eq.i) prslbl = 0 !Then the name will be printed, don't print the symbol
                    end do
                 end if
                 if(prslbl.eq.1) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              end if
           end if
	   spld = spld+1
	end if
     end do
     if(plstar.ge.3) then !Plot star proper names
        call pgsch(schfac*schlbl)
        do i=1,nsn
           if(vm(snr(i)).lt.max(snlim,1.4)) then  !Regulus (1.35) will still be plotted, for conjunction maps
              !call eq2xy(ra(snr(i)),dec(snr(i)),l0,b0,x,y)
              x = real(ra(snr(i)))
              y = real(dec(snr(i)))
              if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
              call pgtext(x+0.02*sz1,y-0.02*sz1,snam(i))
           end if
        end do
     end if !if(plstar.eq.3) then
  end if !if(plstar.gt.0) then
  
  !Restore colours
  call pgscr(1,r1,g1,b1)
  call pgscr(4,r4,g4,b4)
  
end subroutine plotthesky
!************************************************************************************************************************************

!************************************************************************
function getsname(name)               !Get star name from bsc info
  implicit none
  character :: getsname*10,name*10,num*3,grk*3,gn*1
  num = name(1:3)
  grk = name(4:6)
  gn  = name(7:7)
  !      gn = ' '
  
  getsname = '          '
  if(grk.ne.'   ') then  !Greek letter
     if(grk.eq.'Alp') getsname = '\(2127)\u'//gn
     if(grk.eq.'Bet') getsname = '\(2128)\u'//gn
     if(grk.eq.'Gam') getsname = '\(2129)\u'//gn
     if(grk.eq.'Del') getsname = '\(2130)\u'//gn
     if(grk.eq.'Eps') getsname = '\(2131)\u'//gn
     if(grk.eq.'Zet') getsname = '\(2132)\u'//gn
     if(grk.eq.'Eta') getsname = '\(2133)\u'//gn
     if(grk.eq.'The') getsname = '\(2134)\u'//gn
     if(grk.eq.'Iot') getsname = '\(2135)\u'//gn
     if(grk.eq.'Kap') getsname = '\(2136)\u'//gn
     if(grk.eq.'Lam') getsname = '\(2137)\u'//gn
     if(grk.eq.'Mu ') getsname = '\(2138)\u'//gn
     if(grk.eq.'Nu ') getsname = '\(2139)\u'//gn
     if(grk.eq.'Xi ') getsname = '\(2140)\u'//gn
     if(grk.eq.'Omi') getsname = '\(2141)\u'//gn
     if(grk.eq.'Pi ') getsname = '\(2142)\u'//gn
     if(grk.eq.'Rho') getsname = '\(2143)\u'//gn
     if(grk.eq.'Sig') getsname = '\(2144)\u'//gn
     if(grk.eq.'Tau') getsname = '\(2145)\u'//gn
     if(grk.eq.'Ups') getsname = '\(2146)\u'//gn
     if(grk.eq.'Phi') getsname = '\(2147)\u'//gn
     if(grk.eq.'Chi') getsname = '\(2148)\u'//gn
     if(grk.eq.'Psi') getsname = '\(2149)\u'//gn
     if(grk.eq.'Ome') getsname = '\(2150)\u'//gn
  else  !Then number
     if(num(1:1).eq.' ') num = num(2:3)//' '
     if(num(1:1).eq.' ') num = num(2:3)//' '
     getsname = num//'       '
  end if
  return
end function getsname
!************************************************************************

!************************************************************************
function getmag(m,mlim)  !Determine size of stellar 'disk'
  real :: getmag,m,m1,mlim
  m1 = m
  !      if(m1.lt.0.) m1 = m1*0.5  !Less excessive grow in diameter for the brightest objects
  if(m1.lt.-1.e-3) m1 = -sqrt(-m1)  !Less excessive grow in diameter for the brightest objects
  !getmag = max(mlim-m1+0.5,0.)
  getmag = max(mlim-m1+0.5,0.5) !Make sure the weakest stars are still plotted
  !getmag = max(mlim-m1+0.5,0.)+0.5
  return
end function getmag
!************************************************************************


!************************************************************************
function rev(x)        !Returns angle in radians between 0 and 2pi
  real*8 :: x,rev,pi
  pi = 4*datan(1.d0)
  rev = x-floor(x/(2*pi))*2*pi
  return
end function rev
!************************************************************************

!************************************************************************
function rev360(x)        !Returns angle in degrees between 0 and 360
  real :: x,rev360
  rev360 = x-floor(x/(360.))*360.
  return
end function rev360
!************************************************************************

!************************************************************************
function rev24(x)        !Returns angle in hours between 0 and 24
  real :: x,rev24
  rev24 = x-floor(x/(24.))*24.
  return
end function rev24
!************************************************************************

!************************************************************************
function rev2pi(x)        !Returns angle in radians betwee 0 and 2pi
  real :: x,rev2pi,pi
  pi = 4*atan(1.)
  rev2pi = x-floor(x/(2.0*pi))*2.0*pi
  return
end function rev2pi
!************************************************************************


!************************************************************************
subroutine kstwo(data1,n1,data2,n2,d,prob)  !Needs probks(), sort()
  integer :: n1,n2,j1,j2
  real*8 :: d,prob,data1(n1),data2(n2)
  real*8 :: d1,d2,dt,en1,en2,en,fn1,fn2,probks
  call sort(n1,data1)
  call sort(n2,data2)
  en1=n1
  en2=n2
  j1=1
  j2=1
  fn1=0.d0
  fn2=0.d0
  d=0.d0
1 if(j1.le.n1.and.j2.le.n2)then
     d1=data1(j1)
     d2=data2(j2)
     if(d1.le.d2)then
        fn1=j1/en1
        j1=j1+1
     end if
     if(d2.le.d1)then
        fn2=j2/en2
        j2=j2+1
     end if
     dt=dabs(fn2-fn1)
     if(dt.gt.d)d=dt
     goto 1
  end if
  en=dsqrt(en1*en2/(en1+en2))
  prob=probks((en+0.12d0+0.11d0/en)*d)
  return
end subroutine kstwo
!************************************************************************


!************************************************************************
function probks(alam)
  real*8 :: probks,alam,EPS1,EPS2
  PARAMETER (EPS1=1.d-3, EPS2=1.d-8)
  integer :: j
  real*8 :: a2,fac,term,termbf
  a2=-2.d0*alam**2
  fac=2.d0
  probks=0.d0
  termbf=0.d0
  do j=1,100
     term=fac*dexp(a2*j**2)
     probks=probks+term
     if(dabs(term).le.EPS1*termbf.or.dabs(term).le.EPS2*probks)return
     fac=-fac
     termbf=dabs(term)
  end do
  probks=1.d0
  return
end function probks
!************************************************************************

!************************************************************************
subroutine sort(n,arr)
  integer :: n,M,NSTACK
  real*8 :: arr(n)
  PARAMETER (M=7,NSTACK=50)
  integer :: i,ir,j,jstack,k,l,istack(NSTACK)
  real*8 :: a,temp
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        a=arr(j)
        do i=j-1,l,-1
           if(arr(i).le.a)goto 2
           arr(i+1)=arr(i)
        end do
        i=l-1
2       arr(i+1)=a
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     temp=arr(k)
     arr(k)=arr(l+1)
     arr(l+1)=temp
     if(arr(l).gt.arr(ir))then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l+1).gt.arr(ir))then
        temp=arr(l+1)
        arr(l+1)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l).gt.arr(l+1))then
        temp=arr(l)
        arr(l)=arr(l+1)
        arr(l+1)=temp
     end if
     i=l+1
     j=ir
     a=arr(l+1)
3    continue
     i=i+1
     if(arr(i).lt.a)goto 3
4    continue
     j=j-1
     if(arr(j).gt.a)goto 4
     if(j.lt.i)goto 5
     temp=arr(i)
     arr(i)=arr(j)
     arr(j)=temp
     goto 3
5    arr(l+1)=arr(j)
     arr(j)=a
     jstack=jstack+2
     if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
end subroutine sort
!************************************************************************
