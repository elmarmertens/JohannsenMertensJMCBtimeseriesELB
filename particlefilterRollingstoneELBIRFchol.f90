PROGRAM main

  ! particle filter for rollingstone benchmark model 
  ! WITH ELB sampling
  ! AND IRF
  ! simulated IRF: extra dimension to store plus/munus shocks separately

  USE embox, only : hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec,  storeEstimatesTranspose, quantilevec, median, variance, loft, timestampstr, es30d16, int2str
  USE blaspack, only : vech, ivech, ivechU ! recall: ivech is a subroutine, ivechU is a pure function
  USE gibbsbox, only : drawNDXpdf ! , drawNDXsysresample

  USE vslbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  ! exe parameters
  logical :: doIRFactualObserver = .false., doIRFunitscale = .true.
  logical, parameter :: doELB = .true.
  logical, parameter :: doTimestamp = .false., doGains = .true.
  logical, parameter :: doNormalizedGains = .false.
  integer, parameter :: parameterColumns = 12, meanCol = 2 ! meanCol is used to pick parameters from input file; meanCol = 2 chooses median

  ! integer, parameter :: Nshadowproposals = 100


  character(LEN=100), parameter :: mcmcdir = 'datMCMC/'

  ! model parameters
  INTEGER, PARAMETER :: Nyield = 3, Ny = 3 + Nyield, Nbar = 1 + 1 + Nyield, Ngap = Ny

  INTEGER :: p = 2, Nf

  INTEGER, DIMENSION(Ny) :: ordervec = (/ 1,2,3,4,5,6 /)
  character(LEN=100) :: orderCode = '123456'
  INTEGER :: yndxPolicyrate = 3

  ! there are Ngap SV terms plus trend inflation 
  INTEGER, PARAMETER :: Nsv = Ngap + 1, Nw = 1 + Nsv, Ngapshockslopes = Ngap * (Ngap - 1) / 2 
  INTEGER, PARAMETER :: NhSigma = Nsv * (Nsv + 1) / 2, NhGapSigma = Ngap * (Ngap + 1) / 2 ! includes also trend SV which is however orthogonal to gapSV

  INTEGER :: Nx, NsigmaX

  INTEGER :: Nparticles,  Nmixturedraws, Nmcmcdraws, Nsim
  INTEGER :: T,i,j,k,status 

  ! filter particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARTICLEweights
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: ELBtwist
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWsvol, DRAWxhat, DRAWsqrtSigmaX
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: DRAWshadowratehat
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWytilde
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxgain, DRAWxgainstandardized 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sqrtSigmaX

  ! parameters and priors
  DOUBLE PRECISION, ALLOCATABLE :: f(:)
  DOUBLE PRECISION :: gapshockslopes(Ngapshockslopes)
  DOUBLE PRECISION :: hbar(Nsv), hrho(Nsv), sqrtVh(Nsv,Nsv)
  DOUBLE PRECISION :: rbarvol

  ! IRF
  integer, parameter :: irfNhorizons = 48  + 5 * 4
  integer, parameter :: irfNimpulse  = Ny ! report always all choleski shocks 
  integer, parameter :: irfNy = Ny + 2 ! all observables plus pibar, rbar
  INTEGER, PARAMETER :: irfNysim = irfNy + 1 + Nyield ! number of variables for which imulated IRF get computed, all plus the actual rate
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxirf, DRAWyirf
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRAWinnovol
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWysimbaseline
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: DRAWysimirf
  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRAWsvirf ! only 2D, so far nly impact response and in response to ELB shock
  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: IRFweights ! note: IRFweights just tracked for diagnostic purposes, DRAWxirf will already be aggregated across particles inside particlefilter


  ! data
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN
  LOGICAL, ALLOCATABLE, DIMENSION(:)   :: atELB
  double precision, parameter :: elbound = 0.25d0 ! note: this parameter is a late addition to the code, not having thoroughly replaced all instances of hardcoding the elb bound 

  ! logMDD contributions
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: loglike 

  ! other
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Xdraws
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta2, shadowratedraws
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: zdraw
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: shadowrateELBmask
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: theta1

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ndx

  CHARACTER (LEN=200) :: filename, datafile, nandatafile, fileXT, datalabel, parameterlabel, parameterXT

  ! VSL Random Stuff"
  type (vsl_stream_state) :: VSLstream
  type (vsl_stream_state), allocatable, dimension(:) :: VSLarray

  integer :: seed
  ! integer :: brng
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS !, TID

  ! state vector indices
  INTEGER, PARAMETER :: ndxPIbar = 1, ndxRbar = 2, ndxGapStart = ndxRbar + Nyield + 1,  ndxGapStop   = ndxGapStart + Ny - 1, ndxUgap = ndxGapStart + 1, ndxPIgap = ndxGapStart, ndxINTgap=ndxGapStart + 2

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  ! thorough
  Nparticles    = 5 * (10 ** 4)
  Nmixturedraws = 10 ** 4
  ! Nparticles    = 10 ** 5
  Nsim = 250

  ! solid
  Nparticles    = 10 ** 4
  Nmixturedraws = 10 ** 3
  Nsim = 100

  ! ! quick
  ! Nparticles    = 10 ** 3
  ! Nmixturedraws = 10 ** 2
  ! Nsim = 10

  ! quicker
  ! Nparticles    = 10 ** 2
  ! Nmixturedraws = 10 ** 2
  ! Nsim = 10

  datalabel        = 'spectreTB3MSGS020510OutputGapHeadline2018Q4'
  parameterlabel   = 'DEFAULT'


  call getarguments(Nparticles, Nsim, doIRFunitscale, doIRFactualObserver, ordercode, p, datalabel, parameterlabel) 
  IF (parameterlabel == 'DEFAULT')  parameterlabel = datalabel
  ! call getsettings(parameterlabel,Ny)

  ! read out order code
  ! - into ordervec
  ordervec = 0
  do i=1,Ny
     read(orderCode(i:i), '(i1)') ordervec(i)
  end do
  ! - find location of the policyrate, i.e. "3"
  yndxpolicyrate = minloc((ordervec - 3) ** 2,1)


  Nx = Nbar + Ngap * p
  NsigmaX = Nx * (Nx + 1) / 2

  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  call hrulefill


  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  print *, "Number of Threads:", NTHREADS

  ! VSL
  seed    = 0
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203, seed)  

  WRITE(*,'(a25, i20, i20)') 'LAUNCHING VSLSTREAM ', VSLstream%descriptor1, VSLstream%descriptor2
  print *, 'vsl_brng', vsl_brng_mt2203

  allocate(VSLarray(0:NTHREADS-1))
  ! VSL array
  DO i=0,NTHREADS-1
     errcode = vslnewstream(VSLarray(i), vsl_brng_mt2203 + 1 + i, 0)
     if (errcode /= 0) then
        print *, 'error init VSLarray', i
        stop 1
     else
        WRITE(*,'(a25, i4, i20, i20)') 'LAUNCHING VSLARRAY ', i, VSLarray(i)%descriptor1, VSLstream%descriptor2
        print *, 'vsl_brng', vsl_brng_mt2203 + 1 + i
     end if
  END DO
  ! runtime parameters :end: 

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = '.particles.VARlags' // trim(int2str(p)) // '.' // trim(datalabel) 
  if (doELB) fileXT = trim(fileXT) // '.ELB2'
  if (doIRFunitscale) filext = trim(filext) // '.irfUnitSV' 
  if (doIRFactualObserver) filext = trim(filext) // '.irfActualObs'
  fileXT = trim(fileXT) // '.order' // trim(orderCode) // '.dat'
  filext = '.APF' //  trim(filext)
  if (doTimeStamp) filext = '.' // timestampstr() //  trim(filext)

  parameterXT   = '.nopinoise.zlb.SPECTREVAR.' // trim(parameterlabel) // '.VARlags' // trim(int2str(p)) // '.order1234.ar1svcor.T236.Tjumpoff0.dat'

  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  Nf = Ngap * Ngap * p
  allocate(f(Nf))

  ! read data
  T = loft(datafile) 
  IF (T < 10) THEN
     print *, 'Less than 10 observations in input file!', datafile
     STOP 1
  END IF

  ALLOCATE (y(Ny,T), yNaN(Ny,T), atELB(0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Y)'
  END IF

  ! print *, 'trying to read', T, 'obs from', datafile
  CALL readdata(y,datafile,Ny,T, ordervec)
  CALL readnandata(yNaN,nandatafile,Ny,T, ordervec)

  if (doELB) then
     atELB(0) = .false.
     atELB(1:T) = yNaN(yndxpolicyrate,:)
  else
     atELB = .false.
  end if

  ! validate yNaN and y
  DO k=1,T
     DO i = 1, Ny
        if (yNaN(i,k) .AND. y(i,k) /= 0.0d0 ) then
           write (*,*) 'T=', k, ': YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  ! ---------------------------------
  ! LOAD PARAMETERS
  ! ---------------------------------

  ! 1) parameter files with quantiles

  ! f
  allocate(theta2(Nf, parameterColumns))
  filename = trim(adjustl(mcmcdir)) // 'F' // trim(adjustl(parameterXT))
  call loadmatrix(theta2, filename, Nf, parameterColumns)
  f = theta2(:,meanCol)
  deallocate(theta2)

  print *, 'f'
  print *, f
  call hrulefill


  ! shockslopes
  filename = trim(adjustl(mcmcdir)) // 'SHOCKSLOPES' // trim(adjustl(parameterXT))
  allocate(theta2(Ngapshockslopes, parameterColumns))
  call loadmatrix(theta2, filename, Ngapshockslopes, parameterColumns)
  gapshockslopes = theta2(:,meanCol) 
  deallocate(theta2)
  print *, 'gapshockslopes'
  print *, gapshockslopes
  call hrulefill

  ! sqrtVh
  sqrtVh  = 0.0d0
  filename = trim(adjustl(mcmcdir)) // 'HSIGMA' // trim(adjustl(parameterXT))
  allocate(theta2(NhgapSigma,parameterColumns))
  call loadmatrix(theta2, filename, NhgapSigma, parameterColumns)
  call ivech(sqrtVh(2:Nsv,2:Nsv), theta2(:,meanCol))
  ! recall: sqrtVhgap is upper-triangular draw from inverse wishart, i.e. Vhgap = sqrtVhgap * sqrtVhgap'
  deallocate(theta2)
  print *, 'sqrtVh'
  print *, sqrtVh
  call hrulefill


  ! 2) parameters stored as draws

  ! rbarvol

  ! horthvar
  filename =  trim(adjustl(mcmcdir)) // 'RBARVAR' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta1(Nmcmcdraws))
  call loadarray1(theta1, filename, Nmcmcdraws)
  if (meanCol == 1) then
     rbarvol = sum(sqrt(theta1)) / dble(Nmcmcdraws)
  else ! do median
     CALL dlasrt('I', Nmcmcdraws, theta1, status)
     rbarvol = sqrt(theta1(floor(real(Nmcmcdraws) * 0.5)))
  end if
  deallocate(theta1)
  print *, 'rbarvol'
  print *, rbarvol
  call hrulefill

  ! 2a) HORTH (just taking pibar column)

  ! horthvar
  filename =  trim(adjustl(mcmcdir)) // 'HORTHVAR' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta2(Nmcmcdraws,2))
  call loadmatrix(theta2, filename, Nmcmcdraws, 2)
  if (meanCol == 1) then
     sqrtVh(1,1) = sum(sqrt(theta2(:,1))) / dble(Nmcmcdraws)
  else ! do median
     CALL dlasrt('I', Nmcmcdraws, theta2(:,1), status)
     sqrtVh(1,1) = sqrt(theta2(floor(real(Nmcmcdraws) * 0.5),1))
  end if
  deallocate(theta2)
  print *, 'HVARpibar'
  print *, sqrtVh(1,1)
  call hrulefill

  ! horthrho
  filename =  trim(adjustl(mcmcdir)) // 'HORTHRHO' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta2(Nmcmcdraws,2))
  call loadmatrix(theta2, filename, Nmcmcdraws, 2)
  if (meanCol == 1) then
     hrho(1) = sum(theta2(:,1)) / dble(Nmcmcdraws)
  else ! do median
     CALL dlasrt('I', Nmcmcdraws, theta2(:,1), status)
     hrho(1) = theta2(floor(real(Nmcmcdraws) * 0.5),1)
  end if
  deallocate(theta2)
  print *, 'HRHOpibar'
  print *, hrho(1)
  call hrulefill

  ! horthbar
  filename =  trim(adjustl(mcmcdir)) // 'HORTHBAR' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta2(Nmcmcdraws,2))
  call loadmatrix(theta2, filename, Nmcmcdraws, 2)
  if (meanCol == 1) then
     hbar(1) = sum(theta2(:,1)) / dble(Nmcmcdraws)
  else ! do median
     CALL dlasrt('I', Nmcmcdraws, theta2(:,1), status)
     hbar(1) = theta2(floor(real(Nmcmcdraws) * 0.5),1)
  end if
  deallocate(theta2)
  print *, 'HBARpibar'
  print *, hbar(1)
  call hrulefill

  ! 2b) other SV
  ! hgapbar
  filename =  trim(adjustl(mcmcdir)) // 'HBAR' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta2(Nmcmcdraws, Ngap))
  call loadmatrix(theta2, filename, Nmcmcdraws, Ngap)
  if (meanCol == 1) then
     hbar(2:Nsv) = sum(theta2, 1) / dble(Nmcmcdraws)
  else ! do median
     do j=1,Ngap
        CALL dlasrt('I', Nmcmcdraws, theta2(:,j), status)
     end do
     hbar(2:Nsv) = theta2(floor(real(Nmcmcdraws) * 0.5),:)
  end if
  deallocate(theta2)
  print *, 'HGAPBAR'
  print *, hbar(2:Nsv)
  call hrulefill

  ! hgaprho
  filename =  trim(adjustl(mcmcdir)) // 'HRHO' // trim(adjustl(parameterXT))
  Nmcmcdraws = loft(filename)
  allocate(theta2(Nmcmcdraws, Ngap))
  call loadmatrix(theta2, filename, Nmcmcdraws, Ngap)
  if (meanCol == 1) then
     hrho(2:Nsv) = sum(theta2, 1) / dble(Nmcmcdraws)
  else ! do median
     do j=1,Ngap
        CALL dlasrt('I', Nmcmcdraws, theta2(:,j), status)
     end do
     hrho(2:Nsv) = theta2(floor(real(Nmcmcdraws) * 0.5),:)
  end if
  deallocate(theta2)
  print *, 'HGAPRHO'
  print *, hrho(2:Nsv)
  call hrulefill

  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data           = ' // trim(datalabel)
  print *, 'parameters     = ' // trim(parameterlabel)
  print *, 'Ny             = ', Ny
  print *, 'T              = ', T
  print *, 'Nparticles     = ', Nparticles
  print *, 'IRF Nsim       = ', Nsim
  print *, 'Nmixturedraws  = ', Nsim
  print *, 'p              = ', p
  print *, 'ordervec       = ', ordervec
  print *, 'yndxpolicyrate = ', yndxpolicyrate
  print *, 'with APF sampling'
  if (doIRFunitscale) print *, 'IRF: unit scale'
  if (doIRFactualobserver) print *, 'IRF: with actual observer'
  CALL HRULEFILL

  ! allocate memory for draws
  ALLOCATE (ELBtwist(Nparticles,1:T)) 
  ALLOCATE (DRAWytilde(Ny,Nparticles,1:T))
  ALLOCATE (DRAWshadowratehat(Nparticles,1:T))
  ALLOCATE (DRAWxgain(Nx,Ny,Nparticles,T), DRAWxgainstandardized(Nx,Ny,Nparticles,T), STAT=status)
  ALLOCATE (PARTICLEweights(Nparticles,0:T), DRAWxhat(Nx,Nparticles,0:T), DRAWsvol(Nparticles,Nsv,0:T), DRAWsqrtSigmaX(NsigmaX,Nparticles,0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws)'
  END IF
  ALLOCATE (DRAWxirf(Nx,0:irfNhorizons,irfNimpulse,1:T),DRAWyirf(irfNy,0:irfNhorizons,irfNimpulse,1:T), STAT=status) ! IRFweights(irfNimpulse,Nparticles,1:T)
  ALLOCATE (DRAWysimirf(irfNysim,0:irfNhorizons,2,irfNimpulse,1:T),  DRAWysimbaseline(irfNysim,0:irfNhorizons,1:T), STAT=status)
  ALLOCATE (DRAWinnovol(irfNimpulse,1:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (irf draws)'
  END IF

  ALLOCATE (loglike(T), STAT=status)
  loglike = 0.0d0

  PARTICLEweights = 1.0d0 / dble(Nparticles)
  DRAWxhat        = 0.0d0
  DRAWytilde      = 0.0d0
  DRAWshadowratehat  = 0.0d0
  DRAWsvol        = 0.0d0
  DRAWxgain       = 0.0d0
  DRAWxgainstandardized  = 0.0d0

  DRAWytilde      = 0.0d0

  ! IRFweights       = 1.0d0 / dble(Nparticles)
  DRAWxirf         = 0.0d0
  DRAWyirf         = 0.0d0
  DRAWinnovol      = 0.0d0
  DRAWysimirf      = 0.0d0
  DRAWysimbaseline      = 0.0d0
  ! DRAWsvirf        = 0.0d0

  call particlefilter(doIRFunitscale,doIRFactualObserver, elbound, Nyield, Nsim, T, Ny, y, yNaN, ordervec, atELB, loglike, Nparticles, PARTICLEweights, ELBtwist, irfNhorizons, irfNimpulse, DRAWinnovol, DRAWxirf, DRAWyirf, irfNy, DRAWysimirf, DRAWysimbaseline, irfNysim, DRAWytilde, DRAWshadowratehat, DRAWxhat, DRAWsqrtSigmaX, DRAWxgain,  DRAWxgainstandardized, Nx, NsigmaX, Nw, DRAWsvol, Nsv, p, f, Nf, gapshockslopes, Ngapshockslopes, hbar, hrho, sqrtVh, rbarvol, VSLstream, NTHREADS, VSLarray)


  CALL HRULEFILL
  WRITE (*,*) 'PARTICLE FILTER IS DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a60)') 'TIME: ', trim(adjustl(timestampstr()))
  WRITE(4,'(a20,a60)') 'Data: ', datalabel
  WRITE(4,'(a20,a60)') 'Parameters: ', parameterlabel
  WRITE(4,'(a80)') repeat('-',80)
  WRITE(4,'(a20,I10)') 'Nparticles: ', Nparticles
  WRITE(4,'(a20,I10)') 'Nsim: ', Nsim
  WRITE(4,'(a20,I10)') 'p: ', p
  CLOSE(UNIT=4)
  CALL HRULEFILL

  deallocate(f)

  ! ----------------------------------------------------------------------------
  ! STORE
  ! ----------------------------------------------------------------------------

  CALL HRULEFILL
  WRITE (*,*) 'STARTING W/STORAGE !!!'
  CALL HRULEFILL

  ! STORE ESTIMATES
  ! Note: manual reshape avoids segmentation faults

  filename = 'YDATA' // filext
  call savemat(y, filename)

  filename = 'YNAN' // filext
  call savematlogical(yNaN, filename)


  ! LIKELIHOOD
  filename = 'LOGLIKE' // filext
  call savevec(loglike, filename)
  call hrulefill
  WRITE (*,*) 'STORED LFF'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T
  call hrulefill
  ! DEALLOCATE (loglike) -- occurs later, output is repeated toward end of file

  ! store ytilde
  ALLOCATE(theta2(T,Ny))
  filename  = 'YTILDE' // filext
  forall (j=1:T,i=1:Ny) theta2(j,i) = sum(PARTICLEweights(:,j) * DRAWytilde(i,:,j))
  call savemat(theta2, filename)
  WRITE (*,*) 'STORED YTILDE'
  DEALLOCATE(theta2)

  ALLOCATE (theta1(T), STAT=status)

  ! ESS
  filename  = 'ESS' // filext
  theta1    = 1 / sum(PARTICLEweights(:,1:T) ** 2, 1) / Nparticles 
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ESS'

  ! ELBtwist
  filename  = 'ELBTWIST' // filext
  ! call savemat(ELBTWIST, filename)
  theta1    = sum(ELBTWIST, 1)
  forall (j=1:T,i=1:Nparticles) ELBTWIST(i,j) = ELBTWIST(i,j) / theta1(j)
  theta1 = 1 / sum(ELBTWIST(:,1:T) ** 2, 1) / Nparticles 
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ELBTWIST'

  ! store trends
  DO i=1,Nx
     filename  = 'STATEHAT' // trim(int2str(i)) // filext
     theta1 = sum(PARTICLEweights(:,1:T) * DRAWxhat(i,:,1:T), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED STATEHAT', i
  END DO

  ! store shadowrate
  filename  = 'SHADOWRATEHAT' // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWshadowratehat(:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED SHADOWRATEHAT'

  DEALLOCATE (theta1)

  ! draw distribution for linear states
  ALLOCATE (ndx(Nmixturedraws,T),Xdraws(Nx,Nmixturedraws,T))

  ! sample particle indices (note: these will also be used later for the SV particles)
  print *, 'Drawing particle indices  ...'
  DO j=1,T
     call drawNDXpdf(ndx(:,j), Nmixturedraws, PARTICLEweights(:,j), Nparticles, VSLstream)
  END DO
  print *, 'Done drawing particle indices.'

  print *, 'Drawing Xdraws normals ...'
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nx * Nmixturedraws, Xdraws, 0.0d0, 1.0d0)
  print *, 'done drawing Xdraws normals.'

  allocate(shadowratedraws(Nmixturedraws,T))

  ! draw states from their joint distribution
  ! note: these are first proposals, need to apply rejection sampling at ELB in next step
  !$OMP PARALLEL DO SHARED(ndx,Nx,DRAWxhat,DRAWsqrtSigmaX,Xdraws,shadowratedraws,Nmixturedraws) PRIVATE(k,i,sqrtSigmaX)
  DO j=1,T

     allocate(sqrtSigmaX(Nx,Nx))

     DO k=1,Nmixturedraws
        sqrtSigmaX = ivechU(DRAWsqrtSigmaX(:,ndx(k,j),j),Nx)
        call DTRMV('u','t','n',Nx,sqrtSigmaX,Nx,Xdraws(:,k,j),1)
     END DO ! k
     forall (i=1:Nx,k=1:Nmixturedraws) Xdraws(i,k,j) = DRAWxhat(i,ndx(k,j),j) + Xdraws(i,k,j)

     forall (k=1:Nmixturedraws) shadowratedraws(k,j) = Xdraws(ndxPIbar,k,j) + Xdraws(ndxRbar,k,j) + Xdraws(ndxINTgap,k,j)

     deallocate(sqrtSigmaX)
  END DO ! j
  !$OMP END PARALLEL DO

  ! check for ELB and apply rejection sampling as necessary
  allocate(zdraw(Nx))
  DO j=1,T

     if (atELB(j)) then

        ! ! reject draws that do not conform to ELB
        do while (any(shadowratedraws(:,j) > elbound))

           ! loop over all proposals and resample if necessary
           ! note: could do this with a where statement, but the sampling would be less efficient
           do k=1,Nmixturedraws
              if (shadowratedraws(k,j) > elbound) then

                 ! todo: drawing just Nx zdraw at a time is not very efficient (though the most direct way to implement the algorithm)
                 errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx, zdraw, 0.0d0, 1.0d0)

                 Xdraws(:,k,j) = constructXdraw(Nx, DRAWxhat(:,ndx(k,j),j), DRAWsqrtSigmaX(:,ndx(k,j),j), zdraw)

                 shadowratedraws(k,j) = Xdraws(ndxPIbar,k,j) + Xdraws(ndxRbar,k,j) + Xdraws(ndxINTgap,k,j)

              end if
           end do

        end do


     end if

  END DO
  deallocate(zdraw)

  ! TRUNCATED SHADOWRATE DRAWS AT ELB
  filename  = 'SHADOWRATE' // filext
  CALL storeEstimatesTranspose(shadowratedraws,T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED SHADOWRATE'
  ! filename  = 'SHADOWRATE.DRAWS' // filext
  ! call savemat(shadowratedraws, filename)
  ! WRITE (*,*) 'STORED SHADOWRATE DRAWS'

  ! just checking: are all shadowrate draws admissible?
  allocate(shadowrateELBmask(Nmixturedraws,T))
  forall(j=1:T,k=1:Nmixturedraws) shadowrateELBmask(k,j) = (.NOT. atELB(j)) .OR. (shadowratedraws(k,j) < elbound) 
  if (.not. all(shadowrateELBmask)) then
     ! todo: report rejection rate
     do j=1,T
        print *, j, count(shadowrateELBmask(:,j))
     end do
     print *, 'check the shadow rate draws! they do not seem to satisfy the ELB data'
     stop 1
  end if
  deallocate(shadowrateELBmask)

  ! store other states
  allocate(theta2(Nmixturedraws,T))
  ! trend nominal shadow rate
  filename  = 'NOMINALRATETREND' // filext
  forall (j=1:T,k=1:Nmixturedraws) theta2(k,j) = Xdraws(ndxPIbar,k,j) + Xdraws(ndxRbar,k,j) 
  CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED NOMINALRATETREND'

  ! yeild trends
  DO i=1,Nyield
     filename  = 'NOMINALLONGRATE' // trim(int2str(i)) // 'TREND' // filext
     forall (j=1:T,k=1:Nmixturedraws) theta2(k,j) = Xdraws(ndxPIbar,k,j) + Xdraws(ndxRbar+i,k,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED NOMINALLONGRATETREND', i
  END DO

  ! avg term premia
  DO i=1,Nyield
     filename  = 'AVGTERMPREMIA' // trim(int2str(i))  // filext
     forall (j=1:T,k=1:Nmixturedraws) theta2(k,j) = Xdraws(ndxRbar+i,k,j) - Xdraws(ndxRbar,k,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED AVGTERMPREMIA', i
  END DO


  ! Store Quantiles of State Draws
  ! NOTE: do this only after storing linear combinations of the states, since storeEstimatesTranspose sorts the draws! (MKL side effects ...)
  DO i=1,Nx
     filename  = 'STATE' // trim(int2str(i)) // filext
     CALL storeEstimatesTranspose(Xdraws(i,:,:),T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED STATE', i
  END DO


  deallocate(theta2)

  deallocate(DRAWytilde)

  DEALLOCATE (Xdraws, DRAWxhat)
  DEALLOCATE (DRAWshadowratehat)
  DEALLOCATE (DRAWsqrtSigmaX)


  ! TODO: OK *not* to mask the draws for Kgain and SV ?
  ! implicitly, the rejection sampling is thus applies to the normal shocks to the mixture, not the draws from the mixture grid itself

  ! 2D Gain Matrices
  if (doGains) then
     ! store analytical moments of gain
     ALLOCATE (theta2(T,Ny))
     ! trend gains 
     DO i=1,Nx
        ! filename  = 'GAINSTATEHAT' // trim(int2str(i)) 
        ! if (doNormalizedGains) then
        !    filename = trim(filename) // '.normalized'
        ! end if
        ! filename  = trim(filename) // filext
        filename  = 'GAINSTATEHAT' // trim(int2str(i)) // filext
        FORALL (j=1:Ny) theta2(:,j) = sum(PARTICLEweights(:,1:T) * DRAWxgain(i,j,:,:), 1)
        call savemat(theta2, filename)
        filename  = 'GAINSTATEHAT' // trim(int2str(i)) // '.standardized' // filext
        FORALL (j=1:Ny) theta2(:,j) = sum(PARTICLEweights(:,1:T) * DRAWxgainstandardized(i,j,:,:), 1)
        call savemat(theta2, filename)
        WRITE (*,*) 'STORED GAINSTATEHAT', i
     END DO
     DEALLOCATE (theta2)
  end if ! doGains
  DEALLOCATE(DRAWxgain,DRAWxgainstandardized)


  ! 2) Nparticle draws for the other particles
  ALLOCATE (theta1(T), STAT=status)
  DO i=1,Nsv
     filename  = 'SVHAT' // trim(int2str(i)) // filext
     theta1 = sum(PARTICLEweights(:,1:T) * DRAWsvol(:,i,1:T), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED SVHAT', i
  END DO


  ! draw distribution
  ALLOCATE (theta2(Nmixturedraws,T))
  DO i=1,Nsv
     filename  = 'SV' // trim(int2str(i)) // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWsvol(ndx(k,j),i,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED SV', i
  END DO

  DEALLOCATE(theta1)
  DEALLOCATE(theta2)
  DEALLOCATE(ndx)
  DEALLOCATE(DRAWsvol)

  ! ----------------------------------------------------------------------------
  ! STORE IRF
  ! ----------------------------------------------------------------------------

  ! output IRF in separate files per state and per shock
  ! each file reports T x irfNhorizons matrix of IRF coefficients

  DO i=1,Nx
     DO k=1,irfNimpulse
        filename = 'IRFx' // trim(int2str(i)) // 'impulse' // trim(int2str(k)) 
        filename = trim(filename) // filext

        call savemat(transpose(DRAWxirf(i,:,k,:)), filename)
        WRITE (*,*) 'STORED IRF of X', i, 'TO IMPULSE ', k
     END DO
  END DO


  DO i=1,irfNy
     DO k=1,irfNimpulse
        filename = 'IRFy' // trim(int2str(i)) // 'impulse' // trim(int2str(k))
        filename = trim(filename) // filext

        call savemat(transpose(DRAWyirf(i,:,k,:)), filename)
        WRITE (*,*) 'STORED IRF of Y', i, 'TO IMPULSE ', k
     END DO
  END DO

  ! simulated IRF 
  DO i=1,irfNysim

     filename = 'IRFysim' // trim(int2str(i)) // 'baseline'  // filext
     call savemat(transpose(DRAWysimbaseline(i,:,:)), filename)
     WRITE (*,*) 'STORED BASELINE PATH of SIM-Y', i

     DO k=1,irfNimpulse
        filename = 'IRFysim' // trim(int2str(i)) // 'impulse' // trim(int2str(k)) // '.plus'
        filename = trim(filename) // filext
        call savemat(transpose(DRAWysimirf(i,:,1,k,:)), filename)

        filename = 'IRFysim' // trim(int2str(i)) // 'impulse' // trim(int2str(k)) // '.minus'
        filename = trim(filename) // filext
        call savemat(transpose(DRAWysimirf(i,:,2,k,:)), filename)

        WRITE (*,*) 'STORED IRF of SIM-Y', i, 'TO IMPULSE ', k
     END DO
  END DO

  filename = 'IRFinnovols' // filext
  call savemat(transpose(DRAWinnovol), filename)

  ! i = irfNimpulse
  ! filename = 'IRFsv' // trim(int2str(i)) // 'impactELB'
  ! filename = trim(filename) // filext
  ! call savemat(transpose(DRAWsvirf), filename)


  ! DO k=1,irfNimpulse
  !    filename = 'IRFweight' // trim(int2str(k)) // filext
  !    call savemat(IRFweights(i,:,:), filename)
  !    WRITE (*,*) 'STORED IRFWEIGHT', k
  ! END DO

  deallocate(DRAWxirf)
  ! print *, 'done xirf'
  deallocate(DRAWyirf)
  deallocate(DRAWinnovol)
  ! print *, 'done yirf'
  deallocate(DRAWysimirf)
  ! print *, 'done ysimirf'
  deallocate(DRAWysimbaseline)
  ! print *, 'done ysimbaseline'
  ! deallocate(IRFweights)
  ! print *, 'done irfweights'

  ! deallocate(DRAWsvirf)

  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE FILTER
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  ! CLEANUP FILTER
  ! ----------------------------------------------------------------------------

  call hrulefill
  WRITE (*,*) 'LFF:'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T
  call hrulefill
  DEALLOCATE (loglike)

  ! ----------------------------------------------------------------------------
  ! FINAL CLEANUP
  ! ----------------------------------------------------------------------------


  DEALLOCATE (PARTICLEweights)
  DEALLOCATE (ELBtwist)
  DEALLOCATE (y, yNaN)

  ! VSLstreams
  errcode = vsldeletestream(VSLstream)     
  do i=0,NTHREADS-1
     errcode = vsldeletestream(VSLarray(i))
  end do
  deallocate(VSLarray)

  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  FUNCTION constructXdraw(Nx, xhat, sqrtSigmaXvec, zdraw) result(thisXdraw)

    integer, intent(in) :: Nx
    double precision, dimension(:), intent(in) :: xhat, zdraw, sqrtSigmaXvec
    double precision, dimension(Nx) :: thisXdraw

    double precision, dimension(Nx,Nx) :: sqrtSigmaX

    thisXdraw = zdraw
    sqrtSigmaX = ivechU(sqrtSigmaXvec,Nx)
    call DTRMV('u','t','n',Nx,sqrtSigmaX,Nx,thisXdraw,1)
    forall (i=1:Nx) thisXdraw(i) = xhat(i) + thisXdraw(i)

  END FUNCTION constructXdraw

  SUBROUTINE getarguments(Nparticles, Nsim, doIRFunitscale, doIRFactualObserver, orderCode, p, datalabel,parameterlabel)

    CHARACTER (LEN=100), intent(inout) :: orderCode
    CHARACTER (LEN=100), intent(inout) :: datalabel,parameterlabel
    logical, intent(inout) :: doIRFactualobserver,doIRFunitscale
    INTEGER, intent(inout) :: p, Nparticles, Nsim
    INTEGER :: counter, dummy
    CHARACTER(len=32) :: arg

    counter = 0
    ! IF (command_argument_count() == 0) THEN
    !    print *, 'WARNING. No Datalabel specified!'
    !    print *, 'Using default: ' // datalabel
    ! END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Nparticles
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Nsim
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doIRFunitscale = .true.
       else
          doIRFunitscale = .false.
       end if
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doIRFactualObserver = .true.
       else
          doIRFactualObserver = .false.
       end if
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, orderCode) 
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') p
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, datalabel) 
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, parameterlabel)
    END IF

  END SUBROUTINE getarguments


  SUBROUTINE readdata(y,filename,Ny,T,ordervec)
    IMPLICIT NONE

    INTENT(IN) :: filename,Ny,T,ordervec
    INTENT(INOUT) :: y
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    DOUBLE PRECISION, DIMENSION(:,:) :: y
    INTEGER :: i, T, Ny
    INTEGER :: ordervec(Ny)
    double precision :: tmp(Ny)

    fmtstr = es30d16(Ny)
    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')
    DO i=1,T
       READ(4,fmtstr) tmp
       y(:,i) = tmp(ordervec)
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readdata

  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  SUBROUTINE readnandata(nanny,filename,Ny,T,ordervec)
    IMPLICIT NONE

    INTENT(IN) :: filename,T,Ny,ordervec
    INTENT(INOUT) :: nanny
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    LOGICAL, DIMENSION(:,:) :: nanny
    INTEGER :: work(Ny), rework(Ny), ordervec(Ny)

    INTEGER :: i, j, T, Ny

    fmtstr = '(I2' // repeat(',I2', Ny-1) // ')'

    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    DO i=1,T
       READ(4,fmtstr) (work(j), j=1,Ny)
       rework = work(ordervec)
       WHERE (rework == 1) 
          nanny(:,i) = .TRUE.
       ELSEWHERE
          nanny(:,i) = .FALSE.
       END WHERE
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readnandata
  ! -----------------------------------------------------------------


END PROGRAM main
! -----------------------------------------------------------------

! @\newpage\subsection{particlefilter}@
SUBROUTINE particlefilter(doIRFunitscale, doIRFactualObserver, elbound, Nyield, Nsim, T, Ny, y, yNaN, ordervec, atELB, loglike, Nparticles, PARTICLEweights, ELBtwist, irfNhorizons, irfNimpulse, DRAWinnovol, DRAWxirf, DRAWyirf, irfNy, DRAWysimirf, DRAWysimbaseline, irfNysim, DRAWytilde, DRAWshadowratehat, DRAWxhat, DRAWsqrtSigmaX, DRAWxgain, DRAWxgainstandardized, Nx, NsigmaX, Nw, DRAWsvol, Nsv, p, f, Nf, gapshockslopes, Ngapshockslopes, hbar, hrho, sqrtVh, rbarvol, VSLstream, NTHREADS, VSLarray)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery
  use densitybox, only : simySVCORdraws, meanvarianceTruncNormalUB ! meanTruncNormalUB

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  logical, intent(in) :: doIRFactualObserver, doIRFunitscale
  double precision, intent(in) :: elbound
  integer, intent(in) :: Nsim

  INTENT(OUT) :: PARTICLEweights, DRAWxhat, DRAWsqrtSigmaX, DRAWxgain, DRAWxgainstandardized, DRAWsvol
  INTENT(INOUT) :: DRAWytilde
  INTENT(OUT) :: loglike
  INTENT(OUT) :: ELBtwist  
  INTENT(IN) :: T,Ny,y,yNaN, atELB, Nx, NsigmaX, Nw, Nsv, Nparticles
  INTENT(IN) ::  p, f, Nf, gapshockslopes, Ngapshockslopes, hbar, hrho, sqrtVh
  INTENT(IN) ::  rbarvol
  INTENT(OUT) :: DRAWxirf, DRAWyirf
  INTENT(OUT) :: DRAWshadowratehat
  INTENT(OUT) :: DRAWysimirf, DRAWysimbaseline
  INTENT(OUT) :: DRAWinnovol
  INTENT(IN) :: irfNhorizons, irfNimpulse, irfNysim, irfNy


  INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw
  INTEGER :: ii, nn, lag


  integer, intent(in), dimension(Ny) :: ordervec

  ! OPEN MP
  INTEGER :: TID = 0


  type(progresstimer) :: timer

  ! DOUBLE PRECISION, DIMENSION(Ny-1) :: horizons

  double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  LOGICAL, DIMENSION(0:T) :: atELB
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nparticles,1:T) :: ELBtwist ! note: allocation only per time 1
  DOUBLE PRECISION, DIMENSION(Ny,Nparticles,1:T)  :: DRAWytilde
  DOUBLE PRECISION, DIMENSION(Nparticles,1:T)  :: DRAWshadowratehat
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nx,Ny,Nparticles,T) :: DRAWxgain
  DOUBLE PRECISION, DIMENSION(Nx,Ny,Nparticles,T) :: DRAWxgainstandardized
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,0:T) :: DRAWsvol ! NOTE:Nparticles in 1st dimension !!!

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), llf(Nparticles)
  ! DOUBLE PRECISION :: xdraws(Nx,Nparticles)

  DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
  DOUBLE PRECISION :: kernelsum, loglikemax
  DOUBLE PRECISION, DIMENSION(T) :: loglike

  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), logdetSigmaY

  DOUBLE PRECISION :: h(Nsv, Nparticles), SVol(Nsv,Nparticles), hshock(Nsv,Nparticles)

  ! IRF
  INTEGER :: irfNhorizons, irfNimpulse, irfNysim, irfNy
  DOUBLE PRECISION, DIMENSION(Nx,0:irfNhorizons,irfNimpulse,1:T)  :: DRAWxirf
  DOUBLE PRECISION, DIMENSION(irfNy,0:irfNhorizons,irfNimpulse,1:T)  :: DRAWyirf
  DOUBLE PRECISION, DIMENSION(irfNysim,0:irfNhorizons,irfNimpulse,2,1:T)  :: DRAWysimirf
  DOUBLE PRECISION, DIMENSION(irfNysim,0:irfNhorizons,1:T)  :: DRAWysimbaseline
  DOUBLE PRECISION, DIMENSION(irfNimpulse, 1:T)  :: DRAWinnovol
  ! DOUBLE PRECISION, DIMENSION(Nsv,1:T) :: DRAWsvirf
  ! DOUBLE PRECISION, DIMENSION(irfNimpulse,Nparticles,1:T)  :: IRFweights
  DOUBLE PRECISION, DIMENSION(Nx,0:irfNhorizons,Nparticles)  :: Xbaseline
  DOUBLE PRECISION, DIMENSION(Nx,0:irfNhorizons)  :: XbaselineHat ! integral over particles
  ! DOUBLE PRECISION, DIMENSION(Nx,0:irfNhorizons,Ny,Nparticles)  :: Xresponse
  DOUBLE PRECISION, DIMENSION(Nx,0:irfNhorizons,irfNimpulse)  :: XresponseHat ! integral over particles
  DOUBLE PRECISION :: Cobserver(Ny,Nx), Csim(irfNysim,Nx), Cirf(irfNy,Nx)
  LOGICAL, DIMENSION(Ny) :: thisYnan
  LOGICAL, DIMENSION(irfNysim) :: ytrunc

  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny) ! , Kgain(Nx,Ny,Nparticles)
  DOUBLE PRECISION :: Kgain(Nx,Ny) ! new code: Kgain shall now be particle specific, stoarge direclty in DRAWxgain
  DOUBLE PRECISION, DIMENSION(NsigmaX, Nparticles) :: vecSqrtSigmaX

  ! qr
  DOUBLE PRECISION :: qrR(Ny+Nx+Nw,Ny+Nx), qrRelb(1+Nx,1+Nx), qrRlyap(Nx + 1,Nx)
  INTEGER :: qrLwork, qrLelb, qrLlyap

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), ytilde(Ny) ! , sqrtR(Ny,Ny)

  ! SV parameters
  DOUBLE PRECISION, DIMENSION(Nsv) :: hbar, hrho, hintercept
  DOUBLE PRECISION, DIMENSION(Nsv,Nsv) :: sqrtVh !, Ah, sqrtVh0
  DOUBLE PRECISION :: minSVh(Nsv)

  ! Other parameters
  integer :: Nf
  DOUBLE PRECISION :: f(Nf)
  integer :: Ngapshockslopes
  DOUBLE PRECISION :: gapshockslopes(Ngapshockslopes)
  DOUBLE PRECISION :: rbarvol

  INTEGER :: Nynonan

  ! CHARACTER (LEN=200) :: filename

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  integer, intent(in) :: NTHREADS 
  type (vsl_stream_state), intent(inout) :: VSLstream
  type (vsl_stream_state), intent(inout), dimension(0:NTHREADS-1) :: VSLarray


  double precision :: uniformdraws(2,T) ! two draws: one for each APF step

  ! index variables for state space
  INTEGER, INTENT(IN) :: Nyield
  INTEGER :: Ngap
  INTEGER :: ndxPIbar, ndxRbar, ndxGapStart, ndxGapStop, ndxGapCompanionStop, ndxUgap, ndxPIgap, ndxINTgap
  INTEGER :: shockndxPIbar, shockndxGapStart, shockndxGapStop, shockndxRBAR 
  ! INTEGER, PARAMETER :: yndxPolicyrate = 3
  double precision, dimension(Nx) :: Cshadowrate, shadowvolvec
  double precision, dimension(Nparticles) :: critvalues, ELBlike
  double precision, dimension(Nparticles) :: shadowrateMean, shadowrateMeanAtELB, shadowrateVol, shadowrateVolAtELB
  double precision, dimension(Nx) :: shadowrategain

  ! double precision :: critval, elbweight

  integer :: offset, these

  Ngap         = Ny
  ndxPIbar     = 1
  ndxRbar      = 2
  ndxGapStart  = ndxRbar + Nyield + 1
  ndxGapStop   = ndxGapStart + Ny - 1
  ndxGapCompanionStop = ndxGapStart - 1 + Ngap * p
  if (ndxGapCompanionStop /= Nx) then
     print *, 'state space indices are off'
     stop 1
  end if
  ndxPIgap  = ndxGapStart
  ndxUgap   = ndxGapStart + 1
  ndxINTgap = ndxGapStart + 2 

  shockndxPIbar    = 1
  shockndxGapStart = 2
  shockndxGapStop  = shockndxGapStart - 1 + Ngap
  shockndxrbar     = Nw


  minSVh   = log(0.001d0 ** 2)

  ! prior over initial states
  Ex0           = 0.0d0
  call eye(sqrtVx0, 1.0d1)

  ! joint prior over pibar, and ibar; transformed into prior over pibar and rbar
  Ex0(ndxPIbar)     = 2.0d0 ! rough level for pibar
  Ex0(ndxRbar)      = 2.0d0 ! rough level for rbar
  forall (i = 1 : Nyield) Ex0(ndxRbar+i)  = 2.0d0 + dble(i) * 0.5d0 ! rough level for rbar


  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxPIbar,ndxPIbar) = 1.0d0
  A(ndxRbar,ndxRbar) = 1.0d0
  forall (j=1:Nyield) A(ndxRbar+j,ndxRbar+j) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:Ngap*(p-1)) A(ndxGapStop+j,ndxGapStart-1+j) = 1.0d0
  END IF
  ! fill in transition coefficients
  A(ndxGapStart:ndxGapStop,ndxGapStart:ndxGapCompanionStop) = transpose(reshape(f, (/ Ngap * p, Ngap /)))

  ! B
  B                           = 0.0d0
  ! unit loading on SV
  B(ndxPIbar,shockndxPIbar)   = 1.0d0
  ! constant rbarvol parameter for trend real rates
  B(ndxRbar,shockndxRbar)     = rbarvol
  forall (j=1:Nyield) B(ndxRbar+j,shockndxRbar) = rbarvol

  ! copy over to Bsv (for non SV shock columns)
  Bsv = B

  ! unit loadings on SV for gap shocks
  forall (j=0:Ngap-1) B(ndxGapStart+j,shockndxGapStart+j) = 1.0d0
  ! slopes
  offset = 0
  DO i = 2, Ngap 
     these = i-1 
     ! slopes in row i have index offsetslopes+1:offsetslopes + these
     B(ndxGapStart-1+i,shockndxGapStart:shockndxGapStart-1+these) = gapshockslopes(offset+1:offset+these)
     offset = offset + these
  END DO



  ! first setup Cobserver (2D),then reorder it, then copy over into C (3D)

  ! Cobserver -- setup via the usual ordering
  Cobserver         = 0.0d0
  ! ugap
  Cobserver(1,ndxUgap)       = 1.0d0
  ! inflation
  Cobserver(2,ndxPIBAR)      = 1.0d0
  Cobserver(2,ndxPIgap)      = 1.0d0
  ! policy rate
  Cobserver(3,ndxPIBAR)      = 1.0d0
  Cobserver(3,ndxRbar)       = 1.0d0
  Cobserver(3,ndxINTgap)     = 1.0d0
  ! copy into Cshadowrate vector
  Cshadowrate = Cobserver(3,:)

  ! long rates
  forall (j=1:Nyield) Cobserver(3+j,ndxPIBAR)    = 1.0d0
  forall (j=1:Nyield) Cobserver(3+j,ndxRbar+j)   = 1.0d0
  forall (j=1:Nyield) Cobserver(3+j,ndxINTgap+j) = 1.0d0

  ! setup Cirf
  ! for long-rates: neglect average term premia by ripping off the short-rate rbar
  Cirf          = 0.0d0
  Cirf(1:3,:)   = Cobserver(1:3,:)
  ! long rates
  ! forall (j=1:Nyield) Cirf(3+j,:)  = Cobserver(3+j,:)

  forall (j=1:Nyield) Cirf(3+j,ndxPIBAR)         = 1.0d0
  forall (j=1:Nyield) Cirf(3+j,ndxRbar)          = 1.0d0  ! this line is different from Cobserver
  forall (j=1:Nyield) Cirf(3+j,ndxINTgap+j)      = 1.0d0

  Cirf(Ny+1,ndxPIBAR)  = 1.0d0 
  Cirf(Ny+2,ndxRBAR)   = 1.0d0

  ! setup Csim
  Csim               = 0.0d0
  Csim(1:irfNy,:)    = Cirf
  Csim(irfNy+1,:)    = Cshadowrate
  Csim(irfNy+2,:)    = Cirf(4,:)
  Csim(irfNy+3,:)    = Cirf(5,:)
  Csim(irfNy+4,:)    = Cirf(6,:)
  ytrunc             = .false.
  ytrunc(irfNy+1:irfNysim)   = .true.


  ! setup C
  forall (j=1:T) C(:,:,j) = Cobserver
  ! zero out missing obs
  DO j=1,T
     DO i = 1, Ny
        if (yNaN(i,j)) then 
           C(i,:,j)     = 0.0d0 
        end if
     END DO
  END DO

  ! Cobserver -- reorder row as desired
  Cobserver = Cobserver(ordervec,:)

  ! Time 0 particles

  ! SV0 -- set equal to hbar (analogous to MCMC estimation)
  ! errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h, 0.0d0, 1.0d0)
  ! Ah = 0.0d0
  ! forall (i=1:Nsv) Ah(i,i) = hrho(i)
  ! call dlyapsqrt(sqrtVh0, Ah, sqrtVh, Nsv, Nsv, errcode)
  ! if (errcode /= 0) then
  !    print *, 'PF-LIKE: error computing ergodic VCV of log-SV'
  !    call savemat(Ah, 'Ah.debug')
  !    call savemat(sqrtVh, 'sqrtVh.debug')
  !    stop 1
  ! end if

  ! ! h = sqrtVh0' * h
  ! call DTRMM('l','u','t','n',Nsv,Nparticles,1.0d0,sqrtVh0,Nsv,h,Nsv)

  ! forall (k=1:Nparticles,i=1:Nsv) h(i,k) = hbar(i) + h(i,k)

  forall (k=1:Nparticles,i=1:Nsv) h(i,k) = hbar(i)
  SVol = exp(h * 0.5d0)

  ! prepare log-SV intercept
  hintercept = (1 - hrho) * hbar
  ! print *, hintercept
  ! stop 11

  ! RB priors for linear states
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0

  ! prepare prior variance of linear states
  sqrtVx0         = transpose(sqrtVx0) ! needs to be upper triangular RIGHT factor henceforth
  vecSqrtSigmaX   = 0.0d0

  DO k=1,Nparticles
     vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)
  END DO



  FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,0)  = vecSqrtSigmaX(i,k)
  FORALL(i=1:Nx,k=1:Nparticles)  DRAWxhat(i,k,0)            = xposterior(i,k) 
  FORALL (i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,0) = SVol(i,k)  ! note the transpose

  PARTICLEweights = 1.0d0 / dble(Nparticles)
  ELBtwist        = 1.0d0 ! ratio of weights generated by ELB

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR     = 0.0d0
  qrlwork = qrquery(qrR)
  qrRelb  = 0.0d0
  qrlelb  = qrquery(qrRelb)
  qrRlyap = 0.0d0
  qrLlyap = qrquery(qrRlyap)

  CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: IRF at time j
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! 1): Compute baseline expectations (from last period)
     !$OMP PARALLEL DO SHARED(xposterior, Xbaseline, Nx, A, irfNhorizons, Nparticles) SHARED(ndxUGap, ndxPIbar, ndxPIgap, ndxRbar) PRIVATE(lag) DEFAULT(NONE) SCHEDULE(STATIC)
     DO k=1,Nparticles

        lag = 0
        ! recall: time zero is one-step ahead for the baseline
        call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,Xbaseline(:,lag,k),1)

        DO lag=1,irfNhorizons
           call DGEMV('n',Nx,Nx,1.0d0,A,Nx,Xbaseline(:,lag-1,k),1,0.0d0,Xbaseline(:,lag,k),1)
        END DO

     END DO ! k
     !$OMP END PARALLEL DO

     ! integrate baseline across particles (recall: equally weighted particles b/o resampling at every step)
     XbaselineHat = sum(Xbaseline,3) / Nparticles

     ! draw new SVol proposals

     errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hshock, 0.0d0, 1.0d0)
     call DTRMM('l','u','n','n',Nsv,Nparticles,1.0d0,sqrtVh,Nsv,hshock,Nsv) ! recall: sqrtVh is left-upper choleski factor
     forall (i=1:Nsv,k=1:Nparticles) SVol(i,k) = exp((hintercept(i) + hrho(i) * h(i,k) + hshock(i,k)) * 0.5d0)


     ! 2) call IRF routine for each shock
     XresponseHat = 0.0d0

     if (doIRFactualObserver) then
        Cobserver = C(:,:,j)
        thisYnan  = yNaN(:,j)
     else
        thisYnan  = .false.
     end if

     call particleIRF(doIRFunitscale, XresponseHat(:,:,1:Ny), DRAWinnovol(:,j), Nparticles, Nx, Ny, irfNimpulse, irfNhorizons, xposterior, vecSqrtSigmaX, Nw, A, B, Cobserver, Nsv, SVol) 
     ! note: just applied to Ny (not irfNimpulse) impulses

     ! 3a) Compute IRF as difference between response and baseline
     forall (ii=1:Nx,lag=0:irfNhorizons,nn=1:irfNimpulse) DRAWxirf(ii,lag,nn,j) = XresponseHat(ii,lag,nn)  - XbaselineHat(ii,lag)

     ! 4) compute yirf
     !  using loops of DGEMMs Cobserver * X is faster than manual forall construction
     do nn=1,irfNimpulse
        call DGEMM('n','n',irfNy,irfNhorizons+1,Nx,1.0d0,Cirf,irfNy,DRAWxirf(:,:,nn,j),Nx,0.0d0,DRAWyirf(:,:,nn,j),irfNy)
     end do

     ! ALTERNATIVE: call simulated IRF
     ! note: considering only standardized shocks, thus no SV signal and equally weighted IRF
     ! note: IRF are for next-period shocks (xposterior already accounts for jumpoff at j-1)
     call particleIRFsim(doIRFunitscale, Nsim, DRAWysimbaseline(:,:,j), DRAWysimirf(:,:,:,:,j), thisYnan, ytrunc, elbound, Nparticles, Nx, Ny, irfNimpulse, irfNhorizons, irfNysim, xposterior, vecSqrtSigmaX, Nw, A, B, Cobserver, Csim, Nsv, h, hbar, hrho, sqrtVh, VSLstream, NTHREADS, VSLarray)

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! DONE: IRF at time j
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: APF RESAMPLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     Nynonan = count(.not. yNaN(:,j))

     ! 1) Construct proposal for SVol at median (i.e. forecast of AR1 log-SV)
     forall (i=1:Nsv,k=1:Nparticles) SVol(i,k) = exp((hintercept(i) + hrho(i) * h(i,k)) * 0.5d0)
     ! OR: stick with previous SV
     ! forall (i=1:Nsv,k=1:Nparticles) SVol(i,k) = exp(h(i,k) * 0.5d0)



     !$OMP PARALLEL DO SHARED(xposterior, vecSqrtSigmaX, SVol, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nw, Nynonan, h, hrho, hintercept,A,B,C,Cshadowrate,critvalues,elbound,atELB,shadowrateMean,shadowrateVol) FIRSTPRIVATE(qrLWORK) PRIVATE(ytilde, logdetSigmaY, errcode, xprior, sqrtSigmaX, sqrtSigmaY, qrR, shadowvolvec, nn) FIRSTPRIVATE(Bsv) DEFAULT(NONE) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        ! 2) Fill Particles into state space
        ! Bsv
        FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * SVol(i,k) 
        ! Note: the rbarshock column has been preserved by passing Bsv as FIRSTPRIVATE to the OMP loop

        ! 3) Kalman Filter

        ! xprior = A * xposterior(-1)
        xprior      = 0.0d0
        call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)

        sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)

        ! fill directly into qrR
        qrR = 0.0d0
        ! qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular

        ! ytilde = y - C * xprior
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xprior,1,1.0d0,ytilde,1)

        ! singularity fix: insert unit dummies for missing values and set ytilde to zero
        do i=1,Ny
           if (yNaN(i,j)) sqrtSigmaY(i,i) = 1.0d0
           ! ytilde(i) = 0.0d0 -- not needed since y(i)=0 and C(i) * xprior also zero
        end do

        logdetSigmaY = 2.0d0 * sum( (/ (log(abs(sqrtSigmaY(i,i))), i = 1, Ny) /) )

        ! rotate ytilde (up to sign, consistent with rotation of K -- needed for llf computation)
        call dtrsv('U', 'T', 'N', Ny, sqrtSigmaY, Ny, ytilde, 1) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! compute log-likelihood
        ! llf
        llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

        if (atELB(j)) then
           shadowvolvec      = Cshadowrate
           call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,shadowvolvec,1)
           shadowrateVol(k)  =  sqrt(sum(shadowvolvec ** 2))
           shadowrateMean(k) = sum(Cshadowrate * xposterior(:,k))
           critvalues(k)     = (elbound - shadowrateMean(k)) / shadowrateVol(k) 
           ! critvalues will be converted into pdf after particle loop
        end if

     END DO ! k particles
     !$OMP END PARALLEL DO 

     ! compute ELBlike (if at ELB),
     ! note: it is more efficient to call vdcdf over vector of particles rather than individually for each particle across OMP loops
     if (atELB(j)) then
        call vdcdfnorm(Nparticles,critvalues,ELBlike)
     end if

     if (Nynonan > 0) then

        ! Reweight particles for next round

        loglikeAPFmax      = maxval(llf)
        llf                = llf - loglikeAPFmax
        if (atELB(j)) then
           APFlike = exp(llf) * ELBlike          
        else
           APFlike = exp(llf)
        end if
        APFkernelweights     = APFlike * PARTICLEweights(:,j-1)
        APFkernelsum         = sum(APFkernelweights)
        PARTICLEweights(:,j)      = APFkernelweights / APFkernelsum
	   

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(1,j))

        FORALL(k=1:Nparticles) shufflevec(k) = APFlike(ndx(k))
        APFlike = shufflevec

        DO i=1,Nx
           FORALL(k=1:Nparticles) shufflevec(k) = xposterior(i,ndx(k))
           xposterior(i,:) = shufflevec
        END DO

        DO i=1,Nsigmax
           FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaX(i,ndx(k))
           vecSqrtSigmaX(i,:) = shufflevec
        END DO

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = h(i,ndx(k))
           h(i,:) = shufflevec
        END DO

     end if

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: APF STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! RESAMPLE PARTICLES
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! 1) Draw Particles 
     errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hshock, 0.0d0, 1.0d0)

     ! call savemat(hshock, 'zshock.debug')

     call DTRMM('l','u','n','n',Nsv,Nparticles,1.0d0,sqrtVh,Nsv,hshock,Nsv) ! recall: sqrtVh is left-upper choleski factor

     forall (i=1:Nsv,k=1:Nparticles) h(i,k) = hintercept(i) + hrho(i) * h(i,k) + hshock(i,k)
     SVol        = exp(h * 0.5d0)


     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     !$OMP PARALLEL DO SHARED(DRAWytilde,xposterior, vecSqrtSigmaX, SVol, DRAWxgain, DRAWxgainstandardized, llf, Nparticles, Nsv, j, y, yNaN, atELB, Ny, Nx, Nw, Nynonan,A,B,C,Cshadowrate,elbound,critvalues,shadowrateMean, shadowrateVol) FIRSTPRIVATE(qrLWORK) PRIVATE(ytilde, logdetSigmaY, errcode, xprior, sqrtSigmaX, sqrtSigmaY, Kgain, qrR) PRIVATE(shadowvolvec) DEFAULT(NONE) FIRSTPRIVATE(Bsv) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        xprior      = 0.0d0

        ! 2) Fill Particles into state space
        ! Bsv
        FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * SVol(i,k) 

        ! 3) Kalman Filter

        ! xprior = A * xposterior(-1)
        call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)

        ! ------------------------------------------------------------------------
        ! SQRT KALMAN
        ! ------------------------------------------------------------------------
        sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)

        ! fill directly into qrR
        qrR = 0.0d0
        ! qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
        Kgain        = transpose(qrR(1:Ny,Ny+1:Ny+Nx))

        ! ytilde and logdetSigmaY
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xprior,1,1.0d0,ytilde,1)

        ! store ytilde before it gets rotated
        DRAWytilde(:,k,j) = ytilde

        ! singularity fix: insert unit dummies for missing values
        do i=1,Ny
           if (yNaN(i,j)) sqrtSigmaY(i,i) = 1.0d0
           ! ytilde(i) = 0.0d0 -- not needed since y(i)=0 and C(i) * xprior also zero
        end do

        logdetSigmaY = 2.0d0 * sum( (/ (log(abs(sqrtSigmaY(i,i))), i = 1, Ny) /) )

        ! rotate/normalize ytilde (up to sign, consistent with rotation of K)
        call dtrsv('U', 'T', 'N', Ny, sqrtSigmaY, Ny, ytilde, 1) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! xposterior = xprior + K * ytilde
        xposterior(:,k) = xprior
        call DGEMV('N',Nx,Ny,1.0d0,Kgain,Nx,ytilde,1,1.0d0,xposterior(:,k),1)


        ! remove unit dummies -- can be omitted since sqrtSigmaY not used any further
        ! do i=1,Ny
        !    if (yNaN(i,j)) sqrtSigmaY(i,i) = 0.0d0
        ! end do

        ! compute log-likelihood
        ! llf
        llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

        shadowrateMean(k) = sum(Cshadowrate * xposterior(:,k)) ! will be equal to actual rate away from ELB
        if (atELB(j)) then
           shadowvolvec      = Cshadowrate
           call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,shadowvolvec,1)
           shadowrateVol(k)  =  sqrt(sum(shadowvolvec ** 2))
           critvalues(k)     = (elbound - shadowrateMean(k)) / shadowrateVol(k) 
           ! critvalues will be converted into pdf after particle loop
        end if

        vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)


        ! first, store standardized Kgains (just adjusting for sign of innovations)

        FORALL(i=1:Nx,nn=1:Ny) DRAWxgainstandardized(i,nn,k,j) = Kgain(i,nn) 
        do nn=1,Ny
           if (sqrtSigmaY(nn,nn) < 0) then 
              forall (i=1:Nx) DRAWxgainstandardized(i,nn,k,j) = -DRAWxgainstandardized(i,nn,k,j) 
           end if
        end do

        ! second, store rotated Kalman gains, conforming with non-normalized ytilde
        call dtrsm('R', 'U', 'T', 'N', Nx, Ny, 1.0d0, sqrtSigmaY, Ny, Kgain, Nx) ! better to let dttrsm operate on the non-sliced Kgain rather than a slice of DRAWxgain
        FORALL(i=1:Nx,nn=1:Ny) DRAWxgain(i,nn,k,j)  = Kgain(i,nn) 


        ! ------------------------------------------------------------------------
        ! DONE: SQRT KALMAN
        ! ------------------------------------------------------------------------


     END DO ! k particles
     !$OMP END PARALLEL DO 


     ! compute ELBlike (if at ELB),
     ! note: brief testing indicates that it is more efficient to call vdcdf over vector of particles rather than individually for each particle across OMP loops
     if (atELB(j)) then
        ! shadowrateMeanAtELB = meanTruncNormalUB(Nparticles,shadowrateMean,shadowrateVol,ELBound)

        call meanvarianceTruncNormalUB(Nparticles,shadowrateMeanAtELB,shadowrateVolAtELB,shadowrateMean,shadowrateVol,ELBound)

        call vdcdfnorm(Nparticles,critvalues,ELBlike)
        llf = llf + log(ELBlike)

        DRAWshadowratehat(:,j) = shadowrateMeanAtELB

        ! update xposterior and sqrtSigmaX
        !$OMP PARALLEL DO SHARED(Nparticles,Nx,xposterior,vecSqrtSigmaX,Cshadowrate,shadowrateMeanAtELB,shadowrateVolAtELB,shadowrateMean,shadowrateVol) PRIVATE(TID,qrRelb,qrRlyap,shadowrategain,sqrtSigmaX) FIRSTPRIVATE(qrLelb,qrLlyap) DEFAULT(NONE) SCHEDULE(STATIC)
        DO k=1,Nparticles

           !$ TID = OMP_GET_THREAD_NUM()

           sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)
           ! call savemat(sqrtSigmaX, 'sqrtSigmaXprior.debug')

           ! set up QR -- TODO: move this up to the critvals? (second moments independent of realized signal)
           qrRelb                = 0.0d0
           qrRelb(2:Nx+1,2:Nx+1) = sqrtSigmaX
           qrRelb(2:Nx+1,1)      = Cshadowrate
           call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,qrRelb(2:Nx+1,1),1)
           ! QR decomposition
           call qrot(qrRelb, qrLelb)

           sqrtSigmaX         = qrRelb(2:1+Nx,2:1+Nx)
           vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)

           shadowrategain   = qrRelb(1,2:1+Nx) / qrRelb(1,1)

           ! update xposterior
           xposterior(:,k) = xposterior(:,k) + shadowrategain * (shadowrateMeanAtELB(k) - shadowrateMean(k))
           ! recall that shadowrateMean is what would be expected of shadowrate absent ELB; todo: replace previous line with call to DAXPY(N,DA,DX,INCX,DY,INCY) ?

           ! update sqrtSigmaX
           qrRlyap = 0.0d0
           qrRlyap(1:Nx,1:Nx) = sqrtSigmaX
           qrRlyap(1+Nx,1:Nx) = shadowrategain * shadowratevolAtELB(k)

           ! if (TID == 0) call savemat(qrRlyap, 'Rpre.debug')


           call qrot(qrRlyap, qrLlyap)


           sqrtSigmaX         = qrRlyap(1:Nx,1:Nx)
           vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)


        END DO
        !$OMP END PARALLEL DO 
     else
        DRAWshadowratehat(:,j) = shadowrateMean
     end if


     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! log-likelihood
     loglikemax        = maxval(llf)
     llf               = llf - loglikemax
     kernelweights     = exp(llf) / APFlike 
     kernelsum         = sum(kernelweights)
     loglike(j)        = log(kernelsum / Nparticles) + log(APFkernelsum) + loglikemax ! note: loglikeAPFmax correction should be added to APFlike and subtracted from kernelweights (b/o APFlike in denominator), thus cancelling each other


        ! particles weights
        PARTICLEweights(:,j) = kernelweights / kernelsum

     if (atELB(j)) then
        ELBtwist(:,j) = ELBlike
     end if

     ! Store statistics
     FORALL(i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,j)   = SVol(i,k) ! note the transpose

     FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,j)  = vecSqrtSigmaX(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxhat(i,k,j)             = xposterior(i,k)

     CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

  ! cumulate loglike contributions
  do j=T,1,-1
     loglike(j) = sum(loglike(1:j))
  end do
  

CONTAINS

  FUNCTION Cshortratepath(K, Nx, A, cshadow) result(c)

    integer, intent(in) :: Nx, K
    double precision, dimension(Nx,Nx), intent(in) :: A
    double precision, dimension(Nx), intent(in) :: cshadow
    double precision, dimension(Nx) :: c

    double precision, dimension(Nx,Nx) :: sumA, AA
    integer :: i

    sumA = A
    do i=2,K
       AA   = sumA
       sumA = A
       call DGEMM('n','n',Nx,Nx,Nx,1.0d0,A,Nx,AA,Nx,1.0d0,sumA,Nx)
    end do

    ! c = cshadow * sumA / K
    call DGEMV('t',Nx,Nx,1 / dble(K),sumA,Nx,cshadow,1,0.0d0,c,1)

  END FUNCTION Cshortratepath


END SUBROUTINE particlefilter

! @\newpage\subsection{particleIRF}@
SUBROUTINE particleIRF(doIRFunitscale, XresponseHat, IRFinnovol, Nparticles, Nx, Ny, Nimpulse, Nhorizons, xposterior, vecSqrtSigmaX, Nw, A, B, C, Nsv, SVol) !, VSLstream)

  ! Version 0.1:
  ! - without APF step, proposal SVs are assumed to be given
  ! - no ELB shocks, just shadowrate shocks
  ! - todo: sample SV proposals?
  ! = to be seen: in this version, "proposal weights" are equally weighted, coresponding to the case where innovations are standardized (such that the "innovation data" does not provide information about SV)
  ! todo: simplify code by computing IRF directly, instead of difference between XresponseHat and Xbaseline (since same weighting applies)

  USE embox, only : savemat, savevec, int2str

  ! use gibbsbox, only : drawNDXsysresample, igammaDraws
  ! use statespacebox, only : DLYAPsqrt
  use blaspack, only : logtwopi, eye, ivechU, qrot, qrquery

  use vslbox
  use omp_lib


  IMPLICIT NONE

  ! double precision, intent(in) :: ELBound

  logical, intent(in) :: doIRFunitscale 
  integer, intent(in) :: Nx, Ny, Nimpulse, Nparticles, Nhorizons, Nw, Nsv
  integer :: k, i, ii, nn, lag !, NsigmaX

  ! logical, intent(in) :: elb
  ! double precision, intent(inout), dimension(Ny,Nparticles) :: IRFweights
  double precision, intent(out), dimension(Nx,0:Nhorizons,Ny) :: XresponseHat
  double precision, intent(out), dimension(Ny) :: IRFinnovol
  ! parameter inputs
  double precision, intent(in), dimension(Nx,Nparticles) :: xposterior 
  double precision, intent(in), dimension(Nx * (Nx + 1) / 2,Nparticles) :: vecSqrtSigmaX
  double precision, intent(in) :: A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
  double precision, intent(in) :: SVol(Nsv,Nparticles)

  ! work
  double precision, dimension(Nx,0:Nhorizons,Ny,Nparticles) :: Xresponse
  double precision :: ytilde(Ny,Ny) !! matrix of normalized shocks
  double precision, parameter :: ytildeScale = 1.0d0

  double precision, dimension(Ny,Nparticles) :: innovol

  ! state space and Kalman Filtering
  double precision :: Bsv(Nx,Nw)
  double precision :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny), Kgain(Nx,Ny), xprior(Nx)
  ! double precision, dimension(Ny,Nparticles) :: llf ! likelihood per shocks per particle
  ! double precision :: total(Ny)
  double precision, parameter :: minParticleWeight = 1.0d-12

  ! qr
  double precision, dimension(Ny+Nx+Nw,Ny+Nx) ::  qrR
  integer :: qrLwork

  ! logical, dimension(Ny), intent(in) :: yNaN
  ! init stuff
  Xresponse = 0.0d0
  Bsv       = B
  call eye(ytilde,ytildeScale)
  ! todo: zero out ytilde when there are missing obs


  ! workspace query for qr decomposition
  qrR     = 0.0d0
  qrlwork = qrquery(qrR)


  ! ------------------------------------------------------------------------------------------------------------------------------
  ! BEGIN: MAIN IRF PARTICLE STEP
  ! ------------------------------------------------------------------------------------------------------------------------------

  !$OMP PARALLEL DO SHARED(Xresponse,xposterior,vecSqrtSigmaX, SVol, Nparticles, Nsv, Ny, Nimpulse, Nx, Nw, Nhorizons, A,B,C,doIRFunitscale) SHARED(innovol) FIRSTPRIVATE(qrLWORK) PRIVATE(sqrtSigmaX, sqrtSigmaY, Kgain, qrR) PRIVATE(i,xprior) DEFAULT(NONE) FIRSTPRIVATE(Bsv) SCHEDULE(STATIC) ! SHARE(llf)


  DO k = 1,Nparticles

     ! 2) Fill Particles into state space
     ! Bsv
     FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * SVol(i,k) 

     ! ------------------------------------------------------------------------
     ! SQRT KALMAN
     ! ------------------------------------------------------------------------
     sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)

     ! fill directly into qrR
     qrR = 0.0d0
     qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
     ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
     call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
     ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
     call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C,Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

     ! QR decomposition
     call qrot(qrR, qrLWORK)

     ! map qr into Kalman objects
     sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
     sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
     Kgain        = transpose(qrR(1:Ny,Ny+1:Ny+Nx)) ! Gain on normalized shocks

     ! ------------------------------------------------------------------------
     ! DONE: SQRT KALMAN
     ! ------------------------------------------------------------------------

     ! sign / scale shocks: figure out diagonal of sqrtSigmaY = C * K
     do nn=1,Nimpulse
        if (doIRFunitscale) then
           ! scale to 1SV unit
           forall (i=1:Nx) Kgain(i,nn) = Kgain(i,nn) / sqrtSigmaY(nn,nn)
        else 
           ! sign only
           if (sqrtSigmaY(nn,nn) < 0) then
              forall (i=1:Nx) Kgain(i,nn) = -Kgain(i,nn) 
           end if
        end if
     end do

     forall(nn=1:Nimpulse) innovol(nn,k) = abs(sqrtSigmaY(nn,nn))

     ! ------------------------------------------------------------------------
     ! BEGIN: XRESPONSES
     ! ------------------------------------------------------------------------

     ! compute new expectations
     lag = 0
     ! a): form xprior (same for al impulses)
     ! xprior = A * xposterior
     call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)
     ! b) update prior to impact after unit inovation by adding corresponding column of (standardized) Kgain
     forall (nn=1:Nimpulse) Xresponse(:,lag,nn,k) = xprior + Kgain(:,nn)
     ! c) iterate over expectations at higher lags
     do nn=1,Nimpulse
        do lag=1,Nhorizons
           call DGEMV('n',Nx,Nx,1.0d0,A,Nx,Xresponse(:,lag-1,nn,k),1,0.0d0,Xresponse(:,lag,nn,k),1)
        end do
     end do

     ! ------------------------------------------------------------------------
     ! DONE: XRESPONSES
     ! ------------------------------------------------------------------------


  END DO ! k particles
  !$OMP END PARALLEL DO 

  ! ------------------------------------------------------------------------------------------------------------------------------
  ! END: MAIN IRF PARTICLE STEP
  ! ------------------------------------------------------------------------------------------------------------------------------


  !$OMP PARALLEL DO SHARED(Nx,Nhorizons,Nimpulse,XresponseHat,Xresponse,Ny,Nparticles) PRIVATE(ii,lag) DEFAULT(NONE) SCHEDULE(STATIC) 
  ! ingegrate Xresponse across particles
  do nn =1,Nimpulse
     forall (ii=1:Nx,lag=0:Nhorizons) XresponseHat(ii,lag,nn) = sum(Xresponse(ii,lag,nn,:)) / dble(Nparticles)
  end do
  !$OMP END PARALLEL DO 

  IRFinnovol = sum(innovol,2) / dble(Nparticles)

END SUBROUTINE particleIRF


! @\newpage\subsection{particleIRFsim}@
SUBROUTINE particleIRFsim(doIRFunitscale, Nsim, ysimBaseline, yirf, yNaN, ytrunc, elbound, Nparticles, Nx, Ny, Nimpulse, Nhorizons, Nysim, xposterior, vecSqrtSigmaX, Nw, A, B, Cobserver, Csim, Nsv, hprev, hbar, hrho, sqrtVh, VSLstream, NTHREADS, VSLarray)

  ! computes simulated IRF and tracks truncated variables
  USE embox, only : savemat, savevec, int2str
  use blaspack, only : vechU, ivechU, qrot, qrquery
  use timerbox, only : tic, toc

  use vslbox
  use omp_lib


  IMPLICIT NONE

  logical, intent(in) :: doIRFunitscale
  integer, intent(in) :: Nsim 
  double precision :: timermark

  integer, intent(in) :: Nparticles, Nx, Ny, Nimpulse, Nhorizons, Nysim, Nw, Nsv
  logical, dimension(Nysim) :: ytrunc
  logical, dimension(Ny) :: yNaN
  ! logical, intent(in) :: jumpoffAtELB
  double precision, intent(in) :: elbound
  double precision, intent(in) :: xposterior(Nx,Nparticles), vecSqrtSigmaX(Nx * (Nx + 1) / 2, Nparticles)
  ! integer :: NsigmaX
  double precision, intent(in) :: A(Nx,Nx), B(Nx,Nw), Cobserver(Ny,Nx), Csim(Nysim,Nx) ! , Cshadowrate(Nx)
  double precision, intent(in) :: hprev(Nsv,Nparticles), hbar(Nsv), hrho(Nsv), sqrtVh(Nsv,Nsv)
  double precision :: hintercept(Nsv)
  double precision :: Ctilde(Ny,Nx)

  ! double precision, intent(inout), dimension(Nimpulse,Nparticles)  :: IRFweights
  ! NOTE: in response to unit-SV shocks, IRF weights are all equal

  double precision, intent(out), dimension(Nysim,0:Nhorizons,2,Nimpulse)  :: yirf

  ! OMP
  integer, intent(in) :: NTHREADS
  integer :: TID = 0
  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  type (vsl_stream_state), intent(inout) :: VSLstream
  type (vsl_stream_state), intent(inout), dimension(0:NTHREADS-1) :: VSLarray


  integer :: k, ii, jj, nn, mm

  double precision, dimension(Nsv,Nparticles) :: h0
  double precision, dimension(Nsv,Nsim,Nhorizons) :: h

  double precision, dimension(Nw,Nparticles) :: SVol0
  double precision, dimension(Nw,Nsim,Nhorizons) :: SVol
  ! Note: SV is stored for all shocks (set to one in case of fixed variance shocks)

  double precision, dimension(Nw,Nsim,0:Nhorizons) :: wshock ! note: 0 draws only used for baseline simulation
  double precision, dimension(Nx,Nsim) :: x0shock  ! note: using the same shocks for each impulse's jmp-off. todo: better to draw separate ones?
  double precision, dimension(Nx,Nsim) :: xsim, xsimlag
  double precision, dimension(Nysim,Nsim,0:Nhorizons) :: ysim
  double precision, dimension(Nysim,0:Nhorizons,Nparticles) :: ysimBaselineParticles
  double precision, dimension(Nysim,0:Nhorizons), intent(out) :: ysimBaseline
  double precision, dimension(Nysim,0:Nhorizons,2,Nimpulse,Nparticles) :: ysimAlternatives
  integer :: plusminus

  ! work
  double precision, dimension(Nx,2,Nimpulse) :: XpostImpact
  double precision, dimension(Nx * (Nx + 1) / 2, Nimpulse) :: vecSqrtSigmaXpostImpact

  double precision, dimension(Nx,Nparticles) :: Xjumpoff ! , xdraws
  ! double precision, dimension(Nparticles) :: shadowrate
  ! logical :: OK

  ! state space and Kalman Filtering
  double precision :: Bsv(Nx,Nw)
  double precision :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny), Kgain(Nx,Ny), xprior(Nx)
  double precision, parameter :: minParticleWeight = 1.0d-12

  ! prepare 3 qr objects: 
  double precision ::  qrRprior(Nx + Nw, Nx), qrRpostImpact(1+Nx,1+Nx) ! , qrRobserver(Nx,Ny)
  integer :: qrLprior, qrLpostImpact !, qrLobserver

  DOUBLE PRECISION :: qrR(Ny+Nx+Nw,Ny+Nx)
  INTEGER :: qrL

  ! init stuff
  ! XpostImpact = 0.0d0
  Bsv         = B

  ! workspace query for qr decomposition
  qrRprior       = 0.0d0
  qrLprior       = qrquery(qrRprior)
  qrRpostImpact  = 0.0d0
  qrLpostImpact  = qrquery(qrRpostImpact)
  ! qrRobserver    = 0.0d0
  ! qrLobserver    = qrquery(qrRobserver)
  qrR            = 0.0d0
  qrL            = qrquery(qrR)

  ! qrR    = 0.0d0
  ! qrL    = qrquery(qrR)

  ! ------------------------------------------------------------------------
  ! BEGIN: Draw SVol0 and SVol
  ! ------------------------------------------------------------------------

  timermark = tic() 

  ! a: prepare hintercept
  hintercept = (1 - hrho) * hbar

  ! b: time0 -- draw shocks, rotate and propagate AR1 processes
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h0, 0.0d0, 1.0d0)
  call DTRMM('l','u','n','n',Nsv,Nparticles,1.0d0,sqrtVh,Nsv,h0,Nsv) ! recall: sqrtVh is left-upper choleski factor
  forall (ii=1:Nsv,k=1:Nparticles) h0(ii,k) = hintercept(ii) + hrho(ii) * hprev(ii,k) + h0(ii,k)
  SVol0 = 1.0d0
  SVol0(1:Nsv,:) = exp(h0 * 0.5d0)


  ! ------------------------------------------------------------------------
  ! Draw Xjumpoff (for baseline)
  ! ------------------------------------------------------------------------
  ! note:
  ! - need to draw Xjumpoff while observing whether data was at ELB
  ! - ELB already ensured for Xjumpoff due to particle sampler
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Nparticles, xjumpoff, 0.0d0, 1.0d0)

  ! ------------------------------------------------------------------------
  ! BIG K LOOP ACROSS PARTICLES
  ! ------------------------------------------------------------------------
  ! ---- PRIVATE(sqrtSigmaXtilde, sqrtSigmaYtilde, Ktilde)
  !$OMP PARALLEL DO SHARED(VSLarray,Nparticles,Nsv,Nhorizons,Nsim,Nw,Nx,Ny,Nysim,hintercept,hrho,sqrtVh, vecSqrtSigmaX,xjumpoff,xposterior,A,B,Cobserver,Csim,ysimBaselineParticles,yNaN,ytrunc,elbound) SHARED(Svol0, h0) SHARED(ysimAlternatives) SHARED(Nimpulse,doIRFunitscale) PRIVATE(x0shock,SVol,h,wshock,XpostImpact,vecSqrtSigmaXpostimpact) PRIVATE(TID,errcode,ii,jj,nn,ysim,xsim,xsimlag) FIRSTPRIVATE(qrLprior,qrLpostImpact,qrL, Bsv) PRIVATE(qrRprior,qrRpostImpact, qrR, sqrtSigmaX, sqrtSigmaY, Kgain, Ctilde,xprior)  DEFAULT(NONE) SCHEDULE(STATIC)
  do k=1,Nparticles

     TID = 0
     !$ TID = OMP_GET_THREAD_NUM()

     ! print *, 'k', k, 'TID', TID


     ! ------------------------------------------------------------------------
     ! Draw  x0shock
     ! ------------------------------------------------------------------------
     errcode = vdrnggaussian(VSLmethodGaussian, VSLarray(TID), Nx * Nsim, x0shock, 0.0d0, 1.0d0)

     errcode = vdrnggaussian(VSLmethodGaussian, VSLarray(TID), Nsv * Nsim * Nhorizons, h, 0.0d0, 1.0d0)
     call DTRMM('l','u','n','n',Nsv,Nhorizons * Nsim, 1.0d0,sqrtVh,Nsv,h,Nsv)
     ! note: lapack does not care about exact shape of h (as if it performed implicit reshape)

     SVol = 1.0d0
     nn = 1 
     forall (ii=1:Nsv,jj=1:Nsim) h(ii,jj,nn) = hintercept(ii) + hrho(ii) * h0(ii,k) + h(ii,jj,nn)
     forall (ii=1:Nsv,nn=2:Nhorizons,jj=1:Nsim) h(ii,jj,nn) = hintercept(ii) + hrho(ii) * h(ii,jj,nn-1) + h(ii,jj,nn)
     SVol(1:Nsv,:,:) = exp(h * 0.5d0)
     ! ------------------------------------------------------------------------
     ! END: Draw SVol0 and SVol
     ! ------------------------------------------------------------------------


     ! ------------------------------------------------------------------------
     ! Draw and scale wshock
     ! ------------------------------------------------------------------------

     errcode = vdrnggaussian(VSLmethodGaussian, VSLarray(TID), Nw * Nsim * (1 + Nhorizons), wshock, 0.0d0, 1.0d0)

     ! time 0 (special since thre is only one SVol0 per particle)
     nn = 0
     forall (ii=1:Nw,jj=1:Nsim) wshock(ii,jj,nn) = SVol0(ii,k) * wshock(ii,jj,nn)
     ! future periods
     forall (ii=1:Nw,jj=1:Nsim,nn=1:Nhorizons) wshock(ii,jj,nn) = SVol(ii,jj,nn) * wshock(ii,jj,nn)

     ! print *, 'drew shocks', toc(timermark)

     ! ------------------------------------------------------------------------
     ! scale xjumpoff
     ! ------------------------------------------------------------------------
     ! collect posterior variance
     ! (exploit packed storage)
     call DTPMV('u','t','n',Nx,vecSqrtSigmaX(:,k),xjumpoff(:,k),1)
     xJumpoff(:,k) = xposterior(:,k) + xjumpoff(:,k)

     ! ------------------------------------------------------------------------
     ! SIMULATE FORWARD:
     ! ------------------------------------------------------------------------
     ! - Nsimdraws of x0 based on SVol0 (per particle)
     ! - then for the remaining horizons


     ! first copy Xjumpoff into xsimlag
     forall (ii=1:Nx,jj=1:Nsim) xsimlag(ii,jj) = xjumpoff(ii,k)
     ! note: in principle, A * xsimlag would have to be computed only once at nn=1; for the sake of more compact code I ignore this however


     ! simulate over horizons
     do nn=0,Nhorizons
        ! xsim = A * xsimlag
        call DGEMM('n','n',Nx,Nsim,Nx,1.0d0,A,Nx,xsimlag,Nx,0.0d0, xsim, Nx)
        ! xsim = xsim + B * wshock
        call DGEMM('n','n',Nx,Nsim,Nw,1.0d0,B,Nx,wshock(:,:,nn),Nw,1.0d0,xsim,Nx)
        ! ysim = Csim * xsim
        call DGEMM('n','n',Nysim,Nsim,Nx,1.0d0,Csim,Nysim,xsim,Nx,0.0d0,ysim(:,:,nn),Nysim)

        ! prep next horizons
        xsimlag = xsim

     end do ! nn

     ! implement truncation
     forall (jj=1:Nsim,nn=0:Nhorizons)
        where (ytrunc .AND. ysim(:,jj,nn) < elbound) ysim(:,jj,nn) = elbound
     end forall

     ! integrate over sims to compute baseline
     ysimBaselineParticles(:,:,k) = sum(ysim,2) / dble(Nsim) 

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN IRF PARTICLE STEP: delivers posterior mean and variance for XpostImpact
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! collect posterior from t=-1
     sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)
     ! call savemat(sqrtSigmaX, 'sqrtSigmaX0.debug') 

     ! ------------------------------------------------------------------------
     ! Fill Particles into state space
     ! ------------------------------------------------------------------------
     ! Bsv
     FORALL (ii=1:Nsv) Bsv(:,ii)  = B(:,ii) * Svol0(ii,k) 



     ! ------------------------------------------------------------------------
     ! basic/original QR
     ! ------------------------------------------------------------------------
     qrR = 0.0d0
     qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
     ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
     call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
     ! todo: swap previous DGEMM for DTRMM
     ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
     call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,Cobserver,Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

     ! QR decomposition
     call qrot(qrR, qrL)

     ! map qr into Kalman objects
     sqrtSigmaY        = qrR(1:Ny,1:Ny) ! upper triangular
     ! sqrtSigmaXtilde   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
     Kgain            = transpose(qrR(1:Ny,Ny+1:Ny+Nx))
     ! ------------------------------------------------------------------------


     ! impose sign on Kgain
     do nn=1,Ny
        if (doIRFunitscale) then
           ! scale to 1SV unit
           forall (ii=1:Nx) Kgain(ii,nn) = Kgain(ii,nn) / sum(Cobserver(nn,:) * Kgain(:,nn))
        else 
           ! sign only
           if (sum(Cobserver(nn,:) * Kgain(:,nn))  < 0) then
              forall (ii=1:Nx) Kgain(ii,nn) = -Kgain(ii,nn) 
           end if
        end if
     end do

     ! ------------------------------------------------------------------------

     ! ------------------------------------------------------------------------
     ! PRIOR VARIANCE on impact (t=0)
     ! [ sqrtSigmaX' , 0 ] Q = [A * sqrtSigmaX' , Bsv ]
     ! ------------------------------------------------------------------------
     qrRprior = 0.0d0
     qrRprior(Nx+1:Nx+Nw,1:Nx) = transpose(Bsv)

     ! qrR(1:Nx,1:Nx) = sqrtSigmaX * A'
     qrRprior(1:Nx,1:Nx) = transpose(A)
     call DTRMM('l','u','n','n',Nx,Nx,1.0d0,sqrtSigmaX,Nx,qrRprior(1:Nx,1:Nx),Nx)

     call qrot(qrRprior, qrLprior)
     sqrtSigmaX = qrRprior(1:Nx,1:Nx) ! upper triangular

     ! ------------------------------------------------------------------------
     ! VARIANCE of observer 
     ! [ sqrtSigmaY' , 0 ] Q = [ C * sqrtSigmaX' ]
     ! ------------------------------------------------------------------------
     ! qrRobserver = 0.0d0

     ! qrRobserver = sqrtSigmaX * C'
     ! qrRobserver = transpose(Cobserver)
     ! call DTRMM('l','u','n','n',Nx,Ny,1.0d0,sqrtSigmaX,Nx,qrRobserver,Nx)

     ! ! call savemat(qrRobserver, 'M2.debug')
     ! ! call DGEMM('n','t',Nx,Ny,Nx,1.0d0,sqrtSigmaX,Nx,Cobserver,Ny,0.0d0,qrRobserver,Nx)
     ! ! call savemat(qrRobserver, 'M.debug')
     ! ! call savemat(sqrtSigmaX, 'sqrtSigmaX.debug')
     ! ! call savemat(Cobserver, 'C.debug')

     ! call qrot(qrRobserver, qrLobserver)
     ! sqrtSigmaY  = qrRobserver(1:Ny,1:Ny) ! upper triangular

     ! to handle NaN: set diag(sqrtSigmaY) to unity where yNaN
     do nn=1,Ny
        if (yNaN(nn)) sqrtSigmaY(nn,nn) = 1.0d0
     end do

     ! ------------------------------------------------------------------------
     ! Ctilde = sqrtSigmaY' \ C
     ! ------------------------------------------------------------------------
     Ctilde = Cobserver
     call DTRSM('l','u','t','n',Ny,Nx,1.0d0,sqrtSigmaY,Ny,Ctilde,Ny)

     ! ------------------------------------------------------------------------
     ! Iterate over separate Kalman updates for each Column of Ctilde
     ! - note: these are *separate* not sequential updates
     ! - end product: fill Kgain and store posterior VCV
     ! ------------------------------------------------------------------------
     do nn = 1,Nimpulse ! loop over impulses
        ! set up qr update of the Kalman filter
        ! note:
        ! - (Ctilde * sqrtSigmaPrior) (Ctilde * sqrtSigmaPrior)' = I
        ! - Ktilde = SigmaPrior * Ctilde'
        ! - SigmaPosterior = SigmaPrior - Ktilde * Ktilde'
        ! thus [ I, Ktilde'; Ktilde, SigmaPrior ] has a factorization equal to [ I, 0; Ktilde, sqrtSigmaPosterior ]
        ! - sqrtSigmaPosterior is typically singular, the factorization above is thus not necessarily a choleski 
        ! the QR then operates on
        ! [ 0 , Ctilde * sqrtSigmaPrior; 0, sqrtSigmaPrior ] Q =  [ I, 0; Ktilde, sqrtSigmaPosterior ]
        ! we do the above separately for each row of Ctilde

        qrRpostimpact           = 0.0d0
        qrRpostimpact(2:1+Nx,1) = Ctilde(nn,:)   
        call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,qrRpostimpact(2:1+Nx,1),1)
        qrRpostimpact(2:1+Nx,2:1+Nx) = sqrtSigmaX

        ! call savemat(qrRpostImpact, 'A.debug')

        ! perform QR        
        call qrot(qrRpostImpact, qrLpostimpact)

        ! read out qrR
        ! sqrtSigmaXpostimpact(:,:,nn) = qrRpostImpact(2:1+Nx,2:1+Nx)
        vecSqrtSigmaXpostimpact(:,nn) = vechU(qrRpostImpact(2:1+Nx,2:1+Nx),Nx)

        ! Kgain(:,nn) = qrRpostImpact(1,2:1+Nx)

        ! ! sign Kgain such that C * Kgain = is positive (note: that is C not Ctilde since we want to sign the effect on the non-orthogonalized variables!)
        ! if (sum(Cobserver(nn,:) * Kgain(:,nn))  < 0) then
        !    forall (ii=1:Nx) Kgain(ii,nn) = -Kgain(ii,nn) 
        ! end if

     end do

     ! ------------------------------------------------------------------------
     ! DONE: SQRT KALMAN
     ! ------------------------------------------------------------------------


     ! compute posterior for states on impact after each impulse
     ! xprior = A * xposterior
     call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)
     ! update prior to impact after unit inovation by adding corresponding column of (standardized) Kgain

     forall (nn=1:Nimpulse) XpostImpact(:,1,nn) = xprior + Kgain(:,nn)
     forall (nn=1:Nimpulse) XpostImpact(:,2,nn) = xprior - Kgain(:,nn)

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN IRF PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! ------------------------------------------------------------------------
     ! BEGIN: simulate new trajectories
     ! ------------------------------------------------------------------------

     do ii=1,Nimpulse ! loop over impulses
        do plusminus = 1,2 
           ! IMPACT
           nn = 0
           ! draw from impact posterior
           ! xsimlag = sqrtSigmaXpostimpact' * x0shock

           xsimlag = x0shock
           do jj=1,Nsim
              call dtpmv('u','t','n',Nx,vecSqrtSigmaXpostimpact(:,ii),xsimlag(:,jj),1)
           end do ! sim jj
           forall(jj=1:Nsim,mm=1:Nx) xsimlag(mm,jj) = XpostImpact(mm,plusminus,ii) + xsimlag(mm,jj)


           ! ysim = Csim * xsim
           call DGEMM('n','n',Nysim,Nsim,Nx,1.0d0,Csim,Nysim,xsimlag,Nx,0.0d0,ysim(:,:,nn),Nysim)

           ! SUBSEQUENT HORIZONS
           do nn=1,Nhorizons
              ! xsim = A * xsimlag
              call DGEMM('n','n',Nx,Nsim,Nx,1.0d0,A,Nx,xsimlag,Nx,0.0d0, xsim, Nx)
              ! xsim = xsim + B * wshock
              call DGEMM('n','n',Nx,Nsim,Nw,1.0d0,B,Nx,wshock(:,:,nn),Nw,1.0d0,xsim,Nx)
              ! ysim = Csim * xsim
              call DGEMM('n','n',Nysim,Nsim,Nx,1.0d0,Csim,Nysim,xsim,Nx,0.0d0,ysim(:,:,nn),Nysim)
              ! prep next horizons
              xsimlag = xsim
           end do ! horizon nn

           ! implement truncation
           forall (jj=1:Nsim,nn=0:Nhorizons)
              where (ytrunc .AND. ysim(:,jj,nn) < elbound) ysim(:,jj,nn) = elbound
           end forall

           ! integrate over ysim
           ysimAlternatives(:,:,plusminus,ii,k) = sum(ysim,2) / dble(Nsim)

        end do ! plusminus
     end do ! impulse ii

     ! ------------------------------------------------------------------------
     ! DONE: simulate new trajectories
     ! ------------------------------------------------------------------------

  end do ! BIG K
  !$OMP END PARALLEL DO 
  ! ------------------------------------------------------------------------
  ! END BIG K LOOP
  ! ------------------------------------------------------------------------

  ! COMPUTE IRF
  ysimBaseline =  sum(ysimBaselineParticles,3) / dble(Nparticles)

  !$OMP PARALLEL DO SHARED(yirf,ysimAlternatives,ysimBaseline,Nimpulse,Nysim,Nhorizons,Nparticles) PRIVATE(plusminus,nn,jj) DEFAULT(NONE) SCHEDULE(STATIC) 
  do ii=1,Nimpulse
     do plusminus = 1,2
        forall(nn=0:Nhorizons,jj=1:Nysim) yirf(jj,nn,plusminus,ii) =  sum(ysimAlternatives(jj,nn,plusminus,ii,:)) / dble(Nparticles) - ysimBaseline(jj,nn)
     end do
  end do
  !$OMP END PARALLEL DO

  ! print *, 'done IRFsimf', toc(timermark)

END SUBROUTINE particleIRFsim
