PROGRAM main

  ! JAN 2017
  ! VAR model for rollingstone -- observable ugap
  ! QUARTERLY MODEL
  ! handles also missing data
  ! shockslopes
  ! NEW: predictive density, including 4-quarter inflation

  ! added a flag to switch on/off pdf computations

  ! allowing for correlated SV shocks

  ! SETUP:
  ! Nsv SV shocks, but only Ngap of them are correlated (in levels and variances)
  ! i.e. Nsv: hbar, hrho, SVol, but only NGap hSigma and Nsv-Ngap hvar
  ! This version track two SVol blocks: "h" are correlated, "horth" are orthogonal

  USE embox, only: hrulefill, savemat, savevec, storeestimatesOMP, storeestimates, lofT, int2str, timestampstr, es30d16, logtwopi
  USE blaspack, only: eye ! , sandwich ! , sqrtvcv ! , sqrtvcvTR
  USE gibbsbox, only: GelmanTest1, simpriormaxroot, igammadraw
  use densitybox, only : crps, logGaussianScore, logGaussianScoreMV, gewekelogMDD

  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  LOGICAL, PARAMETER :: doTimestamp = .false.

  INTEGER :: Nyield, Ny,  Nbar, Ngap, Nsvorth, Nsv, Nshock,  Nshockslopes

  INTEGER :: gapOrderKey = 1234
  ! gapOrderKey is a 4-digit code that specifies the order of variables within the gap block:
  ! - the first digit (from the left) specifies the position of inflation
  ! - the 2nd digit is for ugap
  ! - 3rd digit is the policy rate gap
  ! - 4th digit is the yield block
  ! - the implementation is handled via a select case statement in thissampler, missing cases can simply be added there

  INTEGER :: p
  INTEGER :: NhSigma 
  INTEGER, PARAMETER :: longrateHorizon = 4 * 5 ! in quarters
  LOGICAL :: ZLBrejection = .true., doPDF = .false., doMDD = .true., doPInoise = .false.
  LOGICAL :: doZLBasZeroes = .false.
  INTEGER, PARAMETER :: ndxSVpinoise = 2 ! hard coded, just needed to whack out Gelman stats and if doPInoise = .false.

  ! forecasting parameters
  integer :: NNyy 
  INTEGER, PARAMETER :: Nforecastdraws = 50
  integer :: NNxx

  ! horizons are consecutive sequence up to maxhorizon
  integer, parameter :: maxhorizons = 8
  integer, dimension(maxhorizons)  :: horizons

  ! forecast obs
  ! double precision, dimension(NNyy,maxhorizons) :: ypred
  ! logical, dimension(NNyy,maxhorizons) :: yNaNpred
  double precision, allocatable, dimension(:,:) :: ypred
  logical, allocatable, dimension(:,:) :: yNaNpred

  ! forecast results
  ! double precision, allocatable, dimension(:,:,:,:) :: DRAWypdf, DRAWyforecast, DRAWymedian, DRAWyprobLB, DRAWycondvar
  double precision, allocatable, dimension(:,:,:,:,:) :: DRAWy

  INTEGER :: Tdata
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: yhatmean, yhatmedian ! , yavgcondvar, yvarcondexp, ycondvar
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: crpsScore, logScore
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: logScoreMV

  INTEGER :: Nx, Nf, Nstates

  INTEGER :: Nsim, Burnin, Nstreams
  INTEGER :: Nparticles = 10
  INTEGER :: T,h,i,j,k,n,status,Ndraws
  INTEGER :: T0 = 0 ! jump off when reading in data
  ! LOGICAL :: OK
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWstates, DRAWsvol, DRAWsvolorth
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWavgtermpremia
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWhsigma, DRAWhbar, DRAWhrho 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWhorthvar, DRAWhorthbar, DRAWhorthrho
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWf, DRAWmaxlambda, DRAWshockslopes

  ! MDD
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRAWMDDlogprior, DRAWMDDloglike
  INTEGER :: Ntheta, ndxTheta
  DOUBLE PRECISION :: gewekeTau(4) = (/ 0.25d0, 0.5d0, 0.75d0, 0.9d0 /), gewekeTruncLogProb
  INTEGER :: gewekeTruncationNsim = 10 ** 5
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ndxf, ndxrho
  DOUBLE PRECISION :: fPriorTruncationProbLog, hrhoPriorTruncationProbLog ! adjustment factor to trunctated MV prior for f

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Ef0, Eshockslopes0
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Eh0, hrho0, minSV, maxSV, SVol0
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: hSigmaT, sqrtVh0, hrhoV0
  INTEGER :: hSigmaDof

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Ehorth0, Vhorth0, horthrho0, horthrhoV0, minSVorth, maxSVorth, horthvarT, SVol0orth
  INTEGER, ALLOCATABLE, DIMENSION(:) :: horthvarDof

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sqrtVf0, sqrtVshockslopes0
  DOUBLE PRECISION :: lambda1Vf, lambda2Vf


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GELMANstates, GELMANsvol, GELMANsvolorth
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta,theta2
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: theta1

  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: thetaDraws
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GELMANf, GELMANshockslopes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GELMANhorthvar, GELMANhorthbar, GELMANhorthrho
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GELMANhbar, GELMANhSigma, GELMANhrho

  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: thetaprob

  INTEGER :: ndxRealrate, ndxPIbar, ndxRbar, ndxLONGRbar, ndxPIgap, ndxPInoise, ndxUgap, ndxINTgap, ndxLONGINTgap, ndxSHADOWRATE !, ndxLONGSHADOW, ndxLONGRATE, ndxZLBPROB1, ndxZLBPROB2, ndxZLBPROB3

  INTEGER :: yndxPolicyRate = 3

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN, zlb

  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, filext, datalabel

  ! VSL Random Stuff
  ! type (vsl_stream_state), allocatable, dimension(:) :: VSLstreams
  type (vsl_stream_state) :: VSLdefaultstream, VSLstream
  integer :: seed
  ! integer :: brng
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS, TID, runnum = 1

  ! Lower Bound
  DOUBLE PRECISION, PARAMETER :: ELBound = 0.25d0

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  Nyield = 3  ! just a default, can be overriden by command line arguments
  p      = 2  ! just a default, can be overriden by command line arguments

  Nsvorth = 3 ! this could actually have been defined as a parameter 


  ! TOC

  ! Nsim       = 10 ** 4
  ! Burnin     = 10 ** 4
  ! Nparticles = 10 ** 4

  ! Nsim       = 10 ** 3
  ! Burnin     = 10 ** 4
  ! Nparticles = 10 ** 5


  Nsim    = (10 ** 3)           
  Burnin  = (10 ** 5)
  Nparticles = 10 ** 4

  runnum = 1
  ! ! QUICK
  ! Nsim    = (10 ** 2)
  ! Burnin  = (10 ** 3)
  ! Nparticles = 10

  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  Nstreams = max(NTHREADS, 4)
  print *, "Number of Threads:", NTHREADS
  print *, "Number of Streams:", Nstreams


  ! read data
  datalabel   = 'spectreTB3MSGS020510Headline2018Q4';

  ! T0 = 120 ! T0 = 120 skips data prior to 1990Q1
  T = Tdata - T0! default
  call getarguments(ZLBrejection, doPDF, doMDD, doPInoise, doZLBasZeroes, T, T0, gapOrderKey, datalabel, Nyield, p, Nparticles, Nsim, burnin, Nstreams, runnum)
  
  ! VSL
  call hrulefill
  write (*,*) 'Allocating VSLstreams ...'

  ! NOTE about initializing vslnewstreams:
  ! This code uses NStreams + 1 random streams, all from the mt2203 family
  ! For each value of runnum, a different set of (Nstreams + 1) mt2203 streams is used
  ! The default stream is set at runnum * (Nstreams + 1), and used outside the parallel MCMC loops
  ! The other Nstreams streams are used for the different MCMC runs, 
  ! ... they are initialized at (runnum - 1) * (Nstreams + 1) + j
  seed    = 0
  errcode = vslnewstream(VSLdefaultstream, vsl_brng_mt2203 + runnum * (Nstreams + 1), seed)  
  WRITE(*,*) "... default VSLstream ", VSLdefaultstream
  write (*,*) '... VSLstreams done.'
  call hrulefill

  ! define forecast horizons
  forall (k=1:maxhorizons) horizons(k) = k

  if (doMDD .AND. (gapOrderKey /= 1234)) then
     WRITE(*,'(a78, i5)') 'MDD calculation currently supported only with gapOrderKey 1234, but you chose ', gapOrderKey
     stop 1
  end if

  ! set parameters that depend on Nyield
  Ny = 3 + Nyield
  Nbar = 1 + 1 + Nyield
  Ngap = Ny
  Nsv = Ngap
  Nshock = 1 + Nsv + Nsvorth
  Nshockslopes = Nsv * (Nsv - 1) / 2
  NhSigma = Nsv * (Nsv + 1) / 2
  NNyy = Ny + 2 ! adding 4q-MA of inflation and the shadowrate
  Nx        = Nbar + (Ngap * p) + 1                      ! Nbar and Ngaps (with lags), pinoise
  Nf        = Ngap * (Ngap * p)


  ! runtime parameters :end: 

  Nstates   = 8 + 2 * Nyield ! see list of ndxXXX below

  ndxRealrate     = 1
  ndxPIbar        = 2
  ndxRbar         = 3
  ndxLONGRbar     = 4
  ndxPIgap        = ndxLONGRbar  + Nyield
  ndxUgap         = ndxPIgap + 1
  ndxINTgap       = ndxUgap + 1
  ndxLONGINTgap   = ndxINTgap + 1
  ndxPInoise      = ndxLONGINTgap + Nyield
  ndxSHADOWRATE   = ndxPInoise + 1


  if (ndxSHADOWrate .ne. Nstates) then
     print *, 'state indices are off'
     stop 1
  end if

  if (doZLBasZeroes) then
     ! doShadowObsGap = .false.
     ZLBrejection   = .false.
  end if

  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  print *, datalabel
  print *, datafile
  Tdata = loft(datafile) 
  IF (Tdata < T) THEN
     print *, 'Data file has less than T obs', datafile, T, Tdata
     STOP 1
  END IF


  ! read in data 
  ALLOCATE (yNaN(Ny,T), zlb(Ny,T), y(Ny,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Y)'
  END IF

  ! print *, 'trying to read', T, 'obs from', datafile
  CALL readdata(y,datafile,Ny,T,T0)
  CALL readnandata(yNaN,nandatafile,Ny,T,T0)

  if (doZLBasZeroes .eqv. .true.) then
     yNaN(yndxPolicyRate,:) = .false.
  end if

  ! call savemat(y, 'y.debug')
  ! call savemat(dble(yNaN), 'yNaN.debug')

  ! define a separate ZLB index to distinguish rates at ZLB from actually missing observations
  zlb = .false.
  ! for now: zlb matters only for poliyc rate not longer-term yields
  ! NOTE: key assumption is that missing obs for policy rate are solely due to ZLB!!
  zlb(yndxPolicyRate,:) = yNaN(yndxPolicyRate,:) 

  ALLOCATE (ypred(NNyy,maxhorizons), ynanpred(NNyy,maxhorizons), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Y)'
  END IF
  ypred    = readpreddata(datafile,Ny,NNyy,maxhorizons,Tdata,T+T0)
  yNaNpred = readprednandata(nandatafile,Ny,NNyy,maxhorizons,Tdata,T+T0)

  ! fix ELB predictions for the policy rate:
  ! they are read in as missing data, but should be set to the ELB
  ! assumptions:
  ! -- there are no missing data for the policy rate in sample
  ! -- ELB obs are encoded as exact zeros
  where (ypred(yndxPolicyRate,:) == 0.0d0)
     ypred(yndxPolicyRate,:)     = ELBound
     yNaNpred(yndxPolicyRate,:)  = .false.
  end where
  ! same for the shadowrate
  where (ypred(NNyy,:) == 0.0d0)
     ypred(NNyy,:)     = ELBound
     yNaNpred(NNyy,:)  = .false.
  end where
  ! call savemat(ypred, 'ypred.debug')
  ! call savemat(dble(ynanpred), 'nanypred.debug')
  ! stop 11

  ! setup number of states in extended state vector (used for forecasting; append 3 inflation lags)
  NNxx          = Nx + 3

  filext = '.SPECTREVAR.' // trim(datalabel) // '.VARlags' // trim(int2str(p)) // '.order' // trim(int2str(gapOrderKey)) // '.ar1svcorRbarSV.T' // trim(int2str(T)) // '.Tjumpoff' // trim(int2str(T0)) // '.run' // trim(int2str(runnum)) // '.dat'
  if (ZLBrejection) then
     filext = '.zlb' // filext
  end if
  if (.not. doPInoise) then
     filext = '.nopinoise' // filext
  end if
  if (doZLBasZeroes) then
     filext = '.zlbaszeroes' // filext
  end if
  if (doTimeStamp) filext = '.' // timestampstr() //  filext


  ! Model parameters and priors
  ALLOCATE (Ef0(Nf), sqrtVf0(Nf,Nf))
  ALLOCATE (SVol0(Nsv), Eh0(Nsv), minSV(Nsv), maxSV(Nsv), sqrtVh0(Nsv,Nsv), hSigmaT(Nsv,Nsv), hrho0(Nsv), hrhoV0(Nsv,Nsv))
  ALLOCATE (SVol0orth(Nsvorth), Ehorth0(Nsvorth), minSVorth(Nsvorth), maxSVorth(Nsvorth), Vhorth0(Nsvorth), horthvarT(Nsvorth), horthvarDof(Nsvorth), horthrho0(Nsvorth), horthrhoV0(Nsvorth))
  ! horthrhoV0 is a vector since prior is independend in case of orthogona;l SV's
  ALLOCATE (Eshockslopes0(Nshockslopes), sqrtVshockslopes0(Nshockslopes,Nshockslopes))

  
  
  ! ----------------------------------------------------------
  ! Sigma
  ! ----------------------------------------------------------


  ! SV
  ! TODO recalibrate?
  SVol0    = 0.1d0
  Eh0      = log(SVol0 ** 2)
  call eye(sqrtVh0, 5.0d0)
  hSigmaDof = Nsv + 1 + 10
  call eye(hSigmaT, (.2d0 ** 2) * dble(hSigmaDof - Nsv - 1))
  hrho0        = 0.8d0
  call eye(hrhoV0, 0.2d0 ** 2)
  minSV    = 0.001d0 ! 0.0001d0 causes failure for t=123
  maxSV    = 100.0d0

  SVol0orth    = 0.1d0
  Vhorth0      = 5.0d0 ** 2
  Ehorth0      = log(SVol0 ** 2)
  horthvarDof = 12
  horthvarT   = (0.2d0 ** 2) * dble(horthvarDof - 2) ! was 0.1 ...
  horthrho0        = 0.8d0
  horthrhoV0       = 0.2d0 ** 2
  minSVorth    = 0.001d0 ! 0.0001d0 causes failure for t=123
  maxSVorth    = 100.0d0



  ! cycle VAR coefficients
  Ef0  = 0.0d0

  lambda1Vf = 0.3d0
  lambda2Vf = 1.0d0
  lambda1Vf = 0.5d0
  lambda2Vf = 0.2d0
  call minnesotaVCVsqrt(sqrtVf0, Ny, p, lambda1Vf, lambda2Vf)

  ! gap-shockslopes
  Eshockslopes0  = 0.0d0
  call eye(sqrtVshockslopes0, 1.0d0) ! NEW: value was 1.0d0

  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data     = ' // trim(datalabel)
  print *, 'T        = ', T
  print *, 'T0       = ', T0
  ! print *, 'longrateHorizon= ', longrateHorizon
  print *, 'Nsim     = ', Nsim
  print *, 'Burnin   = ', Burnin
  print *, 'Nstreams = ', Nstreams
  if (doMDD) print *, 'Nparticles = ', Nparticles
  print *, 'p        = ', p
  print *, 'gapOrderKey= ', gapOrderKey
  print *, 'ELBound  = ', ELBound

  if (doPInoise) then
     print *, 'with PI noise'
  else
     print *, 'without PI noise'
  end if
  if (ZLBrejection) print *, 'with ZLB rejection'
  ! if (doLongRateTruncation) print *, 'with LongRateTruncation'
  if (doPDF) print *, 'with PDF computation'
  if (doMDD) print *, 'with MDD computation'
  if (doZLBasZeroes) print *, 'treating policy rate at ELB as data'
  print *,'runnum', runnum
  CALL HRULEFILL

  ! ----------------------------------------------------------------------------
  ! GENERATE AND STORE DRAWS FROM PRIORS
  ! ----------------------------------------------------------------------------
  Ndraws = max(10000, Nsim * Nstreams)

  ! ! hvar
  ! ALLOCATE (theta(Nsv,Ndraws))
  ! DO j = 1, Nsv
  !    DO k = 1, Ndraws
  !       call igammaDraw(theta(j,k), hvarT(j), hvarDof(j), VSLdefaultstream)
  !    END DO
  ! END DO
  ! filename  = 'hvar' // '.prior' // filext
  ! call savemat(transpose(theta), filename)
  ! WRITE (*,*) 'STORED PRIOR HVAR'
  ! DEALLOCATE (theta)

  ! TODO: hbar

  ! ! hrho
  ! ALLOCATE (theta(Nsv,Ndraws),theta2(Ndraws))
  ! theta = 10.0d0
  ! do j=1,Nsv
  !    OK = .FALSE.
  !    do while (.not. OK) 
  !       errcode = vdrnggaussian(VSLmethodGaussian, VSLdefaultstream, Ndraws, theta2, 0.0d0, 1.0d0)
  !       where (abs(theta(j,:)) > 1.0d0) theta(j,:)  = theta2 * sqrt(hrhoV0(j)) + hrho0(j)
  !       ok = all(abs(theta(j,:)) < 1.0d0)
  !    end do
  ! end do
  ! filename  = 'hrho' // '.prior' // filext
  ! call savemat(transpose(theta), filename)
  ! WRITE (*,*) 'STORED PRIOR hrho'
  ! DEALLOCATE (theta,theta2)

  ! ! SVol0
  ! ALLOCATE (theta(Nsv,Ndraws))
  ! errcode = vdrnggaussian(VSLmethodGaussian, VSLdefaultstream, Ndraws * Nsv, theta, 0.0d0, 1.0d0)
  ! FORALL (j = 1:Nsv) theta(j,:) = theta(j,:) * sqrt(Vh0(j)) + Eh0(j)
  ! theta = exp(theta * 0.5d0)
  ! filename  = 'SVol0' // '.prior' // filext
  ! call savemat(transpose(theta), filename)
  ! WRITE (*,*) 'STORED PRIOR SVol0'
  ! DEALLOCATE (theta)

  ! hrho SV prior truncation probabilities
  ! compute acceptace probability for hrho under prio
  ! ignoring PInoise case
  if (doMDD) then
     if (doPInoise) then
        print *, 'MDD not yet supported with PInoise'
        stop 1
     end if
     ! prepare critical values for upper/lower bound
     allocate (theta(Nsv + 2, 2),theta2(Nsv + 2, 2))

     forall (i=1:Nsv) theta(i,1) = (-1.0d0 - hrho0(i)) / sqrt(hrhoV0(i,i))
     forall (i=1:Nsv) theta(i,2) = (1.0d0 - hrho0(i)) / sqrt(hrhoV0(i,i))
     ! note: the previous line embeds knowledge that hrhoV0 is diagonal

     i = Nsv + 1
     theta(i,1) = (-1.0d0 - horthrho0(1)) / sqrt(horthrhoV0(1))
     theta(i,2) = (1.0d0 - horthrho0(1)) / sqrt(horthrhoV0(1))
     i = Nsv + 2
     theta(i,1) = (-1.0d0 - horthrho0(3)) / sqrt(horthrhoV0(3))
     theta(i,2) = (1.0d0 - horthrho0(3)) / sqrt(horthrhoV0(3))

     call vdcdfnorm((Nsv + 2) * 2,theta,theta2)
     ! call savemat(theta, 'theta.debug')
     ! call savemat(theta2, 'thetacdf.debug')

     hrhoPriorTruncationProbLog = log(1.0d0)

     do i=1,Nsv+2
        hrhoPriorTruncationProbLog =  hrhoPriorTruncationProbLog + log(theta2(i,2) - theta2(i,1))
        ! print *, 'hrhoPriorTruncationProb is', hrhoPriorTruncationProb
     end do

     print *, 'hrhoPriorTruncationProb is', exp(hrhoPriorTruncationProbLog)
     deallocate(theta,theta2)
  end if



  ! maxlambda
  ALLOCATE (theta(Ndraws,1))
  call simPriorMaxroot(theta(:,1), Ndraws, Ef0, sqrtVf0, Ny, p, VSLdefaultstream)
  filename  = 'maxlambda' // '.prior' // filext
  call savevec(theta(:,1), filename)
  WRITE (*,*) 'STORED PRIOR MAXLAMBDA'

  ! COMPUTE TRUNCATION PROBABILITY
  fPriorTruncationProbLog = dble(count(theta(:,1) .le. 1.0d0)) / dble(Ndraws)

  print *, 'f: Truncation Probability is', fPriorTruncationProbLog
  fPriorTruncationProbLog = log(fPriorTruncationProbLog)
  DEALLOCATE (theta)

  ! shockslopes
  ALLOCATE (theta(Nshockslopes,Ndraws))
  errcode = vdrnggaussian(VSLmethodGaussian, VSLdefaultstream, Ndraws * Nshockslopes, theta, 0.0d0, 1.0d0)
  CALL DTRMM('L', 'U', 'N', 'N', Nshockslopes, Ndraws, 1.0d0, sqrtVshockslopes0, Nshockslopes, theta, Nshockslopes)
  filename  = 'shockslopes' // '.prior' // filext
  call savemat(transpose(theta), filename)
  WRITE (*,*) 'STORED PRIOR SHOCKSLOPES'
  DEALLOCATE (theta)

  ! ----------------------------------------------------------------------------
  ! DONE: SIMULATING PRIORS
  ! ----------------------------------------------------------------------------



  ! allocate memory for draws
  ALLOCATE (DRAWmaxlambda(1,Nsim,Nstreams), DRAWf(Nf,Nsim,Nstreams), DRAWshockslopes(Nshockslopes,Nsim,Nstreams), DRAWstates(Nstates,0:T,Nsim,Nstreams), STAT=status)
  ALLOCATE (DRAWavgtermpremia(Nyield,Nsim,Nstreams), STAT=status)
  ALLOCATE (DRAWhorthvar(Nsvorth,Nsim,Nstreams), DRAWhorthbar(Nsvorth,Nsim,Nstreams), DRAWsvolorth(Nsvorth,0:T,Nsim,Nstreams), STAT=status)
  ALLOCATE (DRAWhorthrho(Nsvorth,Nsim,Nstreams), STAT=status)
  ALLOCATE (DRAWhsigma(Nhsigma,Nsim,Nstreams), DRAWhbar(Nsv,Nsim,Nstreams), DRAWsvol(Nsv,0:T,Nsim,Nstreams), STAT=status)
  ALLOCATE (DRAWhrho(Nsv,Nsim,Nstreams), STAT=status)
  ! IF (status /= 0) THEN
  !    WRITE (*,*) 'Allocation problem (draws)'
  ! END IF
  allocate(DRAWy(NNyy,maxhorizons,Nforecastdraws,Nsim,Nstreams))

  allocate(DRAWMDDlogprior(Nsim,Nstreams),DRAWMDDloglike(Nsim,Nstreams))
  DRAWMDDlogprior = 0.0d0
  DRAWMDDloglike  = 0.0d0

  call hrulefill
  WRITE(*,*) 'STARTING SWEEPS'
  call hrulefill

  !$OMP PARALLEL DO SHARED(DRAWMDDlogprior,DRAWMDDloglike,DRAWy,DRAWavgtermpremia,DRAWstates,DRAWsvol,DRAWf,DRAWshockslopes,DRAWhorthvar,DRAWhorthbar,DRAWhorthrho,DRAWsvolorth, DRAWhSigma,DRAWhbar,DRAWhrho,DRAWmaxlambda,runnum, seed), FIRSTPRIVATE(ZLBrejection,doPDF,doMDD,doPInoise, gapOrderKey, Nstreams, Nsim,Nparticles,Burnin,T,y,yNaN,zlb,Nstates,Nx, Nf, Eshockslopes0, sqrtVshockslopes0, Ef0, sqrtVf0, Eh0, sqrtVh0, minSV, maxSV, hSigmaT, hSigmaDof, hrho0, hrhoV0, Ehorth0, Vhorth0, minSVorth, maxSVorth, horthvarT, horthvarDof, horthrho0, horthrhoV0, NNxx), PRIVATE(VSLstream,TID,errcode,timer) SHARED(Nyield,Ny,Nbar,Ngap,Nshock,Nshockslopes,p,Nsv,Nhsigma,Nsvorth,NNyy) DEFAULT(NONE) SHARED(fPriorTruncationProbLog,hrhoPriorTruncationProbLog) SCHEDULE(STATIC)

  DO j=1,Nstreams

     TID = 0
     !$ TID = OMP_GET_THREAD_NUM()

     errcode = vslnewstream(VSLstream, vsl_brng_mt2203  + (runnum - 1) * (Nstreams + 1) + j, seed)  !  could also replace seed by a constant, say, 0
     if (errcode /= 0) then
        print *,'VSL new stream failed'
        stop 1
     end if


     WRITE(*,'(a16, i2, a5, i2, a8, i20, i20)') ' LAUNCHING SWEEP ', j, ' TID ', TID, ' Stream ', VSLstream%descriptor1, VSLstream%descriptor2
     ! print *,' TID ', TID,  VSLstream

     CALL initprogressbar(timer, 15.0d0, j)
     call thissampler(ZLBrejection, ELBound, doPDF, doMDD, doPInoise, gapOrderKey, T, Ny, Nyield, y, yNaN, zlb, DRAWMDDlogprior(:,j), DRAWMDDloglike(:,j), NNyy, NNxx, DRAWy(:,:,:,:,j), Nforecastdraws, maxhorizons, DRAWavgtermpremia(:,:,j), DRAWstates(:,:,:,j), Nstates, Nx, Nshock, Nbar, Ngap, p, DRAWf(:,:,j), Nf, Ef0, sqrtVf0, fPriorTruncationProbLog, DRAWmaxlambda(:,:,j), DRAWsvol(:,:,:,j), Nsv, Eh0, sqrtVh0, minSV, maxSV, DRAWhbar(:,:,j), DRAWhSigma(:,:,j), NhSigma, hSigmaT, hSigmaDof, DRAWhrho(:,:,j), hrho0, hrhoV0, hrhoPriorTruncationProbLog, DRAWsvolorth(:,:,:,j), Nsvorth, Ehorth0, Vhorth0, minSVorth, maxSVorth, DRAWhorthbar(:,:,j), DRAWhorthvar(:,:,j), horthvarT, horthvarDof, DRAWhorthrho(:,:,j), horthrho0, horthrhoV0, DRAWshockslopes(:,:,j), Nshockslopes, Eshockslopes0, sqrtVshockslopes0, Nsim, Burnin, Nparticles, VSLstream, timer)

     ! WRITE(*,*) 'STREAM', j, 'IS DONE.', ' (TID: ', TID, ')'
     errcode = vsldeletestream(VSLstream)     

  END DO
  !$OMP END PARALLEL DO 


  CALL HRULEFILL
  WRITE (*,*) 'ALL STREAMS ARE DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a20)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a20)') 'Data: ', datalabel
  ! WRITE(4,'(a20,I20)') 'longrateHorizon: ', longrateHorizon
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,I20)') 'Sims: ', Nsim
  WRITE(4,'(a20,I20)') 'Burnin: ', Burnin
  WRITE(4,'(a20,I20)') 'Nparticles: ', Nparticles
  WRITE(4,'(a20,I20)') 'Streams: ', Nstreams
  WRITE(4,'(a20,I20)') 'p: ', p
  WRITE(4,'(a20,F4.2)') 'ELBound: ', ELBound
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a60)') repeat('-',60)
  CLOSE(UNIT=4)
  CALL HRULEFILL




  ! ----------------------------------------------------------------------------
  ! STORE
  ! ----------------------------------------------------------------------------

  CALL HRULEFILL
  WRITE (*,*) 'STARTING W/STORAGE'
  CALL HRULEFILL

  ! STORE ESTIMATES
  ! Note: manual reshape avoids segmentation faults
  Ndraws = Nsim * Nstreams

  filename = 'YDATA' // filext
  call savemat(y, filename)

  filename = 'YNAN' // filext
  call savemat(dble(ynan), filename)

  ALLOCATE (theta(T,Ndraws))

  ! filename  = 'INFLATIONTREND.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxRealrate,1:T,k,j)
  ! call savemat(theta, filename)
  filename  = 'REALRATE' // filext
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED REALRATE'

  ! filename  = 'INFLATIONTREND.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxPIbar,1:T,k,j)
  ! call savemat(theta, filename)
  filename  = 'INFLATIONTREND' // filext
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED INFLATIONTREND'


  ! filename  = 'REALRATETREND.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxRbar,1:T,k,j)
  ! call savemat(theta, filename)
  filename  = 'REALRATETREND' // filext
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED REALRATETREND'

  do i = 1,Nyield
     ! filename  = 'LONGREALRATETREND.DRAWS' // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxLONGRbar+i-1,1:T,k,j)
     ! call savemat(theta, filename)
     filename  = 'LONGREALRATETREND'  // trim(int2str(i)) // filext
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED LONGREALRATETREND'

     ! filename  = 'LONGNOMINALRATETREND.DRAWS' // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxPIbar,1:T,k,j) + DRAWstates(ndxLONGRbar+1,1:T,k,j)
     ! call savemat(theta, filename)
     filename  = 'LONGNOMINALRATETREND'  // trim(int2str(i)) // filext
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED LONGNOMINALRATETREND'
  end do

  filename  = 'INFLATIONGAP' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxPIgap,1:T,k,j) 
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED INFLATIONGAP'

  filename  = 'UGAP' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxUgap,1:T,k,j)
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED UGAP'

  filename  = 'INTGAP' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxINTgap,1:T,k,j)
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED INTGAP'

  do i = 1,Nyield
     filename  = 'LONGINTGAP' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxLONGINTgap+i-1,1:T,k,j)
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED LONGINTGAP', i
  end do

  if (doPInoise) then
     filename  = 'INFLATIONNOISE' // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxPINOISE,1:T,k,j)
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED INFLATIONNOISE'
  end if

  ! filename  = 'SHADOWRATE.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxSHADOWRATE,1:T,k,j)
  ! call savemat(theta, filename)
  filename  = 'SHADOWRATE' // filext
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED SHADOWRATE'

  filename  = 'NOMINALRATETREND' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxPIbar,1:T,k,j) + DRAWstates(ndxRbar,1:T,k,j)
  CALL storeEstimatesOMP(theta,T,Ndraws,filename)
  WRITE (*,*) 'STORED NOMINALRATETREND'

  DO i=1,Nsv
     filename  = 'SV' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsvol(i,1:T,k,j)
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED SV', i
  END DO

  DO i=1,Nsvorth
     filename  = 'SVORTH' // trim(int2str(i)) // filext
     FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWsvolorth(i,1:T,k,j)
     CALL storeEstimatesOMP(theta,T,Ndraws,filename)
     WRITE (*,*) 'STORED SVORTH', i
  END DO

  DEALLOCATE(theta)

  ! STORE PARAMETERS

  ALLOCATE(theta(1,Ndraws))
  filename  = 'DELTAREALRATETREND' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWstates(ndxRbar,T,k,j) - DRAWstates(ndxRbar,0,k,j)
  CALL savevec(theta(1,:),filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED DELTAREALRATETREND'


  ALLOCATE(theta(Nf,Ndraws))
  filename  = 'F' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWf(:,k,j)
  CALL storeEstimatesOMP(theta,Nf,Ndraws,filename)
  ! filename  = 'F.draws' // filext
  ! CALL savemat(theta, filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED F'

  ALLOCATE(theta(Nshockslopes,Ndraws))
  filename  = 'SHOCKSLOPES' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWshockslopes(:,k,j)
  ! CALL storeEstimatesOMP(theta,Nshockslopes,Ndraws,filename)
  ! filename  = 'F.draws' // filext
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED SHOCKSLOPES'




  ALLOCATE(theta(NhSigma,Ndraws))
  filename  = 'HSIGMA.DRAWS' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhSigma(:,k,j)
  CALL savemat(transpose(theta), filename)
  filename  = 'HSIGMA' // filext
  CALL storeEstimates(theta,NhSigma,Ndraws,filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HSIGMA'


  ALLOCATE(theta(Nsv,Ndraws))
  filename  = 'HBAR' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhbar(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HBAR'

  ALLOCATE(theta(Nsv,Ndraws))
  filename  = 'HRHO' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhrho(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HRHO'

  ALLOCATE(theta(Nsvorth,Ndraws))
  filename  = 'HORTHVAR' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhorthvar(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HORTHVAR'

  ALLOCATE(theta(Nsvorth,Ndraws))
  filename  = 'HORTHBAR' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhorthbar(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HORTHBAR'

  ALLOCATE(theta(Nsvorth,Ndraws))
  filename  = 'HORTHRHO' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWhorthrho(:,k,j)
  CALL savemat(transpose(theta), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED HORTHRHO'

  ALLOCATE(theta(1,Ndraws))
  filename  = 'MAXLAMBDA' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWmaxlambda(:,k,j)
  ! CALL storeEstimates(theta,1,Ndraws,filename)
  CALL savevec(theta(1,:), filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED MAXLAMBDA'

  ! STORE DRAWS OF AVG TERM PREMIA
  ALLOCATE(theta(Nyield,Ndraws))
  filename  = 'AVGTERMPREMIA' // filext
  FORALL (k=1:Nsim,j=1:Nstreams) theta(:,(j-1) * Nsim + k) = DRAWavgtermpremia(:,k,j)
  CALL savemat(theta, filename)
  DEALLOCATE(theta)
  WRITE (*,*) 'STORED AVG TERM PREMIA'


  if (doPDF) then 
     ! collect 
     call hrulefill
     WRITE (*,*) 'COLLECTING PREDICTIVE DENSITY'
     call hrulefill

     filename = 'YPRED' // filext
     call savemat(ypred(:,horizons), filename)

     filename = 'YNANPRED' // filext
     call savemat(dble(yNaNpred(:,horizons)), filename)

     filename = 'YHORIZONS' // filext
     call savevec(dble(horizons), filename)


     ! DRAWY
     Ndraws = Nsim * Nstreams * Nforecastdraws

     ALLOCATE(theta(maxhorizons,Ndraws))
     ALLOCATE(yhatmedian(NNyy,maxhorizons))

     ! WRITE OUT YDRAW     
     DO i=1,NNyy
        FORALL (h=1:maxhorizons,k=1:Nsim,j=1:Nstreams,n=1:Nforecastdraws) theta(h,(j-1) * (Nsim * Nforecastdraws) + (k - 1) * Nforecastdraws + n) = DRAWy(i,h,n,k,j)

        filename  = 'YDRAW' // trim(int2str(i)) // filext
        CALL storeEstimatesOMP(theta,maxhorizons,Ndraws,filename)
        WRITE (*,*) 'STORED YDRAW', i

        ! note: theta is returned in sorted state
        yhatmedian(i,:) = theta(:,floor(real(Ndraws) * 0.5))
     END DO

     DEALLOCATE(theta)

     ! STORE MEDIAN FORECAST
     filename  = 'YHATMEDIAN' // filext
     CALL savemat(yhatmedian, filename)
     WRITE (*,*) 'STORED YHATMEDIAN'
     filename  = 'YHATMEDIANERROR' // filext
     yhatmedian = ypred - yhatmedian
     where (yNaNpred) yhatmedian = -999.0d0
     CALL savemat(yhatmedian, filename)
     WRITE (*,*) 'STORED YHATMEDIANERROR'
     DEALLOCATE(yhatmedian)

     ! STORE MEAN FORECAST
     ALLOCATE(yhatmean(NNyy,maxhorizons))
     forall (h=1:maxhorizons,i=1:NNyy) yhatmean(i,h) = sum(DRAWy(i,h,:,:,:)) / dble(Ndraws)
     filename  = 'YHATMEAN' // filext
     CALL savemat(yhatmean, filename)
     WRITE (*,*) 'STORED YHATMEAN'
     filename  = 'YHATMEANERROR' // filext
     yhatmean = ypred - yhatmean
     where (yNaNpred) yhatmean = -999.0d0
     CALL savemat(yhatmean, filename)
     WRITE (*,*) 'STORED YHATMEANERROR'
     DEALLOCATE(yhatmean)


     ! write out the scores in a separate loop
     WRITE (*,*) 'COMPUTING SCORES ... '
     ALLOCATE(crpsScore(NNyy,maxhorizons))
     ALLOCATE(logScore(NNyy,maxhorizons))

     !$OMP PARALLEL DO SHARED(theta, DRAWy, ypred, crpsScore, logScore, Nstreams, Nsim, Ndraws, NNyy) PRIVATE(theta1,h,k,j,n) DEFAULT(NONE) SCHEDULE(STATIC)
     DO i=1,NNyy
        ALLOCATE(theta1(Ndraws))
        do h = 1,maxhorizons
           FORALL (k=1:Nsim,j=1:Nstreams,n=1:Nforecastdraws) theta1((j-1) * (Nsim * Nforecastdraws) + (k - 1) * Nforecastdraws + n) = DRAWy(i,h,n,k,j)

           crpsScore(i,h) = crps(ypred(i,h), theta1, Ndraws)
           logScore(i,h)  = logGaussianScore(ypred(i,h), theta1, Ndraws)

        end do
        DEALLOCATE(theta1)
     END DO
     !$OMP END PARALLEL DO 
     ! DEALLOCATE(theta1)

     filename  = 'CRPS' // filext
     CALL savemat(crpsScore, filename)
     WRITE (*,*) 'STORE CRPSSCORE'

     filename  = 'LOGSCORE' // filext
     CALL savemat(logScore, filename)
     WRITE (*,*) 'STORE LOGSCORE'

     deallocate(crpsScore)
     deallocate(logScore)


     ! MV logGaussianScore
     ALLOCATE(logScoreMV(maxhorizons))
     ALLOCATE(theta(Ny,Ndraws))
     do h = 1,maxhorizons
        FORALL (i=1:Ny,k=1:Nsim,j=1:Nstreams,n=1:Nforecastdraws) theta(i,(j-1) * (Nsim * Nforecastdraws) + (k - 1) * Nforecastdraws + n) = DRAWy(i,h,n,k,j)

        logScoreMV(h) = logGaussianScoreMV(ypred(:,h), theta, Ny, Ndraws)
     end do

     filename  = 'MVLOGSCORE' // filext
     CALL savevec(logScoreMV, filename)
     WRITE (*,*) 'STORE MVLOGSCORE'


     DEALLOCATE(theta)
     DEALLOCATE(logScoreMV)

  end if ! doPDF


  ! MDD
  if (doMDD) then

     filename  = 'MDDlogprior' // filext
     CALL savemat(DRAWMDDlogprior, filename)
     WRITE (*,*) 'STORE MDDlogprior'
     filename  = 'MDDloglike' // filext
     CALL savemat(DRAWMDDloglike, filename)
     WRITE (*,*) 'STORE MDDloglike'

     ! collect paramters draws into theta
     ! theta :
     ! - 1 : rbarvar
     ! - Nf : f
     ! - Nshockslopes : shockslopes
     ! - 2 * Nsv : hrho, hbar
     ! - Nhsigma : hSigma ! TODO: does it matter that hSigma is stored as sqrt?
     ! - Nsvorth times 3 (horthrho, horthbar, horthvar)


     ! set up theta
     Ntheta = Nf + Nshockslopes + Nsv * 2 + Nhsigma + 6
     if (doPInoise) Ntheta = Ntheta + 3

     Ndraws = Nsim * Nstreams

     allocate (theta(Ntheta,Ndraws), STAT=status)
     if (status /= 0) then
        write (*,*) 'Allocation problem (MDD theta)'
     end if

     ndxtheta = 0


     ! f
     FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta+1:ndxtheta+Nf,(j-1) * Nsim + k) = DRAWf(:,k,j)
     ndxtheta = ndxtheta + Nf

     ! shockslopes
     FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta+1:ndxtheta+Nshockslopes,(j-1) * Nsim + k) = DRAWshockslopes(:,k,j)
     ndxtheta = ndxtheta + Nshockslopes

     ! hrho
     FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta+1:ndxtheta+Nsv,(j-1) * Nsim + k) = DRAWhrho(:,k,j)
     ndxtheta = ndxtheta + Nsv

     ! hbar
     FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta+1:ndxtheta+Nsv,(j-1) * Nsim + k) = DRAWhbar(:,k,j)
     ndxtheta = ndxtheta + Nsv

     ! hSigma
     FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta+1:ndxtheta+Nhsigma,(j-1) * Nsim + k) = DRAWhsigma(:,k,j)
     ndxtheta = ndxtheta + Nhsigma

     ! horth: hrho
     do i = 1,Nsvorth
        if (i == ndxSVpinoise .and. (.not. doPInoise)) cycle
        ndxtheta = ndxtheta + 1
        FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta,(j-1) * Nsim + k) =  DRAWhorthrho(i,k,j)
     end do

     ! horth: hbar
     do i = 1,Nsvorth
        if (i == ndxSVpinoise .and. (.not. doPInoise)) cycle
        ndxtheta = ndxtheta + 1
        FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta,(j-1) * Nsim + k) =  DRAWhorthbar(i,k,j)
     end do

     ! horth: hvar
     do i = 1,Nsvorth
        if (i == ndxSVpinoise .and. (.not. doPInoise)) cycle
        ndxtheta = ndxtheta + 1
        FORALL (k=1:Nsim,j=1:Nstreams) theta(ndxtheta,(j-1) * Nsim + k) =  DRAWhorthvar(i,k,j)
     end do

     filename  = 'MDDtheta' // filext
     CALL savemat(theta, filename)

     ! simulate truncation probablity for theta when drawn from MV normal
     ! thetabar = sum(theta,2) / dble(Ndraws)
     ! Vtheta   = sqrtvcvTR(Ntheta,Ndraws,theta)
     allocate (ndxf(Nf),ndxrho(Nsv+2))
     forall (k=1:Nf)  ndxf(k)   = k
     forall (k=1:Nsv) ndxrho(k) = Nf + Nshockslopes + k
     ndxrho(Nsv+1) =  Nf + Nshockslopes + Nsv * 2 + Nhsigma + 1
     ndxrho(Nsv+2) =  Nf + Nshockslopes + Nsv * 2 + Nhsigma + 2
     call thisMDDthetaTruncationProb(gewekeTruncLogProb, Ny, p, ndxf, Nsv + 2, ndxrho, gewekeTruncationNsim, Ndraws, Ntheta, theta, VSLdefaultstream)
     deallocate (ndxf,ndxrho)

     ! store various truncation probabilities
     ! - f prior
     ! - hrho prio
     ! - f and rho posterior
     allocate (theta1(3))
     theta1(1) = fPriorTruncationProbLog
     theta1(2) = hrhoPriorTruncationProbLog
     theta1(3) = gewekeTruncLogProb
     filename  = 'MDDtrunclogprobs' // filext
     call savevec(theta1, filename)
     deallocate(theta1)


     ! note: calling gewekeMDD with 2D DRAWMDD* objects performs implicit reshape; note that the draws in theta have been collected in the same order as a reshape of DRAWMDD* would produce
     allocate (theta1(size(gewekeTau)))

     call hrulefill

     do i=1,size(gewekeTau)
        theta1(i) = gewekelogMDD(Ndraws, Ntheta, theta, DRAWMDDlogprior, DRAWMDDloglike, gewekeTau(i), gewekeTruncLogProb)
        print *, 'GEWEKE estimate of logMDD with tau =', gewekeTau(i), ' is: ', theta1(i)
     end do

     call hrulefill

     filename  = 'MDDstats' // filext
     CALL savevec(theta1, filename)

     deallocate(theta,theta1)
  end if
  deallocate(DRAWMDDlogprior,DRAWMDDloglike)




  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
  ! GELMAN
  ! ----------------------------------------------------------------------------

  ALLOCATE (GELMANstates(Nstates,T), GELMANf(Nf), GELMANshockslopes(Nshockslopes), STAT=status)
  ALLOCATE (GELMANsvol(Nsv,T), GELMANhbar(Nsv), GELMANhsigma(NhSigma),  STAT=status)
  ALLOCATE (GELMANhrho(Nsv),  STAT=status)
  ALLOCATE (GELMANsvolorth(Nsvorth,T), GELMANhorthbar(Nsvorth), GELMANhorthvar(Nsvorth),  STAT=status) 
  ALLOCATE (GELMANhorthrho(Nsvorth),  STAT=status)
  ! IF (status /= 0) THEN
  !   WRITE (*,*) 'Allocation problem (Gelman statistics)'
  ! END IF

  CALL HRULEFILL
  WRITE (*,*) 'GELMAN STORAGE ALLOCATED'
  CALL HRULEFILL


  !$OMP PARALLEL SHARED(DRAWstates,GELMANstates,DRAWf,GELMANF,GELMANshockslopes, DRAWshockslopes, DRAWsvol,GELMANsvol,DRAWsvolorth,GELMANsvolorth, DRAWhorthvar, GELMANhorthvar, DRAWhSigma, GELMANhSigma, DRAWhbar, GELMANhbar, DRAWhorthbar, GELMANhorthbar, DRAWhorthrho, GELMANhorthrho, DRAWhrho, GELMANhrho)

  !$OMP DO 
  DO j = 1, Nstates
     if (j == ndxUgap) then
        GELMANstates(j,:) = 1.0d0
        forall (i=1:Nstreams,k=1:Nsim) DRAWstates(j,1:T,k,i) = DRAWstates(j,1:T,k,i) - y(1,:)
        print *,'max ugap delta:', maxval(abs(DRAWstates(j,1:T,:,:)))
     else
        DO k = 1,T
           call GelmanTest1(GELMANstates(j,k), Drawstates(j,k,:,:), Nsim, Nstreams)
        END DO
     end if

  END DO
  !$OMP END DO

  if (.not. doPInoise) GELMANstates(ndxPInoise,:) = 1.0d0

  !$OMP DO 
  DO k = 1,T
     DO j = 1, Nsv
        call GelmanTest1(GELMANsvol(j,k), DRAWsvol(j,k,:,:), Nsim, Nstreams)
     END DO
  END DO
  !$OMP END DO

  !$OMP DO 
  DO k = 1,T
     DO j = 1, Nsvorth
        call GelmanTest1(GELMANsvolorth(j,k), DRAWsvolorth(j,k,:,:), Nsim, Nstreams)
     END DO
  END DO
  !$OMP END DO
  if (.not. doPInoise) GELMANsvolorth(ndxSVpinoise,:) = 1.0d0

  !$OMP DO 
  DO j = 1, Nsvorth
     call GelmanTest1(GELMANhorthvar(j), DRAWhorthvar(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO
  if (.not. doPInoise) GELMANhorthvar(ndxSVpinoise) = 1.0d0
  !$OMP DO 
  DO j = 1, Nsvorth
     call GelmanTest1(GELMANhorthbar(j), DRAWhorthbar(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO
  if (.not. doPInoise) GELMANhorthbar(ndxSVpinoise) = 1.0d0

  !$OMP DO 
  DO j = 1, NhSigma
     call GelmanTest1(GELMANhSigma(j), DRAWhSigma(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nsv
     call GelmanTest1(GELMANhbar(j), DRAWhbar(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nsv
     call GelmanTest1(GELMANhrho(j), DRAWhrho(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nsvorth
     call GelmanTest1(GELMANhorthrho(j), DRAWhorthrho(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO
  if (.not. doPInoise) GELMANhorthrho(ndxSVpinoise) = 1.0d0

  !$OMP DO 
  DO j = 1, Nf
     call GelmanTest1(GELMANf(j), DRAWf(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP DO 
  DO j = 1, Nshockslopes
     call GelmanTest1(GELMANshockslopes(j), DRAWshockslopes(j,:,:), Nsim, Nstreams)
  END DO
  !$OMP END DO

  !$OMP END PARALLEL 

  CALL HRULEFILL
  WRITE (*,*) 'GELMAN STATISTICS ARE DONE!'
  CALL HRULEFILL

  IF (ALL(ABS(GELMANstates - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for the STATES'
  ELSE
     WRITE(*,*) 'STATES: GELMAN FAILURE, Max SRstat=', maxval(GELMANstates)
  END IF
  filename = 'GELMAN.states' // filext
  call savemat(GELMANstates, filename)

  IF (ALL(ABS(GELMANsvol - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for the SVOL'
  ELSE
     WRITE(*,*) 'SVOL: GELMAN FAILURE, Max SRstat=', maxval(GELMANsvol)
  END IF
  filename = 'GELMAN.svol' // filext
  call savemat(GELMANsvol, filename)

  IF (ALL(ABS(GELMANsvolorth - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for the SVOLORTH'
  ELSE
     WRITE(*,*) 'SVOLORTH: GELMAN FAILURE, Max SRstat=', maxval(GELMANsvolorth)
  END IF
  filename = 'GELMAN.svolorth' // filext
  call savemat(GELMANsvolorth, filename)

  IF (ALL(ABS(GELMANhSigma - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HSIGMA'
  ELSE
     WRITE(*,*) 'HSIGMA: GELMAN FAILURE, Max SRstat=', maxval(GELMANhsigma)
  END IF
  filename = 'GELMAN.hsigma' // filext
  call savevec(GELMANhSigma, filename)

  IF (ALL(ABS(GELMANhorthvar - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HORTHVAR'
  ELSE
     WRITE(*,*) 'HORTHVAR: GELMAN FAILURE, Max SRstat=', maxval(GELMANhorthvar)
  END IF
  filename = 'GELMAN.horthvar' // filext
  call savevec(GELMANhorthvar, filename)

  IF (ALL(ABS(GELMANhbar - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HBAR'
  ELSE
     WRITE(*,*) 'HBAR: GELMAN FAILURE, Max SRstat=', maxval(GELMANhbar)
  END IF
  filename = 'GELMAN.hbar' // filext
  call savevec(GELMANhbar, filename)

  IF (ALL(ABS(GELMANhorthbar - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HORTHBAR'
  ELSE
     WRITE(*,*) 'HORTHBAR: GELMAN FAILURE, Max SRstat=', maxval(GELMANhorthbar)
  END IF
  filename = 'GELMAN.horthbar' // filext
  call savevec(GELMANhorthbar, filename)

  IF (ALL(ABS(GELMANhrho - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HRHO'
  ELSE
     WRITE(*,*) 'HRHO: GELMAN FAILURE, Max SRstat=', maxval(GELMANhrho)
  END IF
  filename = 'GELMAN.hrho' // filext
  call savevec(GELMANhrho, filename)

  IF (ALL(ABS(GELMANhorthrho - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for HORTHRHO'
  ELSE
     WRITE(*,*) 'HORTHRHO: GELMAN FAILURE, Max SRstat=', maxval(GELMANhorthrho)
  END IF
  filename = 'GELMAN.horthrho' // filext
  call savevec(GELMANhorthrho, filename)

  IF (ALL(ABS(GELMANshockslopes - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for SHOCKSLOPES'
  ELSE
     WRITE(*,*) 'SHOCKSLOPES: GELMAN FAILURE, Max SRstat=', maxval(GELMANshockslopes)
  END IF
  filename = 'GELMAN.shockslopes' // filext
  call savevec(GELMANshockslopes, filename)

  IF (ALL(ABS(GELMANf - 1) < 0.2)) THEN
     WRITE(*,*) 'Gelman found decent convergence for F'
  ELSE
     WRITE(*,*) 'F: GELMAN FAILURE, Max SRstat=', maxval(GELMANf)
  END IF
  filename = 'GELMAN.f' // filext
  call savevec(GELMANf, filename)

  DEALLOCATE (GELMANstates, GELMANsvol, GELMANf, GELMANshockslopes, GELMANhorthvar, GELMANhorthbar, GELMANsvolorth, GELMANhSigma, GELMANhbar)
  DEALLOCATE (GELMANhrho, GELMANhorthrho)

  ! ----------------------------------------------------------------------------
  ! FINISHED: GELMAN
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! CLEANUP 1
  ! ----------------------------------------------------------------------------

  DEALLOCATE (DRAWmaxlambda, DRAWstates, DRAWsvol, DRAWf, DRAWshockslopes, DRAWhSigma, DRAWhbar)
  DEALLOCATE (DRAWavgtermpremia)
  DEALLOCATE (DRAWsvolorth, DRAWhorthvar, DRAWhorthbar)

  ! DEALLOCATE (DRAWypdf, DRAWyforecast, DRAWymedian, DRAWyproblb, DRAWycondvar)
  DEALLOCATE (DRAWy)
  DEALLOCATE (y, yNaN,zlb)
  DEALLOCATE (ypred, yNaNpred)
  DEALLOCATE (SVol0,Eh0, sqrtVh0, minSV, maxSV, hSigmaT)
  DEALLOCATE (SVol0orth, Ehorth0, Vhorth0, minSVorth, maxSVorth, horthvarT, horthvarDof)
  DEALLOCATE (hrho0, hrhoV0)
  DEALLOCATE (horthrho0, horthrhoV0)
  DEALLOCATE (Ef0,sqrtVf0)
  DEALLOCATE (Eshockslopes0,sqrtVshockslopes0)

  ! ----------------------------------------------------------------------------
  ! CLEANUP 2
  ! ----------------------------------------------------------------------------
  ! VSLstreams
  errcode = vsldeletestream(VSLdefaultstream)     

  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  SUBROUTINE minnesotaVCVsqrt(sqrtVf0, N, p, lambda1, lambda2)

    integer, intent(in) :: N, p
    double precision, intent(inout), dimension(N*N*p,N*N*p) :: sqrtVf0
    double precision, intent(in) :: lambda1, lambda2

    integer :: ndxVec, ndxLHS = 1, ndxLag = 1, ndxRHS = 1
    integer :: Nf

    Nf = N * N * p
    sqrtVf0 = 0.0d0

    ndxVec = 0
    do ndxLHS = 1,N
       do ndxLag = 1,p
          do ndxRHS = 1,N

             ndxVec = ndxVec + 1

             if (ndxLHS == ndxRHS) then
                sqrtVf0(ndxVec,ndxVec) = lambda1 / dble(ndxLag)
             else
                sqrtVf0(ndxVec,ndxVec) = lambda1 * lambda2 / dble(ndxLag)
             end if

          end do
       end do
    end do

  END SUBROUTINE minnesotaVCVsqrt

  SUBROUTINE getarguments(ZLBrejection,doPDF, doMDD, doPInoise, doZLBasZeroes, T, T0, gapOrderKey, datalabel, Nyield, p, Nparticles, Nsim, burnin, Nstreams, runnum)

    INTENT(INOUT) ZLBrejection, doPDF, doMDD, doPInoise, doZLBasZeroes, T, T0, datalabel, Nyield, p, Nsim, Nparticles, burnin, Nstreams

    integer, intent(inout) :: gapOrderKey

    INTEGER :: T, T0, Nparticles, Nsim, Nyield, p, burnin, Nstreams, dummy, runnum
    LOGICAL :: ZLBrejection,doPDF, doMDD, doPInoise, doZLBasZeroes
    INTEGER :: counter
    CHARACTER (LEN=100) :: datalabel
    CHARACTER(len=32) :: arg

    counter = 0

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          ZLBrejection = .true.
       else
          ZLBrejection = .false.
       end if
    END IF


    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doPDF = .true.
       else
          doPDF = .false.
       end if
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doMDD = .true.
       else
          doMDD = .false.
       end if
    END IF


    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doPInoise = .true.
       else
          doPInoise = .false.
       end if
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doZLBasZeroes = .true.
       else
          doZLBasZeroes = .false.
       end if
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') T
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') T0
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') gapOrderKey
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, datalabel) 
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') Nyield
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') p
    END IF


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
       READ(arg, '(i20)') Burnin
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN    
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') Nstreams
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN    
       CALL get_command_argument(counter, arg)
       READ(arg, '(i20)') runnum
    END IF
    
  END SUBROUTINE getarguments



  ! -----------------------------------------------------------------
  SUBROUTINE readdata(y,filename,Ny,T,T0)
    IMPLICIT NONE

    INTENT(IN) :: filename,Ny,T,T0
    INTENT(INOUT) :: y
    CHARACTER (LEN=200) :: filename
    CHARACTER (LEN=500) :: fmtstr

    INTEGER i, T, T0, Ny
    DOUBLE PRECISION :: y(Ny,T), data(Ny,T)


    fmtstr = es30d16(Ny)
    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')
    ! skio first T0 lines
    DO i=1,T0
       READ(4,*) 
    END DO
    DO i=1,T
       READ(4,fmtstr) data(:,i)
       ! print *, i, data(:,i)
    END DO
    CLOSE(UNIT=4)

    ! transform data
    y = data ! trivial, since already done in matlab

  END SUBROUTINE readdata
  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  SUBROUTINE readnandata(nanny,filename,Ny,T,T0)
    IMPLICIT NONE

    INTENT(IN) :: filename,T,T0,Ny
    INTENT(INOUT) :: nanny
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    LOGICAL, DIMENSION(:,:) :: nanny
    INTEGER :: work(Ny)

    INTEGER i, j, T, T0, Ny

    fmtstr = '(I2' // repeat(',I2', Ny-1) // ')'

    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')
    DO i=1,T0
       READ(4,*) 
    END DO
    DO i=1,T
       READ(4,fmtstr) (work(j), j=1,Ny)
       WHERE (work == 1) 
          nanny(:,i) = .TRUE.
       ELSEWHERE
          nanny(:,i) = .FALSE.
       END WHERE
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readnandata
  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  function readpreddata(filename,Ny,NNyy,maxhorizons,Tdata,T) result(ypred)
    IMPLICIT NONE

    INTENT(IN) :: filename,Ny,NNyy,maxhorizons,Tdata,T

    CHARACTER (LEN=200) :: filename
    CHARACTER (LEN=500) :: fmtstr

    INTEGER row, T, Ny, NNyy, maxhorizons, Tdata
    DOUBLE PRECISION :: ypred(NNyy,maxhorizons), data(Ny,Tdata)

    ypred = dble(-999) 

    fmtstr = es30d16(Ny)

    ! Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    ! read all data
    do row=1,Tdata
       READ(4,fmtstr) (data(k,row), k=1,Ny) 
    end do
    CLOSE(UNIT=4)

    ! extract ypred
    row = 1
    DO while (T + row <= min(Tdata,T+maxhorizons)) 
       ypred(1:Ny,row) = data(:,T+row)
       ! construct 4q-MA of inflation
       ypred(Ny+1,row) = sum(data(2,T+row-3:T+row)) / 4.0d0
       ! add shadow rate
       ypred(Ny+2,row) = data(3,T+row)
       row = row + 1
    END DO




  end function readpreddata
  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  function readprednandata(filename,Ny,NNyy,maxhorizons,Tdata,T) result(nanny)

    IMPLICIT NONE

    INTENT(IN) :: filename,Ny, NNyy,maxhorizons,Tdata,T

    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    INTEGER row, k, Ny, NNyy, maxhorizons, Tdata, T
    LOGICAL, DIMENSION(NNyy,maxhorizons) :: nanny
    INTEGER :: data(Ny,Tdata)

    fmtstr = '(I2' // repeat(',I2', Ny-1) // ')'

    nanny = .TRUE.

    ! Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    ! read all data
    do row=1,Tdata
       READ(4,fmtstr) (data(k,row), k=1,Ny) 
    end do
    CLOSE(UNIT=4)

    ! extract yNaNpred
    row = 1
    DO while (T + row <= min(Tdata,T+maxhorizons)) 
       WHERE (data(:,T+row) == 1) 
          nanny(1:Ny,row) = .TRUE.
       ELSEWHERE
          nanny(1:Ny,row) = .FALSE.
       END WHERE
       ! 4q MA of inflation
       if (any(data(2,T+row-3:T+row) == 1)) then
          nanny(Ny+1,row) = .TRUE.
       else
          nanny(Ny+1,row) =- .FALSE.
       end if
       ! shadow rate
       nanny(Ny+2,row) = nanny(3,row)

       row = row + 1
    END DO

  end function readprednandata
  ! -----------------------------------------------------------------


END PROGRAM main

! -----------------------------------------------------------------

! @\newpage\subsection{thissampler}@
SUBROUTINE thissampler(ZLBrejection, ELBound, doPDF, doMDD, doPInoise, gapOrderKey, T, Ny, Nyield, y, yNaN, zlb, DRAWMDDlogprior, DRAWMDDloglike, NNyy, NNxx, DRAWy, Nforecastdraws, maxhorizons, DRAWavgtermpremia, DRAWstates, Nstates, Nx, Nw, Nbar, Ngap, p, DRAWf, Nf, Ef0, sqrtVf0, fPriorTruncationProbLog, DRAWmaxlambda, DRAWsvol, Nsv, Eh0, sqrtVh0, minSV, maxSV, DRAWhbar, DRAWhSigma, NhSigma, hSigmaT, hSigmaDof, DRAWhrho, hrho0, hrhoV0, hrhoPriorTruncationProbLog, DRAWsvolorth, Nsvorth, Ehorth0, Vhorth0, minSVorth, maxSVorth, DRAWhorthbar, DRAWhorthvar, horthvarT, horthvarDof, DRAWhorthrho, horthrho0, horthrhoV0,  DRAWshockslopes, Nshockslopes, Eshockslopes0, sqrtVshockslopes0, Nsim, Burnin, Nparticles, VSLstream,timer)

  use embox, only: hrulefill, savemat, savevec, int2str
  use blaspack, only: eye, vectortimesmatrix, predictstate, vech
  use gibbsbox, only: VARmaxroot, bayesVARSV, bayesdiffuseVARSV, igammadraw, variancedraw, iwishcholdraw, vcvcholDrawTR, bayesregcholeskidiffuseNorth, drawAR1, drawAR1correlated, stochvolKSCar1plus, SVHcholKSCAR1cor, bayesUnivariateRegressionSlope, bayesAR1SUR ! , bayesVARSVbarshock, bayesDiffuseRegressionSlopeSV, bayesRegressionSlopeSV, SVHcholeskiAR1diffuseslopesplus, svhcholeskiKSCar1
  use statespacebox, only: DLYAP, samplerA3B3C3nanscalar
  use densitybox, only : predictivedensitySVtruncLB, predictivedensitySVCORdraws, igammalogpdf, MVnormlogpdf, normlogpdf, iwishlogpdf, logmvgamma, logdetSym

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: DRAWstates, DRAWf, DRAWshockslopes, VSLstream, timer
  INTENT(INOUT) :: DRAWavgtermpremia
  INTENT(INOUT) :: DRAWMDDlogprior, DRAWMDDloglike
  INTENT(INOUT) :: DRAWsvol, DRAWhbar, DRAWhSigma, DRAWhrho
  INTENT(INOUT) :: DRAWsvolorth, DRAWhorthbar, DRAWhorthvar, DRAWhorthrho
  ! INTENT(INOUT) :: DRAWypdf, DRAWyforecast, DRAWymedian, DRAWyprobLB, DRAWycondvar
  INTENT(INOUT) :: DRAWy
  INTENT(IN)    :: ZLBrejection, doPDF, doMDD, doPInoise, T, Ny, Nyield, y, yNaN, zlb, Nforecastdraws, maxhorizons, NNyy, NNxx, Nstates, Nx, Nw, Nbar, Ngap, p, Nf, Nsim, Burnin, Ef0, sqrtVf0, Nshockslopes, Eshockslopes0, sqrtVshockslopes0
  INTENT(IN)    :: Nsv, Eh0, sqrtVh0, minSV, maxSV, NhSigma, hSigmaT, hSigmaDof, hrho0, hrhoV0
  INTENT(IN)    :: Nsvorth, Ehorth0, Vhorth0, minSVorth, maxSVorth, horthvarT, horthvarDof, horthrho0, horthrhoV0

  double precision, intent(in) :: ELBound
  integer, intent(in) :: gapOrderKey

  double precision, intent(in) :: fPriorTruncationProbLog, hrhoPriorTruncationProbLog


  LOGICAL :: ZLBrejection, doPDF, doMDD, doPInoise

  INTEGER :: J, I, K, T, Nsim, thisdraw, joflastdraw, Burnin, TotalSim, Nstates, Nx, Nw, Nbar, Ngap, Ngaplags, Ny, Nf, Nshockslopes, p, status

  INTEGER :: Nsv, Nhsigma
  INTEGER :: Nsvorth

  double precision :: sharedone

  INTEGER :: ndxPIbar, ndxRbar, ndxLONGRbar, ndxUgap, ndxINTgap, ndxLONGINTgap, Nyield, ndxPIgap, ndxPInoise, ndxGapStart, ndxGapStop, ndxGapLagStop
  INTEGER :: ndxShockPIbar, ndxShockRbar, ndxShockGapStart, ndxShockPIgap, ndxShockUgap, ndxShockINTgap, ndxShockLONGINTgap, ndxShockPInoise
  INTEGER :: ndxSVpibar, ndxSVpigap, ndxSVugap, ndxSVintgap, ndxSVlongintgap, ndxSVpinoise, ndxSVrbar
  INTEGER :: yndxLONGINT, yndxInflation, yndxPolicyrate, yndxUgap

  ! forecast stuff
  ! integer :: Nhorizons, horizons(Nhorizons)
  integer :: maxhorizons
  ! double precision, dimension(NNyy, maxhorizons) :: ypred
  ! logical :: yNaNpred(Ny, maxhorizons)
  double precision, dimension(NNyy) :: ytrunclb
  logical, dimension(NNyy) :: ytruncated
  ! double precision, dimension(NNyy, Nhorizons, Nsim) :: DRAWypdf, DRAWyforecast, DRAWymedian, DRAWyprobLB, DRAWycondvar
  integer :: Nforecastdraws
  double precision, dimension(NNyy, maxhorizons, Nforecastdraws, Nsim) :: DRAWy

  integer :: NNyy, NNxx 
  double precision :: AA(NNxx,NNxx), BB(NNxx,Nw), CC(NNyy,NNxx), xx(NNxx)

  ! particlelike
  integer, intent(in) :: Nparticles 

  ! OTHER
  type (vsl_stream_state) :: VSLstream
  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN, zlb
  LOGICAL, DIMENSION(0:T) :: atELB ! for particlelike, considers only whether policyrate is binding
  DOUBLE PRECISION, DIMENSION(Nsim) :: DRAWMDDlogprior, DRAWMDDloglike
  DOUBLE PRECISION, DIMENSION(Nstates,0:T,Nsim) :: DRAWstates
  DOUBLE PRECISION, DIMENSION(Nyield,Nsim) :: DRAWavgtermpremia
  DOUBLE PRECISION, DIMENSION(Nf,Nsim) :: DRAWf
  DOUBLE PRECISION, DIMENSION(Nshockslopes,Nsim) :: DRAWshockslopes
  DOUBLE PRECISION, DIMENSION(1,Nsim) :: DRAWmaxlambda

  DOUBLE PRECISION, DIMENSION(Nsv,0:T,Nsim) :: DRAWsvol
  DOUBLE PRECISION, DIMENSION(Nsv,Nsim) :: DRAWhbar, DRAWhrho
  DOUBLE PRECISION, DIMENSION(NhSigma,Nsim) :: DRAWhSigma

  DOUBLE PRECISION, DIMENSION(Nsvorth,0:T,Nsim) :: DRAWsvolorth
  DOUBLE PRECISION, DIMENSION(Nsvorth,Nsim) :: DRAWhorthbar, DRAWhorthrho
  DOUBLE PRECISION, DIMENSION(Nsvorth,Nsim) :: DRAWhorthvar

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nw), Bsv(Nx,Nw,T), C(Ny,Nx,T), Ci(Nx) ! ,gap0variance(Ngap * p, Ngap * p)
  DOUBLE PRECISION :: invGapTransposeB(Ngap,Ngap), invGapTransposeBsv(Ngap,Ngap,T)

  DOUBLE PRECISION :: f(Nf), Ef0(Nf), sqrtVf0(Nf,Nf), iVf0(Nf,Nf)
  DOUBLE PRECISION :: shockslopes(Nshockslopes), Eshockslopes0(Nshockslopes), sqrtVshockslopes0(Nshockslopes,Nshockslopes), iVshockslopes0(Nshockslopes,Nshockslopes)
  INTEGER :: offsetslopes, these

  DOUBLE PRECISION :: maxlambda, SigmaStarT(Nx,Nx)
  DOUBLE PRECISION :: istar(0:T), rr(0:T), xhat(Nx)

  DOUBLE PRECISION :: gapshock(T,Ngap), gap(-(p-1):T,Ngap), gapVCV(Ngap,Ngap,T)

  DOUBLE PRECISION, DIMENSION (Nsv) :: Eh0, minSV, maxSV, hbar, hrho0, hrho
  DOUBLE PRECISION, DIMENSION(Nsv,Nsv) :: sqrtVh0, sqrtVhshock, hSigmaT, hrhoV0, hrhosqrtV0, iVhrho0
  DOUBLE PRECISION, DIMENSION(Nsv,Nsv) :: sqrthSigmaT
  DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: h, hgap, SVol
  DOUBLE PRECISION, DIMENSION(Nsv,T) :: SVinno, hshock
  DOUBLE PRECISION :: zdraw(Nsv)
  INTEGER :: hSigmaDof

  DOUBLE PRECISION, DIMENSION (Nsvorth) :: Ehorth0, minSVorth, maxSVorth, horthbar, horthrho0, horthrho
  DOUBLE PRECISION, DIMENSION (Nsvorth) :: Vhorth0(Nsvorth), horthvar, horthvarT, horthrhoV0
  DOUBLE PRECISION, DIMENSION(Nsvorth,0:T) :: horth, horthgap, SVolorth
  DOUBLE PRECISION, DIMENSION(Nsvorth,T) :: SVorthinno, horthshock
  INTEGER :: horthvarDof(Nsvorth)
  DOUBLE PRECISION :: horthvarLogGammaDof(Nsvorth)

  ! INTEGER :: longrateHorizon 


  DOUBLE PRECISION, DIMENSION(Nx, 0:T) :: x
  DOUBLE PRECISION, DIMENSION(Nx, T)   :: xshock 


  ! stack management
  DOUBLE PRECISION, PARAMETER :: unity = .999
  integer, parameter :: stacksize = 10, maxShakes = 10, zlbShakes = 100
  integer :: maxRejections 
  double precision, dimension(Nf,stacksize) :: PREV_F
  double precision, dimension(Nshockslopes,stacksize) :: PREV_SHOCKSLOPES

  double precision, dimension(Nsv,0:T,stacksize) :: PREV_SVOL
  double precision, dimension(Nsv,stacksize) :: PREV_HBAR, PREV_HRHO
  double precision, dimension(Nsv,Nsv,stacksize) :: PREV_HSIGMA

  double precision, dimension(Nsvorth,0:T,stacksize) :: PREV_SVOLORTH
  double precision, dimension(Nsvorth,stacksize) :: PREV_HORTHBAR, PREV_HORTHRHO
  double precision, dimension(Nsvorth,stacksize) :: PREV_HORTHVAR

  ! variables to store SV jump-off for density draws
  double precision, dimension(Nsv+Nsvorth) :: h_sim0, hbar_sim0, hrho_sim0
  double precision, dimension(Nsv+Nsvorth,Nsv+Nsvorth) :: sqrtVhshock_sim0

  integer, dimension(stacksize) :: stackRegister
  integer :: lastInStack, shakes, stackResetCount, rejectionCount, nominaltrendELBcount
  logical :: OK 
  CHARACTER (LEN=200) :: resetmsg

  ! MDD
  double precision :: hSigmalogGammaDof, hSigmaTlogdet
  ! double precision :: dummy
  double precision, dimension(Nsv + Nsvorth - 1) :: ABCSVhbar, ABCSVhrho
  double precision, dimension(Nsv + Nsvorth - 1, Nsv + Nsvorth - 1) :: ABCSVsqrtVh

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! state space bounds
  DOUBLE PRECISION, PARAMETER :: maxabsRbar = 1000.0d0

  ! CHARACTER (LEN=200) :: filename
  CHARACTER (LEN=100) :: progresscomment

  stackRegister = -1 ! dummy value for unregistered stack
  stackResetCount = -1

  Ngaplags = Ngap * p
  TotalSim = Burnin + Nsim ! only used for updating the progress bar
  maxRejections = TotalSim / 4 * maxShakes

  ! prepare state space

  ! y: ugap, inflation, fed funds, long-bond yields
  yndxUgap        = 1
  yndxInflation   = 2
  yndxPolicyRate  = 3
  yndxLONGINT     = 4


  ! x: pibar, rbar, gaps, ... gaps(lags), pinoise
  ndxPIbar        = 1
  ndxRbar         = 2
  ndxLONGRbar     = 3
  ndxPIgap        = ndxLONGRbar + Nyield
  ndxUgap         = ndxPIgap + 1
  ndxINTgap       = ndxUgap + 1
  ndxLONGINTgap   = ndxINTgap + 1

  ndxGapStart     = ndxPIgap
  ndxGapStop      = ndxLONGINTgap + (Nyield - 1)
  ndxGapLagStop   = (ndxPIgap - 1) + (Ngap * p)

  ndxPInoise      = ndxGapLagStop + 1

  if (ndxPInoise .ne. Nx) then
     print *, 'indices are off'
     print *, 'Nx:', Nx
     print *, 'ndxPInoise:', ndxPInoise
     stop 1
  end if

  ! ndxshock (order SV shocks first)
  ndxShockPIbar      = 1
  ndxShockGapStart   = ndxShockPibar + 1

  ndxShockPInoise    = ndxShockPIbar + Ngap + 1
  ndxShockRbar       = ndxShockPInoise + 1

  ! ordering gap shocks
  select case (gapOrderKey)
  case (1234)
     ! x
     ndxPIgap        = (ndxGapStart - 1) + 1
     ndxUgap         = (ndxGapStart - 1) + 2
     ndxINTgap       = (ndxGapStart - 1) + 3
     ndxLONGINTgap   = (ndxGapStart - 1) + 4
     ! shock     
     ndxShockPIgap      = (ndxShockGapStart - 1) + 1
     ndxShockUgap       = (ndxShockGapStart - 1) + 2
     ndxShockINTgap     = (ndxShockGapStart - 1) + 3
     ndxShockLONGINTgap = (ndxShockGapStart - 1) + 4

  case (2134)
     ! x
     ndxPIgap        = (ndxGapStart - 1) + 2
     ndxUgap         = (ndxGapStart - 1) + 1
     ndxINTgap       = (ndxGapStart - 1) + 3
     ndxLONGINTgap   = (ndxGapStart - 1) + 4
     ! shock     
     ndxShockPIgap      = (ndxShockGapStart - 1) + 2
     ndxShockUgap       = (ndxShockGapStart - 1) + 1
     ndxShockINTgap     = (ndxShockGapStart - 1) + 3
     ndxShockLONGINTgap = (ndxShockGapStart - 1) + 4

  case (3214)
     ! x
     ndxPIgap        = (ndxGapStart - 1) + 3
     ndxUgap         = (ndxGapStart - 1) + 2
     ndxINTgap       = (ndxGapStart - 1) + 1
     ndxLONGINTgap   = (ndxGapStart - 1) + 4
     ! shock     
     ndxShockPIgap      = (ndxShockGapStart - 1) + 3
     ndxShockUgap       = (ndxShockGapStart - 1) + 2
     ndxShockINTgap     = (ndxShockGapStart - 1) + 1
     ndxShockLONGINTgap = (ndxShockGapStart - 1) + 4

  case (4321)
     ! x
     ndxPIgap        = (ndxGapStart - 1) + 4 + (Nyield - 1)
     ndxUgap         = (ndxGapStart - 1) + 3 + (Nyield - 1)
     ndxINTgap       = (ndxGapStart - 1) + 2 + (Nyield - 1)
     ndxLONGINTgap   = (ndxGapStart - 1) + 1
     ! shock     
     ndxShockPIgap      = (ndxShockGapStart - 1) + 4 + (Nyield - 1)
     ndxShockUgap       = (ndxShockGapStart - 1) + 3 + (Nyield - 1)
     ndxShockINTgap     = (ndxShockGapStart - 1) + 2 + (Nyield - 1)
     ndxShockLONGINTgap = (ndxShockGapStart - 1) + 1 

  case (1324)
     ! x
     ndxPIgap        = (ndxGapStart - 1) + 1
     ndxUgap         = (ndxGapStart - 1) + 3
     ndxINTgap       = (ndxGapStart - 1) + 2
     ndxLONGINTgap   = (ndxGapStart - 1) + 4
     ! shock     
     ndxShockPIgap      = (ndxShockGapStart - 1) + 1
     ndxShockUgap       = (ndxShockGapStart - 1) + 3
     ndxShockINTgap     = (ndxShockGapStart - 1) + 2
     ndxShockLONGINTgap = (ndxShockGapStart - 1) + 4


  case default
     print *,'gapOrderKey= ', gapOrderKey, ' not recognized'
     stop 1
  end select
  ndxSVpibar        = ndxShockPIbar
  ndxSVpigap        = ndxShockPIgap
  ndxSVugap         = ndxShockUgap
  ndxSVintgap       = ndxShockINTgap
  ndxSVlongintgap   = ndxShockLONGINTgap + Nyield - 1 ! NOTE: needs to point at LAST yield
  ndxSVpinoise      = ndxShockPInoise
  ndxSVrbar         = ndxShockRbar

  A = 0.0d0
  ! unit roots for bar  
  FORALL(i=1:T,k=1:Nbar) A(k,k,i) = 1.0d0

  ! kompanion for gap
  FORALL (i = 1 : (Ngap * (p - 1)))
     A(ndxGapStop+i,ndxGapStart-1+i,:) = 1.0d0
  END FORALL

  ! prior over initial states
  Ex0           = 0.0d0
  call eye(sqrtVx0, 1.0d1)

  ! joint prior over pibar, and ibar; transformed into prior over pibar and rbar
  Ex0(ndxPIbar)     = 2.0d0 ! rough level for pibar
  Ex0(ndxRbar)      = 2.0d0 ! rough level for rbar
  forall (i = 0 : Nyield - 1) Ex0(ndxLONGRbar+i)  = 2.0d0 + dble(i+1) * 0.5d0 ! rough level for rbar
  ! sqrtVx0(ndxPIbar,ndxPIbar)     = 4.0d0
  ! sqrtVx0(ndxRbar,ndxPIbar)      = -  sqrtVx0(ndxPIbar,ndxPIbar)
  ! forall (i = 0 : Nyield - 1) sqrtVx0(ndxLONGRbar+i,ndxPIbar)  = -  sqrtVx0(ndxPIbar,ndxPIbar)

  Bsv = 0.0d0
  B   = 0.0d0
  ! unit loadings on SV shocks
  B(ndxPIbar,ndxshockPIBar)     = 1.0d0
  B(ndxPIgap,ndxshockPIgap)     = 1.0d0
  B(ndxRbar,ndxshockRBar)       = 1.0d0
  forall (i = 0 : Nyield - 1) B(ndxLONGRbar+i, ndxShockRbar) =  1.0d0

  B(ndxUgap,ndxshockUgap)               = 1.0d0
  B(ndxINTgap,ndxshockINTgap)           = 1.0d0
  forall (i = 0 : Nyield - 1) B(ndxLONGINTgap+i,ndxshockLONGINTgap+i)   = 1.0d0

  if (doPInoise) B(ndxPInoise,ndxshockPInoise) = 1.0d0

  ! Measurement
  C = 0.0d0

  ! u
  C(yndxUgap,ndxUgap,:)   = 1.0d0 

  ! pi
  C(yndxInflation,ndxPIbar,:)   = 1.0d0    
  C(yndxInflation,ndxPIgap,:)   = 1.0d0
  if (doPInoise) C(yndxInflation,ndxPInoise,:) = 1.0d0    

  ! i
  C(yndxPolicyrate,ndxRbar,:)   = 1.0d0    
  C(yndxPolicyrate,ndxPIbar,:)  = 1.0d0    
  C(yndxPolicyrate,ndxINTgap,:) = 1.0d0    

  ! longi
  forall (i = 0 : Nyield - 1) C(yndxLONGINT+i,ndxLONGRbar+i,:)   = 1.0d0    
  forall (i = 0 : Nyield - 1) C(yndxLONGINT+i,ndxPIbar,:)        = 1.0d0    
  forall (i = 0 : Nyield - 1) C(yndxLONGINT+i,ndxLONGINTgap+i,:) = 1.0d0    

  ! i for Ct
  Ci = 0.0d0
  Ci(ndxRbar)   = 1.0d0    
  Ci(ndxPIbar)  = 1.0d0    
  Ci(ndxINTgap) = 1.0d0    


  ! prediction state space
  ! - init
  AA = 0.0d0
  BB = 0.0d0
  CC = 0.0d0
  ! - copy original state space block
  AA(1:Nx,1:Nx) = A(:,:,T)
  CC(1:Ny,1:Nx) = C(:,:,T)
  BB(1:Nx,:)    = B
  ! - expand state transition
  AA(Nx+1,1:Nx) = C(yndxInflation,:,T)
  forall (j=2:3) AA(Nx+j,Nx+j-1) = 1.0d0
  ! - define 4q inflation 
  CC(Ny+1,1:Nx)      = C(yndxInflation,:,T) / 4.0d0
  CC(Ny+1,Nx+1:Nx+3) = 1.0d0 / 4.0d0
  ! CC(Ny+1,Nx+1) = 1.0d0

  ! - define shadow rate
  CC(Ny+2,1:Nx) = C(yndxPolicyRate,:,T)

  ! - plug data on lagged inflation into the xx state
  xx            = 0.0d0

  ytruncated    = .false.
  ytruncated(yndxPolicyrate) = .true.
  forall (i=0:Nyield-1) ytruncated(yndxLONGINT+i)    = .true.


  ytrunclb = ELBound

  ! prepare C for missing values
  DO k=1,T
     DO i = 1, Ny
        if (yNaN(i,k)) C(i,:,k) = 0.0d0
        if (yNaN(i,k) .AND. y(i,k) /= 0.0d0 ) then
           write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  if (doMDD) then
     sqrthSigmaT = hSigmaT
     call DPOTRF('L', Nsv, sqrthSigmaT, Nsv, status)
     if (status /= 0) then
        write(*,*) 'DPOTRF ERROR, INFO: ', status, ' [sqrthSigmaT]'
        stop 1
     end if

     hSigmaTlogdet = logdetSym(sqrthSigmaT)

  end if

  hrhosqrtV0 = hrhoV0
  call DPOTRF('L', Nsv, hrhosqrtV0, Nsv, status)
  if (status /= 0) then
     write(*,*) 'DPOTRF ERROR, INFO: ', status, ' [hrhosqrtV0]'
     stop 1
  end if

  iVhrho0 = hrhosqrtV0
  call DPOTRI('U', Nsv, iVhrho0, Nsv, status)
  if (status /= 0) then
     write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [iVhrho0]'
     stop 1
  end if

  iVf0 = sqrtVf0
  call DPOTRI('U', Nf, iVf0, Nf, status)
  if (status /= 0) then
     write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [iVf0]'
     stop 1
  end if

  iVshockslopes0 = sqrtVshockslopes0
  call DPOTRI('U', Nshockslopes, iVshockslopes0, Nshockslopes, status)
  if (status /= 0) then
     write(*,*) 'DPOTRI ERROR, INFO: ', status, ' [iVshockslopes0]'
     stop 1
  end if

  ! -------------------------------------------------------------------------
  ! MDD PREP
  ! -------------------------------------------------------------------------
  ! compute log(Gamma(dof)); needed as parameter for MDD computations later
  if (doMDD) then
     call vdlgamma(Nsvorth, dble(horthvarDof), horthvarlogGammaDof)
     hsigmaLogGammaDof = logmvgamma(hSigmaDof * 0.5d0, Nsv)
     ! print *, 'hSigmaLogGammaDof:', hSigmaLogGammaDof
  end if
  atELB(0)   = .false.
  atELB(1:T) = zlb(yndxpolicyrate,:)

  ! -------------------------------------------------------------------------
  ! END: MDD PREP
  ! -------------------------------------------------------------------------  

  shakes = 1
  j = 1
  lastInStack = 1 ! redundant, but avoids compiler warning
  resetmsg = 'DRAW 0'
  rejectionCount = 0
  thisdraw = 1

  DO WHILE (thisdraw <= Nsim)

     ! WRITE (*,*) 'STEP ', j

     IF (j == 1) THEN

        call initprogressbar(timer,timer%wait,timer%rank)
        stackRegister   = -1
        stackResetCount = stackResetCount + 1

        IF (rejectionCount > 1) THEN
           call hrulefill
           WRITE (*, '(" RE-INIT STREAM ", i2, " (reset count = ", i4, ", rejection count = ", i8, ", ", a10, ")")')  timer%rank, stackResetCount, rejectionCount, trim(resetmsg)
           call hrulefill
        END IF

        rejectionCount = 0
        nominaltrendELBcount = 0

        resetmsg = 'NIL'

        ! DRAW SV
        ! draw the Nsv-Ngap non-correlated SV 
        DO i = 1,Nsvorth
           ! Draw horthrho
           OK = .false.
           DO 
              errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, 1, horthrho(i), 0.0d0, 1.0d0) 
              horthrho(i) = horthrho0(i) + horthrho(i) * sqrt(horthrhoV0(i)) 
              OK = abs(horthrho(i)) < unity
              if (OK) exit
           END DO
           ! Draw the rest
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, 1, horthbar(i), 0.0d0, 1.0d0)
           horthbar(i) = Ehorth0(i) + sqrt(Vhorth0(i)) * horthbar(i)
           call igammaDraw(horthvar(i), horthvarT(i), horthvarDof(i), VSLstream)
           call drawAR1(horth(i,:),horthrho(i),T,sqrt(horthvar(i)),VSLstream)
           forall (k=0:T) horth(i,k) = horth(i,k) + horthbar(i)
           SVolorth(i,:) = exp(0.5d0 * horth(i,:))
           WHERE (SVolorth(i,:) .gt. maxSVorth(i)) SVolorth(i,:) = maxSVorth(i)
           WHERE (SVolorth(i,:) .lt. minSVorth(i)) SVolorth(i,:) = minSVorth(i)
        END DO
        ! TODO: reconsider rejection sampling ?

        if (.not. doPInoise) SVolorth(2, :) = minSV(2) ! todo: remove hard-coded 2
        if (.not. doPInoise) horthbar(2)    = 0.0d0


        ! draw Ngap correlated SV's
        ! Draw hrho
        OK = .false.
        hrho = 10.0d0
        DO  ! brute force: draw all jointly until all inside unit circle
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv, zdraw, 0.0d0, 1.0d0)
           hrho = hrho0
           call DGEMV('N',Nsv,Nsv,1.0d0,hrhosqrtV0,Nsv,zdraw,1,1.0d0,hrho,1)
           OK = all(abs(hrho) < unity)
           if (OK) exit
        END DO

        ! 1) hbar 
        errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv, zdraw, 0.0d0, 1.0d0)
        hbar = Eh0
        call DGEMV('N',Nsv,Nsv,1.0d0,sqrtVh0,Nsv,zdraw,1,1.0d0,hbar,1)
        ! hSigma
        call iwishcholDraw(sqrtVhshock, hSigmaT, hSigmaDof, Nsv, VSLstream)
        ! b) draw hgap
        CALL drawAR1correlated (h,Nsv,T,hrho,sqrtVhshock,VSLstream)
        ! c) add together
        forall (i=1:Nsv,k=0:T) h(i,k) = h(i,k) + hbar(i)
        SVol = exp(0.5d0 * h)
        ! TODO: reconsider rejection sampling ?
        DO i=1,Nsv
           WHERE (SVol(i,:) .gt. maxSV(i)) SVol(i,:) = maxSV(i)
           WHERE (SVol(i,:) .lt. minSV(i)) SVol(i,:) = minSV(i)
        END DO


        ! init VAR coefficients f
        DO
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nf, f, 0.0d0, 1.0d0)
           CALL DTRMV('U', 'T', 'N', Nf, sqrtVf0, Nf, f, 1)
           f = f + Ef0
           call VARmaxroot(maxlambda, f, Ngap, p)
           OK = maxlambda < unity
           IF (OK) EXIT
        END DO

        ! init shockslopes
        errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nshockslopes, shockslopes, 0.0d0, 1.0d0)
        CALL DTRMV('U', 'T', 'N', Nshockslopes, sqrtVshockslopes0, Nshockslopes, shockslopes, 1)
        shockslopes = shockslopes + Eshockslopes0

        ! init stack
        lastInStack = 1

        PREV_SVOL(:,:,lastInStack)      = SVol
        PREV_HBAR(:,lastInStack)        = hbar
        PREV_HSIGMA(:,:,lastInStack)    = sqrtVhshock
        PREV_HRHO(:,lastInStack)        = hrho

        PREV_SVOLORTH(:,:,lastInStack)      = SVolorth
        PREV_HORTHBAR(:,lastInStack)        = horthbar
        PREV_HORTHVAR(:,lastInStack)        = horthvar
        PREV_HORTHRHO(:,lastInStack)        = horthrho

        PREV_F(:,lastInStack)           = f
        PREV_SHOCKSLOPES(:,lastInStack) = shockslopes

        stackRegister(lastInStack)   = j

        ! WRITE (*,*) 'STREAM ', timer%rank, 'STEP ', j, '[ DONE INIT ]'
        ! WRITE (*,*) 'STREAM ', timer%rank, 'LAST RESET WAS ', resetmsg
        ! resetmsg = 'NIL'


     ELSE ! j > 1

        SVol          = PREV_SVOL(:,:,lastInStack) 
        hbar          = PREV_HBAR(:,lastInStack)
        sqrtVhshock   = PREV_HSIGMA(:,:,lastInStack) 
        hrho          = PREV_HRHO(:,lastInStack) 

        SVolorth      = PREV_SVOLORTH(:,:,lastInStack) 
        horthbar      = PREV_HORTHBAR(:,lastInStack)
        horthvar      = PREV_HORTHVAR(:,lastInStack) 
	horthrho      = PREV_HORTHRHO(:,lastInStack) 

        f             = PREV_F(:,lastInStack) 
        shockslopes   = PREV_SHOCKSLOPES(:,lastInStack) 

        if (j == BURNIN) call initprogressbar(timer,timer%wait,timer%rank)

     END IF ! j == 1

     ! ---------------------------------------------------------------------------
     ! SMOOTHING SAMPLER
     ! ---------------------------------------------------------------------------


     ! inflation gap loading inserted below 

     ! Fill in VAR coefficients
     FORALL (i=1:T) A(ndxGapStart:ndxGapStop,ndxGapStart:ndxGapLagStop,i) = transpose(reshape(f, (/ Ngap * p, Ngap /)))


     ! Fill in shockslopes
     offsetslopes = 0
     DO i = 2, Nsv 
        these = i-1 
        ! slopes in row i have index offsetslopes+1:offsetslopes + these
        ! note: assumes that gap shocks start at the 2nd column of B
        B(Nbar+i,2:1+these) = shockslopes(offsetslopes+1:offsetslopes+these)
        offsetslopes = offsetslopes + these
     END DO
     invGapTransposeB = 0.0d0
     invGapTransposeB = transpose(B(Nbar+1:Nbar+Ngap,1+1:1+Ngap))
     ! call savemat(invGapTransposeB, 'GapTransposeB.debug')
     call DTRTRI('u', 'u', Ngap, invGapTransposeB, Ngap, errcode )
     if (errcode /= 0) then
        write (*,*) 'DTRTRI error (invGapTransposeB)', errcode
        stop 1
     end if
     ! call savemat(invGapTransposeB, 'invGapTransposeB.debug')
     FORALL (i=1:T,k=1:Ngap) invGapTransposeBsv(:,k,i) = invGapTransposeB(:,k) / SVol(k,i) 

     ! multiply with SV loadings
     Bsv    = 0.0d0
     k = ndxShockPIbar
     FORALL (i=1:T) Bsv(:,k,i)    = B(:,k) * SVolorth(1,i)
     k = ndxShockPInoise
     FORALL (i=1:T) Bsv(:,k,i)    = B(:,k) * SVolorth(2,i)
     k = ndxShockRbar
     FORALL (i=1:T) Bsv(:,k,i)    = B(:,k) * SVolorth(3,i)

     FORALL (i=1:T,k=1:Nsv) Bsv(:,1+k,i)   = B(:,1+k) * SVol(k,i) ! note the offset to account for PIbar in the first column

     ! there are no nonSV shocks anymore; hence the next line can be commented out
     ! FORALL (i=1:T, k=ndxShockPInoise+1:Nw) Bsv(:,k,i) = B(:,k) 

     ! ! Fill in unconditional variance of stationary states
     ! CALL DLYAP(gap0variance, A(ndxGapStart:ndxGapLagStop,ndxGapStart:ndxGapLagStop,1), Bsv(ndxGapStart:ndxGapLagStop,:,1), Ngaplags, Nw, errcode) 
     ! if (errcode /= 0) then
     !    write (*,*) 'DLYAP error (gap0variance)', errcode
     !    call savemat(Bsv(ndxGapStart:ndxGapLagStop,:,1), 'B.debug')
     !    call savemat(A(ndxGapStart:ndxGapLagStop,ndxGapStart:ndxGapLagStop,1), 'A.debug')
     !    stop 1
     ! end if
     ! ! Factorize the unconditional variance
     ! CALL DPOTRF('L', Ngaplags, gap0variance, Ngaplags, errcode)
     ! if (errcode /= 0) then
     !    write (*,*) 'DPOTRF error (gap0variance)', errcode
     !    call savemat(Bsv(ndxGapStart:ndxGapLagStop,:,1), 'B.debug')
     !    call savemat(A(ndxGapStart:ndxGapLagStop,ndxGapStart:ndxGapLagStop,1), 'A.debug')
     !    stop 1
     ! end if

     ! ! zero out the upper triangular
     ! FORALL (i=2:Ngaplags) gap0variance(1:i-1,i) = 0.0d0
     ! ! fill it in
     ! sqrtVx0(ndxGapStart:ndxGapLagStop,ndxGapStart:ndxGapLagStop) = gap0variance

     OK = .FALSE.
     i = 0
     DO WHILE (.NOT. OK .AND. (i < zlbShakes))

        i = i + 1
        ! print *, 'entering abc' 
        call samplerA3B3C3nanscalar(x,xshock,SigmaStarT,y,yNaN,T,Ny,Nx,Nw,A,Bsv,C,Ex0,sqrtVx0,VSLstream,errcode)
        ! print *, 'exiting abc'

        ! errcode = 1
        if (errcode /= 0) then

           OK = .FALSE.
           print *, 'oops', i

           call savemat(A(:,:,1), 'A.debug')
           call savemat(C(:,:,1), 'C.debug')
           call savemat(A(:,:,2), 'A2.debug')
           call savemat(B(:,:), 'B.debug')
           call savemat(Bsv(:,:,1), 'B1.debug')
           call savemat(Bsv(:,:,2), 'B2.debug')
           call savemat(C(:,:,2), 'C2.debug')
           call savevec(Ex0, 'Ex0.debug')
           call savemat(sqrtVx0, 'sqrtVx0.debug')
           call savemat(x, 'x.debug')
           call savemat(xshock, 'xshock.debug')
           call savemat(y, 'y.debug')
           ! call savemat(gap0variance, 'gapvariance')
           call savemat(SVol, 'SVol.debug')
           call savemat(SVolorth, 'SVolorth.debug')
           call savevec(shockslopes, 'shockslopes.debug')

           ! print *, 'ndxUgap', ndxugap
           stop 1 

        else
           OK = .true.
        end if

        ! if (.not. OK) exit
        ! if (.not. OK) print *,'what am i doing here?'

        ! construct real rate
        rr = 0.0d0
        IF (j > BURNIN) THEN

           ! real rate
           xhat  = predictstate(Nx, A(:,:,1), x(:,0), 1)
           rr(0) = x(ndxRbar,0) + x(ndxINTgap,0) - xhat(ndxPIgap)
           DO k=1,T
              xhat  = predictstate(Nx, A(:,:,k), x(:,k), 1)
              if (zlb(yndxPolicyrate,k)) then
                 rr(k) = ELBound - x(ndxPIbar,k) - xhat(ndxPIgap) 
              else
                 rr(k) = x(ndxRbar,k) + x(ndxINTgap,k) - xhat(ndxPIgap)
              end if
           END DO
        END IF


        istar = x(ndxRbar,:) + x(ndxPIbar,:) + x(ndxINTgap,:) 

        ! implement lower bound on nominal rate trend
        OK = .true.
        ! if (ZLBrejection) then ! TODO: check this condition
        !    OK = ALL(x(ndxRbar,:) + x(ndxPIbar,:) .ge. ELBound) 
        ! end if

        OK = OK .AND. ALL(abs(x(ndxRbar,:)) .lt. maxabsRbar)

        if (.NOT. OK) then
           resetmsg = 'RBAR'
           nominaltrendELBcount = nominaltrendELBcount + 1
           ! WRITE (*,*) 'STREAM ', timer%rank, 'hitting ', resetmsg
        else if (ZLBrejection) then 
           OK = ALL((.NOT. zlb(yndxPolicyrate,:)) .OR. (istar(1:T) < ELBound))
           IF (.NOT. OK) resetmsg = 'ZLB'
        end if

     END DO


     IF (OK) THEN

        ! ---------------------------------------------------------------------------
        ! GAP VAR
        ! ---------------------------------------------------------------------------

        ! gapVCV
        gapVCV = 0.0d0
        DO i = 1,T

           call DSYRK('u','n',Ngap,Ngap,1.0d0,invGapTransposeBsv(:,:,i),Ngap,0.0d0,gapVCV(:,:,i),Ngap)

        END DO

        ! collect gap vector
        FORALL (i=1:p) gap(1-i,:) = x(ndxGapStart+(i-1)*Ngap:ndxGapStop+(i-1)*Ngap,0)
        gap(1:T,:)                = transpose(x(ndxGapStart:ndxGapStop,1:T))

        ! d: draw VAR coefficients (recall: gapVCV contains inverse of gapVCV)
        call bayesVARSV(f, gapshock, p, gap, Ngap, T, gapVCV, Ef0, iVf0, VSLstream)
        ! call bayesdiffuseVARSV(f, gapshock, p, gap, Ngap, T, gapVCV, VSLstream)

        if (any(isnan(f))) then
           write(*,*) 'NaN draws from bayesVAR'
           stop 1
        end if

        ! d) check stability
        call VARmaxroot(maxlambda, f, Ngap, p)

        OK = maxlambda < unity
        IF (.NOT. OK) resetmsg = 'VAR'

        xshock(ndxGapStart:ndxGapStop,:) = transpose(gapshock)

        ! print *, 'done VAR'

     END IF


     IF (OK) THEN

        ! ---------------------------------------------------------------------------
        ! Orthogonal SV shocks
        ! ---------------------------------------------------------------------------

        SVorthinno(1,:) = xshock(ndxShockPIbar,:)
        SVorthinno(2,:) = xshock(ndxShockPInoise,:)
        SVorthinno(3,:) = xshock(ndxShockRbar,:)

        ! call hrulefill
        ! print *, j
        ! print *, Ehorth0, Vhorth0, horthrho, horthbar, horthvar
        ! call hrulefill

        i = 1
        CALL stochvolKSCar1plus(SVolorth(i,:), horth(i,:), horthbar(i), horthshock(i,:), SVorthinno(i,:), T, sqrt(horthvar(i)), horthrho(i), Ehorth0(i), Vhorth0(i), VSLstream)
        i = 3
        CALL stochvolKSCar1plus(SVolorth(i,:), horth(i,:), horthbar(i), horthshock(i,:), SVorthinno(i,:), T, sqrt(horthvar(i)), horthrho(i), Ehorth0(i), Vhorth0(i), VSLstream)
        i = 2
        if (doPInoise) then
           CALL stochvolKSCar1plus(SVolorth(i,:), horth(i,:), horthbar(i), horthshock(i,:), SVorthinno(i,:), T, sqrt(horthvar(i)), horthrho(i), Ehorth0(i), Vhorth0(i), VSLstream)
        else
           SVolorth(i,:) = 1.0d0
           horth(i,:) = 0.0d0
           horthbar(i) = 0.0d0
           horthshock(i,:) = 0.0d0
        end if

        do i=1,Nsvorth
           OK = OK .AND. ALL(SVolorth(i,:) >= minSVorth(i)) .AND. ALL(SVolorth(i,:) <= maxSVorth(i))
        end do
        IF (.NOT. OK) resetmsg = 'SVOLORTH'
        ! IF (.NOT. OK) then
        !    print *, 'RESET: ', j, resetmsg
        !    call savemat(SVolorth, 'SVolorth.debug')
        !    call savemat(horth, 'horth.debug')
        !    call savemat(SVorthinno, 'SVorthinno.debug')
        !    call savemat(xshock, 'xshock.debug')
        !    print *, Ehorth0, Vhorth0, horthrho, horthbar, horthvar
        !    stop 11
        ! end if


     END IF

     IF (OK) THEN
        ! ---------------------------------------------------------------------------
        ! Gap SV block
        ! ---------------------------------------------------------------------------

        SVinno = xshock(ndxGapStart:ndxGapStop,:)

        call SVHcholKSCAR1cor(T, Nsv, SVol, h, hshock, hbar, hrho, shockslopes, Nshockslopes, SVinno, Eshockslopes0, iVshockslopes0, sqrtVhshock, Eh0, sqrtVh0, VSLstream)

        do i=1,Nsv
           OK = OK .AND. ALL(SVol(i,:) >= minSV(i)) .AND. ALL(SVol(i,:) <= maxSV(i))
        end do
        IF (.NOT. OK) resetmsg = 'SVOL'

     END IF


     ! draw AR(1)-SUR system
     IF (OK) THEN
        ! note: demeaning h here; do not use it afterwards!
        forall (i=1:Nsv,k=0:T) hgap(i,k) = h(i,k) - hbar(i)
        hshock = hgap(:,1:T)
        call bayesAR1SUR(Nsv,T,hrho,hshock,hgap(:,0:T-1),sqrtVhshock,hrho0,iVhrho0,VSLstream)
        OK = ALL(ABS(hrho) < unity)
        ! print *, hrho
        IF (.NOT. OK) resetmsg = 'HRHO'
     END IF


     ! draw AR(1)-SV parameters horthrho
     ! individual draws, since assumed to be independent
     IF (OK) THEN
        ! note: demeaning h here; do not use it afterwards!
        forall (i=1:Nsvorth,k=0:T) horthgap(i,k) = horth(i,k) - horthbar(i)
        horthshock = horthgap(:,1:T)
        do i=1,Nsvorth
           call bayesUnivariateRegressionSlope(horthrho(i), horthshock(i,:), horthgap(i,0:T-1), T, 1 / horthvar(i), horthrho0(i), 1 / horthrhoV0(i), VSLstream)
        end do
        if (.not. doPInoise) horthrho(2) = 0.0d0
        OK = ALL(ABS(horthrho) < unity)
        IF (.NOT. OK) resetmsg = 'HORTHRHO'
     END IF


     IF (OK) THEN ! move forward

        call vcvcholDrawTR(sqrtVhshock, hSigmaT, hSigmaDof, transpose(hshock), T, Nsv, VSLstream)
        ! if (j == 100) then
        !    call savemat(sqrtVhshock, 'sqrtVhshock.debug')
        !    stop 11
        ! end if

        i = 1
        call varianceDraw(horthvar(i), horthvarT(i), horthvarDof(i), horthshock(i,:), T, VSLstream)
        i = 3
        call varianceDraw(horthvar(i), horthvarT(i), horthvarDof(i), horthshock(i,:), T, VSLstream)

        i = 2
        if (doPInoise) then
           call varianceDraw(horthvar(i), horthvarT(i), horthvarDof(i), horthshock(i,:), T, VSLstream)
        else
           horthvar(i) = 0.0d0
        end if

        ! ---------------------------------------------------------------------------
        ! PREPARE NEXT ITERATION
        ! ---------------------------------------------------------------------------

        ! a) store Draws after Burnin
        IF (j > BURNIN) THEN

           !  WRITE (*,*) 'STORING THISDRAW ', thisdraw

           DRAWstates(1,:,thisdraw)         = rr        
           DRAWstates(2,:,thisdraw)         = x(ndxPIbar,:)
           DRAWstates(3,:,thisdraw)         = x(ndxRbar,:)
           forall (i=1:Nyield) DRAWstates(3+i,:,thisdraw)         = x(ndxLONGRbar-1+i,:)
           DRAWstates(3+Nyield+1,:,thisdraw)         = x(ndxPIgap,:)
           DRAWstates(3+Nyield+2,:,thisdraw)         = x(ndxUgap,:)
           DRAWstates(3+Nyield+3,:,thisdraw)         = x(ndxINTgap,:)
           forall (i = 0 : Nyield - 1) DRAWstates(3+Nyield+4+i,:,thisdraw)  = x(ndxLONGINTgap+i,:)
           DRAWstates(6+2*Nyield+1,:,thisdraw)  = x(ndxPInoise,:)
           DRAWstates(6+2*Nyield+2,:,thisdraw) = istar

           DRAWsvol(:,:,thisdraw)  = SVol
           DRAWhbar(:,thisdraw)    = hbar
	   DRAWhrho(:,thisdraw)    = hrho
           call vech(DRAWhsigma(:,thisdraw), sqrtVhshock)

           DRAWsvolorth(:,:,thisdraw)  = SVolorth
           DRAWhorthbar(:,thisdraw)    = horthbar
	   DRAWhorthrho(:,thisdraw)    = horthrho
           DRAWhorthvar(:,thisdraw)    = horthvar

           DRAWf(:,thisdraw)           = f
           DRAWshockslopes(:,thisdraw) = shockslopes
           DRAWmaxlambda(:,thisdraw)   = maxlambda

           forall (i=1:Nyield) DRAWavgtermpremia(i,thisdraw) = x(ndxLONGRbar+i-1,0) - x(ndxRbar,0) 

           if (doPDF) then
              ! predictive density at end of sample

              ! update BB with rbarvar draw (note: X is on top of expanded state vector, ndx* variables still work
              BB(1:Nx,:)     = B
              AA(1:Nx,1:Nx)  = A(:,:,T)

              ! update (extended) state vector
              xx(1:Nx)                = x(:,T)
              forall (i=1:3) xx(Nx+i) = x(ndxPIbar,T-i) + x(ndxPIgap,T-i) + x(ndxPInoise,T-i) 

              ! SV objects are evaluated at PREV_ values, just to be synchronous with what went into the ABC and xx matrices
              ! a) h_sim0
              h_sim0 = 0.0d0
              h_sim0(ndxSVpibar)                 = PREV_SVOLORTH(1,T,lastInstack)
              h_sim0(ndxSVpinoise)               = PREV_SVOLORTH(2,T,lastInstack)
              h_sim0(ndxSVrbar)                  = PREV_SVOLORTH(3,T,lastInstack)
              h_sim0(ndxSVpigap:ndxSVlongintgap) = PREV_SVOL(:,T,lastInstack)

              h_sim0 = 2.0d0 * log(h_sim0)
              ! b) hbar_sim0
              hbar_sim0 = 0.0d0
              hbar_sim0(ndxSVpibar)                 = PREV_horthbar(1,lastInstack)
              hbar_sim0(ndxSVpinoise)               = PREV_horthbar(2,lastInstack)
              hbar_sim0(ndxSVrbar)               = PREV_horthbar(3,lastInstack)
              hbar_sim0(ndxSVpigap:ndxSVlongintgap) = PREV_hbar(:,lastInstack)
              ! c) hrho_sim0
              hrho_sim0(ndxSVpibar)                 = PREV_horthrho(1,lastInstack)
              hrho_sim0(ndxSVpinoise)               = PREV_horthrho(2,lastInstack)
              hrho_sim0(ndxSVrbar)                  = PREV_horthrho(3,lastInstack)
              hrho_sim0(ndxSVpigap:ndxSVlongintgap) = PREV_hrho(:,lastInstack)
              ! d) sqrtVhshock_sim0
              sqrtVhshock_sim0 = 0.0d0
              sqrtVhshock_sim0(ndxSVpibar,ndxSVpibar)     = sqrt(PREV_horthvar(1,lastInstack))
              sqrtVhshock_sim0(ndxSVpinoise,ndxSVpinoise) = sqrt(PREV_horthvar(2,lastInstack))
              sqrtVhshock_sim0(ndxSVrbar,ndxSVrbar)       = sqrt(PREV_horthvar(3,lastInstack))
              sqrtVhshock_sim0(ndxSVpigap:ndxSVlongintgap,ndxSVpigap:ndxSVlongintgap) = PREV_HSIGMA(:,:,lastInstack)



              call predictiveDensitySVCORdraws(DRAWy(:,:,:,thisdraw), Nforecastdraws, maxhorizons, NNyy, NNxx, Nw, xx, AA, BB, CC, Nsv, h_sim0, hbar_sim0, hrho_sim0, sqrtVhshock_sim0, VSLstream)

              ! truncation 
              forall (i=1:maxhorizons,k=1:Nforecastdraws)
                 where (ytruncated .and. DRAWy(:,i,k,thisdraw) < ytrunclb)  DRAWy(:,i,k,thisdraw) = ytrunclb
              end forall
           end if ! do PDF


           if (doMDD) then

              ! --------------------------------------------------------------
              ! evaluate parameters at prior
              ! --------------------------------------------------------------

              DRAWMDDlogprior(thisdraw) = 0.0d0

              ! general note on the normal priors:
              ! many are independent normals, used block-wise
              ! but no harm in computing the priors jointly (numerically identical)

              ! f
              DRAWMDDlogprior(thisdraw) = MVnormlogpdf(Nf, PREV_f(:,lastInStack), Ef0, sqrtVf0) - fPriorTruncationProbLog + DRAWMDDlogprior(thisdraw)

              ! shockslopes
              DRAWMDDlogprior(thisdraw) = MVnormlogpdf(Nshockslopes, PREV_shockslopes(:,lastInStack), Eshockslopes0, sqrtVshockslopes0) + DRAWMDDlogprior(thisdraw)

              ! hrho (recall: hrhoSqrtV0 is left-lower-triangular factor)
              DRAWMDDlogprior(thisdraw) = MVnormlogpdf(Nsv, PREV_hrho(:,lastInStack), hrho0, hrhosqrtV0, .true.) + DRAWMDDlogprior(thisdraw) - hrhoPriorTruncationProbLog

              ! hbar
              DRAWMDDlogprior(thisdraw) = MVnormlogpdf(Nsv, PREV_hbar(:,lastInStack), Eh0, sqrtVh0) + DRAWMDDlogprior(thisdraw)

              ! hSigma
              DRAWMDDlogprior(thisdraw) = iwishlogpdf(Nsv, PREV_hSigma(:,:,lastInStack), sqrthSigmaT, hSigmaDof, hSigmaTlogdet, hsigmaLogGammaDof) + DRAWMDDlogprior(thisdraw)


              ! horth stuff
              i = 1
              ! hortrho
              DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthrho(i,lastInStack), horthrho0(i), horthrhoV0(i)) + DRAWMDDlogprior(thisdraw)
              ! horthbar
              DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthbar(i,lastInStack), Ehorth0(i),  Vhorth0(i)) + DRAWMDDlogprior(thisdraw)
              ! horthvar
              DRAWMDDlogprior(thisdraw) = igammalogpdf(PREV_horthvar(i,lastInStack),  horthvarT(i), horthvarDof(i), horthvarLogGammaDof(i)) + DRAWMDDlogprior(thisdraw) 

              if (doPInoise) then
                 i = 2
                 DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthbar(i,lastInStack), Ehorth0(i),  Vhorth0(i)) + DRAWMDDlogprior(thisdraw)
                 DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthrho(i,lastInStack), horthrho0(i), horthrhoV0(i)) + DRAWMDDlogprior(thisdraw)
                 DRAWMDDlogprior(thisdraw) = igammalogpdf(PREV_horthvar(i,lastInStack),  horthvarT(i), horthvarDof(i), horthvarLogGammaDof(i)) + DRAWMDDlogprior(thisdraw) 
              end if

              i = 3
              ! hortrho
              DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthrho(i,lastInStack), horthrho0(i), horthrhoV0(i)) + DRAWMDDlogprior(thisdraw)
              ! horthbar
              DRAWMDDlogprior(thisdraw) = normlogpdf(PREV_horthbar(i,lastInStack), Ehorth0(i),  Vhorth0(i)) + DRAWMDDlogprior(thisdraw)
              ! horthvar
              DRAWMDDlogprior(thisdraw) = igammalogpdf(PREV_horthvar(i,lastInStack),  horthvarT(i), horthvarDof(i), horthvarLogGammaDof(i)) + DRAWMDDlogprior(thisdraw) 

              ! --------------------------------------------------------------
              ! evaluate likelihood with PF
              ! --------------------------------------------------------------

              DRAWMDDloglike(thisdraw) = 0.0d0

              ! todo:
              ! - collect h stuff
              ! - pass prior Ex0 Eh0 etc to particlelike
              ! - pass state space matrices (not collection of coefficients) to particlelike

              ! collect SV parameters
              ! note: this structure mirros quite a bit what has been done above with the _sim0 Variables; however, the latter include also PInoise stuff which is not supporte dby particlelike
              ! todo: clean this up!
              ! todo: pass Ex0 etc
              ABCSVhrho = 0.0d0
              ABCSVhrho(1)         = PREV_horthrho(1,lastInstack)
              ABCSVhrho(2)         = PREV_horthrho(3,lastInstack)
              ABCSVhrho(2+1:2+Nsv) = PREV_hrho(:,lastInstack)

              ABCSVhbar = 0.0d0
              ABCSVhbar(1)         = PREV_horthbar(1,lastInstack)
              ABCSVhbar(2)         = PREV_horthbar(3,lastInstack)
              ABCSVhbar(2+1:2+Nsv) = PREV_hbar(:,lastInstack)

              ABCSVsqrtVh = 0.0d0
              ABCSVsqrtVh(1,1) = sqrt(PREV_horthvar(1,lastInstack))
              ABCSVsqrtVh(2,2) = sqrt(PREV_horthvar(3,lastInstack))
              ABCSVsqrtVh(2+1:2+Nsv,2+1:2+Nsv) = PREV_hSigma(:,:,lastInStack)

              if (doPInoise) then
                 print *, 'MDD not yet supported when PInoise true'
                 stop 1
              else
                 call particlelike(DRAWMDDloglike(thisdraw), ELBound, Nyield, T, Ny, y, yNaN, atELB, Nparticles, Nx-1, Nw-1, Nsv+Nsvorth-1, p, PREV_f(:,lastInStack), Nf, PREV_shockslopes(:,lastInStack), Nshockslopes, ABCSVhbar, ABCSVhrho, ABCSVsqrtVh, VSLstream)
              end if

           end if

           joflastdraw = j
           thisdraw = thisdraw + 1

        END IF ! J > BURNIN

        ! b) update stack
        lastInStack = minloc(stackRegister,1) ! find a free element on the stack

        stackRegister(lastInStack) = j + 1 ! notice the +1

        PREV_SVOL(:,:,lastInStack)      = SVol
        PREV_HBAR(:,lastInStack)        = hbar
        PREV_HSIGMA(:,:,lastInStack)    = sqrtVhshock
	PREV_HRHO(:,lastInStack)        = hrho

        PREV_SVOLORTH(:,:,lastInStack)      = SVolorth
        PREV_HORTHBAR(:,lastInStack)        = horthbar
        PREV_HORTHVAR(:,lastInStack)        = horthvar
        PREV_HORTHRHO(:,lastInStack)        = horthrho

        PREV_F(:,lastInStack)           = f
        PREV_SHOCKSLOPES(:,lastInStack) = shockslopes

        ! WRITE (*,'("Stream ", i2, " at step ", i5, " updated lastInStack to ", i5)') timer%rank, j, lastInStack


        j = j + 1



     ELSE ! NOT(OK): shake or move back

        OK = .TRUE. ! a little superfluous, but good for the sake of clarity

        rejectionCount = rejectionCount + 1

        IF (rejectionCount > maxRejections) THEN
           j = 1
        ELSE IF (shakes < maxshakes) THEN
           shakes = shakes + 1  ! WRITE (*,'("Stream ", i2, ": shaking at step ", i5)') timer%rank, j
        ELSE ! move back
           shakes = 1
           IF (j > 1) THEN
              ! WRITE (*,'("Stream ", i2, ": moving back at step ", i5, " b/o ", a30)') timer%rank, j, resetmsg
              stackRegister(lastInStack) = -1
              lastInStack = maxloc(stackRegister,1)
              ! WRITE(*,*) "stackRegister=", stackRegister(lastInStack), 'j=', j

              if ((thisdraw > 1) .and. (j == joflastdraw)) then
                 WRITE (*,'("Stream ", i2, ": moving back at thisdraw ", i5)') timer%rank, thisdraw
                 thisdraw = thisdraw - 1
              end if

              j = j - 1
              if (stackRegister(lastInStack) < 1) then
                 j = 1
              else
                 if (stackRegister(lastInStack) .ne. j) then
                    print *, 'stack problem', ' j=', j, ' stackregister=', stackregister(lastinstack)
                 end if
              end if
           END IF ! j > 1

        END IF ! shakes 

     END IF ! OK


     ! Report Progress
     progresscomment = 'RC = ' // trim(int2str(rejectioncount))
     if (j > BURNIN) then
        sharedone = dble(j - Burnin) / dble(Nsim)
     else
        progresscomment = 'Burnin, ' // trim(progresscomment)
        sharedone = dble(j) / dble(BURNIN)
     end if
     call progressbarcomment(sharedone, timer, trim(progresscomment))

  END DO ! while thisdraw < Nsim

  call hrulefill
  WRITE (*, '(" STREAM ", i2, " IS DONE [reset count = ", i4, ", rejection count = ", i8, ", nominalrate trend ELB rejections = ", i8, ", Ndraws =  ", i8, "]")')  timer%rank, stackResetCount, rejectionCount, nominaltrendELBcount, Nsim + Burnin
  call hrulefill


END SUBROUTINE thissampler

! @\newpage\subsection{particlelike}@
SUBROUTINE particlelike(loglike, elbound, Nyield, T, Ny, y, yNaN, atELB, Nparticles, Nx, Nw, Nsv, p, f, Nf, gapshockslopes, Ngapshockslopes, hbar, hrho, sqrtVh, VSLstream)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery
  use densitybox, only : simySVCORdraws, drawXelb

  use vslbox
  ! use timerbox

  IMPLICIT NONE


  logical, parameter :: doSecondResamplingStep = .false.

  double precision, intent(in) :: elbound

  INTENT(IN) :: T,Ny,y,yNaN, atELB, Nx, Nw, Nsv, Nparticles
  INTENT(IN) ::  p, f, Nf, gapshockslopes, Ngapshockslopes, hbar, hrho, sqrtVh



  INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw

  integer, dimension(Ny) :: ordervec
  double precision, intent(out) :: loglike

  ! double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  LOGICAL, DIMENSION(0:T) :: atELB

  DOUBLE PRECISION, DIMENSION(Nparticles)     :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: prevPARTICLEweights
  DOUBLE PRECISION, DIMENSION(T)              :: loglikecontrib

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), llf(Nparticles)
  DOUBLE PRECISION :: xdraws(Nx,Nparticles)

  DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
  DOUBLE PRECISION :: kernelsum, loglikemax

  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), logdetSigmaY

  DOUBLE PRECISION :: h(Nsv, Nparticles), SVol(Nsv,Nparticles), hshock(Nsv,Nparticles)

  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny) 
  DOUBLE PRECISION :: Kgain(Nx,Ny) 
  DOUBLE PRECISION, DIMENSION(Nx * (Nx + 1) / 2, Nparticles) :: vecSqrtSigmaX

  ! qr
  DOUBLE PRECISION :: qrR(Ny+Nx+Nw,Ny+Nx)
  INTEGER :: qrLwork

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), ytilde(Ny), Cobserver(Ny,Nx)

  ! DOUBLE PRECISION :: gapB0(Ny*p,Ny)

  ! SV parameters
  DOUBLE PRECISION, DIMENSION(Nsv) :: hbar, hrho, hintercept
  DOUBLE PRECISION, DIMENSION(Nsv,Nsv) :: sqrtVh ! , Ah, sqrtVh0
  DOUBLE PRECISION :: minSVh(Nsv)

  ! Other parameters
  integer :: Nf
  DOUBLE PRECISION :: f(Nf)
  integer :: Ngapshockslopes
  DOUBLE PRECISION :: gapshockslopes(Ngapshockslopes)


  INTEGER :: Nynonan

  ! CHARACTER (LEN=200) :: filename

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  type (vsl_stream_state), intent(inout) :: VSLstream
  double precision :: uniformdraws(2,T) ! two draws: one for each APF step

  ! index variables for state space
  INTEGER, INTENT(IN) :: Nyield
  INTEGER :: Ngap
  INTEGER :: ndxPIbar, ndxRbar, ndxGapStart, ndxGapStop, ndxGapCompanionStop, ndxUgap, ndxPIgap, ndxINTgap
  INTEGER :: shockndxPIbar, shockndxGapStart, shockndxGapStop, shockndxRBAR 
  double precision, dimension(Nx) :: Cshadowrate, zvec
  double precision, dimension(Nparticles) :: critvalues, ELBlike

  integer :: offset, these


  ordervec = (/ ( i, i=1,Ny) /)
  NsigmaX = Nx * (Nx + 1) / 2


  Ngap         = Ny
  ndxPIbar     = 1
  ndxRbar      = 2
  ndxGapStart  = ndxRbar + Nyield + 1
  ndxGapStop   = ndxGapStart + Ny - 1
  ndxGapCompanionStop = ndxGapStart - 1 + Ngap * p
  if (ndxGapCompanionStop /= Nx) then
     print *, 'particlelike: state space indices are off'
     stop 1
  end if
  ndxPIgap  = ndxGapStart
  ndxUgap   = ndxGapStart + 1
  ndxINTgap = ndxGapStart + 2 

  shockndxPIbar    = 1
  shockndxrbar     = 2
  shockndxGapStart = 3
  shockndxGapStop  = shockndxGapStart - 1 + Ngap


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
  B(ndxRbar,shockndxRbar)     = 1.0d0
  forall (j=1:Nyield) B(ndxRbar+j,shockndxRbar) = 1.0d0

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



  ! first set up Cobserver (2D),then reorder it, then copy over into C (3D)

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
  ! Cobserver -- reorder row as desired
  Cobserver = Cobserver(ordervec,:)

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
  sqrtVx0         = transpose(sqrtVx0) ! needs to be upper triangular RIGHT factor henceforth; redundant as long as prior is spherical
  ! gapB0 = B(ndxGapStart:ndxGapCompanionStop,shockndxGapStart:shockndxGapStop)
  ! forall(i=1:Ngap) gapB0(:,i) = gapB0(:,i) * exp(hbar(i) * 0.5d0) ! todo: could choose a fancier value based on the mean, not the median, of the log-normal 
  ! call dlyapsqrt(sqrtVx0(ndxGapStart:ndxGapCompanionStop,ndxGapStart:ndxGapCompanionStop), A(ndxGapStart:ndxGapCompanionStop,ndxGapStart:ndxGapCompanionStop), gapB0, Ngap*p, Ngap, errcode)
  ! if (errcode /= 0) then
  !    print *, 'PF-LIE: error computing ergodic VCV of PF state space'
  !    call savemat(A(ndxGapStart:ndxGapCompanionStop,ndxGapStart:ndxGapCompanionStop), 'A.debug')
  !    call savemat(A, 'AA.debug')
  !    call savemat(gapB0, 'B.debug')
  !    stop 1
  ! end if
  vecSqrtSigmaX   = 0.0d0

  DO k=1,Nparticles
     vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)
  END DO


  PARTICLEweights = 1.0d0 / dble(Nparticles)

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR = 0.0d0
  qrlwork = qrquery(qrR)

  ! CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     prevPARTICLEweights = PARTICLEweights

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! ELB adjustment for last period's linear states (if needed)
     ! ------------------------------------------------------------------------------------------------------------------------------
     if (atELB(j-1)) then

        ! print *, 'ELB Xdraw at ', j, ' ...'
        call drawXelb(Nx, Nparticles, xdraws, elbound, xposterior, vecSqrtSigmaX, Cshadowrate, VSLstream)
        ! print *, 'done'

        xposterior    = xdraws
        vecSqrtSigmaX = 0.0d0
     end if



     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: APF RESAMPLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     Nynonan = count(.not. yNaN(:,j))

     ! 1) Construct proposal for SVol at median (i.e. forecast of AR1 log-SV)
     forall (i=1:Nsv,k=1:Nparticles) SVol(i,k) = exp((hintercept(i) + hrho(i) * h(i,k)) * 0.5d0)
     ! OR: stick with previous SV
     ! forall (i=1:Nsv,k=1:Nparticles) SVol(i,k) = exp(h(i,k) * 0.5d0)

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
           zvec    = Cshadowrate
           call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,zvec,1)
           critvalues(k) = (elbound - sum(Cshadowrate * xposterior(:,k))) / sqrt(sum(zvec ** 2))
           ! critvalues will be converted into pdf after particle loop
        end if

     END DO ! k particles
     ! --- $OMP END PARALLEL DO 

     ! compute ELBwlike (if at ELB),
     ! note: it is more efficient to call vdcdf over vector of particles rather than individually for each particle across OMP loops
     if (atELB(j)) then
        call vdcdfnorm(Nparticles,critvalues,ELBlike)
     end if

     ! todo: handle elblike also when Nynonan == 0
     if (Nynonan > 0) then

        ! Reweight particles for next round

        loglikeAPFmax      = maxval(llf)
        llf                = llf - loglikeAPFmax
        if (atELB(j)) then
           APFlike = exp(llf) * ELBlike          
        else
           APFlike = exp(llf)
        end if
        if (doSecondResamplingStep) then
           APFkernelweights     = APFlike / dble(Nparticles) 
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights      = APFkernelweights / APFkernelsum
        else
           APFkernelweights     = APFlike * prevPARTICLEweights
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights      = APFkernelweights / APFkernelsum
        end if

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights, Nparticles, uniformdraws(1,j))

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

     else ! i.e. Nynonan == 0
        ! do nothing since no information received
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
        ! DRAWytilde(:,k,j) = ytilde

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

        if (atELB(j)) then
           zvec    = Cshadowrate
           call DTRMV('u','n','n',Nx,sqrtSigmaX,Nx,zvec,1)
           critvalues(k) = (elbound - sum(Cshadowrate * xposterior(:,k))) / sqrt(sum(zvec ** 2))
           ! critvalues will be converted into pdf after particle loop
        end if

        vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)


        ! ------------------------------------------------------------------------
        ! DONE: SQRT KALMAN
        ! ------------------------------------------------------------------------


     END DO ! k particles



     ! compute ELBlike (if at ELB),
     ! note: brief testing indicates that it is more efficient to call vdcdf over vector of particles rather than individually for each particle across OMP loops
     if (atELB(j)) then
        call vdcdfnorm(Nparticles,critvalues,ELBlike)
        llf = llf + log(ELBlike)
     end if


     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! log-likelihood
     loglikemax        = maxval(llf)
     llf               = llf - loglikemax
     kernelweights     = exp(llf) / APFlike 
     kernelsum         = sum(kernelweights)
     loglikecontrib(j) = log(kernelsum / Nparticles) + log(APFkernelsum) + loglikemax ! note: loglikeAPFmax correction should be added to APFlike and subtracted from kernelweights (b/o APFlike in denominator), thus cancelling each other


     ! Store NON-reweighted statistics
     if (Nynonan > 0) then

        ! particles weights
        PARTICLEweights = kernelweights / kernelsum


        ! resample
        if (doSecondResamplingStep) then
           call drawNDXsysresample(ndx, Nparticles, PARTICLEweights, Nparticles, uniformdraws(2,j))


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

        end if ! doSecondResamplingStep


     end if ! Nynonan > 0

  END DO ! j=1,T

  loglike = sum(loglikecontrib)

END SUBROUTINE particlelike


! @\newpage\subsection{MDDthetaTruncationProb}@
SUBROUTINE thisMDDthetaTruncationProb(prob, Ny, p, ndxf, Nsv, ndxrho, Nsim, Ndraws, Ntheta, theta, VSLstream) 

  use vslbox
  use embox, only : savevec, savemat
  use blaspack, only : sqrtvcvTR
  use densitybox, only :  MVNrnd
  use gibbsbox, only :  maxrootVAR

  IMPLICIT NONE

  integer, intent(in) :: Ndraws, Ntheta, Nsim
  integer, intent(in) :: Ny, p, Nsv
  integer, dimension(Ny * p * Ny), intent(in) :: ndxf
  integer, dimension(Nsv), intent(in) :: ndxrho
  double precision, dimension(Ntheta, Ndraws), intent(in) :: theta
  double precision, intent(out) :: prob
  logical, dimension(Nsim) :: acceptdraw

  double precision :: thetabar(Ntheta), sqrtVtheta(Ntheta,Ntheta), thetadraws(Ntheta,Nsim)

  integer :: i

  ! VSL
  type (vsl_stream_state), intent(inout) :: VSLstream

  ! compute linear moments
  thetabar     = sum(theta,2) / dble(Ndraws)
  sqrtVtheta   = sqrtvcvTR(Ntheta,Ndraws,theta)

  ! sample
  thetadraws = MVNrnd(Nsim, Ntheta, thetabar, sqrtVtheta, VSLstream)

  ! evaluate rejection
  acceptdraw = .true.

  ! check that VAR roots are all stable

  !$OMP PARALLEL DO SHARED(thetadraws,ndxf,Ny,p,acceptdraw,Nsim) DEFAULT(NONE) SCHEDULE(STATIC)
  do i=1,Nsim
     acceptdraw(i) = maxrootVAR(thetadraws(ndxf,i), Ny, p) .lt. 1.0d0
  end do
  !$OMP END PARALLEL DO

  ! call savemat(thetadraws(ndxf,:), 'fdraws.debug')
  ! call savevec(dble(acceptdraw), 'accepts.debug')
  ! stop 2
  ! check that SV-rho are all stable
  ! do i=1,Nsim
  !    ! acceptdraw(i) = acceptdraw(i) .AND. all(abs(thetadraws(ndxrho,i)) .lt. 1.0d0)
  ! end do
  ! call savevec(dble(acceptdraw), 'accepts.debug')

  ! check that SV-rho are all stable
  acceptdraw = acceptdraw .AND. all(abs(thetadraws(ndxrho,:)) .lt. 1.0d0, 1)
  ! call savevec(dble(all(abs(thetadraws(ndxrho,:)) .lt. 1.0d0, 1)), 'accepts.debug')
  ! stop 2
  ! call savevec(dble(acceptdraw), 'accepts2.debug')

  ! call savemat(thetadraws(ndxrho,:), 'rhodraws.debug')
  ! call savemat(theta(ndxrho,:), 'rho.debug')

  prob = log(dble(count(acceptdraw))) - log(dble(Nsim))
  print *, 'acceptance prob', exp(prob)

END SUBROUTINE thisMDDthetaTruncationProb

