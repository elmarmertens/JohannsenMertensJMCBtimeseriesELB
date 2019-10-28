MODULE densitybox

  USE embox, only : sqrttwopi, logpi, logtwopi, mean, variance, normpdf, savemat, savevec, savearray3
  USE blaspack, only: sandwich, qrecon, sqrtvcvTR, qrot, qrquery, ivechu !, sandwich2
  USE statespacebox, only : simyABCsvol0, simyABC
  USE timerbox, only : tidomp
  USE vslbox


  IMPLICIT NONE

CONTAINS

  ! @\newpage\subsection{meanTruncNormal}@
  subroutine meanvarianceTruncNormalUB(N,mubar,volbar,mu,vol,ub) 
    ! mean of normal subject to a upper bound

    integer, intent(in) :: N
    double precision, dimension(N), intent(out) :: mubar, volbar
    double precision, dimension(N), intent(in) :: mu, vol
    double precision, intent(in)   :: ub
    double precision, dimension(N) :: beta, phi, probbelowUB

    beta    = (ub - mu) / vol
    call vdcdfnorm(N, beta, probbelowUB)

    ! phi     = normpdf(beta)
    ! mubar   = mu - phi / probbelowUB * vol
    ! volbar  = vol * sqrt(1 - phi / probbelowUB * (beta + phi / probbelowUB))

    phi     = normpdf(beta) / probbelowUB
    mubar   = mu - phi * vol
    volbar  = vol * sqrt(1 - phi * (beta + phi))

    ! call savevec(mu, 'mu.debug')
    ! call savevec(mubar, 'mubar.debug')
    ! call savevec(vol, 'vol.debug')
    ! call savevec(volbar, 'volbar.debug')
    ! call savevec(beta, 'beta.debug')
    ! call savevec(phi * probbelowUB, 'phipdf.debug')
    ! call savevec(probbelowUB, 'phicdf.debug')
    ! stop 12

  end subroutine meanvarianceTruncNormalUB

  ! @\newpage\subsection{meanTruncNormal}@
  function meanTruncNormal(N,mu,vol,lb,ub) result(mubar)
    ! mean of normal subject to a upper bound

    integer, intent(in) :: N
    double precision, dimension(N), intent(in) :: mu, vol
    double precision, intent(in)   :: lb,ub
    double precision, dimension(N) :: mubar
    double precision, dimension(N) :: alpha, beta, phialpha, phibeta, probbelowLB, probbelowUB

    alpha    = (lb - mu) / vol
    call vdcdfnorm(N, alpha, probbelowLB)
    phialpha = normpdf(alpha)

    beta    = (ub - mu) / vol
    call vdcdfnorm(N, beta, probbelowUB)
    phibeta = normpdf(beta)

    mubar   = mu + (phialpha - phibeta) / (probbelowUB - probbelowLB) * vol

  end function meanTruncNormal

  ! @\newpage\subsection{meanTruncNormalUB}@
  function meanTruncNormalUB(N,mu,vol,ub) result(mubar)
    ! mean of normal subject to a upper bound

    integer, intent(in) :: N
    double precision, dimension(N), intent(in) :: mu, vol
    double precision, intent(in)   :: ub
    double precision, dimension(N) :: mubar
    double precision, dimension(N) :: beta, phi, probbelowUB

    beta    = (ub - mu) / vol
    call vdcdfnorm(N, beta, probbelowUB)
    phi     = normpdf(beta)
    mubar   = mu - phi / probbelowUB * vol

    ! call savevec(mu, 'mu.debug')
    ! call savevec(mubar, 'mubar.debug')
    ! call savevec(vol, 'vol.debug')
    ! call savevec(beta, 'beta.debug')
    ! call savevec(phi, 'phipdf.debug')
    ! call savevec(probbelowUB, 'phicdf.debug')
    ! stop 12

  end function meanTruncNormalUB

  ! @\newpage\subsection{meanTruncNormalLB}@
  function meanTruncNormalLB(N,mu,vol,lb) result(mubar)
    ! mean of normal subject to a upper bound

    integer, intent(in) :: N
    double precision, dimension(N), intent(in) :: mu, vol
    double precision, intent(in)   :: lb
    double precision, dimension(N) :: mubar
    double precision, dimension(N) :: alpha, phi, probbelowLB

    alpha    = (lb - mu) / vol
    call vdcdfnorm(N, alpha, probbelowLB)
    phi     = normpdf(alpha)
    mubar   = mu + phi / (1.0d0 - probbelowLB) * vol

  end function meanTruncNormalLB

  ! @\newpage\subsection{predictiveDensity}@
  SUBROUTINE predictiveDensity(ydraws, Ndraws, Nhorizons, Ny, Nx, Nw, x0, A, B, C, VSLstream)
    
    ! simulate predictive density when all shocks are constant variance
    ! x0 is given

    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ydraws
    INTENT(IN)    :: Nhorizons, Ny, Nx, Nw, x0, A, B, C

    INTEGER :: Ndraws
    INTEGER :: Nhorizons
    INTEGER :: Ny, Nx, Nw

    double precision :: x0(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)

    double precision, dimension(Ny,Nhorizons,Ndraws) :: ydraws


    INTEGER :: ii

    type (vsl_stream_state) :: VSLstream
    ! INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ydraws           = 0.0d0
    do ii=1,Ndraws
       ydraws(:,:,ii) = simyABC(Nhorizons,Ny,Nx,Nw,A,B,C,x0,VSLstream) 
    end do


  END SUBROUTINE predictiveDensity

  ! @\newpage\subsection{predictiveDensitySVdraws}@
  SUBROUTINE predictiveDensitySVdraws(ydraws, Ndraws, Nhorizons, Ny, Nx, Nw, x0, A, B, C, Nsv, h0, hbar, hrho, hinno, VSLstream)

    ! simulate predictive density when first Nsv shcoks are SV (rest is constant variance)
    ! x0 is given

    ! same as predictiveDensitySVdrawsOLD; except using a different ordering of the loops; seems to be a little swifter

    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ydraws
    INTENT(IN)    :: Nhorizons, Ny, Nx, Nw, x0, A, B, C, Nsv, h0, hbar, hrho, hinno

    INTEGER :: Ndraws
    INTEGER :: Nhorizons
    INTEGER :: Ny, Nx, Nw, Nsv

    double precision :: x0(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
    double precision, dimension(Nsv) :: h0, hbar, hrho, hinno
    DOUBLE PRECISION, DIMENSION(Nsv,Nhorizons,Ndraws) :: h, hshock
    DOUBLE PRECISION, DIMENSION(Nw,Nhorizons) :: ShockScale

    double precision, dimension(Ny,Nhorizons,Ndraws) :: ydraws


    INTEGER :: hh, ii, jj

    type (vsl_stream_state) :: VSLstream
    integer :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ydraws           = 0.0d0

    ! SIMULATE SVs
    ! draw random numbers
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nhorizons * Ndraws, hshock, 0.0d0, 1.0d0)
    ! compute transition from jump off (demeaned)
    DO ii=1,Ndraws
       hh = 1
       FORALL (jj=1:Nsv) h(jj,hh,ii)    = hrho(jj) * h0(jj)        + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii)
       DO hh = 2,Nhorizons
          FORALL (jj=1:Nsv) h(jj,hh,ii) = hrho(jj) * h(jj,hh-1,ii) + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii) 
       END DO
    END DO

    ShockScale = 1.0d0 ! note: Nw can be greater than Nsv

    ! simulate state space; TODO: do this inline and optimized?
    do ii=1,Ndraws
       ShockScale(1:Nsv,:)  = exp(0.5d0 * h(:,:,ii))
       ydraws(:,:,ii) = simyABCsvol0(Nhorizons,Ny,Nx,Nw,A,B,C,ShockScale,x0,VSLstream) 
    end do

    ! ShockScale(1:Nsv,:) = 0.5d0
    ! call savemat(ShockScale, 'ShockScale.debug')

    ! print *, 'Nsv', Nsv
    ! print *, 'Nw', Nw

    ! call savevec(x0, 'x0.debug')
    ! call savemat(A, 'A.debug')
    ! call savemat(C, 'C.debug')
    ! call savemat(B, 'B.debug')
    ! call savemat(ydraws(1,:,:), 'y1.debug')
    ! call savemat(ydraws(2,:,:), 'y2.debug')
    ! call savemat(ydraws(3,:,:), 'y3.debug')
    ! call savemat(ydraws(4,:,:), 'y4.debug')
    ! call savemat(ydraws(5,:,:), 'y5.debug')
    ! call savemat(ydraws(6,:,:), 'y6.debug')
    ! stop 11


  END SUBROUTINE predictiveDensitySVdraws

  ! @\newpage\subsection{predictiveDensitySVCORdraws}@
  SUBROUTINE predictiveDensitySVCORdraws(ydraws, Ndraws, Nhorizons, Ny, Nx, Nw, x0, A, B, C, Nsv, h0, hbar, hrho, sqrtVhshock, VSLstream)


    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ydraws
    INTENT(IN)    :: Nhorizons, Ny, Nx, Nw, x0, A, B, C, Nsv, h0, hbar, hrho, sqrtVhshock

    INTEGER :: Ndraws
    INTEGER :: Nhorizons
    INTEGER :: Ny, Nx, Nw, Nsv

    double precision :: x0(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
    double precision, dimension(Nsv) :: h0, hbar, hrho, sqrtVhshock
    DOUBLE PRECISION, DIMENSION(Nsv,Nhorizons,Ndraws) :: h, hshock
    DOUBLE PRECISION, DIMENSION(Nw,Nhorizons) :: ShockScale

    double precision, dimension(Ny,Nhorizons,Ndraws) :: ydraws


    INTEGER :: hh, ii , jj

    type (vsl_stream_state) :: VSLstream
    integer :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ydraws           = 0.0d0

    ! SIMULATE SVs
    ! draw random numbers
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nhorizons * Ndraws, hshock, 0.0d0, 1.0d0)
    ! compute transition from jump off (demeaned)
    DO ii=1,Ndraws
       hh = 1
       h(:,hh,ii)    = hrho * (h0 - hbar) 
       call DGEMV('n',Nsv,Nsv,1.0d0,sqrtVhshock,Nsv,hshock(:,hh,ii),1,1.0d0,h(:,hh,ii),1)
       DO hh = 2,Nhorizons
          h(:,hh,ii)    = hrho * h(:,hh-1,ii)
          call DGEMV('n',Nsv,Nsv,1.0d0,sqrtVhshock,Nsv,hshock(:,hh,ii),1,1.0d0,h(:,hh,ii),1)
       END DO
    END DO
    ! add means to h
    forall (jj=1:Nsv,hh=1:Nhorizons,ii=1:Ndraws) h(jj,hh,ii) = h(jj,hh,ii) + hbar(jj)

    ShockScale = 1.0d0 ! note: Nw can be greater than Nsv

    ! simulate state space; TODO: do this inline and optimized?
    do ii=1,Ndraws
       ShockScale(1:Nsv,:)  = exp(0.5d0 * h(:,:,ii))
       ydraws(:,:,ii) = simyABCsvol0(Nhorizons,Ny,Nx,Nw,A,B,C,ShockScale,x0,VSLstream) 
    end do


  END SUBROUTINE predictiveDensitySVCORdraws

  ! @\newpage\subsection{simySVCORdraws}@
  SUBROUTINE simySVCORdraws(ydraws, Ndraws, Nhorizons, Ny, Nx, Nw, Ex0, sqrtVx0, A, B, C, Nsv, h0, hbar, hrho, sqrtVhshock, VSLstream)

    ! this routine is a slight variant of predictiveDensitySVCORdraws:
    ! - instead of jumping off a gixed x0, x0 is drawn from an initial distribution, characterized by Ex0, sqrtVx0
    ! - sqrtVx0 likely comes from a qr-posterior of the Kalman filter, assumed to be right-upper triangular!!!
    ! - (sqrtVhshock continues to be left-choleski factor!)

    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ydraws
    INTENT(IN)    :: Nhorizons, Ny, Nx, Nw, Ex0, sqrtVx0, A, B, C, Nsv, h0, hbar, hrho, sqrtVhshock

    INTEGER :: Ndraws
    INTEGER :: Nhorizons
    INTEGER :: Ny, Nx, Nw, Nsv

    double precision :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
    double precision, dimension(Nsv) :: h0, hbar, hrho, sqrtVhshock
    double precision, dimension(Nsv,Nhorizons,Ndraws) :: h, hshock
    double precision, dimension(Nw,Nhorizons) :: ShockScale

    double precision, dimension(Nx,Ndraws)  :: x0, draw0

    double precision, dimension(Ny,Nhorizons,Ndraws) :: ydraws


    INTEGER :: hh, ii , jj

    type (vsl_stream_state) :: VSLstream
    integer :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ydraws           = 0.0d0

    ! SIMULATE SVs
    ! draw random numbers
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nhorizons * Ndraws, hshock, 0.0d0, 1.0d0)
    ! compute transition from jump off (demeaned)
    DO ii=1,Ndraws
       hh = 1
       h(:,hh,ii)    = hrho * (h0 - hbar) 
       call DGEMV('n',Nsv,Nsv,1.0d0,sqrtVhshock,Nsv,hshock(:,hh,ii),1,1.0d0,h(:,hh,ii),1)
       DO hh = 2,Nhorizons
          h(:,hh,ii)    = hrho * h(:,hh-1,ii)
          call DGEMV('n',Nsv,Nsv,1.0d0,sqrtVhshock,Nsv,hshock(:,hh,ii),1,1.0d0,h(:,hh,ii),1)
          ! todo: replace DGEMV with DTRMV ?
       END DO
    END DO
    ! add means to h
    forall (jj=1:Nsv,hh=1:Nhorizons,ii=1:Ndraws) h(jj,hh,ii) = h(jj,hh,ii) + hbar(jj)

    ShockScale = 1.0d0 ! note: Nw can be greater than Nsv

    ! draw x0
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Ndraws, draw0, 0.0d0, 1.0d0)
    ! option 1: DGEMV
    ! forall (ii=1:Ndraws) x0(:,ii) = Ex0
    ! DO ii=1,Ndraws
    !    ! todo check transpose
    !    ! todo: call dgemm
    !    call DGEMV('t',Nx,Nx,1.0d0,sqrtVx0,Nx,draw0(:,ii),1,1.0d0,x0(:,ii),1)
    ! END DO
    ! call savemat(x0, 'x0mv.debug')

    ! option 2: DGEMM
    forall (ii=1:Ndraws) x0(:,ii) = Ex0
    call DGEMM('t','n',Nx,Ndraws,Nx,1.0d0,sqrtVx0,Nx,draw0,Nx,1.0d0,x0,Nx)
    ! todo: replace DGEMM with DTRMM ?

    ! call savemat(x0, 'x0mm.debug')
    ! stop 11

    ! simulate state space; TODO: do this inline and optimized?
    do ii=1,Ndraws
       ShockScale(1:Nsv,:)  = exp(0.5d0 * h(:,:,ii))
       ydraws(:,:,ii) = simyABCsvol0(Nhorizons,Ny,Nx,Nw,A,B,C,ShockScale,x0(:,ii),VSLstream) 
    end do

    ! ShockScale(1:Nsv,:) = 0.5d0
    ! call savemat(ShockScale, 'ShockScale.debug')

    ! print *, 'Nsv', Nsv
    ! print *, 'Nw', Nw

    ! call savevec(x0, 'x0.debug')
    ! call savemat(A, 'A.debug')
    ! call savemat(C, 'C.debug')
    ! call savemat(B, 'B.debug')
    ! call savemat(ydraws(1,:,:), 'y1.debug')
    ! call savemat(ydraws(2,:,:), 'y2.debug')
    ! call savemat(ydraws(3,:,:), 'y3.debug')
    ! call savemat(ydraws(4,:,:), 'y4.debug')
    ! call savemat(ydraws(5,:,:), 'y5.debug')
    ! call savemat(ydraws(6,:,:), 'y6.debug')
    ! stop 11


  END SUBROUTINE simySVCORdraws

  ! @\newpage\subsection{predictiveDensitySV}@
  SUBROUTINE predictiveDensitySV(ypdf, yforecast, ycondvar, Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno, VSLstream)


    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ypdf, yforecast, ycondvar
    INTENT(IN)    :: Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno

    INTEGER, PARAMETER :: NsimSV = 100
    INTEGER :: Nhorizons, horizons(Nhorizons), maxhorizons
    INTEGER :: Ny, Nx, Nw, Nsv

    double precision :: x(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
    double precision, dimension(Nsv) :: h0, hbar, hrho, hinno

    double precision, dimension(Ny, Nhorizons) :: ypdf, yforecast, ycondvar
    double precision, dimension(Ny, maxhorizons) :: ypred
    logical, dimension(Ny, maxhorizons) :: ynanpred

    DOUBLE PRECISION, dimension(Ny,maxhorizons) :: forecastMean  
    DOUBLE PRECISION, dimension(NsimSV,Ny,maxhorizons) :: forecastCondVar

    INTEGER :: hh, ii, jj, kk, nn

    ! variables for generating mean forecasts
    double precision, dimension(Nx) :: xlag, xhat


    ! helper
    DOUBLE PRECISION, DIMENSION(Nsv,maxhorizons,NsimSV) :: h, hshock
    DOUBLE PRECISION :: SigmaY(Ny,Ny,maxhorizons), VMA(Nx,Nw,0:maxhorizons-1)

    type (vsl_stream_state) :: VSLstream
    integer :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ypdf             = 0.0d0
    yforecast        = 0.0d0
    ycondvar         = 0.0d0

    forecastMean     = 0.0d0
    forecastCondVar  = 0.0d0

    ! FORECAST MEAN
    xhat=x ! necessary copy since x is intent(in)

    DO hh = 1, maxhorizons

       xlag = xhat
       CALL DGEMV('N', Nx, Nx, 1.0d0, A,  Nx, xlag, 1, 0.0d0, xhat, 1)
       CALL DGEMV('N', Ny, Nx, 1.0d0, C,  Ny, xhat, 1, 0.0d0, forecastMean(:,hh), 1)

    END DO ! hh 


    ! TODO: rewrite a separate RW-SV version (code below works also for rho=1)

    ! SIMULATE SVs (needed for variance computations) 
    ! draw random numbers
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * maxhorizons * NsimSV, hshock, 0.0d0, 1.0d0)
    ! compute transition from jump off (demeaned)
    hh = 1
    FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii)    = hrho(jj) * h0(jj)        + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii)
    DO hh = 2,maxhorizons
       FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii) = hrho(jj) * h(jj,hh-1,ii) + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii) ! for these loops, hh should ideally be in third dimension, but next set of loops works better with ii in last dimension
    END DO

    VMA        = 0.0d0
    VMA(:,:,0) = B
    DO hh = 1,maxhorizons-1
       call dgemm('N', 'N', Nx, Nw, Nx, 1.0d0, A, Nx, VMA(:,:,hh-1), Nx, 0.0d0, VMA(:,:,hh), Nx)
    END DO


    ! FORECAST VARIANCES (unconstrained) FOR EACH SIMULATED SV
    DO ii =1,NsimSV
       ! SigmaY = condvarABCSV(A,B,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
       SigmaY = condvarVMACSV(VMA,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
       forall (nn=1:Ny,hh=1:maxhorizons) forecastCondVar(ii,nn,hh) = SigmaY(nn,nn,hh) ! not quite efficient to have ii as the first-column argument in *this* loop, but helps with mean() further below
    END DO


    ! call savemat(Sigma, 'Sigma.debug')
    ! call savemat(SigmaY, 'SigmaY.debug')
    ! call savemat(B, 'B.debug')
    ! call savemat(A, 'A.debug')
    ! call savemat(C, 'C.debug')
    ! call savemat(forecastcondvar(1,:,:), 'condvar.debug')
    ! call savearray3(VMA, 'VMA', 'debug')
    ! call savemat(h(:,:,1), 'h.debug')
    ! stop 1

    ! yforecast
    forall (nn=1:Ny,kk=1:Nhorizons) yforecast(nn,kk) = forecastMean(nn,horizons(kk))
    ! ycondvar
    forall (nn=1:Ny,kk=1:Nhorizons) ycondvar(nn,kk)  = mean(forecastCondVar(:,nn,horizons(kk))) ! Note: there is no variance across conditional expectations here

    ! compute llf 
    do kk=1,Nhorizons
       hh = horizons(kk)
       do nn = 1, Ny ! could replace this with a where constuct, but might make code less obvious
          ! pdf (if predicted obesrvations available)
          if (.not. ynanpred(nn,hh)) then
             ! normal pdf
             ypdf(nn,kk) = mean(exp(-0.5d0 * (((ypred(nn,hh) - forecastMean(nn,hh)) ** 2) / forecastCondVar(:,nn,hh))) /  (sqrttwopi * sqrt(forecastCondVar(:,nn,hh))))
          end if

       end do ! nn 
    end do ! kk


  END SUBROUTINE predictiveDensitySV

  ! @\newpage\subsection{predictiveDensitySVtruncLB}@
  SUBROUTINE predictiveDensitySVtruncLB(ypdf, yforecast, ymedian, yprobLB, ycondvar, Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, ytruncated, ytrunclb, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno, VSLstream)

    ! returns a vector of *univariate* llfs for realized values yplus at horizon 
    ! note: forecasts of truncated variables assume that lagged *truncated* variables do *not* influence the forecast (TODO)

    ! note: ymedian is niot actually median, but rather mean of conditional medians!!

    INTENT(INOUT) :: VSLstream
    INTENT(OUT)   :: ypdf, yforecast, ymedian, yproblb, ycondvar
    INTENT(IN)    :: Nhorizons, horizons, maxhorizons, Ny, ypred, ynanpred, ytruncated, ytrunclb, Nx, Nw, x, A, B, C, Nsv, h0, hbar, hrho, hinno

    INTEGER, PARAMETER :: NsimSV = 100
    INTEGER :: Nhorizons, horizons(Nhorizons), maxhorizons
    INTEGER :: Ny, Nx, Nw, Nsv

    double precision :: x(Nx), A(Nx,Nx), B(Nx,Nw), C(Ny,Nx)
    double precision, dimension(Nsv) :: h0, hbar, hrho, hinno

    double precision, dimension(Ny, Nhorizons) :: ypdf, yforecast, ymedian, yproblb, ycondvar
    double precision, dimension(Ny, maxhorizons) :: ypred
    logical, dimension(Ny, maxhorizons) :: ynanpred
    double precision, dimension(Ny) :: ytrunclb
    logical, dimension(Ny) :: ytruncated

    DOUBLE PRECISION, dimension(Ny,maxhorizons) :: forecastMean  
    DOUBLE PRECISION, dimension(NsimSV,Ny,maxhorizons) :: forecastCondVar

    INTEGER :: hh, ii, jj, kk, nn

    ! variables for generating mean forecasts
    double precision, dimension(Nx) :: xlag, xhat


    ! helper
    DOUBLE PRECISION, DIMENSION(Nsv,maxhorizons,NsimSV) :: h, hshock
    DOUBLE PRECISION :: SigmaY(Ny,Ny,maxhorizons), VMA(Nx,Nw,0:maxhorizons-1)
    DOUBLE PRECISION, DIMENSION(NsimSV) ::  truncAlpha, PROBaboveLB, PROBatLB, truncPDFalpha, truncMean, censMean, truncVar, censVar, sqrtCondVar, simpdf ! , debug
    DOUBLE PRECISION, PARAMETER :: alphaThreshold = 6.0d0




    type (vsl_stream_state) :: VSLstream
    integer :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    ! init
    ypdf             = 0.0d0
    yforecast        = 0.0d0
    ymedian          = 0.0d0
    yproblb          = 0.0d0
    ycondvar         = 0.0d0

    forecastMean     = 0.0d0
    forecastCondVar  = 0.0d0

    ! FORECAST MEAN
    xhat=x ! necessary copy since x is intent(in)

    DO hh = 1, maxhorizons

       xlag = xhat
       CALL DGEMV('N', Nx, Nx, 1.0d0, A,  Nx, xlag, 1, 0.0d0, xhat, 1)
       CALL DGEMV('N', Ny, Nx, 1.0d0, C,  Ny, xhat, 1, 0.0d0, forecastMean(:,hh), 1)

    END DO ! hh 

    ! call savevec(x, 'x.debug')
    ! call savemat(forecastMean, 'y.debug')
    ! call savemat(A, 'A.debug')
    ! call savemat(C, 'C.debug')
    ! call savemat(B, 'B.debug')
    ! stop 11

    ! SIMULATE SVs (needed for variance computations) 
    ! draw random numbers
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * maxhorizons * NsimSV, hshock, 0.0d0, 1.0d0)
    ! compute transition from jump off (demeaned)
    hh = 1
    FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii)    = hrho(jj) * h0(jj)        + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii)
    DO hh = 2,maxhorizons
       FORALL (jj=1:Nsv,ii=1:NsimSV) h(jj,hh,ii) = hrho(jj) * h(jj,hh-1,ii) + (1 - hrho(jj)) * hbar(jj) + hinno(jj) * hshock(jj,hh,ii) ! for these loops, hh should ideally be in third dimension, but next set of loops works better with ii in last dimension
    END DO


    VMA        = 0.0d0
    VMA(:,:,0) = B
    DO hh = 1,maxhorizons-1
       call dgemm('N', 'N', Nx, Nw, Nx, 1.0d0, A, Nx, VMA(:,:,hh-1), Nx, 0.0d0, VMA(:,:,hh), Nx)
    END DO

    ! FORECAST VARIANCES (unconstrained) FOR EACH SIMULATED SV
    DO ii =1,NsimSV
       ! SigmaY = condvarABCSV(A,B,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
       SigmaY = condvarVMACSV(VMA,C,h(:,:,ii),Ny,Nx,Nw,Nsv,maxhorizons) 
       forall (nn=1:Ny,hh=1:maxhorizons) forecastCondVar(ii,nn,hh) = SigmaY(nn,nn,hh) ! not quite efficient to have ii as the first-column argument in *this* loop, but helps with mean() further below
    END DO


    ! compute llf and truncated moments (if needed) etc
    do kk=1,Nhorizons
       hh = horizons(kk)
       do nn = 1, Ny ! could replace this with a where constuct, but might make cde less obvious

          ! init dummy values
          truncMean  = 0.0d0
          truncVar   = 1.0d0

          ! store the cond vol (for ease of reference)
          sqrtCondVar      = sqrt(forecastCondVar(:,nn,hh))

          if (ytruncated(nn)) then
             ! moments of truncated normal


             ! compute avg-median
             if (forecastMean(nn,hh) > ytrunclb(nn)) then
                ymedian(nn,kk) = forecastMean(nn,hh)
             else
                ymedian(nn,kk) = ytrunclb(nn)
             end if


             ! normalized lower bound value, its cdf and pdf
             truncAlpha       = (ytrunclb(nn) - forecastMean(nn,hh)) / sqrtCondVar
             call vdcdfnorm(NsimSV, truncAlpha, PROBatLB)
             PROBaboveLB      = 1 - PROBatLB ! corresponds to "Z" in wikipedia article on truncated normal

             ! CATCH CASES WHERE UNCENSORED DISTRIBUTION HAS NOT MUCH SUPPORT ABOVE LB
             where (truncAlpha < alphaThreshold)  ! alternatively: PROBatLB a little below unity

                ! good case; nondegenerate support
                truncPDFalpha    = normpdf(truncAlpha)

                ! mean of truncated normal
                truncMean        = forecastMean(nn,hh) + truncPDFalpha / PROBaboveLB * sqrtCondVar

                ! variance of truncated normal
                truncVar         = forecastCondVar(:,nn,hh) * (1 + truncAlpha * truncPDFalpha / PROBaboveLB - (truncPDFalpha / PROBaboveLB) ** 2)

             elsewhere ! untruncated density has barely any support above truncation point

                truncMean = ytrunclb(nn)
                truncVar  = 0.0d0 ! TODO: or some nominal value?

             end where

             if (any(truncMean < ytrunclb(nn))) then
                print *, 'mean below support!', minval(truncMean)
                print *, 'mu', forecastMean(nn,hh)
                print *, 'ytrunclb', ytrunclb(nn)
                call savevec(truncAlpha, 'truncAlpha.debug')
                call savevec(truncMean, 'truncMean.debug')
                call savevec(truncVar, 'truncVar.debug')
                call savevec(truncPDFalpha, 'truncPDFalpha.debug')
                call savevec(PROBaboveLB, 'PROBaboveLB.debug')
                call savevec(sqrtCondVar, 'sqrtCondVar.debug')
                stop 11
             end if


             ! integrate over SV-sims

             ! prob of being at LB
             yproblb(nn,kk) = mean(PROBatLB)


             ! adjust truncMean and truncVar for pointMass
             censMean = PROBaboveLB * truncMean + PROBatLB * ytruncLB(nn)
             censVar  = PROBaboveLB * truncVar  + PROBaboveLB * ((truncMean - censMean) ** 2) + PROBatLB * ((ytruncLB(nn) - censMean) ** 2)


             ! forecast
             yforecast(nn,kk) = mean(censMean)

             ! conditional variance
             ycondvar(nn,kk)  = variance(censMean) + mean(censVar)

             ! pdf (if predicted observations available)
             if (.not. ynanpred(nn,hh)) then


                if (ypred(nn,hh) > ytrunclb(nn)) then

                   where (truncAlpha < alphaThreshold)
                      ! pdf cond on SV draw
                      simpdf = PROBaboveLB * exp(-0.5d0 * ((ypred(nn,hh) - forecastMean(nn,hh)) ** 2) / forecastCondVar(:,nn,hh))  / (sqrttwopi * PROBaboveLB * forecastCondVar(:,nn,hh)) 
                   elsewhere
                      simpdf = 0.0d0 
                   end where

                else ! realized values are at the censoring constraint

                   simpdf = PROBatLB

                   ! where (truncAlpha < alphaThreshold)
                   !    ! pdf cond on SV draw
                   !    simpdf = exp(-0.5d0 * ((ypred(nn,hh) - forecastMean(nn,hh)) ** 2) / forecastCondVar(:,nn,hh))  / (sqrttwopi * PROBaboveLB * forecastCondVar(:,nn,hh)) 
                   ! elsewhere
                   !    simpdf = 1.0d0 
                   ! end where
                end if

                ! integrate over simulated pdfs
                ypdf(nn,kk) = mean(simpdf)
             end if


          else ! UNTRUNCATED PREDICTIONS

             ! forecast
             yforecast(nn,kk) = forecastMean(nn,hh)
             ymedian(nn,kk)   = forecastMean(nn,hh)

             ! conditional variance
             ycondvar(nn,kk)  = mean(forecastCondVar(:,nn,hh)) ! Note: there is no variance across conditional expectations here

             ! pdf (if predicted observations available)
             if (.not. ynanpred(nn,hh)) then

                ! normal pdf
                simpdf = mean(exp(-0.5d0 * (((ypred(nn,hh) - forecastMean(nn,hh)) ** 2) / forecastCondVar(:,nn,hh))) /  (sqrttwopi * sqrtCondVar))

                ypdf(nn,kk) = mean(simpdf)
             end if
          end if ! ytruncated
       end do ! nn 
    end do ! kk



  END SUBROUTINE predictiveDensitySVtruncLB

  ! @\newpage\subsection{condvarABCSV}@
  FUNCTION condvarABCSV(A,B,C,hSV,Ny,Nx,Nw,Nsv,horizons) result(SigmaY)

    INTENT(IN) :: A,B,C,hSV,Ny,Nx,Nw,Nsv,horizons

    integer :: Ny,Nx,Nw,Nsv,horizons
    double precision :: A(Nx,Nx), B(Nx,Nw), C(Ny,Nx), hSV(Nsv,horizons)
    double precision :: SigmaY(Ny,Ny,horizons)
    double precision :: VMAsv(Nx,Nw), VMA(Nx,Nw,0:horizons-1), SigmaX(Nx,Nx)

    integer jj,hh,ll,pit

    VMA        = 0.0d0
    VMA(:,:,0) = B
    DO hh = 1,horizons-1
       call dgemm('N', 'N', Nx, Nw, Nx, 1.0d0, A, Nx, VMA(:,:,hh-1), Nx, 0.0d0, VMA(:,:,hh), Nx)
    END DO

    SigmaY = 0.0d0
    DO hh=1,horizons
       SigmaX  = 0.0d0

       DO ll=0,hh-1
          ! rescale VMA to account for SV components
          VMAsv  = VMA(:,:,ll)
          pit    = hh-ll 
          ! scale SV columns
          FORALL (jj=1:Nsv) VMAsv(:,jj) = VMAsv(:,jj) * exp(0.5d0 * hSV(jj,pit)) 
          ! Sigma = VMAsv * VMAsv' + Sigma
          call dsyrk('U', 'N', Nx, Nw, 1.0d0, VMAsv, Nx, 1.0d0, SigmaX, Nx)
       END DO
       call sandwich(SigmaY(:,:,hh), C, Ny, SigmaX, Nx)

    END DO


  END FUNCTION condvarABCSV

  ! @\newpage\subsection{condvarA3BC3SV}@
  FUNCTION condvarA3BC3SV(A,B,C,hSV,Ny,Nx,Nw,Nsv,horizons) result(SigmaY)

    INTENT(IN) :: A,B,C,hSV,Ny,Nx,Nw,Nsv,horizons

    integer :: Ny,Nx,Nw,Nsv,horizons
    double precision :: A(Nx,Nx,horizons), B(Nx,Nw), C(Ny,Nx,horizons), hSV(Nsv,horizons)
    double precision :: SigmaY(Ny,Ny,horizons)
    double precision :: Bsv(Nx,Nw), SigmaX(Nx,Nx), ASigmaXA(Nx,Nx)

    integer jj,hh


    SigmaY = 0.0d0
    SigmaX = 0.0d0

    hh = 1
    FORALL (jj=1:Nsv) Bsv(:,jj) = B(:,jj) * exp(0.5d0 * hSV(jj,hh)) 
    Bsv(:,Nsv+1:Nw) = B(:,Nsv+1:Nw)
    call dsyrk('U', 'N', Nx, Nw, 1.0d0, Bsv, Nx, 1.0d0, SigmaX, Nx)
    call sandwich(SigmaY(:,:,hh), C(:,:,hh), Ny, SigmaX, Nx)

    DO hh=2,horizons

       call sandwich(ASigmaXA, A(:,:,hh), Nx, SigmaX, Nx)
       SigmaX = ASigmaXa
       FORALL (jj=1:Nsv) Bsv(:,jj) = B(:,jj) * exp(0.5d0 * hSV(jj,hh)) 
       ! Bsv(:,Nsv+1:Nw) = B(:,Nsv+1:Nw)
       call dsyrk('U', 'N', Nx, Nw, 1.0d0, Bsv, Nx, 1.0d0, SigmaX, Nx)
       call sandwich(SigmaY(:,:,hh), C(:,:,hh), Ny, SigmaX, Nx)

    END DO

    ! call savearray3(SigmaY, 'SigmaY', '.debug')
    ! call savearray3(C, 'C', '.debug')
    ! call savearray3(A, 'A', '.debug')
    ! call savemat(hSV, 'h.debug')
    ! call savemat(B, 'B.debug')
    ! stop 1

  END FUNCTION condvarA3BC3SV

  ! @\newpage\subsection{condvarVMACSV}@
  FUNCTION condvarVMACSV(VMA,C,hSV,Ny,Nx,Nw,Nsv,horizons) result(SigmaY)

    INTENT(IN) :: VMA,C,hSV,Ny,Nx,Nw,Nsv,horizons

    integer :: Ny,Nx,Nw,Nsv,horizons
    double precision :: VMA(Nx,Nw,0:horizons-1), C(Ny,Nx), hSV(Nsv,horizons)
    double precision :: SigmaY(Ny,Ny,horizons)
    double precision :: VMAsv(Nx,Nw), SigmaX(Nx,Nx)

    integer jj,hh,ll,pit

    SigmaY = 0.0d0
    DO hh=1,horizons
       SigmaX  = 0.0d0

       DO ll=0,hh-1
          ! rescale VMA to account for SV components
          VMAsv  = VMA(:,:,ll)
          pit    = hh-ll 
          ! scale SV columns
          FORALL (jj=1:Nsv) VMAsv(:,jj) = VMAsv(:,jj) * exp(0.5d0 * hSV(jj,pit)) 
          ! Sigma = VMAsv * VMAsv' + Sigma
          call dsyrk('U', 'N', Nx, Nw, 1.0d0, VMAsv, Nx, 1.0d0, SigmaX, Nx)
       END DO
       call sandwich(SigmaY(:,:,hh), C, Ny, SigmaX, Nx)

    END DO


  END FUNCTION condvarVMACSV

  ! @\newpage\subsection{crps}@
  FUNCTION crps(y,Xdraws,Ndraws) result(score)

    INTENT(IN) :: y,Xdraws,Ndraws

    integer :: Ndraws, status
    double precision :: y, score, Xdraws(Ndraws), X(Ndraws) 
    ! todo: optional flag that signals if inputs are already sorted

    integer :: i
    double precision, dimension(1:Ndraws-1) :: alpha, beta
    double precision :: iom ! "i over m"

    X = Xdraws ! NB:keeping Xdraws and X separate to avoid side effects (though FORTRAN does not complain about side effects from dlasrt dspite the INTENT(IN) attribute of forecastDraws

    ! order draws
    CALL dlasrt('I', Ndraws, X, status)
    IF (status /= 0) THEN
       write (*,*), 'DLASORT ERROR ', status, ' [COMPUTE CRPS]'
       stop 1
    END IF

    ! todo: replace with nested where statements? (note that all three conditions are mutually exclusive (and exhaustive)

    ! alpha
    where (y .lt. X(1:Ndraws-1))
       alpha = 0.0d0
       beta  = X(2:Ndraws) - X(1:Ndraws-1)
    end where

    where ((X(1:Ndraws-1) .le. y) .AND. (y .le. X(2:Ndraws)))
       alpha = y - X(1:Ndraws-1)
       beta  = X(2:Ndraws) - y
    end where

    where (X(2:Ndraws) .lt. y)
       alpha = X(2:Ndraws) - X(1:Ndraws-1)
       beta = 0.0d0
    end where

    ! sum it all together
    score = 0.0d0
    do i=1,Ndraws-1
       iom = dble(i)/dble(Ndraws)
       score = score + alpha(i) * (iom ** 2) + beta(i) * ((1.0d0 - iom) ** 2)
    end do

    ! add cases of i=0 and i=Ndraws
    if (y .gt. X(Ndraws)) then
       ! i = 0 ! NOTE: this step is void since iom = 0
       ! iom = dble(i)/dble(Ndraws)
       ! score = score + (y - X(Ndraws)) * (iom ** 2)

       ! add alpha term for i=Ndraws
       ! i = Ndraws
       ! iom = dble(i)/dble(Ndraws)
       score = score + (y - X(Ndraws)) ! * (iom ** 2)

    elseif (y .lt. X(1)) then

       ! add beta term i=0
       ! i = 0 
       ! iom = dble(i)/dble(Ndraws)
       score = score + (X(1) - y) ! * ((1.0d0 - iom) ** 2)

       ! NOTE: the following step is void since 1-iom=0
       ! i = Ndraws 
       ! iom = dble(i)/dble(Ndraws)
       ! score = score + (y - X(Ndraws)) * (iom ** 2)

    end if

  END FUNCTION crps

  ! @\newpage\subsection{logGaussianScore}@
  FUNCTION logGaussianScore(y,X,Ndraws) result(score)

    INTENT(IN) :: y,X,Ndraws

    integer :: Ndraws ! , status
    double precision :: y, score, X(Ndraws)

    double precision :: mu, sigma2

    mu           = sum(X) / dble(Ndraws)
    sigma2       = sum((X - mu) ** 2) / dble(Ndraws)

    score = -0.5d0 * (logtwopi + log(sigma2) + ((y - mu) ** 2) / sigma2)

  END FUNCTION logGaussianScore

  ! @\newpage\subsection{logGaussianScoreMV}@
  FUNCTION logGaussianScoreMV(y,X,Ny,Ndraws) result(score)

    INTENT(IN) :: y,X,Ndraws,Ny

    double precision :: score
    integer :: Ndraws, Ny, jj, ii

    double precision, dimension(Ny) :: y, ydev
    double precision, dimension(Ny,Ndraws) :: X, Xdev

    double precision :: mu(Ny), sqrtSigma(Ny,Ny), logdetSigmaHalf

    ! compute mean and deviations from mean
    mu  = sum(X,2) / dble(Ndraws)
    forall(jj=1:Ndraws,ii=1:Ny) Xdev(ii,jj) = X(ii,jj) - mu(ii)

    ! a) perform QR on resid
    sqrtSigma = qrecon(transpose(Xdev), Ndraws, Ny) / sqrt(dble(Ndraws))
    ! note: sqrtSigma may not have positive elements on its diagonals (b/o QR)
    ! but for the computations below, only squared/absolute values matter

    ! b) compute 0.5 * log(det(Sigma)) 
    logdetSigmaHalf = 0.0d0
    do ii =1,Ny
       logdetSigmaHalf = logdetSigmaHalf + abs(sqrtSigma(ii,ii))
    end do


    ! c) standardize residuals
    ydev = y - mu
    ! call savevec(ydev, 'ydev.debug')
    call dtrsv('u', 't', 'n', Ny, sqrtSigma, Ny, ydev, 1)

    ! call savemat(sqrtSigma, 'sqrtSigma.debug')
    ! call savevec(ydev, 'zdev.debug')
    ! stop 12

    ! d) compute MV normal score 
    score = -0.5d0 * Ny * logtwopi - logdetSigmaHalf - 0.5 * sum(ydev ** 2)

  END FUNCTION logGaussianScoreMV

  ! @\newpage\subsection{normlogpdf}@
  ELEMENTAL FUNCTION normlogpdf(x,mu,sigma2) 

    double precision, intent(in) :: x, mu, sigma2
    double precision :: normlogpdf

    normlogpdf = -0.5d0 * (logtwopi + log(sigma2) + (x - mu) ** 2 / sigma2)

  END FUNCTION normlogpdf

  ! @\newpage\subsection{MVnormpdf}@
  FUNCTION MVnormlogpdf(N,x,mu,sqrtVCV, flag) 

    integer, intent(in) :: N
    integer :: i
    double precision, intent(in), dimension(N)  :: x, mu
    logical, optional, intent(in) :: flag
    logical :: isLowerLeft
    double precision, dimension(N)  :: e
    double precision, intent(in), dimension(N,N)  :: sqrtVCV
    double precision :: MVnormlogpdf
    double precision :: logdetVCV

    ! sqrtVCV is supposed to be stored as right-upper-triangualr choleski

    ! parse optional flag
    if (present(flag)) then
       isLowerLeft = flag
    else
       isLowerLeft = .false.
    end if


    ! standardize e = sqrtVCV \ (x - mu)
    e = x - mu
    if (isLowerLeft) then
       call DTRSV('u','t','n',N,sqrtVCV,N,e,1)
    else
       call DTRSV('l','n','n',N,sqrtVCV,N,e,1)
    end if


    ! log of determinant of VCV 
    logdetVCV = 2 * sum( (/ (log(abs(sqrtVCV(i,i))), i = 1,N) /) )

    ! pdf
    MVnormlogpdf = -0.5d0 * (logtwopi + N + logdetVCV + sum(e ** 2))

  END FUNCTION MVnormlogpdf

  ! @\newpage\subsection{logdetSym}@
  FUNCTION logdetSym(sqrtVCV) RESULT(logdet)

    ! computes log of determinant of (strictly) p.d. symmetric matrix, given its square root
    ! uses only diagonal elements of sqrtVCV, works for either L/R factorization 

    double precision, dimension(:,:), intent(in) :: sqrtVCV
    integer :: N, i
    double precision :: logdet

    N = size(sqrtVCV,1)
    logdet = 2 * sum( (/ (log(abs(sqrtVCV(i,i))), i = 1,N) /) )


  END FUNCTION logdetSym

  ! @\newpage\subsection{trace}@
  FUNCTION trace(X) RESULT(tr)

    ! computes trace of matrix

    double precision, dimension(:,:), intent(in) :: X
    integer :: N, i
    double precision :: tr

    N = size(X,1)
    tr = sum( (/ (X(i,i), i = 1,N) /) )


  END FUNCTION trace

  ! @\newpage\subsection{igammapdf}@
  ELEMENTAL FUNCTION igammalogpdf(x, beta, dof, logGammaDof) RESULT(logp)

    double precision, intent(in) :: x, beta, logGammaDof
    integer, intent(in) :: dof
    double precision :: logp, alpha

    ! scale is "beta"
    ! dof is "alpha"
    alpha = dble(dof)

    ! first compute logs
    logp = alpha * log(beta) - (alpha + 1) * log(x) - beta / x - logGammaDof 

    ! p = exp(logp)

  END FUNCTION igammalogpdf

  ! @\newpage\subsection{logmvgamma}@
  FUNCTION logmvgamma(a,p) RESULT(logg)

    double precision, intent(in) :: a
    integer, intent(in) :: p
    double precision :: logg
    double precision, dimension(p) :: loggamma
    double precision, dimension(p) :: avec
    integer :: i

    forall (i=1:p) avec(i) = a + (1.0d0 - dble(i)) * 0.5d0

    call vdlgamma(p, avec, loggamma)

    logg = sum(loggamma) + dble(p) * dble(p - 1) * 0.25d0 * logpi

  END FUNCTION logmvgamma

  ! @\newpage\subsection{iwishlogpdf}@
  FUNCTION iwishlogpdf(N, sqrtX, sqrtPsi, dof, OPTARG_logdetPsi, OPTARG_logGammaDof) RESULT(logp)

    integer, intent(in) :: N
    double precision, dimension(N,N), intent(in) :: sqrtX, sqrtPsi
    double precision, optional, intent(in) :: OPTARG_logdetPsi, OPTARG_logGammaDof
    double precision :: logdetPsi, logGammaDof
    integer, intent(in) :: dof
    double precision :: logp
    double precision, dimension(N,N) :: sqrtInvXPsi, InvXPsi

    ! parse optional arguments
    if (present(OPTARG_logdetPsi)) then
       logdetPsi = OPTARG_logdetPsi
    else
       logdetPsi = logdetSym(sqrtPsi)
    end if
    if (present(OPTARG_logGammaDof)) then
       logGammaDof = OPTARG_logGammaDof
    else
       logGammaDof = logmvgamma(dof * 0.5d0, N)
    end if

    ! compute: trace(Psi * invX) = trace(sqrt(Psi) * sqrt(Psi)' * inv(sqrtX)' * inv(sqrtX)) = ...
    ! ... trace(sqrt(Psi)' * inv(sqrtX)' * inv(sqrtX) * sqrt(Psi))
    ! using two steps:
    ! - first, inv(sqrtX) * sqrt(Psi) vis dtsrm, hence need sqrtPsi lower-triangular
    ! - second, dsyrk 

    ! both sqrtX and sqrtPsi are assumed to be upper-triangular right factor
    ! note: must make sure that sqrtPsi has all upper-diagonal elements set to zero
    InvXPsi = 0.0d0
    sqrtInvXPsi = transpose(sqrtPsi)
    call DTRSM('l','u','n','n',N,N,1.0d0,sqrtX,N,sqrtInvXPsi,N)
    call DSYRK('u','t',N,N,1.0d0,sqrtInvXPsi,N,0.0d0,InvXPsi,N)

    logp = 0.5d0 * (dof * logdetPsi - (dof + N + 1) * logdetSym(sqrtX) - trace(InvXPsi) - dof * N * log(2.0d0))
    logp = logp - logGammaDof

  END FUNCTION iwishlogpdf

  ! @\newpage\subsection{gewekelogMDD}@
  FUNCTION gewekelogMDD(Ndraws, Ntheta, theta, priorloglike, loglike, tau, logpstar) RESULT(mdd)

    integer, intent(in) :: Ndraws, Ntheta
    double precision :: mdd
    double precision, intent(in) :: tau,logpstar
    double precision, dimension(Ntheta,Ndraws), intent(in) :: theta
    double precision, dimension(Ntheta,Ndraws) :: resid
    double precision :: Vtheta(Ntheta,Ntheta), thetabar(Ntheta)
    double precision, dimension(Ndraws), intent(in) :: priorloglike, loglike
    double precision, dimension(Ndraws) :: loglikekernel
    double precision, dimension(Ndraws) :: gewekeF, harmonicScore, thetastat
    double precision :: chi2CriticalValue
    integer :: i,j,status

    double precision :: offsetllf


    loglikekernel = priorloglike + loglike

    offsetllf = maxval(loglikekernel)

    loglikekernel = loglikekernel- offsetllf

    ! critical value for chi2
    ! note: select case works only with integers
    select case (Ntheta)
    case (95)
       if (tau == 0.5d0) then
          chi2CriticalValue = 0.943341714463404d2
          ! print *, chi2CriticalValue
       elseif  (tau == 0.25d0) then
          chi2CriticalValue = 0.853756812673199d2
       elseif  (tau == 0.75d0) then
          chi2CriticalValue = 1.038988437863585d2
       elseif  (tau == 0.9d0) then
          chi2CriticalValue = 1.130376862860290d2
       else 
          print *, 'chi2INVCDF currently not implemented for tau = ', tau
          stop 1
       end if
    case (124)
       if (tau == 0.5d0) then
          chi2CriticalValue = 1.233339742870299d2
          ! print *, chi2CriticalValue
       elseif  (tau == 0.25d0) then
          chi2CriticalValue = 1.130464068421386d2
       elseif  (tau == 0.75d0) then
          chi2CriticalValue = 1.342278181202454d2
       elseif  (tau == 0.9d0) then
          chi2CriticalValue = 1.445615572535341d2
       else 
          print *, 'chi2INVCDF currently not implemented for tau = ', tau
          stop 1
       end if
    case (126)
       if (tau == 0.5d0) then
          chi2CriticalValue = 1.253339640543887d2
          ! print *, chi2CriticalValue
       elseif  (tau == 0.25d0) then
          chi2CriticalValue = 1.149608327005966d2
          ! print *, chi2CriticalValue
       elseif  (tau == 0.75d0) then
          chi2CriticalValue =  1.363133766836778d2
       elseif  (tau == 0.9d0) then
          chi2CriticalValue = 1.467240517072355d2
       else 
          print *, 'chi2INVCDF currently not implemented for tau = ', tau
          stop 1
       end if
    case default
       print *, 'chi2INVCDF currently only implemented for Ntheta = 124, but you have Ntheta =', Ntheta
       stop 1
    end select



    ! moments of theta
    thetabar = sum(theta,2) / dble(Ndraws)
    forall (i=1:Ntheta,j=1:Ndraws) resid(i,j) = theta(i,j) - thetabar(i)
    Vtheta = 0.0d0
    call DSYRK('l','n',Ntheta,Ndraws,1.0d0 / dble(Ndraws),resid,Ntheta,0.0d0,Vtheta,Ntheta)
    ! call savemat(Vtheta, 'Vtheta.debug')
    ! call savemat(resid, 'resid.debug')

    ! sqrtVtheta
    call DPOTRF('l', Ntheta, Vtheta, Ntheta, status)
    if (status /= 0) then
       write (*,*), 'DPOTRF ERROR ', status, ' [COMPUTE GEWEKE MDD]'
       stop 1
    END IF

    ! call savemat(Vtheta, 'sqrtVtheta.debug')

    ! normalize resid
    call DTRSM('l','l','n','n',Ntheta,Ndraws,1.0d0,Vtheta,Ntheta,resid,Ntheta)
    ! call savemat(resid, 'zresid.debug')

    ! quadraticform
    thetastat = sum(resid ** 2, 1)
    ! call savevec(thetastat, 'thetastat.debug')

    ! call savevec(priorloglike, 'priorloglike.debug')
    ! call savevec(loglike, 'loglike.debug')

    ! mark 1 code:
    ! - compute ingredients of harmonic mean for all draws and in logs
    ! - then evaluate indicator function
    gewekeF    = - log(tau) - 0.5d0 * (Ntheta * logtwopi - logdetsym(Vtheta) - thetastat) - logpstar

    ! print *, 'logpstar', logpstar

    harmonicScore = gewekeF - loglikekernel
    ! call savevec(harmonicscore, 'logharmonicscore.debug')
    ! call savevec(gewekeF, 'loggewekeF.debug')

    where (thetastat .le. chi2CriticalValue)
       ! gewekeF = exp(gewekeF)
       harmonicScore = exp(harmonicScore)
    elsewhere
       ! gewekeF       = 0.0d0
       harmonicscore = 0.0d0
    end where


    ! mdd = 1.0d0 / (sum(harmonicscore) / dble(Ndraws))

    ! print *, log(mdd),  log(dble(Ndraws)) - log(sum(harmonicScore))

    ! mdd = log(mdd)

    ! ! note: computing *log* of harmonic mean in case of indicator = 1; level otherwise
    ! where (thetastat .le. chi2CriticalValue)  ! when indicator equals one
    !    gewekeF    = - log(tau) - 0.5d0 * (Ntheta * logtwopi - logdetsym(Vtheta) - thetastat)
    !    harmonicScore = exp(gewekeF - priorloglike - loglike)
    !    ! gewekeF = exp(gewekeF)
    ! elsewhere !  when indicator equals zero
    !    gewekeF       = 0.0d0
    !    harmonicscore = 0.0d0
    ! end where



    mdd = log(dble(Ndraws)) - log(sum(harmonicScore)) + offsetllf

  END FUNCTION gewekelogMDD

  ! @\newpage\subsection{MVNrnd}@
  FUNCTION MVNrnd(Ndraws, N, mu, sqrtV, VSLstream) result(x) 


    IMPLICIT NONE

    integer, intent(in) :: Ndraws, N
    double precision, dimension(N), intent(in) :: mu
    double precision, dimension(N,N), intent(in) :: sqrtV

    double precision, dimension(N,Ndraws) :: x

    integer :: i,j

    ! VSL
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state), intent(inout) :: VSLstream

    ! sample
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N * Ndraws, x, 0.0d0, 1.0d0)
    ! call savemat(x, 'z.debug')

    call DTRMM('l','u','t','n',N,Ndraws,1.0d0,sqrtV,N,x,N)

    forall(i=1:N,j=1:Ndraws) x(i,j) = x(i,j) + mu(i)

  END FUNCTION MVNrnd

  ! ! @\newpage\subsection{drawXelb}@
  ! SUBROUTINE drawXelb(Nx, Nparticles, xdraws, elbound, xprior, vecSqrtSigmaX, Cshadowrate, VSLstream)

  !   IMPLICIT NONE

  !   integer, intent(in) :: Nx, Nparticles
  !   double precision, intent(in) :: elbound
  !   ! logical, intent(in) :: elb
  !   double precision, intent(in) :: xprior(Nx,Nparticles), vecSqrtSigmaX(Nx * (Nx + 1) / 2, Nparticles)
  !   double precision, intent(in) :: Cshadowrate(Nx)

  !   double precision, dimension(Nx, Nparticles), intent(out) :: xdraws
  !   double precision, dimension(Nx, Nparticles):: zdraws
  !   double precision, dimension(Nparticles) :: shadowmean, shadowvol, zstat
  !   double precision, parameter :: critval = -3.0d0
  !   double precision :: beta(Nx), sqrtSigmaX(Nx,Nx)

  !   ! work
  !   double precision :: shadowrate(Nparticles)
  !   logical, dimension(Nparticles) :: ok
  !   integer :: k
  !   integer :: counter
  !   double precision :: ddot
  !   ! qr
  !   DOUBLE PRECISION :: qrR(1+Nx,1+Nx)
  !   INTEGER :: qrLwork

  !   ! variables for drawing from conditional posterior
  !   ! double precision :: thisposterior(Nx)

  !   ! VSL
  !   INTEGER :: errcode
  !   INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  !   type (vsl_stream_state), intent(inout) :: VSLstream

  !   ! init
  !   ok = .false.

  !   ! workspace query for QR
  !   qrR = 0.0d0
  !   qrlwork = qrquery(qrR)

  !   ! compute:
  !   ! - mean and vol of shadowrate 
  !   ! - as well as zstat for prob of elb binding
  !   ! shadowmean
  !   call DGEMV('t',Nx,Nparticles,1.0d0,xprior,Nx,Cshadowrate,1,0.0d0,shadowmean,1)
  !   ! shadowvol
  !   DO k = 1,Nparticles 
  !      beta = Cshadowrate ! helper variable
  !      call DTPMV('u','n','n',Nx,vecSqrtSigmaX(:,k),beta,1)
  !      shadowvol(k) = sqrt(sum(beta ** 2))
  !   END DO
  !   ! zstat
  !   zstat = (elbound - shadowmean) / shadowvol

  !   counter = 0
  !   do while (.not. all(ok))

  !      counter = counter + 1
  !      ! draw standard normals (always for all Nparticles, even though not all will be used after the first iteration of the while loops
  !      errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Nparticles, zdraws, 0.0d0, 1.0d0)


  !      DO k = 1,Nparticles 

  !         if (.not. ok(k)) then 
  !            ! check whether prob of elb binding too high
  !            if (zstat(k) .le. critval) then
  !               ! set shadowrate equal to elbound and draw from conditional
  !               shadowrate(k) = elbound
  !               ! compute conditional distribution of x given shadowrate using QR
  !               sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)
  !               ! call savemat(sqrtSigmaX, 'sqrtSigmaXprior.debug')
  !               qrR = 0.0d0
  !               qrR(2:1+Nx,2:1+Nx) = sqrtSigmaX
  !               qrR(2:1+Nx,1) = Cshadowrate
  !               call DTPMV('u','n','n',Nx,vecSqrtSigmaX(:,k),qrR(2:1+Nx,1),1)
  !               ! QR decomposition
  !               call qrot(qrR, qrLWORK)
  !               ! map qr into Kalman objects
  !               beta        = qrR(1,2:1+Nx) /  qrR(1,1)
  !               sqrtSigmaX  = qrR(2:1+Nx,2:1+Nx) ! upper triangular

  !               ! scale zdraws
  !               call DTRMV('u','t','n',Nx,sqrtSigmaX,Nx,zdraws(:,k),1)
  !               ! add mean
  !               xdraws(:,k) = xprior(:,k) + beta * (elbound - shadowmean(k))  + zdraws(:,k)

  !               ok(k) = .true.

  !            else ! regular draw

  !               ! scale standard normal draw by variance
  !               call DTPMV('u','t','n',Nx,vecSqrtSigmaX(:,k),zdraws(:,k),1)
  !               ! add mean
  !               xdraws(:,k) = xprior(:,k) + zdraws(:,k)
  !               ! shadowrate(k) = sum(Cshadowrate * xdraws(:,k))
  !               shadowrate(k) = ddot(Nx, Cshadowrate, 1, xdraws(:,k), 1)
  !               ok(k) = shadowrate(k) .le. elbound
  !            end if ! if zstat
  !         end if ! if shadowrate > elbound
  !      end do ! k particles

  !      ! print *, counter, 'draws', count(ok), 'fits'
  !      if (counter .gt. 10000) then
  !         call savevec(shadowrate, 'shadowrate.debug')
  !         call savemat(xdraws, 'xdraws.debug')
  !         call savevec(zstat, 'zstat.debug')
  !         call savevec(Cshadowrate, 'Cshadowrate.debug')
  !         stop 11
  !      end if


  !   end do

  ! END SUBROUTINE drawXelb

  ! @\newpage\subsection{drawXelb}@
  SUBROUTINE drawXelb(Nx, Nparticles, xdraws, elbound, xprior, vecSqrtSigmaX, Cshadowrate, VSLstream)

    ! alternative algorithm: first, always draws a conforming shadowrate, then draws rest from the conditional MVN (given the shadow rate)

    IMPLICIT NONE

    integer, intent(in) :: Nx, Nparticles
    double precision, intent(in) :: elbound
    double precision, intent(in) :: xprior(Nx,Nparticles), vecSqrtSigmaX(Nx * (Nx + 1) / 2, Nparticles)
    double precision, intent(in) :: Cshadowrate(Nx)

    double precision, dimension(Nx, Nparticles), intent(out) :: xdraws
    double precision, dimension(Nx, Nparticles):: zdraws
    double precision, dimension(Nparticles) :: shadowmean, shadowvol, zstat
    double precision, parameter :: critval = -3.0d0
    double precision :: beta(Nx), sqrtSigmaX(Nx,Nx)

    ! proposal draws for rejection sampling the shadowrate
    ! problem: we won't know how many proposals are needed to sample a conforming shadowrate
    ! idea: 
    ! - create a buffer "propdraws" of standard-normal proposal draws
    ! - move through the buffer as long as there are unused elements
    ! - once all elements have been used, draw a new set of propdraws
    integer, parameter :: Nproposals = 1000
    integer :: propcount
    double precision, dimension(Nproposals) :: propdraws

    ! work
    double precision :: shadowrate(Nparticles)
    logical, dimension(Nparticles) :: ok
    integer :: k
    ! double precision :: ddot
    ! qr
    DOUBLE PRECISION :: qrR(1+Nx,1+Nx)
    INTEGER :: qrLwork

    ! VSL
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state), intent(inout) :: VSLstream

    ! logical :: debugflag 

    ! init
    ok = .false.

    ! workspace query for QR
    qrR = 0.0d0
    qrlwork = qrquery(qrR)

    ! compute:
    ! - mean and vol of shadowrate 
    ! - as well as zstat for prob of elb binding
    ! shadowmean
    call DGEMV('t',Nx,Nparticles,1.0d0,xprior,Nx,Cshadowrate,1,0.0d0,shadowmean,1)
    ! shadowvol
    DO k = 1,Nparticles 
       beta = Cshadowrate ! helper variable
       call DTPMV('u','n','n',Nx,vecSqrtSigmaX(:,k),beta,1)
       shadowvol(k) = sqrt(sum(beta ** 2))
    END DO
    ! zstat
    zstat = (elbound - shadowmean) / shadowvol

    ! STEP 1: generate shadowrate that conforms with elb
    ! 1a: set equal to elbound when only negligible mass below elbound
    where (zstat .le. critval) 
       shadowrate = elbound
       ok = .true.
    end where
    ! if (any(ok)) then
    !    print *, 'found some zstat cases'
    !    debugflag = .true.
    ! else
    !    debugflag = .false.
    ! end if

    ! now loop over all particles and sample where not .ok.
    propcount = Nproposals + 1
    do k=1,Nparticles
       do while (.not. ok(k)) 
          if (propcount .gt. Nproposals) then
             errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nproposals, propdraws, 0.0d0, 1.0d0)
             propcount = 1
             ! print *, 'regenerating propdraws', tidomp()
          end if
          shadowrate(k) = shadowmean(k) + shadowvol(k) * propdraws(propcount)
          ok(k) = shadowrate(k) .le. elbound
          propcount = propcount + 1
       end do
    end do
    
    ! if (debugflag) then
    !    call savevec(shadowrate, 'shadowrate.debug')
    !    call savevec(zstat, 'zstat.debug')
    !    stop 12
    ! end if

    ! STEP 2: draw xdraws conditional on the shadowrate
    ! draw standard normals
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Nparticles, zdraws, 0.0d0, 1.0d0)
    DO k = 1,Nparticles 

       ! compute conditional distribution of x given shadowrate using QR
       sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)
       ! call savemat(sqrtSigmaX, 'sqrtSigmaXprior.debug')
       qrR = 0.0d0
       qrR(2:1+Nx,2:1+Nx) = sqrtSigmaX
       qrR(2:1+Nx,1) = Cshadowrate
       call DTPMV('u','n','n',Nx,vecSqrtSigmaX(:,k),qrR(2:1+Nx,1),1)
       ! QR decomposition
       call qrot(qrR, qrLWORK)
       ! map qr into Kalman objects
       beta        = qrR(1,2:1+Nx) /  qrR(1,1)
       sqrtSigmaX  = qrR(2:1+Nx,2:1+Nx) ! upper triangular
       ! scale zdraws
       call DTRMV('u','t','n',Nx,sqrtSigmaX,Nx,zdraws(:,k),1)
       ! add mean
       xdraws(:,k) = xprior(:,k) + beta * (shadowrate(k) - shadowmean(k))  + zdraws(:,k)
       
    end do ! k particles

  END SUBROUTINE drawXelb

END MODULE densitybox


