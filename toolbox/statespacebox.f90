MODULE statespacebox

  use embox, only : savemat, savevec, hrulefill, display ! i/o for debugging
  use blaspack, only : eye, triu, sandwichplus, sandwich, pi
  use vslbox
  use omp_lib


  IMPLICIT NONE

CONTAINS


  ! @\newpage\subsection{DLYAP}@
  SUBROUTINE DLYAP(Sigma, A, B, Nx, Nw, status)
    ! doubling, calling DSYMM and DGEMM
    ! Sigma = A * Sigma * A' + B * B'
    ! output Sigma is symmetric
    IMPLICIT NONE
    INTENT(INOUT) :: Sigma
    INTENT(IN) :: A, B, Nx, Nw
    INTENT(OUT) status

    INTEGER :: iter, i, status, Nx, Nw
    INTEGER, PARAMETER :: maxiter = 100
    DOUBLE PRECISION, PARAMETER :: tol = 1.0d-8

    LOGICAL :: converged
    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: Sigma0, Sigma, A, AA, AAA, AASigma
    DOUBLE PRECISION, DIMENSION(Nx,Nw) :: B

    ! Sigma0 = B B'
    Sigma0 = 0.0d0
    call DSYRK('U','N',Nx,Nw,1.0d0,B,Nx,0.0d0,Sigma0,Nx)
    ! fill up lower triangular -- necessary for DGEMM below 
    FORALL (i=2:Nx) Sigma0(i,1:i-1) = Sigma0(1:i-1,i)

    ! call savemat(B, 'B.debug')
    ! call savemat(Sigma0, 'S0.debug')

    converged = .false.
    iter = 0

    AA = A
    DO 

       iter = iter + 1


       ! some debugging
       ! call savemat(AA, 'A' // trim(int2str(iter)) // '.debug')
       ! call savemat(Sigma0, 'Sigma' // trim(int2str(iter)) // '.debug')

       ! call sandwichplus(Sigma, AA, Nx, Sigma0, Nx)
       ! MANUAL SANDWICHPLUS: Sigma = AA * Sigma0 * AA' + Sigma
       call DSYMM('R','U',Nx,Nx,1.0d0,Sigma0,Nx,AA,Nx,0.0d0,AASigma,Nx)
       Sigma = Sigma0
       call DGEMM('N','T',Nx,Nx,Nx,1.0d0,AASigma,Nx,AA,Nx,1.0d0,Sigma,Nx)

       ! balance for symmetry
       Sigma = 0.5d0 * (Sigma + transpose(Sigma)) 

       IF (abs(maxval(Sigma - Sigma0)) < tol) converged = .true.       

       ! print *, iter, abs(maxval(Sigma - Sigma0)), tol
       ! Sigma = (Sigma + transpose(Sigma)) / dble(2)

       IF (converged .OR. (iter > maxiter)) EXIT

       ! AAA = AA * AA
       call DGEMM('N','N',Nx,Nx,Nx,1.0d0,AA,Nx,AA,Nx,0.0d0,AAA,Nx)
       AA     = AAA
       Sigma0 = Sigma

    END DO

    ! call savemat(Sigma, 'Sigma.debug')

    IF (converged) THEN
       status = 0
    ELSE
       status = -1
    END IF


  END SUBROUTINE DLYAP

  ! @\newpage\subsection{DLYAPsqrt}@
  SUBROUTINE DLYAPsqrt(Sigma, A, B, Nx, Nw, status)
    ! doubling, calling DSYMM and DGEMM, then factorization
    ! Sigma = A * Sigma * A' + B * B'
    ! output Sigma is symmetric
    IMPLICIT NONE
    INTENT(INOUT) :: Sigma
    INTENT(IN) :: A, B, Nx, Nw
    INTENT(OUT) status

    INTEGER :: iter, i, status, Nx, Nw
    INTEGER, PARAMETER :: maxiter = 100
    DOUBLE PRECISION, PARAMETER :: tol = 1.0d-8

    LOGICAL :: converged
    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: Sigma0, Sigma, A, AA, AAA, AASigma
    DOUBLE PRECISION, DIMENSION(Nx,Nw) :: B

    ! Sigma0 = B B'
    Sigma0 = 0.0d0
    call DSYRK('U','N',Nx,Nw,1.0d0,B,Nx,0.0d0,Sigma0,Nx)
    ! fill up lower triangular -- necessary for DGEMM below 
    FORALL (i=2:Nx) Sigma0(i,1:i-1) = Sigma0(1:i-1,i)

    ! call savemat(B, 'B.debug')
    ! call savemat(Sigma0, 'S0.debug')

    converged = .false.
    iter = 0

    AA = A
    DO 

       iter = iter + 1


       ! some debugging
       ! call savemat(AA, 'A' // trim(int2str(iter)) // '.debug')
       ! call savemat(Sigma0, 'Sigma' // trim(int2str(iter)) // '.debug')

       ! call sandwichplus(Sigma, AA, Nx, Sigma0, Nx)
       ! MANUAL SANDWICHPLUS: Sigma = AA * Sigma0 * AA' + Sigma
       call DSYMM('R','U',Nx,Nx,1.0d0,Sigma0,Nx,AA,Nx,0.0d0,AASigma,Nx)
       Sigma = Sigma0
       call DGEMM('N','T',Nx,Nx,Nx,1.0d0,AASigma,Nx,AA,Nx,1.0d0,Sigma,Nx)

       ! balance for symmetry
       Sigma = 0.5d0 * (Sigma + transpose(Sigma)) 

       IF (abs(maxval(Sigma - Sigma0)) < tol) converged = .true.       

       ! print *, iter, abs(maxval(Sigma - Sigma0)), tol
       ! Sigma = (Sigma + transpose(Sigma)) / dble(2)

       IF (converged .OR. (iter > maxiter)) EXIT

       ! AAA = AA * AA
       call DGEMM('N','N',Nx,Nx,Nx,1.0d0,AA,Nx,AA,Nx,0.0d0,AAA,Nx)
       AA     = AAA
       Sigma0 = Sigma

    END DO

    ! call savemat(Sigma, 'Sigma.debug')

    IF (converged) THEN
       status = 0
       call DPOTRF('U', Nx, Sigma, Nx, status)
       if (status == 0) call triu(Sigma)
    ELSE
       status = -1
    END IF


  END SUBROUTINE DLYAPsqrt

  ! @\newpage\subsection{simyABC}@
  FUNCTION simyABC(T,Ny,Nx,Nw,A,B,C,x0,VSLstream) result(y) 

    ! INTENT(OUT) :: y
    INTENT(IN) :: Ny,Nx,Nw,T,A,B,C,x0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: j,errcode,T,Ny,Nx,Nw
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx) :: C

    DOUBLE PRECISION, DIMENSION(Nx)     :: x0

    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: x
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y
    DOUBLE PRECISION, DIMENSION(Nw,1:T) :: w

    ! draw random numbers
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nw, w, 0.0d0, 1.0d0)

    x = 0.0d0
    y = 0.0d0

    ! scale shocks
    DO j=1,T
       call DGEMV('N',Nx,Nw,1.0d0,B,Nx,w(:,j),1,0.0d0,x(:,j),1)
    END DO

    ! simulate state
    j = 1
    call DGEMV('N',Nx,Nx,1.0d0,A,Nx,x0,1,1.0d0,x(:,j),1)
    DO j=2,T
       call DGEMV('N',Nx,Nx,1.0d0,A,Nx,x(:,j-1),1,1.0d0,x(:,j),1)
    END DO


    ! simulate observer
    ! TODO: replace with a single DGEMM?
    DO j=1,T
       call DGEMV('N',Ny,Nx,1.0d0,C,Ny,x(:,j),1,0.0d0,y(:,j),1)
    END DO

  END FUNCTION simyABC

  ! @\newpage\subsection{simyABCsvol0}@
  FUNCTION simyABCsvol0(T,Ny,Nx,Nw,A,B,C,SVol,x0,VSLstream) result(y) 

    ! INTENT(OUT) :: y
    INTENT(IN) :: Ny,Nx,Nw,T,A,B,C,SVol,x0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: j,i, errcode,T,Ny,Nx,Nw
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx) :: C

    DOUBLE PRECISION, DIMENSION(Nx)     :: x0

    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: x
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y
    DOUBLE PRECISION, DIMENSION(Nw,1:T) :: w, SVol

    ! draw random numbers
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nw, w, 0.0d0, 1.0d0)

    x = 0.0d0
    y = 0.0d0

    ! apply SVol to shocks
    FORALL (j=1:T,i=1:Nw) w(i,j) = SVol(i,j) * w(i,j) 

    ! scale shocks
    DO j=1,T
       call DGEMV('N',Nx,Nw,1.0d0,B,Nx,w(:,j),1,0.0d0,x(:,j),1)
    END DO

    ! simulate state
    j = 1
    call DGEMV('N',Nx,Nx,1.0d0,A,Nx,x0,1,1.0d0,x(:,j),1)
    DO j=2,T
       call DGEMV('N',Nx,Nx,1.0d0,A,Nx,x(:,j-1),1,1.0d0,x(:,j),1)
    END DO


    ! simulate observer
    ! TODO: replace wiht a single DGEMM?
    DO j=1,T
       call DGEMV('N',Ny,Nx,1.0d0,C,Ny,x(:,j),1,0.0d0,y(:,j),1)
    END DO

  END FUNCTION simyABCsvol0

  ! @\newpage\subsection{simA3B3C3noise}@
  SUBROUTINE simA3B3C3noise(y,x,ynoise,xshock,T,Ny,Nx,Nw,A,B,C,noiseVol,Ex0,sqrtVx0,VSLstream)

    ! same as A3B3C3 but with noise shocks attached to each measurement
    ! noise is mutually uncorrelated (but possibly time-varying) iid, noiseVol is thus Ny x T

    INTENT(OUT) :: y,x,xshock,ynoise
    INTENT(IN) :: Ny,Nx,Nw,T,A,B,C,noiseVol,Ex0,sqrtVx0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: j,k,errcode,T,Ny,Nx,Nw
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C
    DOUBLE PRECISION, DIMENSION(Ny,T) :: noiseVol

    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, x0
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: x
    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: xshock
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y, ynoise
    DOUBLE PRECISION, DIMENSION(Nw,1:T) :: w

    ! draw random numbers
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nw, w, 0.0d0, 1.0d0)
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Ny, ynoise, 0.0d0, 1.0d0)
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx, x0, 0.0d0, 1.0d0)

    ! construct x0
    x(:,0) = Ex0
    call DGEMV('N',Nx,Nx,1.0d0,sqrtVx0,Nx,x0,1,1.0d0,x(:,0),1)

    ! scale shocks
    DO j=1,T
       call DGEMV('N',Nx,Nw,1.0d0,B(:,:,j),Nx,w(:,j),1,0.0d0,xshock(:,j),1)
    END DO

    ! scale ynoise
    FORALL (j=1:T,k=1:Ny) ynoise(k,j) = noiseVol(k,j) * ynoise(k,j) 
    ! should also work: ynoise = noiseVol * ynoise 

    ! copy shocks into x and y
    x(:,1:T) = xshock
    y = ynoise

    ! simulate state and observer
    DO j=1,T
       call DGEMV('N',Nx,Nx,1.0d0,A(:,:,j),Nx,x(:,j-1),1,1.0d0,x(:,j),1)
       call DGEMV('N',Ny,Nx,1.0d0,C(:,:,j),Ny,x(:,j),  1,1.0d0,y(:,j),1)
    END DO

  END SUBROUTINE simA3B3C3noise

  ! @\newpage\subsection{simA3B3C3}@
  SUBROUTINE simA3B3C3(y,x,xshock,T,Ny,Nx,Nw,A,B,C,Ex0,sqrtVx0,VSLstream)

    INTENT(OUT) :: y,x,xshock
    INTENT(IN) :: Ny,Nx,Nw,T,A,B,C,Ex0,sqrtVx0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: j, errcode,T,Ny,Nx,Nw
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C

    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, x0
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: x
    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: xshock
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y
    DOUBLE PRECISION, DIMENSION(Nw,1:T) :: w

    ! draw random numbers
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nw, w, 0.0d0, 1.0d0)
    errcode    = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx, x0, 0.0d0, 1.0d0)

    ! construct x0
    x(:,0) = Ex0
    call DGEMV('N',Nx,Nx,1.0d0,sqrtVx0,Nx,x0,1,1.0d0,x(:,0),1)

    ! scale shocks
    DO j=1,T
       call DGEMV('N',Nx,Nw,1.0d0,B(:,:,j),Nx,w(:,j),1,0.0d0,xshock(:,j),1)
    END DO

    ! copy shocks into x
    x(:,1:T) = xshock

    ! simulate state and observer
    DO j=1,T
       call DGEMV('N',Nx,Nx,1.0d0,A(:,:,j),Nx,x(:,j-1),1,1.0d0,x(:,j),1)
       call DGEMV('N',Ny,Nx,1.0d0,C(:,:,j),Ny,x(:,j),1,0.0d0,y(:,j),1)
    END DO


  END SUBROUTINE simA3B3C3


  ! @\newpage\subsection{disturbancesmootherA3B3C3scalar}@
  SUBROUTINE disturbancesmootherA3B3C3scalar(xhat,xshockhat,ynoisehat,y,T,Ny,Nx,Nw,A,B,C,noiseVariance,Ex0,Vx0,status)

    ! same as A3B3C3 but with noise shocks attached to each measurement
    ! noise is mutually uncorrelated (but possibly time-varying) iid, noiseVariance is thus Ny x T

    INTENT(OUT) :: xhat,xshockhat,ynoisehat,status
    INTENT(IN) :: y,Ny,Nx,Nw,T,A,B,C,noiseVariance,Ex0,Vx0

    INTEGER :: j,i,ii,T,Ny,Nx,Nw,status

    DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, ZERO = 0.0d0

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C
    DOUBLE PRECISION, DIMENSION(Ny,T) :: noiseVariance

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xhat
    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y,ytilde,ynoisehat

    DOUBLE PRECISION :: invSigmaY, K(Nx,Ny,T), ImKC(Nx,Nx), Ex0(Nx), Vx0(Nx,Nx), ddot

    DOUBLE PRECISION, DIMENSION(Nx,Nx,0:T) :: Sigma
    DOUBLE PRECISION, DIMENSION(Nx,Nx,1:T) :: BB

    DOUBLE PRECISION, DIMENSION(Nx)     :: S, Sprior

    status = 0

    ! init memory
    xhat  = zero
    Sigma = zero
    BB    = zero
    status = 0

    ! initial parameters
    xhat(:,0) = Ex0
    Sigma(:,:,0) = Vx0
    ytilde = 0.0d0

    ! NOTES:
    ! - ytilde is scaled by invSigmaY

    ! forward loop
    DO j=1,T

       ! prepare priors
       ! xhat(:,:,j) = A(:,:,j) * xhat(:,:,j-1)
       call DGEMV('n',Nx,Nx,1.0d0,A(:,:,j),Nx,xhat(:,j-1),1,0.0d0,xhat(:,j),1)

       ! VCV: Sigma(j) = A * Sigma(j-1) * A' + BB
       call DSYRK('U','N',Nx,Nw,1.0d0,B(:,:,j),Nx,0.0d0,BB(:,:,j),Nx)
       Sigma(:,:,j) = BB(:,:,j)
       call sandwichplus(Sigma(:,:,j), A(:,:,j), Nx, Sigma(:,:,j-1), Nx)

       ! loop over the individual measurements
       DO i=1,Ny

          ! Kalman Gain (not yet scaled by invSigmaY)
          CALL DSYMV('U', Nx, 1.0d0, Sigma(:,:,j), Nx, C(i,:,j), 1, 0.0d0, K(:,i,j), 1)
          ! invSigmaY
          invSigmaY = ddot(Nx, C(i,:,j), 1, K(:,i,j), 1) + noiseVariance(i,j)
          IF (invSigmaY .lt. 0.0d0) THEN
             call hrulefill
             print *, 'Error in Scalar Smoother. Negative SigmaY. [STATESPACEBOX]'
             print *, 'SigmaY=', invSigmaY
             print *, 'j=', j, 'i=', i
             call hrulefill
             status = 1
             return
          ELSE
             invSigmaY = 1.0d0 / invSigmaY
          END IF

          ! innovation ytilde (scaled by invSigmaY)
          ytilde(i,j) = (y(i,j) - ddot(Nx, C(i,:,j), 1, xhat(:,j), 1)) * invSigmaY

          ! update estimate
          xhat(:,j) = xhat(:,j) + K(:,i,j) * ytilde(i,j)

          ! Posterior VCV: Sigma(:,:,j) = Sigma(:,:,j) - K * invSigmaY * K'
          CALL DSYR('U', Nx, -invSigmaY, K(:,i,j), 1, Sigma(:,:,j), Nx)

          ! scale K (needed for smoother)
          K(:,i,j) = K(:,i,j) * invSigmaY

       END DO

    END DO

    ! START SMOOTHING

    j = T
    S = 0.0d0
    DO i=Ny,1,-1

       Sprior = S
       ! I - K(i)* C(i)
       ImKC = 0.0d0
       FORALL (ii=1:Nx) ImKC(ii,ii) = 1.0d0
       call DGER(Nx,Nx,-1.0d0,K(:,i,j),1,C(i,:,j),1,ImKC,Nx)

       ! S(i) = C(i)' ytilde(i) + (I-K(i)C(i))'S(i+1) 
       S = C(i,:,j) * ytilde(i,j)
       call DGEMV('T',Nx,Nx,1.0d0,ImKC,Nx,Sprior,1,1.0d0,S,1)

    END DO
    ! xshockhat = BB * S
    call DSYMV('U',Nx,1.0d0,BB(:,:,j),Nx,S,1,0.0d0,xshockhat(:,j),1)
    ! Sprior = A(j)' S
    call DGEMV('T',Nx,Nx,1.0d0,A(:,:,j),Nx,S,1,0.0d0,Sprior,1)

    ! backward loop
    DO j=T-1,1,-1

       ! xhat(j) = xhat(j) + Sigma(j) * Sprior
       call DSYMV('U',Nx,1.0d0,Sigma(:,:,j),Nx,Sprior,1,1.0d0,xhat(:,j),1)

       DO i=Ny,1,-1
          ! I - K(i)* C(i)
          ImKC = 0.0d0
          FORALL (ii=1:Nx) ImKC(ii,ii) = 1.0d0
          call DGER(Nx,Nx,-1.0d0,K(:,i,j),1,C(i,:,j),1,ImKC,Nx)
          ! S(i) = C(i)' ytilde(i) + (I-K(i)C(i))'S(i+1), where S(i) is S, S(i+1) is Sprior
          S = C(i,:,j) * ytilde(i,j)
          call DGEMV('T',Nx,Nx,1.0d0,ImKC,Nx,Sprior,1,1.0d0,S,1)
          Sprior = S
       END DO

       ! xshockhat(j) = BB(j) * S
       call DSYMV('U',Nx,1.0d0,BB(:,:,j),Nx,S,1,0.0d0,xshockhat(:,j),1)

       ! Sprior = A(j)' S
       call DGEMV('T',Nx,Nx,1.0d0,A(:,:,j),Nx,S,1,0.0d0,Sprior,1)

    END DO

    ! x(0|0) = x(0) + Sigma(0) * Sprior
    ! call DGEMV('N',Nx,Nx,1.0d0,Sigma(:,:,0),Nx,Sprior,1,1.0d0,xhat(:,0),1)
    call DSYMV('U',Nx,1.0d0,Sigma(:,:,0),Nx,Sprior,1,1.0d0,xhat(:,0),1)

    ! DONE SMOOTHING

    ! ynoisehat = y - C * xhat
    ynoisehat = y
    DO j=1,T
       call DGEMV('N', Ny, Nx, -1.0d0, C(:,:,j), Ny, xhat(:,j), 1, 1.0d0, ynoisehat(:,j), 1)
    END DO

  END SUBROUTINE disturbancesmootherA3B3C3scalar

  ! @\newpage\subsection{disturbancesmootherA3B3C3NANscalar}@
  SUBROUTINE disturbancesmootherA3B3C3NANscalar(xhat,xshockhat,SigmaStarT,ynoisehat,y,yNaN,T,Ny,Nx,Nw,A,B,C,noiseVariance,Ex0,Vx0,status)

    INTENT(OUT) :: xhat,xshockhat,SigmaStarT,ynoisehat,status
    INTENT(IN) :: y,yNaN,Ny,Nx,Nw,T,A,B,C,noiseVariance,Ex0,Vx0

    INTEGER :: j,i,ii,T,Ny,Nx,Nw,status

    DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, ZERO = 0.0d0

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C
    DOUBLE PRECISION, DIMENSION(Ny,T) :: noiseVariance

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xhat
    DOUBLE PRECISION, DIMENSION(Nx,1:T) :: xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y,ytilde,ynoisehat
    LOGICAL, DIMENSION(Ny,1:T) :: yNaN

    DOUBLE PRECISION :: invSigmaY, K(Nx,Ny,T), ImKC(Nx,Nx), Ex0(Nx), Vx0(Nx,Nx), SigmaStarT(Nx,Nx), ddot

    DOUBLE PRECISION, DIMENSION(Nx,Nx,0:T) :: Sigma
    DOUBLE PRECISION, DIMENSION(Nx,Nx,1:T) :: BB

    DOUBLE PRECISION, DIMENSION(Nx)     :: S, Sprior

    status = 0

    ! init memory
    xhat  = zero
    Sigma = zero
    BB    = zero
    status = 0

    ! initial parameters
    xhat(:,0) = Ex0
    Sigma(:,:,0) = Vx0
    ytilde = 0.0d0

    ! NOTES:
    ! - ytilde is scaled by invSigmaY

    ! forward loop
    DO j=1,T

       ! prepare priors
       ! xhat(:,:,j) = A(:,:,j) * xhat(:,:,j-1)
       call DGEMV('n',Nx,Nx,1.0d0,A(:,:,j),Nx,xhat(:,j-1),1,0.0d0,xhat(:,j),1)

       ! VCV: Sigma(j) = A * Sigma(j-1) * A' + BB
       call DSYRK('U','N',Nx,Nw,1.0d0,B(:,:,j),Nx,0.0d0,BB(:,:,j),Nx)
       Sigma(:,:,j) = BB(:,:,j)
       call sandwichplus(Sigma(:,:,j), A(:,:,j), Nx, Sigma(:,:,j-1), Nx)

       ! loop over the individual measurements
       DO i=1,Ny

          IF (.NOT. yNaN(i,j)) THEN

             ! Kalman Gain (not yet scaled by invSigmaY)
             CALL DSYMV('U', Nx, 1.0d0, Sigma(:,:,j), Nx, C(i,:,j), 1, 0.0d0, K(:,i,j), 1)
             ! invSigmaY
             invSigmaY = ddot(Nx, C(i,:,j), 1, K(:,i,j), 1) + noiseVariance(i,j)
             IF (invSigmaY .lt. 0.0d0) THEN
                call hrulefill
                print *, 'Error in Scalar Smoother (A3B3C3nanscalar). Negative SigmaY. [STATESPACEBOX]'
                print *, 'SigmaY=', invSigmaY
                print *, 'j=', j, 'i=', i
                call hrulefill
                status = 1
                return
             ELSE
                invSigmaY = 1.0d0 / invSigmaY
             END IF

             ! innovation ytilde (scaled by invSigmaY)
             ytilde(i,j) = (y(i,j) - ddot(Nx, C(i,:,j), 1, xhat(:,j), 1)) * invSigmaY

             ! update estimate
             xhat(:,j) = xhat(:,j) + K(:,i,j) * ytilde(i,j)

             ! Posterior VCV: Sigma(:,:,j) = Sigma(:,:,j) - K * invSigmaY * K'
             CALL DSYR('U', Nx, -invSigmaY, K(:,i,j), 1, Sigma(:,:,j), Nx)

             ! scale K (needed for smoother)
             K(:,i,j) = K(:,i,j) * invSigmaY

          ELSE

             K(:,i,j)    = 0.0d0 ! just to avoid bad effects, C(i,:,j) should be zero anyway
             ytilde(i,j) = 0

          END IF


       END DO

    END DO

    ! START SMOOTHING

    j = T
    S      = 0.0d0
    Sprior = 0.0d0 ! needed if all yNaN(:,t) are true (to support S = Sprior below)
    DO i=Ny,1,-1
       IF (.NOT. yNaN(i,j)) THEN
          Sprior = S
          ! I - K(i)* C(i)
          ImKC = 0.0d0
          FORALL (ii=1:Nx) ImKC(ii,ii) = 1.0d0
          call DGER(Nx,Nx,-1.0d0,K(:,i,j),1,C(i,:,j),1,ImKC,Nx)

          ! S(i) = C(i)' ytilde(i) + (I-K(i)C(i))'S(i+1) 
          S = C(i,:,j) * ytilde(i,j)
          call DGEMV('T',Nx,Nx,1.0d0,ImKC,Nx,Sprior,1,1.0d0,S,1)
       ELSE
          S = Sprior ! will be used after the DO loop for xshockhat and Sprior
       END IF
    END DO
    ! xshockhat = BB * S
    call DSYMV('U',Nx,1.0d0,BB(:,:,j),Nx,S,1,0.0d0,xshockhat(:,j),1)
    ! Sprior = A(j)' S
    call DGEMV('T',Nx,Nx,1.0d0,A(:,:,j),Nx,S,1,0.0d0,Sprior,1)


    SigmaStarT = Sigma(:,:,T)

    ! backward loop
    DO j=T-1,1,-1

       ! xhat(j) = xhat(j) + Sigma(j) * Sprior
       call DSYMV('U',Nx,1.0d0,Sigma(:,:,j),Nx,Sprior,1,1.0d0,xhat(:,j),1)

       DO i=Ny,1,-1
          IF (.NOT. yNaN(i,j)) THEN
             ! I - K(i)* C(i)
             ImKC = 0.0d0
             FORALL (ii=1:Nx) ImKC(ii,ii) = 1.0d0
             call DGER(Nx,Nx,-1.0d0,K(:,i,j),1,C(i,:,j),1,ImKC,Nx)
             ! S(i) = C(i)' ytilde(i) + (I-K(i)C(i))'S(i+1), where S(i) is S, S(i+1) is Sprior
             S = C(i,:,j) * ytilde(i,j)
             call DGEMV('T',Nx,Nx,1.0d0,ImKC,Nx,Sprior,1,1.0d0,S,1)
             Sprior = S
          ELSE
             S = Sprior ! will be used after the DO loop for xshockhat and Sprior
          END IF

       END DO

       ! xshockhat(j) = BB(j) * S
       call DSYMV('U',Nx,1.0d0,BB(:,:,j),Nx,S,1,0.0d0,xshockhat(:,j),1)

       ! Sprior = A(j)' S
       call DGEMV('T',Nx,Nx,1.0d0,A(:,:,j),Nx,S,1,0.0d0,Sprior,1)

    END DO

    ! x(0|0) = x(0) + Sigma(0) * Sprior
    ! call DGEMV('N',Nx,Nx,1.0d0,Sigma(:,:,0),Nx,Sprior,1,1.0d0,xhat(:,0),1)
    call DSYMV('U',Nx,1.0d0,Sigma(:,:,0),Nx,Sprior,1,1.0d0,xhat(:,0),1)

    ! DONE SMOOTHING

    ! ynoisehat = y - C * xhat
    ynoisehat = y
    DO j=1,T
       call DGEMV('N', Ny, Nx, -1.0d0, C(:,:,j), Ny, xhat(:,j), 1, 1.0d0, ynoisehat(:,j), 1)
    END DO

  END SUBROUTINE disturbancesmootherA3B3C3NANscalar

  ! @\newpage\subsection{samplerA3B3C3nanscalar}@
  SUBROUTINE samplerA3B3C3nanscalar(xdraw,xshockdraw,SigmaStarT,y,yNaN,T,Ny,Nx,Nw,A,B,C,Ex0,sqrtVx0,VSLstream,status)

    INTENT(OUT) :: xdraw,xshockdraw,SigmaStarT,status
    INTENT(IN) :: y,Ny,Nx,Nw,T,A,B,C,Ex0,sqrtVx0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: T,Ny,Nx,Nw,status

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C

    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, xNull
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0, Vx0, SigmaStarT

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xdraw, xhat
    DOUBLE PRECISION, DIMENSION(Nx,T)   :: xshockdraw, xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y, ydraw, dummynoise, dummyvols
    LOGICAL, DIMENSION(Ny,1:T) :: yNaN


    ! 1) generate plus data
    call simA3B3C3(ydraw,xdraw,xshockdraw,T,Ny,Nx,Nw,A,B,C,Ex0,sqrtVx0,VSLstream)
    ydraw = ydraw - y

    ! 2) filter the difference
    Vx0 = 0.0d0
    call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,Vx0,Nx)

    ! mean adjustment since projecting onto ydraw - y
    xNull = 0.0d0
    dummyvols = 0.0d0
    call disturbancesmootherA3B3C3nanscalar(xhat,xshockhat,SigmaStarT,dummynoise,ydraw,yNaN,T,Ny,Nx,Nw,A,B,C,dummyvols,xNull,Vx0,status)

    if (status /= 0) return

    ! 3) adjust draws
    xdraw = xdraw - xhat
    xshockdraw = xshockdraw - xshockhat

  END SUBROUTINE samplerA3B3C3nanscalar
  ! @\newpage\subsection{samplerA3B3C3noisenan}@
  SUBROUTINE samplerA3B3C3noisenan(xdraw,xshockdraw,ynoisedraw,y,yNaN,T,Ny,Nx,Nw,A,B,C,noisevol,Ex0,sqrtVx0,VSLstream,status)

    INTENT(OUT) :: xdraw,xshockdraw,ynoisedraw,status
    INTENT(IN) :: y,Ny,Nx,Nw,T,A,B,C,Ex0,sqrtVx0
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: T,Ny,Nx,Nw,status

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C
    DOUBLE PRECISION, DIMENSION(Ny,T)    :: noisevol

    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, xNull
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0, Vx0, SigmaStarT

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xdraw, xhat
    DOUBLE PRECISION, DIMENSION(Nx,T)   :: xshockdraw, xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y, ydraw, ynoisedraw, ynoisehat
    LOGICAL, DIMENSION(Ny,1:T) :: yNaN


    ! 1) generate plus data
    call simA3B3C3noise(ydraw,xdraw,ynoisedraw,xshockdraw,T,Ny,Nx,Nw,A,B,C,noisevol,Ex0,sqrtVx0,VSLstream)
    ydraw = ydraw - y

    ! 2) filter the difference
    Vx0 = 0.0d0
    call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,Vx0,Nx)

    ! mean adjustment since projecting onto ydraw - y
    xNull = 0.0d0
    call disturbancesmootherA3B3C3nanscalar(xhat,xshockhat,SigmaStarT,ynoisehat,ydraw,yNaN,T,Ny,Nx,Nw,A,B,C,noisevol,xNull,Vx0,status)

    if (status /= 0) return

    ! 3) adjust draws
    xdraw = xdraw - xhat
    xshockdraw = xshockdraw - xshockhat
    ynoisedraw = ynoisedraw - ynoisehat

  END SUBROUTINE samplerA3B3C3noisenan

  ! @\newpage\subsection{samplerA3BSVC3nan}@
  SUBROUTINE samplerA3BSVC3nan(xdraw,xshockdraw,y,yNaN,T,Ny,Nx,Nw,A,B,C,SVol,Ex0,sqrtVx0,VSLstream,status)

    ! same as beta, just calling scalar smoother (with zero noise)

    INTENT(OUT) :: xdraw,xshockdraw,status
    INTENT(IN) :: y,yNaN,Ny,Nx,Nw,T,A,B,C,Ex0,sqrtVx0,SVol
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: T,Ny,Nx,Nw,j,k,status

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C

    DOUBLE PRECISION :: BB(Nx,Nw,T), dummySigma(Nx,Nx)


    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, xNull
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0, Vx0

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xdraw, xhat
    DOUBLE PRECISION, DIMENSION(Nw,T)   :: SVol
    DOUBLE PRECISION, DIMENSION(Nx,T)   :: xshockdraw, xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y, ydraw, dummyNoise, dummyNoiseVar
    LOGICAL, DIMENSION(Ny,1:T) :: yNaN

    FORALL (k=1:T,j=1:Nw) BB(:,j,k) = B(:,j) * SVol(j,k)

    ! 1) generate plus data
    call simA3B3C3(ydraw,xdraw,xshockdraw,T,Ny,Nx,Nw,A,BB,C,Ex0,sqrtVx0,VSLstream)
    ydraw = ydraw - y

    ! 2) filter the difference
    Vx0 = 0.0d0
    call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,Vx0,Nx)
    ! mean adjustment since projecting onto ydraw - y
    xNull = 0.0d0

    dummyNoiseVar = 0.0d0
    call disturbancesmootherA3B3C3nanscalar(xhat,xshockhat,dummySigma, dummynoise,ydraw,yNaN,T,Ny,Nx,Nw,A,BB,C,dummyNoiseVar,xNull,Vx0,status)

    if (status == 0) then
       ! 3) adjust draws
       xdraw = xdraw - xhat
       xshockdraw = xshockdraw - xshockhat
    end if

  END SUBROUTINE samplerA3BSVC3nan

  ! @\newpage\subsection{samplerA3B3C3noise}@
  SUBROUTINE samplerA3B3C3noise(xdraw,xshockdraw,ynoisedraw,y,T,Ny,Nx,Nw,A,B,C,noiseVol,Ex0,sqrtVx0,VSLstream,status)

    ! same as A3B3C3 but with noise shocks attached to each measurement
    ! noise is mutually uncorrelated (but possibly time-varying) iid, noiseVol is thus Ny x T

    INTENT(OUT) :: xdraw,xshockdraw,ynoisedraw,status
    INTENT(IN) :: y,Ny,Nx,Nw,T,A,B,C,Ex0,sqrtVx0,noiseVol
    INTENT(INOUT) :: VSLstream

    TYPE (vsl_stream_state) :: VSLstream
    INTEGER :: T,Ny,Nx,Nw,status

    DOUBLE PRECISION, DIMENSION(Nx,Nx,T) :: A
    DOUBLE PRECISION, DIMENSION(Nx,Nw,T) :: B
    DOUBLE PRECISION, DIMENSION(Ny,Nx,T) :: C
    DOUBLE PRECISION, DIMENSION(Ny,T) :: noiseVol

    DOUBLE PRECISION, DIMENSION(Nx) :: Ex0, xNull
    DOUBLE PRECISION, DIMENSION(Nx,NX) :: sqrtVx0, Vx0

    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: xdraw, xhat
    DOUBLE PRECISION, DIMENSION(Nx,T)   :: xshockdraw, xshockhat
    DOUBLE PRECISION, DIMENSION(Ny,1:T) :: y, ydraw, ynoisedraw, ynoisehat

    ! debug helper
    ! DOUBLE PRECISION :: shock(Nx,1:T)
    ! INTEGER :: j

    ! 1) generate plus data
    call simA3B3C3noise(ydraw,xdraw,ynoisedraw,xshockdraw,T,Ny,Nx,Nw,A,B,C,noiseVol,Ex0,sqrtVx0,VSLstream)
    ydraw = ydraw - y

    ! 2) filter the difference
    Vx0 = 0.0d0
    call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,Vx0,Nx)
    ! mean adjustment since projecting onto ydraw - y
    xNull = 0.0d0

    ! call disturbancesmootherA3B3C3noise(xhat,xshockhat,ynoisehat,ydraw,T,Ny,Nx,Nw,A,B,C,noiseVol,xNull,Vx0,status)

    call disturbancesmootherA3B3C3scalar(xhat,xshockhat,ynoisehat,ydraw,T,Ny,Nx,Nw,A,B,C,noiseVol ** 2,xNull,Vx0,status)



    ! delta = maxval(abs(xhat - xhat2))
    ! IF (delta .gt. 1.0d-8) PRINT *, 'delta = ', delta, '(xhat)'
    ! IF (delta .gt. 1.0d-7) THEN
    !    call savemat(xhat, 'xhat.debug')
    !    call savemat(xhat2, 'xhat2.debug')
    !    call savemat(ydraw, 'y.debug')
    !    call savemat(A(:,:,1), 'a.debug')
    !    call savemat(B(:,:,1), 'b.debug')
    !    call savemat(C(:,:,1), 'c.debug')
    !    call savemat(noiseVol, 'noiseVol.debug')
    !    call savemat(Vx0, 'Vx0.debug')

    !    stop 9
    ! END IF

    ! delta = maxval(abs(xshockhat - xshockhat2))
    ! IF (delta .gt. 1.0d-8) PRINT *, 'delta = ', delta, '(xshockhat)'

    ! delta = maxval(abs(ynoisehat - ynoisehat2))
    ! IF (delta .gt. 1.0d-8) PRINT *, 'delta = ', delta, '(ynoisehat)'

    if (status == 0) then
       ! 3) adjust draws
       xdraw      = xdraw      - xhat
       xshockdraw = xshockdraw - xshockhat
       ynoisedraw = ynoisedraw - ynoisehat
    end if


    ! check shockdraw
    ! do j=1,T
    !    shock(:,j) = xdraw(:,j)
    !    call dgemv('N', Nx,Nx, -1.0d0, A(:,:,j), Nx, xdraw(:,j-1), 1, 1.0d0, shock(:,j), 1)
    ! end do
    ! call savemat(shock, 'shock.debug')
    ! call savemat(xshockdraw, 'xshockdraw.debug')


  END SUBROUTINE samplerA3B3C3noise

  ! @\newpage\subsection{A3B3C3kalmanlike}@
  SUBROUTINE A3B3C3kalmanlike(loglike, T, Ny, y, yNaN, Nx, Nw, A, B, C, Ex0, sqrtVx0)

    ! sqrtVx0 is expected as lower triangular left factor (as in other routiens above)

    USE embox, only : savemat, savevec, int2str, mean, hrulefill

    use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery

    IMPLICIT NONE

    double precision, intent(out) :: loglike
    INTEGER, INTENT(IN) :: T, Ny, Nx, Nw
    DOUBLE PRECISION, INTENT(IN), DIMENSION(Ny,T) :: y
    LOGICAL, INTENT(IN), DIMENSION(Ny,T) :: yNaN
    DOUBLE PRECISION, INTENT(IN) :: Ex0(Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nw,T), C(Ny,Nx,T)  

    ! indices for looping and other integers
    INTEGER :: J, I
    INTEGER :: Nynonan


    ! Kalman state space objects
    DOUBLE PRECISION, DIMENSION(Nx) :: xposterior, xprior
    DOUBLE PRECISION ::  ytilde(Ny), logdetSigmaY, llf

    ! SQRT objects
    DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny) 
    DOUBLE PRECISION :: Kgain(Nx,Ny) 

    ! qr
    DOUBLE PRECISION :: qrR(Ny+Nx+Nw,Ny+Nx)
    INTEGER :: qrLwork


    ! init xposterior at prior
    xposterior      = Ex0
    sqrtSigmaX      = transpose(sqrtVx0)

    ! workspace query for qr decomposition
    qrR = 0.0d0
    qrlwork = qrquery(qrR)

    ! CALL initprogressbar(timer, 15.0d0)
    loglike = 0
    DO j=1,T

       Nynonan = count(.not. yNaN(:,j))

       if (Nynonan > 0) then

          ! xprior = A * xposterior(-1)
          xprior      = 0.0d0
          call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior,1,0.0d0,xprior,1)

          ! ------------------------------------------------------------------------
          ! SQRT KALMAN
          ! ------------------------------------------------------------------------

          ! fill directly into qrR
          qrR = 0.0d0
          ! qrR(1:Ny,1:Ny) = transpose(sqrtR)
          qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(B(:,:,j))
          ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
          call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A(:,:,j),Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
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

          ! singularity fix: insert unit dummies for missing values
          ! as long as sqrtSigmaY not used any further, need not undo these dummies ...
          do i=1,Ny
             if (yNaN(i,j)) sqrtSigmaY(i,i) = 1.0d0
             ! ytilde(i) = 0.0d0 -- not needed since y(i)=0 and C(i) * xprior also zero
          end do

          logdetSigmaY = 2.0d0 * sum( (/ (log(abs(sqrtSigmaY(i,i))), i = 1, Ny) /) )

          ! rotate/normalize ytilde (up to sign, consistent with rotation of K)
          call dtrsv('U', 'T', 'N', Ny, sqrtSigmaY, Ny, ytilde, 1) ! recall: sqrtSigmaY is returned as upper triangular, right factor

          ! xposterior = xprior + K * ytilde
          xposterior = xprior
          call DGEMV('N',Nx,Ny,1.0d0,Kgain,Nx,ytilde,1,1.0d0,xposterior,1)

          ! compute log-likelihood
          ! llf
          llf  = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

          ! ------------------------------------------------------------------------
          ! DONE: SQRT KALMAN
          ! ------------------------------------------------------------------------

          loglike = loglike + llf

       end if ! Nynonan > 0
    END DO ! j=1,T


  END SUBROUTINE A3B3C3kalmanlike


END MODULE statespacebox


