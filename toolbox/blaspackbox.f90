MODULE blaspack

  use embox , only : savemat, savevec, pi, logtwopi

  IMPLICIT NONE

CONTAINS

  ! @\newpage\subsection{vcv}@
  FUNCTION vcv(Nobs,Nvar,X) result(S)

    integer, intent(in) :: Nobs,Nvar
    double precision, dimension(Nobs,Nvar), intent(in) :: X
    double precision, dimension(Nobs,Nvar) :: resid
    double precision, dimension(Nvar) :: xbar
    double precision, dimension(Nvar,Nvar) :: S
    integer :: i,j
    
    xbar = sum(x,1) / dble(Nobs)
    forall (i=1:Nobs,j=1:Nvar) resid(i,j) = x(i,j) - xbar(j)

    call DSYRK('u','t',Nvar,Nobs,1.0d0 / dble(Nobs),resid,Nobs,0.0d0,S,Nvar)
    
  END FUNCTION vcv

  ! @\newpage\subsection{sqrtvcv}@
  FUNCTION sqrtvcv(Nobs,Nvar,X) result(S)

    integer, intent(in) :: Nobs,Nvar
    double precision, dimension(Nobs,Nvar), intent(in) :: X
    double precision, dimension(Nobs,Nvar) :: resid
    double precision, dimension(Nvar) :: xbar
    double precision, dimension(Nvar,Nvar) :: S
    integer :: i,j,status

    S = 0.0d0
    xbar = sum(x,1) / dble(Nobs)
    forall (i=1:Nobs,j=1:Nvar) resid(i,j) = x(i,j) - xbar(j)

    call DSYRK('u','t',Nvar,Nobs,1.0d0 / dble(Nobs),resid,Nobs,0.0d0,S,Nvar)
    call DPOTRF('u', Nvar, S, Nvar, status)
    if (status /= 0) then
       write (*,*), 'DPOTRF ERROR ', status, ' [SQRTVCV]'
       stop 1
    END IF

  END FUNCTION sqrtvcv

  ! @\newpage\subsection{sqrtvcvTR}@
  FUNCTION sqrtvcvTR(Nvar,Nobs,X) result(S)

    integer, intent(in) :: Nobs,Nvar
    double precision, dimension(Nvar,Nobs), intent(in) :: X
    double precision, dimension(Nvar,Nobs) :: resid
    double precision, dimension(Nvar) :: xbar
    double precision, dimension(Nvar,Nvar) :: S
    integer :: i,j,status

    S = 0.0d0
    xbar = sum(x,2) / dble(Nobs)
    forall (i=1:Nobs,j=1:Nvar) resid(j,i) = x(j,i) - xbar(j)

    call DSYRK('u','n',Nvar,Nobs,1.0d0 / dble(Nobs),resid,Nvar,0.0d0,S,Nvar)
    call DPOTRF('u', Nvar, S, Nvar, status)
    if (status /= 0) then
       write (*,*), 'DPOTRF ERROR ', status, ' [SQRTVCVTR]'
       stop 1
    END IF

  END FUNCTION sqrtvcvTR
  
! @\newpage\subsection{symmetric}@
  PURE SUBROUTINE symmetric(S)
    ! ensures symmetry of S (assumung S has at least upper storage)
    IMPLICIT NONE
    INTENT(INOUT) :: S

    INTEGER :: j,i,N
    DOUBLE PRECISION, DIMENSION(:,:) :: S

    N = size(S,1)

       DO i=1,N
          DO j=1,i-1
             S(i,j) = S(j,i)
          END DO
       END DO
  END SUBROUTINE SYMMETRIC

! @\newpage\subsection{maxroot}@
  FUNCTION maxroot(A,n)
    INTENT(IN) :: A, n

    double precision :: maxroot

    integer :: N, status, lwork
    double precision, dimension(n,n) :: A, Awork
    double precision, dimension(n) :: lambdaR, lambdaI

    double precision :: dummy(1,1)
    double precision, allocatable, dimension(:) :: WORK

    Awork = A ! do not touch A 

    ! workspace query
    ALLOCATE (WORK(1))

    LWORK = -1
    call DGEEV('N', 'N', n , Awork, n, lambdaR, lambdaI, dummy, 1, dummy, 1, WORK, LWORK, status)
    IF (status /= 0) THEN
       WRITE(*,*) 'DGEEV error (LWORK QUERY)'
       STOP 1
    END IF

    LWORK = ceiling(WORK(1))

    DEALLOCATE(WORK)

    ! setup optimal workspace
    ALLOCATE(WORK(LWORK))
    ! compute eigenvalues
    call DGEEV('N', 'N', n , Awork, n, lambdaR, lambdaI, dummy, 1, dummy, 1, WORK, LWORK, status)
    IF (status /= 0) THEN
       WRITE(*,*) 'DGEEV error, status =', status
       print *, 'WORK:', WORK(1:5)
       print *, 'LWORK:', LWORK
       call savemat(Awork, 'debug.A.dat')
       STOP 1
    END IF

    maxroot      = maxval(abs(dcmplx(lambdaR, lambdaI)))
    ! wrap up
    DEALLOCATE(WORK)

  END FUNCTION maxroot
! @\newpage\subsection{predictstate}@
  FUNCTION predictstate(Nx, A, x0, horizon) result(xhat)
    ! xhat = A^horizons * xprior

    INTENT(IN)  :: Nx, A, x0, horizon

    INTEGER :: h,horizon,Nx

    DOUBLE PRECISION :: A(Nx,Nx), xprior(Nx), xhat(Nx), x0(Nx)


    ! FORECAST MEAN
    xhat=x0

    DO h = 1, horizon

       xprior = xhat

       ! update s
       CALL DGEMV('N', Nx, Nx, 1.0d0, A,  Nx, xprior, 1, 0.0d0, xhat, 1)

    END DO ! h


  END FUNCTION predictstate
! @\newpage\subsection{quadraticform}@
  PURE FUNCTION quadraticform(v, S, N) result(q)
    ! q = v' * S * v where S is symmetric (uppper triangular)

    INTENT(IN)  :: v, S, N

    INTEGER :: N,j,i

    DOUBLE PRECISION :: v(N), S(N,N), q

    q = 0.0d0
    ! off diagonals
    DO j = 2,N
       DO i=1,j-1
          q = q + v(i) * v(j) * S(i,j)
       END DO
    END DO
    q = 2.0d0 * q
    ! diagonals
    DO i=1,N 
       q = q + (v(i) ** 2) * S(i,i)
    END DO
    
  END FUNCTION quadraticform
! @\newpage\subsection{sandwich}@
  SUBROUTINE sandwich(ASA, A, LDA, S, LDS)
    ! ASA = A * S * A' 
    ! input S is upper-triangular symmetric 
    ! output ASA is dense-symmetric
    INTENT(OUT) :: ASA
    INTENT(IN)  :: A, LDA, S, LDS

    INTEGER :: LDA, LDS

    DOUBLE PRECISION :: A(LDA, LDS), S(LDS,LDS), ASA(LDA,LDA), AS(LDA,LDS)

    AS  = 0.0d0
    ASA = 0.0d0

    call DSYMM('R','U', LDA, LDS,1.0d0,S,LDS,A,LDA,0.0d0,AS,LDA)
    call DGEMM('N','T',LDA,LDA,LDS,1.0d0,AS,LDA,A,LDA,0.0d0,ASA,LDA)

  END SUBROUTINE sandwich

! @\newpage\subsection{sandwichplus}@
  SUBROUTINE sandwichplus(ASA, A, LDA, S, LDS)
    ! ASA = A * S * A' + ASA where S is symmetric (uppper triangular)
    ! input S is upper-triangular symmetric 
    ! input/output ASA has to be dense-symmetric !!

    INTENT(INOUT) :: ASA
    INTENT(IN)  :: A, LDA, S, LDS

    INTEGER :: LDA, LDS

    DOUBLE PRECISION :: A(LDA, LDS), S(LDS,LDS), ASA(LDA,LDA), AS(LDA,LDS)

    AS  = 0.0d0

    call DSYMM('R','U', LDA, LDS,1.0d0,S,LDS,A,LDA,0.0d0,AS,LDA)
    call DGEMM('N','T',LDA,LDA,LDS,1.0d0,AS,LDA,A,LDA,1.0d0,ASA,LDA)

  END SUBROUTINE sandwichplus

! @\newpage\subsection{vectorTimesMatrix}@
  SUBROUTINE vectorTimesMatrix(rows, cols, A, X, Y)
    ! computes y = x' A

    INTENT(IN) :: rows, cols, A, X
    INTENT(INOUT)  :: Y

    INTEGER :: rows, cols
    DOUBLE PRECISION :: A(rows,cols), x(rows), y(cols)


    call DGEMV('T', rows, cols, 1.0d0, A, rows, X, 1, 0.0d0, Y, 1)

  END SUBROUTINE vectorTimesMatrix

! @\newpage\subsection{XprimeX}@
  SUBROUTINE XprimeX(XX, X)

    INTENT(OUT) :: XX
    INTENT(IN)  :: X

    INTEGER :: N, T

    DOUBLE PRECISION, DIMENSION(:,:) :: XX, X

    N = size(X,2)
    T = size(X,1)
    XX = 0.0d0 ! to clean out lower triangular part of XX
    call DSYRK('U','T',N,T,1.0d0,X,T,0.0d0,XX,N)

  END SUBROUTINE XprimeX

! @\newpage\subsection{XXprime}@
  SUBROUTINE XXprime(XX, X)

    INTENT(OUT) :: XX
    INTENT(IN)  :: X

    INTEGER :: Ncols, Nrows

    DOUBLE PRECISION, DIMENSION(:,:) :: XX, X

    Ncols = size(X,2)
    Nrows = size(X,1)
    XX = 0.0d0 ! to clean out lower triangular part of XX
    call DSYRK('U','N',Nrows,Ncols,1.0d0,X,Nrows,0.0d0,XX,Nrows)

  END SUBROUTINE XXprime


! @\newpage\subsection{invsym}@
  SUBROUTINE invsym(xx)
    ! inverts p.d. symmetric real matrix, assuming upper triangular storage

    INTENT(INOUT) :: XX

    INTEGER :: n, info
    DOUBLE PRECISION, DIMENSION(:,:) :: XX

    n = size(xx,1)
    call DPOTRF('U', n, XX, n, info )
    IF (info /= 0) THEN
       write(*,*) "DPOTRF ERROR:", INFO, "[INVSYM]"
       STOP 1
    END IF
    
    call DPOTRI('U', n, XX, n, info )
    IF (info /= 0) THEN
       write(*,*) "DPOTRI ERROR:", INFO, "[INVSYM]"
       STOP 1
    END IF

  END SUBROUTINE invsym


! @\newpage\subsection{symkronecker}@
  SUBROUTINE symkronecker(alpha,A,Na,B,Nb,beta,C)
    ! C = kron(A,B) + C
    ! assumes symmetry A,B and C 
    ! notice: A can be upper triangular, but B must be full storage

    INTENT(OUT) :: C
    INTENT(IN) :: A,Na,B,Nb,alpha,beta

    INTEGER :: Na,Nb, i, j

    DOUBLE PRECISION :: A(Na,Na), B(Nb,Nb), C(Na * Nb, Na * Nb), alpha, beta

    ! loop over rows and columns of A
    DO j = 1 , Na
       FORALL(i=1:j) C((i-1) * Nb + 1 : i * Nb, (j-1) * Nb + 1 : j * Nb) = alpha * A(i,j) * B + beta * C((i-1) * Nb + 1 : i * Nb, (j-1) * Nb + 1 : j * Nb)
    END DO


  END SUBROUTINE symkronecker

! @\newpage\subsection{eye}@
  SUBROUTINE eye(I,alpha)
    ! identity matrix of order n, scaled by alpha (default = 1.0d0)

    INTEGER :: N
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(:,:) :: I
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: alpha
    DOUBLE PRECISION :: a
    INTEGER :: ii


    IF (.NOT. PRESENT(alpha)) THEN
       a = 1.0d0
    ELSE
       a = alpha 
    END IF

    N = size(I,1)
    I = 0.0d0
    FORALL (ii=1:N) I(ii,ii) = a

  END SUBROUTINE eye

! @\newpage\subsection{vec}@
  SUBROUTINE vec(v,x)
    ! v = vec(x)

    INTENT(OUT) :: v
    INTENT(IN) :: x

    INTEGER :: i,j,rows,cols

    DOUBLE PRECISION, DIMENSION(:,:) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    rows = size(x,1)
    cols = size(x,2)

    FORALL (i=1:rows,j=1:cols) v((j-1) * rows + i) = x(i,j)

  END SUBROUTINE vec

! @\newpage\subsection{ivech}@
  SUBROUTINE ivech(x,v)
    ! v = vech(x)
    ! assuming upper triangular storage 

    INTENT(IN) :: v
    INTENT(INOUT) :: x

    INTEGER :: i,j,n,s

    DOUBLE PRECISION, DIMENSION(:,:) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    n = size(x,1)
    s = 0
    
    x = 0.0d0
    DO j=1,n
       DO i = 1 ,j
          s = s + 1
          x(i,j) = v(s)
          ! x(j,i) = v(s)
       END DO
    END DO

  END SUBROUTINE ivech

  ! @\newpage\subsection{ivechU}@
  pure function ivechU(v,n) result(x)
    ! returns x s,t, v = vech(x)
    ! assuming upper triangular storage 

    INTENT(IN) :: v,n

    INTEGER :: row,col,n,s

    DOUBLE PRECISION, DIMENSION(n,n) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    s = 0
    
    x = 0.0d0
    DO col=1,n
       DO row = 1, col
          s = s + 1
          x(row,col) = v(s)
       END DO
    END DO

  END FUNCTION ivechU

! @\newpage\subsection{vechU}@
  pure function vechU(x,n) result(v)
    ! returns v = vech(x)
    ! assuming upper triangular storage 

    INTENT(IN) :: x,n

    INTEGER :: row,col,n,s

    DOUBLE PRECISION, DIMENSION(n,n) :: x
    DOUBLE PRECISION, DIMENSION(n * (n + 1) / 2) :: v

    s = 0
    
    v = 0.0d0
    DO col=1,n
       DO row = 1, col
          s = s + 1
          v(s) = x(row,col) 
       END DO
    END DO

  END FUNCTION vechU
! @\newpage\subsection{vech}@
  SUBROUTINE vech(v,x)
    ! v = vech(x)
    ! assuming upper triangular storage 

    INTENT(INOUT) :: v
    INTENT(IN) :: x

    INTEGER :: i,j,n,s

    DOUBLE PRECISION, DIMENSION(:,:) :: x
    DOUBLE PRECISION, DIMENSION(:) :: v

    n = size(x,1)
    s = 0
    
    DO j=1,n
       DO i = 1, j
          s = s + 1
          v(s) = x(i,j)
       END DO
    END DO

  END SUBROUTINE vech
! @\newpage\subsection{triu}@
  pure SUBROUTINE triu(x)
    ! zeros out lower triangular elements

    INTENT(INOUT) :: x

    INTEGER :: row,col,n

    DOUBLE PRECISION, DIMENSION(:,:) :: x

    n = size(x,1)

    DO col=1,n-1
       DO row = col+1,n
          x(row,col) = 0.0d0
       END DO
    END DO

  END SUBROUTINE triu

! @\newpage\subsection{choleski}@
  SUBROUTINE choleski(s)
    intent(inout) :: s

    integer :: n, i
    double precision :: s(:,:)

    n = size(s,1)

    ! factorize
    call dpotrf('u', n, s, n, i)

    ! check for errors
    if (i /= 0) then
       write(*,*) 'CHOLESKI ERROR:', i, '[CHOLESKI BLASPACKBOX]'
       stop 1
    end if

    ! zero out lower triangular
    forall (i = 1 : n-1) s(i+1:n,i) = 0.0d0


  END SUBROUTINE choleski

! @\newpage\subsection{qrecon}@
  FUNCTION qrecon(A,Nrows,Ncols) result(R)

    ! assuems that A(Nrows, Ncols) Nrows > Ncols

    integer, intent(in) :: Nrows, Ncols
    double precision, intent(in), dimension(Nrows,Ncols) :: A
    double precision, dimension(Nrows,Ncols) :: RR
    double precision, dimension(Ncols,Ncols) :: R
    integer :: lwork, i, errcode
    double precision, dimension(1) :: workquery
    double precision, allocatable, dimension(:) :: work
    double precision, dimension(Ncols) :: qrReflectors


    ! step 1: workspace query
    RR = A ! just making sure that lapack does not get to overwrite the input A
    call dgeqrf(Nrows, Ncols, RR, Nrows, qrReflectors, workquery, -1, errcode)
    if (errcode .ne. 0) then
       print *,'QRECON workspace query failed errocde=', errcode
       stop 1
    end if
    LWORK = ceiling(workquery(1))

    ! step 2: allocate optimal workspace
    allocate(work(LWORK))
    
    ! step 3: perform QR
    RR = A
    call dgeqrf(Nrows, Ncols, RR, Nrows, qrReflectors, work, lwork, errcode)
    if (errcode .ne. 0) then
       print *,'QRECON decomp failed errocde=', errcode
       stop 1
    end if
    
    deallocate(work)
    
    ! step 4: copy econ sized RR into R
    R = 0.0d0
    forall (i=1:Ncols) R(1:i,i) = RR(1:i,i)

    ! call savemat(R, 'R.debug')
    ! call savemat(RR, 'RR.debug')
    ! call savemat(A, 'A.debug')
    ! stop 11

  END FUNCTION qrecon

  ! @\newpage\subsection{qrot}@
  SUBROUTINE qrot(R,lwork)

    ! assuems that m(Nrows, Ncols) with Nrows > Ncols

    intent(inout) :: R
    intent(in) :: lwork

    double precision, dimension(:,:) :: R
    integer :: lwork, Nrows, Ncols, errcode, i
    double precision, dimension(lwork) :: work
    double precision, dimension(size(R,2)) :: qrReflectors

    Nrows = size(R,1)
    Ncols = size(R,2)
    
    call dgeqrf(Nrows, Ncols, R, Nrows, qrReflectors, work, lwork, errcode)
    if (errcode .ne. 0) then
       print *,'QROT failed errocde=', errcode
       stop 1
    end if

    ! zero out lower triangular part of R
    forall (i=1:Ncols) R(i+1:Nrows,i) = 0.0d0


  END SUBROUTINE qrot
 
! @\newpage\subsection{qrquery}@
  FUNCTION qrquery(m) result(lwork)

    ! assumes that m(Nrows, Ncols) with Nrows > Ncols

    intent(inout) :: m
    
    double precision, dimension(:,:) :: m
    integer :: lwork, Nrows, Ncols, errcode
    double precision, dimension(1) :: work
    double precision, dimension(size(m,2)) :: qrReflectors

    Nrows = size(M,1)
    Ncols = size(M,2)
    
    call dgeqrf(Nrows,Ncols,m,Nrows, qrReflectors, work, -1, errcode)
    if (errcode .ne. 0) then
       print *,'QRQUERY failed errocde=', errcode
       stop 1
    end if

    LWORK = ceiling(work(1))


  END FUNCTION qrquery
END MODULE blaspack
