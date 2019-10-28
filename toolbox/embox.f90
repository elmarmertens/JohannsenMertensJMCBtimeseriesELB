MODULE embox

  USE vslbox
  USE omp_lib
  
  IMPLICIT NONE


  DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793115997963468544d0, logtwopi = log(2.0d0 * pi), logpi = log(pi),  sqrttwopi = sqrt(2.0d0 * pi), invsqrttwopi = 1 / sqrt(2.0d0 * pi)

  INTEGER, parameter :: Nquantiles = 10
  REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)


CONTAINS


  ! @\newpage\subsection{loadmatrix}@
  SUBROUTINE loadmatrix(A,filename,rows,cols)
    ! loads rows x cols matrix A from filename
    ! assumes vectorized storage in filenamae

    INTENT(IN) :: filename,rows,cols
    INTENT(OUT) :: A
    CHARACTER (LEN=200) :: filename

    INTEGER r,rows, cols
    DOUBLE PRECISION, DIMENSION(rows,cols) :: A
    INTEGER :: nunit
    
    !Open File for reading
    OPEN (newunit=nunit, FILE=filename, STATUS='OLD', ACTION='READ')
    DO r=1,rows
       READ(nunit, es30d16(cols)) A(r,:)
    END DO
    CLOSE(UNIT=nunit)

  END SUBROUTINE loadmatrix


  ! @\newpage\subsection{loadarray1}@
  SUBROUTINE loadarray1(x,filename,n)
    ! loads n-element vector x from filename

    INTENT(IN) :: filename,n
    INTENT(OUT) :: x
    CHARACTER (LEN=200) :: filename

    INTEGER j,n
    DOUBLE PRECISION, DIMENSION(n) :: x
    INTEGER :: nunit
    
    !Open File for reading
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='OLD', ACTION='READ')
    DO j=1,n
       READ(nunit,'(ES30.16)') x(j)
    END DO

    CLOSE(UNIT=nunit)

  END SUBROUTINE loadarray1
  
  ! @\newpage\subsection{loft}@
  FUNCTION loft(filename)
    INTENT(IN) :: filename


    character(len=200) :: filename, dummy
    integer :: loft,status 
    INTEGER :: nunit
    
    open(newunit=nunit,file=filename)
    loft = 0
    do
       read(nunit,*,IOSTAT=status) dummy
       IF (status /=0) EXIT
       loft = loft + 1
    end do

    close(unit=nunit)

  END FUNCTION loft

  ! @\newpage\subsection{mean}@
  PURE FUNCTION mean(x) result(m)
    INTENT(IN) :: x

    double precision, dimension(:) :: x
    double precision :: m

    m = sum(x) / size(x)

  END FUNCTION mean

  ! @\newpage\subsection{variance}@
  PURE FUNCTION variance(x) result(v)
    INTENT(IN) :: x

    double precision, dimension(:) :: x
    double precision :: v
    integer :: N

    N = size(x)
    v = sum((x - (sum(x) / N)) ** 2) / N

  END FUNCTION variance

  ! @\newpage\subsection{normpdf}@
  ELEMENTAL FUNCTION normpdf(x) 

    INTENT(IN) :: x

    double precision :: x, normpdf

    normpdf = exp(-0.5d0 * (logtwopi + x ** 2))

  END FUNCTION normpdf

  ! @\newpage\subsection{timestampstr}@
  FUNCTION timestampstr()
    character(18) :: timestampstr
    character(8)  :: date
    character(10) :: time
    call date_and_time(date,time)

    timestampstr = date//'d'//time(1:4)//'h'//time(5:6)//time(8:9)
  END FUNCTION timestampstr

  ! @\newpage\subsection{timestamp}@
  FUNCTION timestamp()
    INTEGER :: timestamp(8)
    call date_and_time(VALUES=timestamp)

  END FUNCTION timestamp

  ! @\newpage\subsection{HRULEFILL}@
  SUBROUTINE HRULEFILL

    WRITE (*,*) '---------------------------------------------------------------------------'

  END SUBROUTINE HRULEFILL

  ! @\newpage\subsection{display}@
  SUBROUTINE display(x,fs)
    IMPLICIT NONE

    INTENT(IN) :: x

    INTEGER rows,cols,j
    CHARACTER (LEN=*) :: fs
    DOUBLE PRECISION, DIMENSION(:,:) :: x

    cols = size(x,2)
    rows = size(x,1)

    WRITE(*,fmtstr(fs,cols)) (x(j,:), j=1,rows)

  END SUBROUTINE display

  ! @\newpage\subsection{displayvec}@
  SUBROUTINE displayvec(x,fs)
    IMPLICIT NONE

    INTENT(IN) :: x

    INTEGER rows,j
    CHARACTER (LEN=*) :: fs
    DOUBLE PRECISION, DIMENSION(:) :: x

    rows = size(x,1)

    WRITE(*,fmtstr(fs,1)) (x(j), j=1,rows)

  END SUBROUTINE displayvec

  ! @\newpage\subsection{printmat}@
  SUBROUTINE printmat(x)
    IMPLICIT NONE

    INTENT(IN) :: x
    INTEGER rows,cols,j
    CHARACTER (LEN=10) :: fs
    DOUBLE PRECISION, DIMENSION(:,:) :: x

    fs ='e12.4'

    cols = size(x,2)
    rows = size(x,1)

    WRITE(*,fmtstr(fs,cols)) (x(j,:), j=1,rows)

  END SUBROUTINE printmat

  ! @\newpage\subsection{printvec}@
  SUBROUTINE printvec(x)
    IMPLICIT NONE

    INTENT(IN) :: x
    INTEGER rows,j
    CHARACTER (LEN=200) :: fs
    DOUBLE PRECISION, DIMENSION(:) :: x

    fs ='f20.10'

    rows = size(x,1)

    WRITE(*,fmtstr(fs,1)) (x(j), j=1,rows)

  END SUBROUTINE printvec

  ! @\newpage\subsection{int2str}@
  FUNCTION int2str(n) result(str)
    INTENT(IN) :: n
    INTEGER :: n
    CHARACTER (LEN=200) :: str
    WRITE(str,'(i10)') n
    str = trim(adjustl(str))
  END FUNCTION int2str

  ! @\newpage\subsection{es30d16}@
  FUNCTION es30d16(N) result(fmtstr)
    INTENT(IN) :: N
    INTEGER :: n
    CHARACTER (LEN=10000) :: fmtstr
    fmtstr = '(ES30.16' // repeat(',ES30.16', N-1) // ')'
  END FUNCTION es30d16

  ! @\newpage\subsection{nummer84}@
  FUNCTION nummer84(N) result(fmtstr)
    INTENT(IN) :: N
    INTEGER :: n
    CHARACTER (LEN=10000) :: fmtstr
    fmtstr = '(f8.4' // repeat(',f8.4', N-1) // ')'
  END FUNCTION nummer84

  ! @\newpage\subsection{fmtstr}@
  FUNCTION fmtstr(str,N) 
    INTENT(IN) :: N,str
    INTEGER :: n
    CHARACTER (LEN=10000) :: fmtstr ! todo: make allocatable as function of len(str) * N plus change
    CHARACTER (LEN=*) :: str
    fmtstr = '(' // trim(str) // repeat(',' // trim(str), N-1) // ')'
  END FUNCTION fmtstr

  ! @\newpage\subsection{savemat}@
  SUBROUTINE savemat (x,filename)
    IMPLICIT NONE

    INTEGER :: rows, cols
    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: x(:,:)
    INTEGER j
    INTEGER :: nunit

    rows = size(x,1)
    cols = size(x,2)

    ! WRITE (*,*) 'THIS IS SAVEMAT WITH FILENAME ', filename

    !Open File for writing
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(nunit,'(ES30.16' // repeat(',ES30.16', cols-1) // ')') (x(j,:), j=1,rows)
    CLOSE(UNIT=nunit)

  END SUBROUTINE savemat

  ! @\newpage\subsection{savearray3}@
  SUBROUTINE savearray3(x,label,filext)
    IMPLICIT NONE

    INTEGER :: rows, cols, obs
    CHARACTER (LEN=*), INTENT(IN) :: label, filext
    DOUBLE PRECISION, INTENT(IN) :: x(:,:,:)
    CHARACTER (LEN=200) :: filename
    INTEGER jj,ii
    INTEGER :: nunit

    rows = size(x,1)
    cols = size(x,2)
    obs  = size(x,3)

    DO ii=1,obs

       filename = trim(label) // trim(int2str(ii)) // "." // trim(filext)

       !Open File for writing
       OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
       WRITE(nunit,'(ES30.16' // repeat(',ES30.16', cols-1) // ')') (x(jj,:,ii), jj=1,rows)
       CLOSE(UNIT=nunit)

    END DO

  END SUBROUTINE savearray3

  ! @\newpage\subsection{savematlogical}@
  SUBROUTINE savematlogical(x,filename)
    IMPLICIT NONE

    INTEGER :: rows, cols
    CHARACTER (LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: x(:,:)
    INTEGER j
    INTEGER :: nunit
    
    rows = size(x,1)
    cols = size(x,2)

    !Open File for writing
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(nunit,'(I5' // repeat(',I5', cols-1) // ')') (merge(1,0,x(j,:)), j=1,rows)
    CLOSE(UNIT=nunit)

  END SUBROUTINE savematlogical

  ! @\newpage\subsection{savevec}@
  SUBROUTINE savevec (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y(:)
    INTEGER j
    INTEGER :: nunit
    
    !Open File for writing
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    DO j=1,size(y)
       WRITE(nunit,'(ES30.16)') y(j)
    END DO
    CLOSE(UNIT=nunit)

  END SUBROUTINE savevec

  ! @\newpage\subsection{savevecX}@
  SUBROUTINE savevecX (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y(:)
    INTEGER j
    INTEGER :: nunit
    
    !Open File for writing
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    DO j=1,size(y)
       WRITE(nunit,'(ES40.16E3)') y(j)
    END DO
    CLOSE(UNIT=nunit)

  END SUBROUTINE savevecX

  ! @\newpage\subsection{savescalar}@
  SUBROUTINE savescalar (y,filename)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    DOUBLE PRECISION, INTENT(IN) :: y
    INTEGER :: nunit
    !Open File for writing
    OPEN (NEWUNIT=nunit, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(nunit,'(ES30.16)') y
    CLOSE(UNIT=nunit)

  END SUBROUTINE savescalar

  ! @\newpage\subsection{quantilevec}@
  FUNCTION quantilevec(x) result(q)
    ! store mean, median and quantiles of draws into file

    IMPLICIT NONE

    INTENT(IN)    :: x ! NOTE: will be INOUT because of sorting

    ! INTEGER, parameter :: Nquantiles = 10
    ! REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nquantiles) :: fracndx 
    INTEGER :: Ndraws, n,  status
    DOUBLE PRECISION, DIMENSION(:) :: x
    DOUBLE PRECISION :: mittel, median
    DOUBLE PRECISION :: q(2 + Nquantiles)


    Ndraws = size(x)
    CALL dlasrt('I', Ndraws, x, status)
    IF (status /= 0) THEN
       write (*,*), 'DLASORT ERROR ', status, ' [GET QUANTILES]'
       stop 1
    END IF
    ! generate index for each fractile
    FORALL (n = 1:Nquantiles)
       fracndx(n)  = floor(real(Ndraws) * quantilegrid(n))
    END FORALL
    where (fracndx == 0) fracndx = 1

    ! compute results
    mittel    = sum(x) / dble(Ndraws)
    median    = x(floor(real(Ndraws) * 0.5))

    ! collect in results vector
    ! write to file
    q(1) = mittel
    q(2) = median
    forall(n=1:Nquantiles) q(2+n) = x(fracndx(n))

  END FUNCTION quantilevec

  ! @\newpage\subsection{median}@
  FUNCTION median(x) result(m)
    ! computes median

    IMPLICIT NONE

    INTENT(IN)    :: x ! NOTE: will be INOUT because of sorting

    INTEGER :: Ndraws, status
    DOUBLE PRECISION, DIMENSION(:) :: x
    DOUBLE PRECISION :: m
    
    Ndraws = size(x)
    CALL dlasrt('I', Ndraws, x, status)
    IF (status /= 0) THEN
       write (*,*), 'DLASORT ERROR ', status, ' [MEDIAN]'
       stop 1
    END IF

    m  = x(floor(real(Ndraws) * 0.5))

  END FUNCTION median
  
  ! @\newpage\subsection{storeEstimates}@
  SUBROUTINE storeEstimates(theta,Ntheta,Ndraws,filename)
    ! store mean, median and quantiles of draws into file

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! has to be INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename

    ! INTEGER, parameter :: Nquantiles = 10
    ! REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nquantiles) :: fracndx 
    INTEGER :: Ntheta, Ndraws, n, k, j, status
    DOUBLE PRECISION, DIMENSION(Ntheta,Ndraws) :: theta
    DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median

    CHARACTER (LEN=200) :: filename

    ! sort draws
    DO n = 1, Ntheta
       CALL dlasrt('I', Ndraws, theta(n,:), status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF
    END DO

    ! generate index for each fractile
    FORALL (n = 1:Nquantiles)
       fracndx(n)  = floor(real(Ndraws) * quantilegrid(n))
    END FORALL
    where (fracndx == 0) fracndx = 1

    ! compute results
    mittel    = sum(theta, 2) / dble(Ndraws)
    median    = theta(:,floor(real(Ndraws) * 0.5))

    ! write to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), median(j), (theta(j,fracndx(k)), k = 1,Nquantiles), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeEstimates

  ! @\newpage\subsection{storeEstimatesTranspose}@
  SUBROUTINE storeEstimatesTranspose(theta,Ntheta,Ndraws,filename)
    ! store mean, median and quantiles of draws into file
    ! transpose version of storeEstimates, theta is Ndraws, Ntheta, but will be stored as Ntheta x 12

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! has to be INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename

    ! INTEGER, parameter :: Nquantiles = 10
    ! REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nquantiles) :: fracndx 
    INTEGER :: Ntheta, Ndraws, n, k, j, status
    DOUBLE PRECISION, DIMENSION(Ndraws,Ntheta) :: theta
    DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median

    CHARACTER (LEN=200) :: filename

    ! sort draws
    DO n = 1, Ntheta
       CALL dlasrt('I', Ndraws, theta(:,n), status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF
    END DO

    ! generate index for each fractile
    FORALL (n = 1:Nquantiles)
       fracndx(n)  = floor(real(Ndraws) * quantilegrid(n))
    END FORALL
    where (fracndx == 0) fracndx = 1

    ! compute results
    mittel    = sum(theta, 1) / dble(Ndraws)
    median    = theta(floor(real(Ndraws) * 0.5),:)

    ! write to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), median(j), (theta(fracndx(k),j), k = 1,Nquantiles), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeEstimatesTranspose
  ! -----------------------------------------------------------------

    ! @\newpage\subsection{storeMaskedEstimatesTranspose}@
  SUBROUTINE storeMaskedEstimatesTranspose(theta,Ntheta,Ndraws,filename,mask)
    ! store mean, median and quantiles of draws into file
    ! transpose version of storeEstimates, theta is Ndraws, Ntheta, but will be stored as Ntheta x 12

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename,mask

    ! INTEGER, parameter :: Nquantiles = 10
    ! REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nquantiles) :: fracndx 
    INTEGER :: Ntheta, Ndraws, Nslicedraws, n, j, status
    DOUBLE PRECISION, DIMENSION(Ndraws,Ntheta) :: theta
    LOGICAL, DIMENSION(Ndraws,Ntheta) :: mask
    ! DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median
    DOUBLE PRECISION, DIMENSION(Ntheta,2+Nquantiles) :: outputtable
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: maskedslice
    
    CHARACTER (LEN=200) :: filename


    ! note: the explicit shape arguments Ndraws and Ntheta could be dropped from this subroutine; I have kept them for better comparability with the other routines (and since it helps catching bugs arising from passing badly shaped inputs to the routine)
    
    ! operate on each column of theta separately, since each column may have a different number of elements that satisfy the mask 

    ! sort draws
    DO j = 1, Ntheta

       ! maskedslice = theta(:,j)
       Nslicedraws = count(mask(:,j))
       allocate (maskedslice(Nslicedraws))
       maskedslice = pack(theta(:,j), mask(:,j)) ! works fine using auto-allocate; to be safe, I allocate explicitly, however
       
       ! Nslicedraws = size(maskedslice)
       ! print *, 'Nscliced: ', Nslicedraws

       ! generate index for each fractile
       FORALL (n = 1:Nquantiles) fracndx(n)  = floor(real(Nslicedraws) * quantilegrid(n))
       where (fracndx == 0) fracndx = 1

       ! sort
       CALL dlasrt('I', Nslicedraws, maskedslice, status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF

       ! collect into output table
       outputtable(j,1) = sum(maskedslice) / dble(Nslicedraws)
       outputtable(j,2) = maskedslice(floor(real(Nslicedraws) * 0.5))
       forall (n=1:Nquantiles) outputtable(j,2+n) = maskedslice(fracndx(n))

       deallocate(maskedslice)
       
    END DO

    ! write output table to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
   WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') ((outputtable(j,n), n = 1,2+Nquantiles), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeMaskedEstimatesTranspose
  ! -----------------------------------------------------------------


  ! -----------------------------------------------------------------
  ! @\newpage\subsection{storeEstimatesOMP}@
  SUBROUTINE storeEstimatesOMP(theta,Ntheta,Ndraws,filename)
    ! store mean, median and quantiles of draws into file

    IMPLICIT NONE

    INTENT(INOUT) :: theta ! has to be INOUT because of sorting
    INTENT(IN)    :: Ntheta,Ndraws,filename

    ! INTEGER, parameter :: Nquantiles = 10
    ! REAL, DIMENSION(Nquantiles), PARAMETER :: quantilegrid = (/ 0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995 /)
    INTEGER, DIMENSION(Nquantiles) :: fracndx 
    INTEGER :: Ntheta, Ndraws, n, k, j, status
    DOUBLE PRECISION, DIMENSION(Ntheta,Ndraws) :: theta
    DOUBLE PRECISION, DIMENSION(Ntheta) :: mittel, median

    CHARACTER (LEN=200) :: filename

    ! sort draws
    !$OMP PARALLEL DO SHARED(Ndraws,theta) PRIVATE(status)
    DO n = 1, Ntheta
       CALL dlasrt('I', Ndraws, theta(n,:), status)
       IF (status /= 0) THEN
          write (*,*), 'DLASORT ERROR ', status, ' [STORE ESTIMATES]'
          stop 1
       END IF
    END DO
    !$OMP END PARALLEL DO

    ! generate index for each fractile
    FORALL (n = 1:Nquantiles)
       fracndx(n)  = floor(real(Ndraws) * quantilegrid(n))
    END FORALL
    where (fracndx == 0) fracndx = 1 

    ! compute results
    mittel    = sum(theta, 2) / dble(Ndraws)
    median    = theta(:,floor(real(Ndraws) * 0.5))

    ! write to file
    OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
    WRITE(4,'(ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16,ES30.16)') (mittel(j), median(j), (theta(j,fracndx(k)), k = 1,Nquantiles), j=1,Ntheta)
    CLOSE(UNIT=4)

  END SUBROUTINE storeEstimatesOMP
  ! -----------------------------------------------------------------
END MODULE embox

