MODULE timerbox

  use omp_lib

  IMPLICIT NONE

  type progresstimer
     double precision :: start, last, wait
     integer :: rank
  end type progresstimer

CONTAINS


  ! @\newpage\subsection{tidomp}@
  FUNCTION tidomp() result(TID)

    integer :: tid
    TID = 0
    !$ TID = OMP_GET_THREAD_NUM()

  END FUNCTION tidomp

  ! @\newpage\subsection{tic}@
  FUNCTION tic() result(walltime)
    ! same as walltime

    double precision :: walltime

    ! CALL CPU_TIME(walltime) 
    walltime = omp_get_wtime()

  END FUNCTION tic

  ! @\newpage\subsection{timespent}@
  FUNCTION toc(t) result(timespent)
    ! like timespent, except for taking double precision as input
    INTENT(IN) :: t

    DOUBLE PRECISION :: now, timespent, t

    now = walltime()
    timespent = now-t

  END FUNCTION toc
  
  ! @\newpage\subsection{walltime}@
  FUNCTION walltime()

    double precision :: walltime

    ! CALL CPU_TIME(walltime) 
    walltime = omp_get_wtime()

  END FUNCTION walltime

  ! @\newpage\subsection{timespent}@
  FUNCTION timespent(timer)
    INTENT(INOUT) :: timer

    DOUBLE PRECISION :: now, timespent
    type(progresstimer) :: timer

    now = walltime()
    timespent = now-timer%last

  END FUNCTION timespent

  ! @\newpage\subsection{progressbar}@
  SUBROUTINE progressbar(fractiondone, timer)
    INTENT(IN) :: fractiondone
    INTENT(INOUT) :: timer

    DOUBLE PRECISION :: now,timeleft, fractiondone
    INTEGER :: daysleft, hoursleft, minutesleft, secondsleft
    type(progresstimer) :: timer

    now = walltime()

    !     call hrulefill
    !     print *, now
    !     print *, now-timer%last 
    !     print *, timer%wait
    !     call hrulefill

    IF (now-timer%last > timer%wait) THEN
       timeleft = (now-timer%start) / fractiondone * (1 - fractiondone)

       daysleft    = floor(timeleft / (3600.0) / 24.0)
       timeleft    = timeleft - 24.0 * 3600.0 * daysleft
       hoursleft   = floor(timeleft / 3600.0)
       timeleft    = timeleft - 3600.0 * hoursleft
       minutesleft = floor(timeleft / 60.0)
       secondsleft = int(timeleft - 60.0 * minutesleft)

       IF (timer%rank < 0) THEN ! negative stream means "no stream"
          IF (daysleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " days ", i3, " hours and ", i3, " minutes remaining")')  fractiondone * 100.0, daysleft, hoursleft, minutesleft
          ELSEIF (hoursleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " hours and ", i3, " minutes remaining")')  fractiondone * 100.0, hoursleft, minutesleft
          ELSEIF (minutesleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " minutes and ", i2, " sec. remaining")')  fractiondone * 100.0, minutesleft, secondsleft
          ELSE 
             WRITE(*,'(f6.2, " % -- ", i3, " sec. remaining")')  fractiondone * 100.0, secondsleft
          END IF
       ELSE
          IF (daysleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " days ", i3, " hours and ", i3, " minutes remaining")')  timer%rank, fractiondone * 100.0, daysleft, hoursleft, minutesleft
          ELSEIF (hoursleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " hours and ", i3, " minutes remaining")')  timer%rank, fractiondone * 100.0, hoursleft, minutesleft
          ELSEIF (minutesleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " minutes and ", i2, " sec. remaining")')  timer%rank, fractiondone * 100.0, minutesleft, secondsleft
          ELSE 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " sec. remaining")')  timer%rank, fractiondone * 100.0, secondsleft
          END IF
       END IF
       timer%last = now
    END IF

  END SUBROUTINE progressbar

  ! @\newpage\subsection{progressbarcomment}@
  SUBROUTINE progressbarcomment(fractiondone, timer, comment)
    INTENT(IN) :: fractiondone, comment
    INTENT(INOUT) :: timer

    DOUBLE PRECISION :: now,timeleft, fractiondone
    INTEGER :: daysleft, hoursleft, minutesleft, secondsleft
    type(progresstimer) :: timer
    character (len=*) :: comment

    now = walltime()

    !     call hrulefill
    !     print *, now
    !     print *, now-timer%last 
    !     print *, timer%wait
    !     call hrulefill

    IF (now-timer%last > timer%wait) THEN
       timeleft = (now-timer%start) / fractiondone * (1 - fractiondone)

       daysleft    = floor(timeleft / (3600.0) / 24.0)
       timeleft    = timeleft - 24.0 * 3600.0 * daysleft
       hoursleft   = floor(timeleft / 3600.0)
       timeleft    = timeleft - 3600.0 * hoursleft
       minutesleft = floor(timeleft / 60.0)
       secondsleft = int(timeleft - 60.0 * minutesleft)

       IF (timer%rank < 0) THEN ! negative stream means "no stream"
          IF (daysleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " days ", i3, " hours and ", i3, " minutes remaining -- ", a)')  fractiondone * 100.0, daysleft, hoursleft, minutesleft, comment
          ELSEIF (hoursleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " hours and ", i3, " minutes remaining -- ", a)')  fractiondone * 100.0, hoursleft, minutesleft, comment
          ELSEIF (minutesleft > 0) THEN 
             WRITE(*,'(f6.2, " % -- ", i3, " minutes and ", i2, " sec. remaining -- ", a)')  fractiondone * 100.0, minutesleft, secondsleft, comment
          ELSE 
             WRITE(*,'(f6.2, " % -- ", i3, " sec. remaining -- ", a)')  fractiondone * 100.0, secondsleft, comment
          END IF
       ELSE
          IF (daysleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " days ", i3, " hours and ", i3, " minutes remaining -- ", a)')  timer%rank, fractiondone * 100.0, daysleft, hoursleft, minutesleft, comment
          ELSEIF (hoursleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " hours and ", i3, " minutes remaining -- ", a)')  timer%rank, fractiondone * 100.0, hoursleft, minutesleft, comment
          ELSEIF (minutesleft > 0) THEN 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " minutes and ", i2, " sec. remaining -- ", a)')  timer%rank, fractiondone * 100.0, minutesleft, secondsleft, comment
          ELSE 
             WRITE(*,'(" Stream ", i2, ": ", f6.2, " % -- ", i3, " sec. remaining -- ", a)')  timer%rank, fractiondone * 100.0, secondsleft, comment
          END IF
       END IF
       timer%last = now
    END IF

  END SUBROUTINE progressbarcomment

  ! @\newpage\subsection{initprogressbar}@
  SUBROUTINE initprogressbar(timer,wait,rank)
    INTENT(INOUT) :: timer
    INTENT(IN) :: rank, wait
    OPTIONAL :: wait, rank

    type(progresstimer) :: timer
    INTEGER :: rank
    DOUBLE PRECISION :: wait

    timer%start = walltime()
    timer%last = timer%start

    IF (PRESENT(rank)) THEN
       timer%rank = rank
    ELSE
       timer%rank = -1
    END IF
    IF (PRESENT(wait)) THEN
       timer%wait = wait
    ELSE
       timer%wait = 10.0
    END IF


  END SUBROUTINE initprogressbar

  ! @\newpage\subsection{timerstart}@
  SUBROUTINE timerstart(timer)
    INTENT(INOUT) :: timer
    type(progresstimer) :: timer
    timer%start = walltime()
  END SUBROUTINE timerstart


  ! @\newpage\subsection{timerstop}@
  SUBROUTINE timerstop(timer, comment)
    INTENT(IN) :: comment
    INTENT(INOUT) :: timer

    DOUBLE PRECISION :: now, timespent
    INTEGER :: daysspent, hoursspent, minutesspent, secondsspent
    type(progresstimer) :: timer
    character (len=*) :: comment

    now = walltime()

    timespent = now-timer%start

    daysspent    = floor(timespent / (3600.0) / 24.0)
    timespent    = timespent - 24.0 * 3600.0 * daysspent
    hoursspent   = floor(timespent / 3600.0)
    timespent    = timespent - 3600.0 * hoursspent
    minutesspent = floor(timespent / 60.0)
    secondsspent = int(timespent - 60.0 * minutesspent)

    IF (daysspent > 0) THEN 
       WRITE(*,'(i3, " days ", i3, " hours and ", i3, " minutes -- ", a)')  daysspent, hoursspent, minutesspent, comment
    ELSEIF (hoursspent > 0) THEN 
       WRITE(*,'(i3, " hours and ", i3, " minutes -- ", a)')  hoursspent, minutesspent, comment
    ELSEIF (minutesspent > 0) THEN 
       WRITE(*,'(i3, " minutes and ", i2, " sec. -- ", a)')  minutesspent, secondsspent, comment
    ELSE 
       WRITE(*,'(i3, " sec. remaining -- ", a)')  secondsspent, comment
    END IF

  END SUBROUTINE timerstop

END MODULE timerbox


