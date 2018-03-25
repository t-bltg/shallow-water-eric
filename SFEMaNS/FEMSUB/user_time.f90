!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
FUNCTION user_time(dummy) RESULT(time)
  IMPLICIT NONE

  !EXTERNAL etime
  !REAL*4 t(2),  etime
  REAL(KIND=8) :: dummy, time 

  !time = etime(t)
  INTEGER :: count, count_rate, count_max

  CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
  time = (1.d0*count)/count_rate

END FUNCTION user_time
