MODULE character_strings
CONTAINS
  SUBROUTINE read_until(unit, string, error)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    LOGICAL, OPTIONAL                  :: error
    IF (PRESENT(error)) error =.FALSE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       IF (trim(adjustl(control))==string) RETURN
    END DO
    RETURN
11  WRITE(*,*) ' Error reading data file '; IF (PRESENT(error)) error=.TRUE.; RETURN
22  WRITE(*,*) ' Data string ',string,' not found '; IF (PRESENT(error)) error=.TRUE.; RETURN
  END SUBROUTINE read_until

  SUBROUTINE find_string(unit, string, okay)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    LOGICAL                            :: okay
    okay = .TRUE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       IF (trim(adjustl(control))==string) RETURN
    END DO
11  WRITE(*,*) ' Erreur de lecture '; STOP
22  okay = .FALSE.; RETURN
  END SUBROUTINE find_string

END MODULE character_strings
