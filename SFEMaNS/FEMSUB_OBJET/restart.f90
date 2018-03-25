!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!
MODULE restart

CONTAINS

  SUBROUTINE write_restart_ns(vv_mesh, pp_mesh, time, list_mode, un, un_m1, pn, pn_m1, incpn, incpn_m1, filename, it, freq_restart)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: vv_mesh,pp_mesh     
    REAL(KIND=8),                                   INTENT(IN) :: time 
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: un, un_m1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: pn, pn_m1, incpn, incpn_m1
    CHARACTER(LEN=64),                              INTENT(IN) :: filename 
    INTEGER,                                        INTENT(IN) :: it, freq_restart

    INTEGER                           :: rang, code, nb_procs, n, i
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO

    DO n = 1, nb_procs
       IF (rang == n-1) THEN
          IF (rang == 0) THEN
             OPEN(UNIT = 10, FILE = 'suite_ns'//'_'//tit//'.'//filename, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace') 
             WRITE(10) time, vv_mesh%np , pp_mesh%np , nb_procs, SIZE(list_mode)
          ELSE
             OPEN(UNIT = 10, FILE = 'suite_ns'//'_'//tit//'.'//filename, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown') 
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) un(:,:,i)
             WRITE(10) un_m1(:,:,i)
             WRITE(10) pn(:,:,i)
             WRITE(10) pn_m1(:,:,i)
             WRITE(10) incpn(:,:,i)
             WRITE(10) incpn_m1(:,:,i)
          END DO
          CLOSE(10)        
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,code)
    END DO

  END SUBROUTINE write_restart_ns


  SUBROUTINE read_restart_ns(vv_mesh, pp_mesh, time, list_mode, un, un_m1, pn, pn_m1, incpn, incpn_m1, &
       filename, val_init, interpol)

    USE def_type_mesh

    IMPLICIT NONE

    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: vv_mesh,pp_mesh     
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: un, un_m1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: pn, pn_m1, incpn, incpn_m1
    CHARACTER(LEN=64),                              INTENT(IN) :: filename 
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol
    INTEGER     :: rang, code, nb_procs, n, i, mode, j, m_max
    INTEGER     :: npr_vv, npr_pp, m_max_cr, nb_procs_r
    INTEGER     :: compt, m_max_c, nb_mode_r, mode_cherche   
    LOGICAL     :: trouve, okay

    WRITE(*,*) 'restart Navier-Stokes'
    OPEN(UNIT = 10, FILE = 'suite_ns.'//filename, FORM = 'unformatted', STATUS = 'unknown')

    READ(10) time, npr_vv , npr_pp, nb_procs_r, m_max_cr
    IF ((npr_vv /= vv_mesh%np) .OR. (npr_pp /= pp_mesh%np)) THEN
       WRITE(*,*) 'READ RESTART INCOMPATIBLE'
       WRITE(*,*) 'np(P1) = ',pp_mesh%np ,' ; np_restart(P1) = ',npr_pp
       WRITE(*,*) 'np(P2) = ',vv_mesh%np ,' ; np_restart(P2) = ',npr_vv
       STOP
    ENDIF
    CLOSE(10)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    IF (rang == 0) THEN
       WRITE(*,*) 'proprietes fichier ', 'suite_ns.'//filename
       WRITE(*,*) 'time =',time
       WRITE(*,*) 'nombre de processeurs = ',nb_procs_r
       WRITE(*,*) 'nombre de modes par processeur = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul 
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs*m_max_c /= nb_mode_r) THEN
       WRITE(*,*) ' BUG '
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN 
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    WRITE(*,*) 'Relecture des modes Navier-Stokes ...'
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = 'suite_ns.'//filename, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premiere ligne du fichier qui contient des donnees
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang*m_max_c+i) THEN
                DO n=1, 6
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) un(:,:,i)
             READ(10) un_m1(:,:,i)
             READ(10) pn(:,:,i)
             READ(10) pn_m1(:,:,i)
             READ(10) incpn(:,:,i)
             READ(10) incpn_m1(:,:,i)
             WRITE(*,*) 'mode ns',mode_cherche,' trouve '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 6 lignes
             DO n=1, 6
                READ(10)
             ENDDO
          ENDIF
       ENDDO
       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN
             un(:,:,i)    = val_init  ; un_m1(:,:,i)    = val_init 
             pn(:,:,i)    = val_init ; pn_m1(:,:,i)    = val_init
             incpn(:,:,i) = val_init ; incpn_m1(:,:,i) =  val_init
             WRITE(*,*) 'mode ns',mode_cherche,' non trouve'
          ELSE
             un(:,:,i)    = 0.d0 ; un_m1(:,:,i)    = 0.d0
             pn(:,:,i)    = 0.d0 ; pn_m1(:,:,i)    = 0.d0
             incpn(:,:,i) = 0.d0 ; incpn_m1(:,:,i) = 0.d0
             WRITE(*,*) 'mode ns',mode_cherche,' non trouve'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite 
    ENDDO

  END SUBROUTINE read_restart_ns

  SUBROUTINE write_restart_maxwell(H_mesh, phi_mesh, time, list_mode, Hn, Hn1, phin, phin1, filename, it, freq_restart)

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh     
    REAL(KIND=8),                                   INTENT(IN) :: time 
    INTEGER,      DIMENSION(:),                     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(IN) :: phin, phin1
    CHARACTER(LEN=64),                              INTENT(IN) :: filename 
    INTEGER,                                        INTENT(IN) :: it, freq_restart

    INTEGER                           :: rang, code, nb_procs, n, i
    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

    WRITE(tit,'(i3)') it/freq_restart
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO

    DO n = 1, nb_procs
       IF (rang == n-1) THEN
          IF (rang == 0) THEN
             OPEN(UNIT = 10, FILE = 'suite_maxwell'//'_'//tit//'.'//filename, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'replace') 
             WRITE(10) time, H_mesh%np , phi_mesh%np , nb_procs, SIZE(list_mode)
          ELSE
             OPEN(UNIT = 10, FILE = 'suite_maxwell'//'_'//tit//'.'//filename, POSITION='append', &
                  FORM = 'unformatted', STATUS = 'unknown') 
          END IF

          DO i= 1, SIZE(list_mode)
             WRITE(10) list_mode(i)
             WRITE(10) Hn(:,:,i)
             WRITE(10) Hn1(:,:,i)
             WRITE(10) phin(:,:,i)
             WRITE(10) phin1(:,:,i)
          END DO
          CLOSE(10)        
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,code)
    END DO

  END SUBROUTINE write_restart_maxwell


  SUBROUTINE read_restart_maxwell(H_mesh, phi_mesh, time, list_mode, Hn, Hn1, phin, phin1, &
       filename, val_init, interpol)

    USE def_type_mesh

    IMPLICIT NONE

    INCLUDE 'mpif.h'
    TYPE(mesh_type), TARGET                                    :: H_mesh,phi_mesh     
    REAL(KIND=8),                                   INTENT(OUT):: time
    INTEGER,      DIMENSION(:)                                 :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:),                 INTENT(OUT):: phin, phin1
    CHARACTER(LEN=64),                              INTENT(IN) :: filename 
    REAL(KIND=8), OPTIONAL,                         INTENT(IN) :: val_init
    LOGICAL     , OPTIONAL,                         INTENT(IN) :: interpol

    INTEGER     :: rang, code, nb_procs, n, i, mode, j
    INTEGER     :: npr_H, npr_phi, m_max_cr, nb_procs_r
    INTEGER     :: m_max_c, nb_mode_r, mode_cherche
    LOGICAL     :: trouve, okay

    WRITE(*,*) 'restart Maxwell'
    OPEN(UNIT = 10, FILE = 'suite_maxwell.'//filename, FORM = 'unformatted', STATUS = 'unknown')
    READ(10) time, npr_H , npr_phi, nb_procs_r, m_max_cr
    IF ((npr_H /= H_mesh%np) .OR. (npr_phi /= phi_mesh%np)) THEN
       WRITE(*,*) 'READ RESTART INCOMPATIBLE'
       WRITE(*,*) 'np phi = ',phi_mesh%np ,' ; np_restart np phi = ',npr_phi
       WRITE(*,*) 'np H   = ',H_mesh%np ,'   ; np_restart np H   = ',npr_H
       STOP
    ENDIF

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CLOSE(10)
    IF (rang == 0) THEN
       WRITE(*,*) 'proprietes fichier ', 'suite_maxwell.'//filename
       WRITE(*,*) 'time =',time
       WRITE(*,*) 'nombre de processeurs = ',nb_procs_r
       WRITE(*,*) 'nombre de modes par processeur = ',m_max_cr
    ENDIF

    m_max_c   = SIZE(list_mode)      !nombre de modes par proc pour le calcul
    nb_mode_r = nb_procs_r*m_max_cr  !nombre total de modes contenus dans le suite

    !June 7 2007, JLG
    IF (nb_procs*m_max_c /= nb_mode_r) THEN
       WRITE(*,*) ' BUG '
       !STOP
    END IF

    okay = .FALSE.
    IF (PRESENT(interpol)) THEN 
       IF (interpol) THEN
          okay =.TRUE.
       END IF
    END IF
    !June 7 2007, JLG

    WRITE(*,*) 'Relecture des modes Maxwell...'
    DO i=1, m_max_c                  !pour tout les modes du processeur courant
       !ouverture du fichier
       OPEN(UNIT = 10, FILE = 'suite_maxwell.'//filename, FORM = 'unformatted', STATUS = 'unknown')
       !on saute la premier ligne du fichier qui contient des donnes
       READ(10)
       mode_cherche = list_mode(i)
       !recherche du bon mode
       trouve = .FALSE.
       DO j=1, nb_mode_r             !pour tout les modes ecris dans le suite.
          !lecture du mode
          READ(10) mode
          !June 7 2007, JLG
          IF (okay) THEN
             IF (j/=rang*m_max_c+i) THEN
                DO n=1, 4
                   READ(10)
                ENDDO
                CYCLE
             ELSE
                list_mode(i) = mode
                mode_cherche = mode
             END IF
          END IF
          !June 7 2007, JLG
          IF (mode == mode_cherche) THEN   !on a trouve le bon mode
             READ(10) Hn(:,:,i)
             READ(10) Hn1(:,:,i)
             READ(10) phin(:,:,i)
             READ(10) phin1(:,:,i)
             WRITE(*,*) 'mode maxwell',mode_cherche,' trouve '
             trouve = .TRUE.
             EXIT                        !car on a trouve le bon mode
          ELSE                             !on passe au mode suivant en sautant 4 lignes
             DO n=1, 4
                READ(10)
             ENDDO
          ENDIF
       ENDDO
       IF (.NOT.trouve) THEN               !mode_cherche non trouve
          IF (PRESENT(val_init)) THEN 
             Hn(:,:,i)   = val_init ; Hn1(:,:,i)   = val_init 
             phin(:,:,i) = val_init ; phin1(:,:,i) = val_init
             WRITE(*,*) 'mode maxwell',mode_cherche,' non trouve'
          ELSE
             Hn(:,:,i)   = 0.d0 ; Hn1(:,:,i)   = 0.d0
             phin(:,:,i) = 0.d0 ; phin1(:,:,i) = 0.d0
             WRITE(*,*) 'mode maxwell',mode_cherche,' non trouve'
          ENDIF
       ENDIF
       CLOSE(10)                          !fermeture du fichier suite
    ENDDO

  END SUBROUTINE read_restart_maxwell

END MODULE restart

