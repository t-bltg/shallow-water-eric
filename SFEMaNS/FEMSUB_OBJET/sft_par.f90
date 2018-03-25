!
!Authors Jean-Luc Guermond, Rapahel Laguerre, Copyrights 2005
!
MODULE sft_parallele

  IMPLICIT NONE
  PUBLIC ::  FFT_PAR_MP, FFT_PAR_PROD_VECT, FFT_PAR_MP_b, FFT_PAR_PROD_VECT_b,&
       smb_ns_gauss_sft_par, FFT_PAR_CROSS_PROD, smb_cross_prod_gauss_sft_par
  PRIVATE
CONTAINS
  SUBROUTINE FFT_PAR_CROSS_PROD(V1_in, V2_in, V_out, temps)
    ! Format: V_1in(1:np,1:6,1:m_max_c)    
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r,   V1t,   V1z,   V2r,   V2t,   V2z
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p
    REAL(KIND=8) :: t
    INTEGER :: m_max_c, np, i, n, nb_bloc, np_loc, np_alloc, reste, &
         code, nb_procs, n_inf, n_sup, paq, n_up
    LOGICAL :: init, clean

    !Temps 1 = Temps de Comm
    !Temps 2 = Temps de calcul
    !Temps 3 = Temps de Changement de format

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

    m_max_c = SIZE(V1_in,3)
    np      = SIZE(V1_in,1)

    nb_bloc =  MAX(MIN(200*nb_procs,np/2),1)

100 np_loc = np/nb_bloc
    reste = np - np_loc*nb_bloc
    np_alloc = np_loc + reste
    reste = np_alloc*nb_bloc - np
    IF (reste>np_alloc) THEN
       nb_bloc = nb_bloc - 1
       GO TO 100
    END IF

    ALLOCATE(V1r  (2*m_max_c,np_alloc),V1t  (2*m_max_c,np_alloc),V1z  (2*m_max_c,np_alloc))
    ALLOCATE(V2r  (2*m_max_c,np_alloc),V2t  (2*m_max_c,np_alloc),V2z  (2*m_max_c,np_alloc))
    ALLOCATE(V1r_p(2*m_max_c,np_alloc),V1t_p(2*m_max_c,np_alloc),V1z_p(2*m_max_c,np_alloc))
    ALLOCATE(V2r_p(2*m_max_c,np_alloc),V2t_p(2*m_max_c,np_alloc),V2z_p(2*m_max_c,np_alloc))

    n_sup = 0
    DO paq = 1, nb_bloc
       !init = .true.
       !clean= .true.
       init = .FALSE.
       clean= .FALSE.
       IF (paq==1) init=.TRUE.
       IF (paq==nb_bloc) clean=.TRUE.

       n_inf = n_sup + 1
       IF (paq == nb_bloc) THEN
          n_up  = np - n_inf + 1
       ELSE
          n_up  = np_alloc
       END IF
       n_sup = n_inf + n_up - 1

       t = MPI_WTIME()
       DO i =1, m_max_c
          V1r(2*i-1,1:n_up) = V1_in(n_inf:n_sup,1,i)
          V1r(2*i  ,1:n_up) = V1_in(n_inf:n_sup,2,i)
          V1t(2*i-1,1:n_up) = V1_in(n_inf:n_sup,3,i)
          V1t(2*i  ,1:n_up) = V1_in(n_inf:n_sup,4,i)
          V1z(2*i-1,1:n_up) = V1_in(n_inf:n_sup,5,i)
          V1z(2*i  ,1:n_up) = V1_in(n_inf:n_sup,6,i)
          V2r(2*i-1,1:n_up) = V2_in(n_inf:n_sup,1,i)
          V2r(2*i  ,1:n_up) = V2_in(n_inf:n_sup,2,i)
          V2t(2*i-1,1:n_up) = V2_in(n_inf:n_sup,3,i)
          V2t(2*i  ,1:n_up) = V2_in(n_inf:n_sup,4,i)
          V2z(2*i-1,1:n_up) = V2_in(n_inf:n_sup,5,i)
          V2z(2*i  ,1:n_up) = V2_in(n_inf:n_sup,6,i)
       END DO

       IF (paq == nb_bloc .AND. n_up==np_alloc) THEN
          V1r(2*i-1,n_up+1:) = 0
          V1r(2*i  ,n_up+1:) = 0
          V1t(2*i-1,n_up+1:) = 0
          V1t(2*i  ,n_up+1:) = 0
          V1z(2*i-1,n_up+1:) = 0
          V1z(2*i  ,n_up+1:) = 0
          V2r(2*i-1,n_up+1:) = 0
          V2r(2*i  ,n_up+1:) = 0
          V2t(2*i-1,n_up+1:) = 0
          V2t(2*i  ,n_up+1:) = 0
          V2z(2*i-1,n_up+1:) = 0
          V2z(2*i  ,n_up+1:) = 0
       END IF
       IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       IF(PRESENT(temps)) THEN

          CALL FFT_PAR_MP(V2r, V2r_p, 1, temps, init=.TRUE.)
          CALL FFT_PAR_MP(V2t, V2t_p, 1, temps)
          CALL FFT_PAR_MP(V2z, V2z_p, 1, temps)

          CALL FFT_PAR_MP(V1r, V1r_p, 1, temps)
          CALL FFT_PAR_MP(V1t, V1t_p, 1, temps)
          CALL FFT_PAR_MP(V1z, V1z_p, 1, temps)
       ELSE
          CALL FFT_PAR_MP(V2r, V2r_p, 1, init=.TRUE.)
          CALL FFT_PAR_MP(V2t, V2t_p, 1)
          CALL FFT_PAR_MP(V2z, V2z_p, 1)

          CALL FFT_PAR_MP(V1r, V1r_p, 1)
          CALL FFT_PAR_MP(V1t, V1t_p, 1)
          CALL FFT_PAR_MP(V1z, V1z_p, 1)
       END IF

       t = MPI_WTIME()
       V1r = V1t_p*V2z_p - V1z_p*V2t_p
       V1t = V1z_p*V2r_p - V1r_p*V2z_p
       V1z = V1r_p*V2t_p - V1t_p*V2r_p
       IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t



       IF(PRESENT(temps)) THEN
          CALL FFT_PAR_MP(V1r, V2r, 2, temps)
          CALL FFT_PAR_MP(V1t, V2t, 2, temps)
          CALL FFT_PAR_MP(V1z, V2z, 2, temps, clean=.TRUE.)
       ELSE
          CALL FFT_PAR_MP(V1r, V2r, 2)
          CALL FFT_PAR_MP(V1t, V2t, 2)
          CALL FFT_PAR_MP(V1z, V2z, 2, clean=.TRUE.)
       END IF

       t = MPI_WTIME()
       DO i =1 , m_max_c
          V_out(n_inf:n_sup,1,i) = V2r(2*i-1,1:n_up)
          V_out(n_inf:n_sup,2,i) = V2r(2*i  ,1:n_up)
          V_out(n_inf:n_sup,3,i) = V2t(2*i-1,1:n_up)
          V_out(n_inf:n_sup,4,i) = V2t(2*i  ,1:n_up)
          V_out(n_inf:n_sup,5,i) = V2z(2*i-1,1:n_up)
          V_out(n_inf:n_sup,6,i) = V2z(2*i  ,1:n_up)
       END DO
       IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t
    END DO

    DEALLOCATE(V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p)
    DEALLOCATE(V1r, V1t, V1z)
    DEALLOCATE(V2r, V2t, V2z)

  END SUBROUTINE FFT_PAR_CROSS_PROD

  SUBROUTINE FFT_PAR_PROD_VECT(V1_in, V2_in, V_out, temps)    
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r,   V1t,   V1z,   V2r,   V2t,   V2z
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p
    REAL(KIND=8) :: t
    INTEGER :: m_max_c, np, i, n, nb_bloc, np_loc, np_alloc, reste, &
         code, nb_procs, n_inf, n_sup, paq, n_up
    LOGICAL :: init, clean

    !Temps 1 = Temps de Comm
    !Temps 2 = Temps de calcul
    !Temps 3 = Temps de Changement de format

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

    m_max_c = SIZE(V1_in,3)
    np      = SIZE(V1_in,2)

    nb_bloc =  MAX(MIN(200*nb_procs,np/10),1)

100 np_loc = np/nb_bloc
    reste = np - np_loc*nb_bloc
    np_alloc = np_loc + reste
    reste = np_alloc*nb_bloc - np
    IF (reste>np_alloc) THEN
       nb_bloc = nb_bloc - 1
       GO TO 100
    END IF

    ALLOCATE(V1r  (2*m_max_c,np_alloc),V1t  (2*m_max_c,np_alloc),V1z  (2*m_max_c,np_alloc))
    ALLOCATE(V2r  (2*m_max_c,np_alloc),V2t  (2*m_max_c,np_alloc),V2z  (2*m_max_c,np_alloc))
    ALLOCATE(V1r_p(2*m_max_c,np_alloc),V1t_p(2*m_max_c,np_alloc),V1z_p(2*m_max_c,np_alloc))
    ALLOCATE(V2r_p(2*m_max_c,np_alloc),V2t_p(2*m_max_c,np_alloc),V2z_p(2*m_max_c,np_alloc))

    n_sup = 0
    DO paq = 1, nb_bloc
       !init = .true.
       !clean= .true.
       init = .FALSE.
       clean= .FALSE.
       IF (paq==1) init=.TRUE.
       IF (paq==nb_bloc) clean=.TRUE.

       n_inf = n_sup + 1
       IF (paq == nb_bloc) THEN
          n_up  = np - n_inf + 1
       ELSE
          n_up  = np_alloc
       END IF
       n_sup = n_inf + n_up - 1

       t = MPI_WTIME()
       DO i =1, m_max_c
          V1r(2*i-1,1:n_up) = V1_in(1,n_inf:n_sup,i)
          V1r(2*i  ,1:n_up) = V1_in(2,n_inf:n_sup,i)
          V1t(2*i-1,1:n_up) = V1_in(3,n_inf:n_sup,i)
          V1t(2*i  ,1:n_up) = V1_in(4,n_inf:n_sup,i)
          V1z(2*i-1,1:n_up) = V1_in(5,n_inf:n_sup,i)
          V1z(2*i  ,1:n_up) = V1_in(6,n_inf:n_sup,i)
          V2r(2*i-1,1:n_up) = V2_in(1,n_inf:n_sup,i)
          V2r(2*i  ,1:n_up) = V2_in(2,n_inf:n_sup,i)
          V2t(2*i-1,1:n_up) = V2_in(3,n_inf:n_sup,i)
          V2t(2*i  ,1:n_up) = V2_in(4,n_inf:n_sup,i)
          V2z(2*i-1,1:n_up) = V2_in(5,n_inf:n_sup,i)
          V2z(2*i  ,1:n_up) = V2_in(6,n_inf:n_sup,i)
       END DO

       IF (paq == nb_bloc .AND. n_up==np_alloc) THEN
          V1r(2*i-1,n_up+1:) = 0
          V1r(2*i  ,n_up+1:) = 0
          V1t(2*i-1,n_up+1:) = 0
          V1t(2*i  ,n_up+1:) = 0
          V1z(2*i-1,n_up+1:) = 0
          V1z(2*i  ,n_up+1:) = 0
          V2r(2*i-1,n_up+1:) = 0
          V2r(2*i  ,n_up+1:) = 0
          V2t(2*i-1,n_up+1:) = 0
          V2t(2*i  ,n_up+1:) = 0
          V2z(2*i-1,n_up+1:) = 0
          V2z(2*i  ,n_up+1:) = 0
       END IF
       IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       IF(PRESENT(temps)) THEN
          CALL FFT_PAR_MP(V2r, V2r_p, 1, temps, init=.TRUE.)
          CALL FFT_PAR_MP(V2t, V2t_p, 1, temps)
          CALL FFT_PAR_MP(V2z, V2z_p, 1, temps)

          CALL FFT_PAR_MP(V1r, V1r_p, 1, temps)
          CALL FFT_PAR_MP(V1t, V1t_p, 1, temps)
          CALL FFT_PAR_MP(V1z, V1z_p, 1, temps)
       ELSE
          CALL FFT_PAR_MP(V2r, V2r_p, 1, init=.TRUE.)
          CALL FFT_PAR_MP(V2t, V2t_p, 1)
          CALL FFT_PAR_MP(V2z, V2z_p, 1)

          CALL FFT_PAR_MP(V1r, V1r_p, 1)
          CALL FFT_PAR_MP(V1t, V1t_p, 1)
          CALL FFT_PAR_MP(V1z, V1z_p, 1)
       END IF

       t = MPI_WTIME()
       V1r = V1t_p*V2z_p - V1z_p*V2t_p
       V1t = V1z_p*V2r_p - V1r_p*V2z_p
       V1z = V1r_p*V2t_p - V1t_p*V2r_p
       IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t



       IF(PRESENT(temps)) THEN
          CALL FFT_PAR_MP(V1r, V2r, 2, temps)
          CALL FFT_PAR_MP(V1t, V2t, 2, temps)
          CALL FFT_PAR_MP(V1z, V2z, 2, temps, clean=.TRUE.)
       ELSE
          CALL FFT_PAR_MP(V1r, V2r, 2)
          CALL FFT_PAR_MP(V1t, V2t, 2)
          CALL FFT_PAR_MP(V1z, V2z, 2, clean=.TRUE.)
       END IF

       t = MPI_WTIME()
       DO i =1 , m_max_c
          V_out(1,n_inf:n_sup,i) = V2r(2*i-1,1:n_up)
          V_out(2,n_inf:n_sup,i) = V2r(2*i  ,1:n_up)
          V_out(3,n_inf:n_sup,i) = V2t(2*i-1,1:n_up)
          V_out(4,n_inf:n_sup,i) = V2t(2*i  ,1:n_up)
          V_out(5,n_inf:n_sup,i) = V2z(2*i-1,1:n_up)
          V_out(6,n_inf:n_sup,i) = V2z(2*i  ,1:n_up)
       END DO
       IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t
    END DO

    DEALLOCATE(V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p)
    DEALLOCATE(V1r, V1t, V1z)
    DEALLOCATE(V2r, V2t, V2z)

  END SUBROUTINE FFT_PAR_PROD_VECT

  SUBROUTINE FFT_PAR_MP(V_in, V_out, choix, temps, init, clean)

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:),  INTENT(OUT) :: V_out  
    INTEGER,                       INTENT(IN)  :: choix
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT)      :: temps
    LOGICAL,                    OPTIONAL, INTENT(IN) :: init, clean

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    LOGICAL, SAVE :: once=.TRUE.
    INTEGER, SAVE :: m_max, m_max_c, rang, nb_procs, etiquette, np
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE  :: col_Asp, col_Aps
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: angle  ! valeurs des angles
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: sous_prod   !sous produits locaux   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: res  
    REAL(KIND=8)                           , SAVE  :: PI, t

    !pour l'envoi
    INTEGER :: proc_recep, proc_cible, b_inf, b_sup, indice, occur, proc, occur1, occur2
    INTEGER :: code, valeur, i, j, k, l, n, l1, l2, i_glob
    INTEGER :: requete0, requete1

    !--------END OF DECLARATION ------------------------------

    !--------INITIALISATION----------------------------------


    IF (PRESENT(init)) once = init

    IF (once) THEN
       once = .FALSE.

       np = SIZE(V_in,2)
       IF (MOD(SIZE(V_in,1),2)/=0) THEN
          WRITE(*,*) ' BUG FFT_PROD'
          STOP
       END IF

       m_max_c = SIZE(V_in,1)/2   !nombre de modes par processeur
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
       m_max = nb_procs*m_max_c !CONVENTION: m_max = mode_max +1
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
       etiquette = 100

       !Pour les colonnes col_cos et col_sin, je considere que le mode 0 se rattache a une partie
       !en cos et en sin pour systematiser les envois, 
       !on a donc une dimension de 2*mode_max +2 pour les colonnes !
       !Pour les angles, il n'y a qu'un angle pour le mode 0, et 2 pour les
       !autres modes, d'ou une dimension 2*mod_max + 1 pour le tableau 'angle'

       ALLOCATE(angle(2*m_max))
       ALLOCATE(col_Asp(2*m_max,2*m_max_c),col_Aps(2*m_max,2*m_max_c))
       ALLOCATE(sous_prod(2*m_max*2*np))
       ALLOCATE(res(2*m_max*2*np))
       !Calcul des angles qui sont les memes sur tous les processeurs
       PI = ACOS(-1.d0)
       !attention, theta = 0 correspond a i=0 !
       angle(1:2) = 0.d0  !angle(2) = angle bidon
       DO i= 3, 2*m_max
          angle(i) = 2.d0*PI/(2*m_max-1)*(i-2)
       ENDDO

       !Calcul des colonnes pour la fft spec -> phys

       DO i= 1, 2*m_max
          DO j= 1, 2*m_max_c
             l = rang*m_max_c + (j-1)/2 
             IF (MOD(j,2)==0) THEN
                col_Asp(i,j) = SIN(l*angle(i))
             ELSE
                col_Asp(i,j) = COS(l*angle(i))
             END IF
          ENDDO
       ENDDO
       DO j= 1, 2*m_max_c
          col_Asp(2,j) = 0.d0
       ENDDO

       !Calcul colonnes pour la fft phys -> spec

       DO i= 1, 2*m_max
          l = (i-1)/2 
          DO j= 1, 2*m_max_c
             IF (i == 1) THEN
                col_Aps(i,j) =  0.5d0
             ELSEIF (i == 2) THEN
                col_Aps(i,j) = 0.d0
             ELSE
                k = rang*2*m_max_c + j
                IF (MOD(i,2)==0) THEN
                   col_Aps(i,j) = SIN(l*angle(k))
                ELSE
                   col_Aps(i,j) = COS(l*angle(k))
                END IF
             ENDIF
          ENDDO
       ENDDO
       IF (rang==0) THEN
          DO i= 1, 2*m_max
             col_Aps(i,2) = 0.d0
          ENDDO
       END IF
       col_Aps = (2.d0/(2*m_max-1))*col_Aps
    END IF
    !------END OF INITIALISATION--------------------

    t = MPI_WTIME()
    IF (choix==1) THEN
       sous_prod = 0
       DO proc = 1, nb_procs
          DO i = 1, 2*m_max_c
             i_glob = (proc-1)*2*m_max_c + i
             DO n = 1, np
                l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i
                DO j = 1, 2*m_max_c
                   !DO j = 1, 2*m_max_c
                   !    DO n = 1, np
                   !      l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                   sous_prod(l1) = sous_prod(l1) + col_Asp(i_glob,j)*V_in(j,n)
                END DO
             END DO
          END DO
       END DO
       !l1  = 0
       !DO proc = 1, nb_procs
       !   l2  = (proc-1)*2*m_max_c
       !   DO n = 1, np
       !      DO i = 1, 2*m_max_c
       !         i_glob = l2 + i
       !         l1 = l1 + 1
       !         DO j = 1, 2*m_max_c
       !            sous_prod(l1) = sous_prod(l1) + col_Asp(i_glob,j)*V_in(j,n)
       !         END DO
       !      END DO
       !   END DO
       !END DO
    ELSE
       sous_prod = 0
       DO proc = 1, nb_procs 
          DO i = 1, 2*m_max_c
             i_glob = (proc-1)*2*m_max_c + i
             DO n = 1, np
                l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                DO j = 1, 2*m_max_c
                   !DO j = 1, 2*m_max_c
                   !   DO n = 1, np
                   !      l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                   sous_prod(l1) = sous_prod(l1) + col_Aps(i_glob,j)*V_in(j,n)
                END DO
             END DO
          END DO
       END DO
       !l1  = 0
       !DO proc = 1, nb_procs
       !   l2  = (proc-1)*2*m_max_c
       !   DO n = 1, np
       !      DO i = 1, 2*m_max_c
       !         i_glob = l2 + i
       !         l1 = l1 + 1
       !         DO j = 1, 2*m_max_c
       !            sous_prod(l1) = sous_prod(l1) + col_Aps(i_glob,j)*V_in(j,n)
       !         END DO
       !      END DO
       !   END DO
       !END DO
    END IF
    IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

    !Rangement
    !proc1i        -> V1_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c) |
    !procx         -> ...... |
    !proc(nb_procs)-> V1_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c)  
    !En esperant que ca marche?

    b_inf = rang*2*m_max_c*np + 1 
    b_sup = b_inf + 2*m_max_c*np - 1 
    res(b_inf:b_sup) = sous_prod(b_inf:b_sup) ! pour V_in

    !solution avec MPI_SENDRECV
    t = MPI_WTIME()
    DO proc =1, nb_procs 
       IF ((proc - 1) /= rang) THEN
          occur = 2*m_max_c*np 
          b_inf = (proc-1)*occur + 1
          b_sup = b_inf + occur - 1
          CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc-1, etiquette, &
               res(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc-1, etiquette, &
               MPI_COMM_WORLD, statut, code)
       ENDIF
    ENDDO
    IF(PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() - t

    t = MPI_WTIME()
    V_out = 0
    occur = 2*m_max_c*np

    DO n = 1, np 
       DO i = 1, 2*m_max_c
          DO proc = 1, nb_procs
             indice = i + (proc-1)*occur + (n-1)*2*m_max_c
             V_out(i,n) = V_out(i,n) + res(indice)
          END DO
       END DO
    ENDDO

    !indice = 0
    !DO proc = 1, nb_procs
    !   DO n = 1, np 
    !      DO i = 1, 2*m_max_c
    !         indice = indice + 1
    !         V_out(i,n) = V_out(i,n) + res(indice)
    !     END DO
    !  END DO
    !ENDDO

    IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

    IF (PRESENT(clean)) THEN
       IF(clean) THEN
          once =.TRUE.
          DEALLOCATE(angle)
          DEALLOCATE(col_Asp,col_Aps)
          DEALLOCATE(sous_prod)
          DEALLOCATE(res)
       END IF
    END IF

  END SUBROUTINE FFT_PAR_MP

  SUBROUTINE FFT_PAR_PROD_VECT_b(V1_in, V2_in, V_out, temps)    
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r,   V1t,   V1z,   V2r,   V2t,   V2z
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE    :: V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p
    REAL(KIND=8) :: t
    INTEGER :: m_max_c, np, i, n

    !Temps 1 = Temps de Comm
    !Temps 2 = Temps de calcul
    !Temps 3 = Temps de Changement de format


    m_max_c = SIZE(V1_in,3)/2
    np      = SIZE(V1_in,2)

    ALLOCATE(V1r(2*m_max_c,np),V1t(2*m_max_c,np),V1z(2*m_max_c,np))
    ALLOCATE(V2r(2*m_max_c,np),V2t(2*m_max_c,np),V2z(2*m_max_c,np))

    t = MPI_WTIME()
    DO i =1, m_max_c
       V1r(2*i-1,:) = V1_in(1,:,i)
       V1r(2*i,:)   = V1_in(2,:,i)
       V1t(2*i-1,:) = V1_in(3,:,i)
       V1t(2*i,:)   = V1_in(4,:,i)
       V1z(2*i-1,:) = V1_in(5,:,i)
       V1z(2*i,:)   = V1_in(6,:,i)
       V2r(2*i-1,:) = V2_in(1,:,i)
       V2r(2*i,:)   = V2_in(2,:,i)
       V2t(2*i-1,:) = V2_in(3,:,i)
       V2t(2*i,:)   = V2_in(4,:,i)
       V2z(2*i-1,:) = V2_in(5,:,i)
       V2z(2*i,:)   = V2_in(6,:,i)
    END DO
    IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ALLOCATE(V1r_p(2*m_max_c,np),V1t_p(2*m_max_c,np),V1z_p(2*m_max_c,np))
    ALLOCATE(V2r_p(2*m_max_c,np),V2t_p(2*m_max_c,np),V2z_p(2*m_max_c,np))

    IF(PRESENT(temps)) THEN
       CALL FFT_PAR_MP_b(V2r, V2r_p, 1, temps)
       CALL FFT_PAR_MP_b(V2t, V2t_p, 1, temps)
       CALL FFT_PAR_MP_b(V2z, V2z_p, 1, temps)

       CALL FFT_PAR_MP_b(V1r, V1r_p, 1, temps)
       CALL FFT_PAR_MP_b(V1t, V1t_p, 1, temps)
       CALL FFT_PAR_MP_b(V1z, V1z_p, 1, temps)
    ELSE
       CALL FFT_PAR_MP_b(V2r, V2r_p, 1)
       CALL FFT_PAR_MP_b(V2t, V2t_p, 1)
       CALL FFT_PAR_MP_b(V2z, V2z_p, 1)

       CALL FFT_PAR_MP_b(V1r, V1r_p, 1)
       CALL FFT_PAR_MP_b(V1t, V1t_p, 1)
       CALL FFT_PAR_MP_b(V1z, V1z_p, 1)
    END IF

    t = MPI_WTIME()
    V1r = V1t_p*V2z_p - V1z_p*V2t_p
    V1t = V1z_p*V2r_p - V1r_p*V2z_p
    V1z = V1r_p*V2t_p - V1t_p*V2r_p
    IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

    DEALLOCATE(V1r_p, V1t_p, V1z_p, V2r_p, V2t_p, V2z_p)

    IF(PRESENT(temps)) THEN
       CALL FFT_PAR_MP_b(V1r, V2r, 2, temps)
       CALL FFT_PAR_MP_b(V1t, V2t, 2, temps)
       CALL FFT_PAR_MP_b(V1z, V2z, 2, temps)
    ELSE
       CALL FFT_PAR_MP_b(V1r, V2r, 2)
       CALL FFT_PAR_MP_b(V1t, V2t, 2)
       CALL FFT_PAR_MP_b(V1z, V2z, 2)
    END IF

    t = MPI_WTIME()
    DO i =1 , m_max_c
       V_out(1,:,i) = V2r(2*i-1,:)
       V_out(2,:,i) = V2r(2*i,:)
       V_out(3,:,i) = V2t(2*i-1,:)
       V_out(4,:,i) = V2t(2*i,:)
       V_out(5,:,i) = V2z(2*i-1,:)
       V_out(6,:,i) = V2z(2*i,:)
    END DO
    IF(PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t

    DEALLOCATE(V1r, V1t, V1z)
    DEALLOCATE(V2r, V2t, V2z)

  END SUBROUTINE FFT_PAR_PROD_VECT_b

  SUBROUTINE FFT_PAR_MP_b(V_in, V_out, choix, temps, init, clean)

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:),  INTENT(OUT) :: V_out  
    INTEGER,                       INTENT(IN)  :: choix
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT)      :: temps
    LOGICAL,                    OPTIONAL, INTENT(IN) :: init, clean

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    LOGICAL, SAVE :: once=.TRUE.
    INTEGER, SAVE :: m_max, m_max_c, rang, nb_procs, etiquette, np
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE  :: col_Asp, col_Aps
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: angle   ! valeurs des angles
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: sous_prod   !sous produits locaux   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: res  
    REAL(KIND=8)                           , SAVE  :: PI, t

    !pour l'envoi
    INTEGER :: proc_recep, proc_cible, b_inf, b_sup, indice, occur, proc
    INTEGER :: code, valeur, i, j, k, l, n, l1, l2, i_glob
    INTEGER :: requete0, requete1

    !--------END OF DECLARATION ------------------------------

    !--------INITIALISATION----------------------------------



    IF (PRESENT(init)) once = init

    IF (once) THEN
       once = .FALSE.

       np = SIZE(V_in,2)
       IF (MOD(SIZE(V_in,1),2)/=0) THEN
          WRITE(*,*) ' BUG FFT_PROD'
          STOP
       END IF

       m_max_c = SIZE(V_in,1)/2   !nombre de modes par processeur
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
       m_max = nb_procs*m_max_c !CONVENTION: m_max = mode_max +1
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
       etiquette = 100

       !Pour les colonnes col_cos et col_sin, je considere que le mode 0 se rattache a une partie
       !en cos et en sin pour systematiser les envois, 
       !on a donc une dimension de 2*mode_max +2 pour les colonnes !
       !Pour les angles, il n'y a qu'un angle pour le mode 0, et 2 pour les
       !autres modes, d'ou une dimension 2*mod_max + 1 pour le tableau 'angle'

       ALLOCATE(angle(2*m_max))
       ALLOCATE(col_Asp(2*m_max,2*m_max_c),col_Aps(2*m_max,2*m_max_c))
       ALLOCATE(sous_prod(2*m_max*2*np))
       ALLOCATE(res(2*m_max*2*np))
       !Calcul des angles qui sont les memes sur tous les processeurs
       PI = ACOS(-1.d0)
       !attention, theta = 0 correspond a i=0 !
       angle(1:2) = 0.d0  !angle(2) = angle bidon
       DO i= 3, 2*m_max
          angle(i) = 2.d0*PI/(2*m_max-1)*(i-2)
       ENDDO

       !Calcul des colonnes pour la fft spec -> phys

       DO i= 1, 2*m_max
          DO j= 1, 2*m_max_c
             l = rang*m_max_c + (j-1)/2 
             IF (MOD(j,2)==0) THEN
                col_Asp(i,j) = SIN(l*angle(i))
             ELSE
                col_Asp(i,j) = COS(l*angle(i))
             END IF
          ENDDO
       ENDDO
       DO j= 1, 2*m_max_c
          col_Asp(2,j) = 0.d0
       ENDDO

       !Calcul colonnes pour la fft phys -> spec

       DO i= 1, 2*m_max
          l = (i-1)/2 
          DO j= 1, 2*m_max_c
             IF (i == 1) THEN
                col_Aps(i,j) =  0.5d0
             ELSEIF (i == 2) THEN
                col_Aps(i,j) = 0.d0
             ELSE
                k = rang*2*m_max_c + j
                IF (MOD(i,2)==0) THEN
                   col_Aps(i,j) = SIN(l*angle(k))
                ELSE
                   col_Aps(i,j) = COS(l*angle(k))
                END IF
             ENDIF
          ENDDO
       ENDDO
       IF (rang==0) THEN
          DO i= 1, 2*m_max
             col_Aps(i,2) = 0.d0
          ENDDO
       END IF
       col_Aps = (2.d0/(2*m_max-1))*col_Aps
    END IF
    !------END OF INITIALISATION--------------------

    t = MPI_WTIME()
    IF (choix==1) THEN
       sous_prod = 0
       DO proc = 1, nb_procs 
          DO i = 1, 2*m_max_c
             i_glob = (proc-1)*2*m_max_c + i
             DO n = 1, np
                DO j = 1, 2*m_max_c
                   l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                   sous_prod(l1) = sous_prod(l1) + col_Asp(i_glob,j)*V_in(j,n)
                END DO
             END DO
          END DO
       END DO
    ELSE
       sous_prod = 0
       DO proc = 1, nb_procs 
          DO i = 1, 2*m_max_c
             i_glob = (proc-1)*2*m_max_c + i
             DO n = 1, np
                DO j = 1, 2*m_max_c
                   l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                   sous_prod(l1) = sous_prod(l1) + col_Aps(i_glob,j)*V_in(j,n)
                END DO
             END DO
          END DO
       END DO
    END IF
    IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

    !Rangement
    !proc1i        -> V1_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c) |
    !procx         -> ...... |
    !proc(nb_procs)-> V1_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c)  
    !En esperant que ca marche?

    b_inf = rang*2*m_max_c*np + 1 
    b_sup = b_inf + 2*m_max_c*np - 1 
    res(b_inf:b_sup) = sous_prod(b_inf:b_sup) ! pour V_in

    !solution avec MPI_SENDRECV
    t = MPI_WTIME()
    DO proc =1, nb_procs 
       IF ((proc - 1) /= rang) THEN
          occur = 2*m_max_c*np 
          b_inf = (proc-1)*occur + 1
          b_sup = b_inf + occur - 1
          CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc-1, etiquette, &
               res(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc-1, etiquette, &
               MPI_COMM_WORLD, statut, code)
       ENDIF
    ENDDO
    IF(PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() - t

    t = MPI_WTIME()
    V_out = 0
    occur = 2*m_max_c*np
    DO n = 1, np 
       DO proc = 1, nb_procs
          DO i = 1, 2*m_max_c
             indice = i + (proc-1)*occur + (n-1)*2*m_max_c
             V_out(i,n) = V_out(i,n) + res(indice)
          END DO
       END DO
    ENDDO
    IF(PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

    IF (PRESENT(clean)) THEN
       IF(clean) THEN
          once =.TRUE.
          DEALLOCATE(angle)
          DEALLOCATE(col_Asp,col_Aps)
          DEALLOCATE(sous_prod)
          DEALLOCATE(res)
       END IF
    END IF

  END SUBROUTINE FFT_PAR_MP_b

  SUBROUTINE FFT_PROD(V1_in, V2_in, V_out)
    ! V1_in(m_max_c,np), idem v2_in, V_out

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:),  INTENT(OUT) :: V_out  

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    LOGICAL, SAVE :: once=.TRUE.
    INTEGER, SAVE :: m_max, m_max_c, rang, nb_procs, etiquette, np
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE  :: col_Asp, col_Aps, V1_p, V2_p
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: angle   ! valeurs des angles
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: sous_prod   !sous produits locaux   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: res  
    REAL(KIND=8)                           , SAVE  :: PI

    !pour l'envoie
    INTEGER :: proc_recep, proc_cible, b_inf, b_sup, indice, occur, proc
    INTEGER :: code, valeur, i, j, k, l, n, l1, l2, i_glob
    INTEGER :: requete0, requete1

    !--------END OF DECLARATION ------------------------------

    !--------INITIALISATION----------------------------------


    IF (once) THEN
       once = .FALSE.
       np = SIZE(V1_in,2)

       IF (MOD(SIZE(V1_in,1),2)/=0) THEN
          WRITE(*,*) ' BUG FFT_PROD'
          STOP
       END IF
       m_max_c = SIZE(V1_in,1)/2   !nombre de modes par processeur
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
       m_max = nb_procs*m_max_c !CONVENTION: m_max = mode_max +1
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
       etiquette = 100

       !Pour les colonnes col_cos et col_sin, je considere que le mode 0 se rattache a une partie
       !en cos et en sin pour systematiser les envois, 
       !on a donc une dimension de 2*mode_max +2 pour les colonnes !
       !Pour les angles, il n'y a qu'un angle pour le mode 0, et 2 pour les
       !autres modes, d'ou une dimension 2*mod_max + 1 pour le tableau 'angle'

       ALLOCATE(angle(2*m_max))
       ALLOCATE(col_Asp(2*m_max,2*m_max_c),col_Aps(2*m_max,2*m_max_c))
       ALLOCATE(sous_prod(2*m_max*2*np))
       ALLOCATE(res(2*m_max*2*np))
       ALLOCATE(V1_p(2*m_max_c,np), V2_p(2*m_max_c,np))
       !Calcul des angles qui sont les memes sur tous les processeurs
       PI = ACOS(-1.d0)
       !attention, theta = 0 correspond a i=0 !
       angle(1:2) = 0.d0  !angle(2) = angle bidon
       DO i= 3, 2*m_max
          angle(i) = 2.d0*PI/(2*m_max-1)*(i-2)
       ENDDO

       !Calcul des colonnes pour la fft spec -> phys

       DO i= 1, 2*m_max
          DO j= 1, 2*m_max_c
             l = rang*m_max_c + (j-1)/2 
             IF (MOD(j,2)==0) THEN
                col_Asp(i,j) = SIN(l*angle(i))
             ELSE
                col_Asp(i,j) = COS(l*angle(i))
             END IF
          ENDDO
       ENDDO
       DO j= 1, 2*m_max_c
          col_Asp(2,j) = 0.d0
       ENDDO

       !Calcul colonnes pour la fft phys -> spec

       DO i= 1, 2*m_max
          l = (i-1)/2 
          DO j= 1, 2*m_max_c
             IF (i == 1) THEN
                col_Aps(i,j) =  0.5d0
             ELSEIF (i == 2) THEN
                col_Aps(i,j) = 0.d0
             ELSE
                k = rang*2*m_max_c + j
                IF (MOD(i,2)==0) THEN
                   col_Aps(i,j) = SIN(l*angle(k))
                ELSE
                   col_Aps(i,j) = COS(l*angle(k))
                END IF
             ENDIF
          ENDDO
       ENDDO
       IF (rang==0) THEN
          DO i= 1, 2*m_max
             col_Aps(i,2) = 0.d0
          ENDDO
       END IF
       col_Aps = (2.d0/(2*m_max-1))*col_Aps

    END IF
    !------END OF INITIALISATION--------------------

    sous_prod = 0
    DO proc = 1, nb_procs 
       DO i = 1, 2*m_max_c
          i_glob = (proc-1)*2*m_max_c + i
          DO n = 1, np
             DO j = 1, 2*m_max_c
                l1 = (proc-1)*2*m_max_c*2*np + 2*(n-1)*2*m_max_c + i 
                sous_prod(l1) = sous_prod(l1) + col_Asp(i_glob,j)*V1_in(j,n)
                l2 = l1 + 2*m_max_c
                sous_prod(l2) = sous_prod(l2) + col_Asp(i_glob,j)*V2_in(j,n)
             END DO
          END DO
       END DO
    END DO

    !Rangement
    !proc1i        -> V1_p1(1:2*m_max_c) V2_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c) V2_pnp(1:2*m_max_c) |
    !procx         -> ...... |
    !proc(nb_procs)-> V1_p1(1:2*m_max_c) V2_p1(1:2*m_max_c) ... V1_pnp(1:2*m_max_c) V2_pnp(1:2*m_max_c)  
    !En esperant que ca marche?

    b_inf = rang*2*m_max_c*2*np + 1 
    b_sup = b_inf + 2*m_max_c*2*np - 1 
    res(b_inf:b_sup) = sous_prod(b_inf:b_sup) ! pour V1_in V2_in


    !solution avec MPI_SENDRECV
    DO proc =1, nb_procs 

       proc_cible = proc - 1 !vrai numero du processeur

       IF (proc_cible /= rang) THEN

          occur = 2*m_max_c*2*np 
          b_inf = proc_cible*occur + 1
          b_sup = b_inf + occur - 1
          CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
               res(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
               MPI_COMM_WORLD, statut, code)

       ENDIF
    ENDDO


    V1_p = 0
    V2_p = 0
    occur = 2*m_max_c*2*np
    DO n = 1, np 
       DO proc = 1, nb_procs
          DO i = 1, 2*m_max_c
             indice = i + (proc-1)*occur + 2*(n-1)*2*m_max_c
             V1_p(i,n) = V1_p(i,n) + res(indice)
             indice = indice + 2*m_max_c
             V2_p(i,n) = V2_p(i,n) + res(indice)
          END DO
       END DO
    ENDDO

    !Produit V1 * V2
    occur = 2*m_max_c*np
    DO n = 1, np 
       DO proc = 1, nb_procs
          DO i = 1, 2*m_max_c
             indice = 2*m_max*np + i + (proc-1)*occur + (n-1)*2*m_max_c
             sous_prod(indice) = V1_p(i,n)* V2_p(i,n) ! On utilise la seconde partie de sous_prod
             ! pour economiser la memoire
          END DO
       END DO
    ENDDO


    sous_prod(1:2*m_max*np) = 0 !Attention a ne pas ecraser la seconde partie du tableau.
    DO proc = 1, nb_procs 
       DO i = 1, 2*m_max_c
          i_glob = (proc-1)*2*m_max_c + i
          DO n = 1, np
             DO j = 1, 2*m_max_c 
                l1 = (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + i 
                indice = 2*m_max*np + (proc-1)*2*m_max_c*np + (n-1)*2*m_max_c + j
                sous_prod(l1) = sous_prod(l1) + col_Aps(i_glob,j)*sous_prod(indice)
             END DO
          END DO
       END DO
    END DO

    b_inf = rang*2*m_max_c*np + 1 
    b_sup = b_inf + 2*m_max_c*np - 1 
    res(b_inf:b_sup) = sous_prod(b_inf:b_sup)


    !solution avec MPI_SENDRECV
    DO proc = 1, nb_procs
       proc_cible = proc - 1 !vrai numero du processeur
       IF (proc_cible /= rang) THEN
          occur = 2*m_max_c*np
          b_inf = proc_cible*occur + 1
          b_sup = b_inf + occur - 1
          CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
               res(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
               MPI_COMM_WORLD, statut, code)
       ENDIF
    ENDDO

    V_out = 0
    occur = 2*m_max_c*np
    DO n = 1, np 
       DO proc = 1, nb_procs
          DO i = 1, 2*m_max_c
             indice = i + (proc-1)*occur + (n-1)*2*m_max_c
             V_out(i,n) = V_out(i,n) + res(indice)
          END DO
       END DO
    ENDDO


  END SUBROUTINE FFT_PROD

  SUBROUTINE FFT_PAR(V_in, V_out, choix_fft)
    !fft parallele spec -> phys si choix_fft = 1 et inverse sinon

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(:),  INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:),  INTENT(OUT) :: V_out  
    INTEGER,                     INTENT(IN)  :: choix_fft

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    LOGICAL, SAVE :: once=.TRUE.
    INTEGER, SAVE :: m_max, m_max_c, rang, nb_procs, etiquette
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE  :: col_Asp, col_Aps
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: angle   ! valeurs des angles
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: sous_prod   !sous produits locaux   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE  :: res  
    REAL(KIND=8)                           , SAVE  :: PI

    !pour l'envoie
    INTEGER :: proc_recep, proc_cible, b_inf, b_sup, indice, occur, proc
    INTEGER :: code, valeur, i, j, k, l, n
    INTEGER :: sol, requete0, requete1

    !--------END OF DECLARATION ------------------------------

    !--------INITIALISATION----------------------------------

    sol = 1

    IF (once) THEN
       once = .FALSE.
       IF (MOD(SIZE(V_in),2)/=0) THEN
          WRITE(*,*) ' BUG FF_PAR'
          STOP
       END IF
       m_max_c = SIZE(V_in)/2   !nombre de modes par processeur
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
       m_max = nb_procs*m_max_c !CONVENTION: m_max = mode_max +1
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
       etiquette = 100

       ALLOCATE(angle(2*m_max))
       ALLOCATE(col_Asp(2*m_max,2*m_max_c),col_Aps(2*m_max,2*m_max_c))
       ALLOCATE(sous_prod(2*m_max))
       ALLOCATE(res(2*m_max))

       !Calcul des angles qui sont les memes sur tous les processeurs
       PI = ACOS(-1.d0)
       !attention, theta = 0 correspond a i=0 !
       angle(1:2) = 0.d0  !angle(2) = angle bidon
       DO i= 3, 2*m_max
          angle(i) = 2.d0*PI/(2*m_max-1)*(i-2)
       ENDDO

       !Calcul des colonnes pour la fft spec -> phys

       DO i= 1, 2*m_max
          DO j= 1, 2*m_max_c
             l = rang*m_max_c + (j-1)/2 
             IF (MOD(j,2)==0) THEN
                col_Asp(i,j) = SIN(l*angle(i))
             ELSE
                col_Asp(i,j) = COS(l*angle(i))
             END IF
          ENDDO
       ENDDO
       DO j= 1, 2*m_max_c
          col_Asp(2,j) = 0.d0
       ENDDO

       !Calcul colonnes pour la fft phys -> spec

       DO i= 1, 2*m_max
          l = (i-1)/2 
          DO j= 1, 2*m_max_c
             IF (i == 1) THEN
                col_Aps(i,j) =  0.5d0
             ELSEIF (i == 2) THEN
                col_Aps(i,j) = 0.d0
             ELSE
                k = rang*2*m_max_c + j
                IF (MOD(i,2)==0) THEN
                   col_Aps(i,j) = SIN(l*angle(k))
                ELSE
                   col_Aps(i,j) = COS(l*angle(k))
                END IF
             ENDIF
          ENDDO
       ENDDO
       IF (rang==0) THEN
          DO i= 1, 2*m_max
             col_Aps(i,2) = 0.d0
          ENDDO
       END IF
       col_Aps = (2.d0/(2*m_max-1))*col_Aps
    END IF
    !------END OF INITIALISATION--------------------

    !calcul des sous produits locaux : 
    IF (choix_fft == 1) THEN  !spectral -> physique
       sous_prod = MATMUL(col_Asp,V_in)
    ELSE                      !physique -> spectral
       sous_prod = MATMUL(col_Aps,V_in)
    ENDIF

    b_inf = rang*2*m_max_c + 1
    b_sup = b_inf + 2*m_max_c - 1
    res(b_inf:b_sup) = sous_prod(b_inf:b_sup)

    !chaque processeur envoie la partie dont il n'a pas besoin pour le calul
    !de la fft et recoit la partie utile : 


    IF (sol == 1) THEN
       !solution avec MPI_SENDRECV
       DO proc =1, nb_procs 

          proc_cible = proc - 1 !vrai numero du processeur

          IF (proc_cible /= rang) THEN

             occur = 2*m_max_c 
             b_inf = proc_cible*occur + 1
             b_sup = b_inf + occur - 1
             CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
                  res(b_inf:b_sup), occur, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
                  MPI_COMM_WORLD, statut, code)

          ENDIF
       ENDDO

    ELSEIF (sol == 2) THEN
       STOP   
    ENDIF


    !Il faut maintenant additionner les valeurs du vecteurs sous_prod pour obtenir 
    !la fft pour les angles locaux a chaque proc
    V_out = 0
    occur = 2*m_max_c 
    DO i = 1, 2*m_max_c
       DO proc = 1, nb_procs 
          indice = i+(proc-1)*occur
          V_out(i) = V_out(i) + res(indice)
       END DO
    ENDDO

  END SUBROUTINE FFT_PAR

  SUBROUTINE smb_ns_gauss_sft_par(mesh,list_mode,V_in,V_out)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points     

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                                  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: V_out

    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G*mesh%me,SIZE(list_mode)) :: RotV, W
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc     
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray
    REAL(KIND=8)   :: user_time, tps, dummy
    REAL(KIND=8), DIMENSION(3)                  :: temps

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me


    tps = user_time(dummy)
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, me
          j_loc = jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !--------On calcul le rayon du point gauss
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(1,index,i) = SUM(Vs(:,1)*ww(:,l))
             W(3,index,i) = SUM(Vs(:,3)*ww(:,l))
             W(5,index,i) = SUM(Vs(:,5)*ww(:,l))

             W(2,index,i) = SUM(Vs(:,2)*ww(:,l))
             W(4,index,i) = SUM(Vs(:,4)*ww(:,l))
             W(6,index,i) = SUM(Vs(:,6)*ww(:,l))

             !-----------------rotational sur les points de Gauss---------------------------
             !coeff sur les cosinus 
             RotV(1,index,i) = mode/ray*W(6,index,i) &
                  -SUM(Vs(:,3)*dw_loc(2,:))
             RotV(4,index,i) =          SUM(Vs(:,2)*dw_loc(2,:)) &
                  -SUM(Vs(:,6)*dw_loc(1,:))
             RotV(5,index,i) =    1/ray*W(3,index,i) &
                  +SUM(Vs(:,3)*dw_loc(1,:)) &
                  -mode/ray*W(2,index,i)
             !coeff sur les sinus       
             RotV(2,index,i) =-mode/ray*W(5,index,i) &
                  -SUM(Vs(:,4)*dw_loc(2,:))
             RotV(3,index,i) =         SUM(Vs(:,1)*dw_loc(2,:)) &
                  -SUM(Vs(:,5)*dw_loc(1,:))
             RotV(6,index,i) =   1/ray*W(4,index,i) &
                  +SUM(Vs(:,4)*dw_loc(1,:))&
                  +mode/ray*W(1,index,i)
          ENDDO
       ENDDO
    END DO
    !tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    !tps = user_time(dummy)
    temps = 0
    CALL FFT_PAR_PROD_VECT(RotV, W, V_out, temps)  
    tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_ns_gauss_sft_par

  SUBROUTINE smb_cross_prod_gauss_sft_par(mesh,list_mode,V_in,V_out)
    !=================================
    !Calcul du second membre pour le probleme de Stokes
    !sans terme non lineaire et source (ff)

    USE Gauss_points     

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                                  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode    
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: V_out

    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: RotV, W
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc     
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray
    REAL(KIND=8)   :: user_time, tps, dummy
    REAL(KIND=8), DIMENSION(3)                  :: temps

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me


    tps = user_time(dummy)
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, me
          j_loc = jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !--------On calcul le rayon du point gauss
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Vs(:,1)*ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*ww(:,l))

             !-----------------rotational sur les points de Gauss---------------------------
             !coeff sur les cosinus 
             RotV(index,1,i) = mode/ray*W(index,6,i) &
                  -SUM(Vs(:,3)*dw_loc(2,:))
             RotV(index,4,i) =          SUM(Vs(:,2)*dw_loc(2,:)) &
                  -SUM(Vs(:,6)*dw_loc(1,:))
             RotV(index,5,i) =    1/ray*W(index,3,i) &
                  +SUM(Vs(:,3)*dw_loc(1,:)) &
                  -mode/ray*W(index,2,i)
             !coeff sur les sinus       
             RotV(index,2,i) =-mode/ray*W(index,5,i) &
                  -SUM(Vs(:,4)*dw_loc(2,:))
             RotV(index,3,i) =         SUM(Vs(:,1)*dw_loc(2,:)) &
                  -SUM(Vs(:,5)*dw_loc(1,:))
             RotV(index,6,i) =    1/ray*W(index,4,i) &
                  +SUM(Vs(:,4)*dw_loc(1,:))&
                  +mode/ray*W(index,1,i)
          ENDDO
       ENDDO
    END DO
    !tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    !tps = user_time(dummy)
    temps = 0
    CALL FFT_PAR_CROSS_PROD(RotV, W, V_out, temps)  
    tps = user_time(dummy) - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_cross_prod_gauss_sft_par



END MODULE sft_parallele
