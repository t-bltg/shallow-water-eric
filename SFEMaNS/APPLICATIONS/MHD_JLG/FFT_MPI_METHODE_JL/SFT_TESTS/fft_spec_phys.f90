PROGRAM fft_parallele
  
  IMPLICIT NONE
  include 'mpif.h'

  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

  INTEGER :: code, nb_procs, rang, n, valeur, i, j, l
  INTEGER :: mode, mode_max, compt
  INTEGER :: etiquette

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: col_cos, col_sin
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: angle   ! valeurs des angles
  REAL(KIND=8)                             :: PI
  REAL(KIND=8), DIMENSION(10)              :: tt

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: vect_spec   
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: vect_phys       
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: sous_prod   !sous produits locaux   
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: res  

  !pour l'envoie
  INTEGER        :: proc_recep, proc_cible, b_inf, b_sup, indice
  
!--------END OF DECLARATION ------------------------------
  
!--------INITIALISATION----------------------------------

  CALL MPI_INIT(code)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  mode_max = nb_procs-1
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  mode = rang 
  etiquette = 100
  
  tt(1) = MPI_WTIME()
!Pour les colonnes col_cos et col_sin, je considere que le mode 0 se rattache a une partie
!en cos et en sin pour systematiser les envois, 
!on a donc une dimension de 2*mode_max +2 pour les colonnes !
!Pour les angles, il n'y a qu'un angle pour le mode 0, et 2 pour les
!autres modes, d'ou une dimension 2*mod_max + 1 pour le tableau 'angle'

  ALLOCATE(angle(2*mode_max+1))
  ALLOCATE(col_cos(2*mode_max+2),col_sin(2*mode_max+2))
  col_cos = 0.d0 ; col_sin = 0.d0 ; angle = 0.d0
  ALLOCATE(sous_prod(2*mode_max+2))
  sous_prod = 0.d0 
  ALLOCATE(res(2*mode_max+2))
  res = 0.d0

  !Calcul des angles qui sont les memes sur tous les processeurs
  PI = ACOS(-1.d0)
 !attention, theta = 0 correspond a i=0 !
  DO i= 1, 2*mode_max+1
     angle(i) = 2.d0*PI/(2*mode_max+1)*(i-1)
  ENDDO
  !Calcul des colonnes pour la fft spec -> phys
  DO i= 1, 2*mode_max+2
     IF (i == 1) THEN
        col_cos(i) = cos(mode*angle(i))
        col_sin(i) = sin(mode*angle(i))
     ELSEIF (i > 2) THEN
        col_cos(i) = cos(mode*angle(i-1))
        col_sin(i) = sin(mode*angle(i-1))
     ENDIF
  ENDDO
!INITIALISATION POUR TEST FFT spec -> phys ---------------------
  ALLOCATE(vect_spec(2))
  ALLOCATE(vect_phys(2))
  IF (mode == 0) THEN
     vect_spec(1) = 0.d0
     vect_spec(2) = 0.d0
  ELSEIF (mode == 1) THEN
     vect_spec(1) = 1.d0
     vect_spec(2) = 1.d0
  ELSEIF (mode == 2) THEN
     vect_spec(1) = 1.d0
     vect_spec(2) = 1.d0 
  ENDIF
 
!pour les tests
 IF (mode == 0) THEN
   ! write(*,*) 'valeurs a obtenir'
    DO i=1, 2*mode_max + 1
   !    write(*,*) cos(angle(i)) + sin(angle(i)) &
   !  + cos(2*angle(i)) +sin(2*angle(i))
    ENDDO
 ENDIF

!------END OF INITIALISATION--------------------

!calcul des sous produits locaux : 
!sous_prod(i) = c_mode * cos(mode*theta_i) + s_mode * sin(mode*theta_i)
!attention au mode 0 qui decale tout ...
 
  DO i=1, 2*mode_max+2
     sous_prod(i) = vect_spec(1)*col_cos(i) + &
                    vect_spec(2)*col_sin(i)   
  ENDDO
  res(2*mode+1) = sous_prod(2*mode+1)
  res(2*mode+2) = sous_prod(2*mode+2)

!chaque processeur envoie la partie dont il n'a pas besoin pour le calul
!de la fft et recoit la partie utile : 
!le processeur qui s'occupe du mode 'mode' garde les valeur 
!sous_prod(2*mode+1) et sous_prod(2*mode+2) et envoie le reste.

!Chaque proc envoie a chacun des autres un doublet et recoit de lui un doublet
!a mettre au meme endroit que le doubler qu'il a envoye
  
  tt(2) = MPI_WTIME() 

  DO i=0, mode_max
     IF (i /= mode) THEN
        proc_cible = i
        b_inf = 2*i+1
        b_sup = 2*i+2
        
        CALL MPI_SENDRECV(sous_prod(b_inf:b_sup), 2, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
                          res(b_inf:b_sup), 2, MPI_DOUBLE_PRECISION, proc_cible, etiquette, &
                          MPI_COMM_WORLD, statut, code)
        !CALL MPI_SENDRECV_REPLACE(sous_prod(b_inf:b_sup), 2, MPI_DOUBLE_PRECISION, proc_cible,&
        !                  etiquette, proc_cible, etiquette, &
        !                  MPI_COMM_WORLD, statut, code)


     ENDIF
  ENDDO   
 
  tt(3) = MPI_WTIME()
  write(*,*) 'temps communication mode',mode,' = ', tt(3) - tt(2) 

!Il faut maintenant additionner les valeurs du vecteurs sous_prod pour obtenir 
!la fft pour les angles locaux a chaque proc
  DO i=0,mode_max 
     indice = 2*i+1
     vect_phys(1) = vect_phys(1) + res(indice)
     !vect_phys(1) = vect_phys(1) + sous_prod(indice) 
     indice = 2*i+2
     vect_phys(2) = vect_phys(2) + res(indice) 
     !vect_phys(2) = vect_phys(2) + sous_prod(indice)
  ENDDO
  !write(*,*) 'mode',mode,': f(theta_1)=',vect_phys(1),'et f(theta_2) =',vect_phys(2)

  tt(4) = MPI_WTIME()
  write(*,*) 'Temps total mode',mode,' = ',tt(4) - tt(1) 
  !write(*,*) 'mode',mode,'temps_comm / temps_tot = ',(tt(3)-tt(2))/(tt(4)-tt(1))

  CALL MPI_FINALIZE(code)

END PROGRAM fft_parallele
   
