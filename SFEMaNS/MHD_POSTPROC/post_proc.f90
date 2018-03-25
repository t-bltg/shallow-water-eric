MODULE post_processing

CONTAINS


  SUBROUTINE gradphi(phi_mesh, phi, list_mode, data_fichier, grphi)

     USE post_processing_vtk
     USE def_type_mesh
     USE chaine_caractere

     IMPLICIT NONE
     TYPE(mesh_type), TARGET                   :: phi_mesh
     REAL(KIND=8), DIMENSION(:,:,:),INTENT(IN) :: phi
     INTEGER,      DIMENSION(:),    INTENT(IN) :: list_mode
     CHARACTER(LEN=*),              INTENT(IN) :: data_fichier
     REAL(KIND=8), DIMENSION(:,:,:),INTENT(OUT):: grphi
     INTEGER,      DIMENSION(4)                :: ipar_sp
     REAL(KIND=8), DIMENSION(4)                :: fpar_sp
     INTEGER :: k

     OPEN (UNIT=21, FILE = data_fichier, FORM='formatted', STATUS='unknown')
     CALL read_until(21, 'data_solver_maxwell')
     DO k=1,4;      READ(21,*) ipar_sp(k)
     END DO
     DO k=1,2;      READ(21,*) fpar_sp(k)
     END DO
     CLOSE(21)
     !------------------------------------------------------------------------------

     CALL calcul_grad_champ_scal(phi_mesh, list_mode, phi, grphi, ipar_sp, fpar_sp)

  END SUBROUTINE gradphi

  SUBROUTINE val_ener(mesh, list_mode, v, e_mode, z_dom)

    USE Gauss_points
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ  
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode !energie par mode
    REAL(KIND=8), OPTIONAL,                     INTENT(IN) :: z_dom

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    REAL(KIND=8) :: x, ray
    INTEGER      :: i, k, j, l, m, ni
    INTEGER      :: m_max_c


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me   
    m_max_c = size(list_mode)
    e_mode = 0.d0

    !calcul de l'energie
    !===================

    DO j=1, m_max_c

       DO m = 1, me
! CN
       IF(PRESENT(z_dom)) THEN
        IF (z_dom.GT.0.d0) THEN  ! z_dom=+1 pour Hemisphere Nord; z_dom=-1 pour Hemisphere Sud
         IF(MINVAL(mesh%rr(2,mesh%jj(:,m))).LT.0.d0) CYCLE

         ELSE

         IF(MAXVAL(mesh%rr(2,mesh%jj(:,m))).GT.0.d0) CYCLE

         ENDIF
       ENDIF ! IF present(z_dom)

          DO l = 1, l_G
             x = 0.d0
             !--------On calcul le rayon du point gauss
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + mesh%rr(1,i)*ww(ni,l)
             END DO

             DO k= 1,size(v,2)
                x = x + SUM(v(jj(:,m),k,j)*ww(:,l))**2
             ENDDO

             e_mode(j)  = e_mode(j) + 0.5d0*x*ray*rj(l,m) 

          ENDDO

       END DO ! boucle m

    ENDDO ! boucle j

    !normalisation
    !=============
    DO i=1, m_max_c
       IF (list_mode(i) /= 0) THEN
          e_mode(i) = e_mode(i)/2.d0
       END IF
    END DO

  END SUBROUTINE val_ener

!  SUBROUTINE val_ener(mesh, list_mode, v, e_mode)
!
!    USE Gauss_points
!    USE def_type_mesh
!
!    IMPLICIT NONE
!    TYPE(mesh_type), TARGET                                :: mesh
!    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
!    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ  
!    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode !energie par mode
!
!    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
!    INTEGER,                      POINTER                  :: me
!    REAL(KIND=8) :: x, ray
!    INTEGER      :: i, k, j, l, m, ni
!    INTEGER      :: m_max_c
!
!
!    CALL gauss(mesh)
!    jj => mesh%jj
!    me => mesh%me   
!    m_max_c = size(list_mode)
!    e_mode = 0.d0
!
!    !calcul de l'energie
!    !===================
!
!    DO j=1, m_max_c
!
!       DO m = 1, me
!
!          DO l = 1, l_G
!
!             x = 0.d0
!
!             !--------On calcul le rayon du point gauss
!             ray = 0
!             DO ni = 1, n_w;  i = jj(ni,m)
!                ray = ray + mesh%rr(1,i)*ww(ni,l)
!             END DO
!
!             DO k= 1,size(v,2)
!                x = x + SUM(v(jj(:,m),k,j)*ww(:,l))**2
!             ENDDO
!
!             e_mode(j)  = e_mode(j) + 0.5d0*x*ray*rj(l,m) 
!
!          ENDDO
!
!       END DO
!
!    ENDDO
!
!    !normalisation
!    !=============
!    DO i=1, m_max_c
!       IF (list_mode(i) /= 0) THEN
!          e_mode(i) = e_mode(i)/2.d0
!       END IF
!    END DO
!
!  END SUBROUTINE val_ener

  SUBROUTINE val_ener_sym(mesh, list_mode, v, e_mode, e_mode_sym, e_mode_anti,type_sym)
  !type_sym = 1 pour un champ pair
  !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ  
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8),                               INTENT(IN) :: type_sym
    INTEGER, ALLOCATABLE, DIMENSION(:),         SAVE       :: point_sym
    LOGICAL,                                    SAVE       :: once =  .true.
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti      !champ  

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    REAL(KIND=8) :: x, ray
    INTEGER      :: i, k, j, l, m, ni
    INTEGER      :: m_max_c


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me   
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0
    
    IF (once) THEN 
       once = .false.
       ALLOCATE(point_sym(mesh%np))
       !determination des points symetriques les uns par rapport aux autres
       DO i= 1, mesh%np
          point_sym(i) = find_point(mesh,mesh%rr(1,i), -mesh%rr(2,i)) 
       ENDDO
       write(*,*) 'symetrie du maillage en r : ', MAXVAL(ABS(mesh%rr(1,:)-mesh%rr(1,point_sym(:))))
       write(*,*) 'symetrie du maillage en z : ', MAXVAL(ABS(mesh%rr(2,:)+mesh%rr(2,point_sym(:))))
    ENDIF
    

    !calcul du champ anti-symetrique
    champ_anti(:,:,:) = 0.5d0*(v(:,:,:) - type_sym*v(point_sym(:),:,:))
   
    !calcul des energies
    CALL val_ener(mesh, list_mode, v, e_mode)
    CALL val_ener(mesh, list_mode, champ_anti, e_mode_anti)
    e_mode_sym = e_mode - e_mode_anti


  END SUBROUTINE val_ener_sym

  

  SUBROUTINE calcul_point_sym(mesh, point_sym)
  !pour tous les points du maillages mesh, 
  !point_sym contient les numeros des points symetriques par
  !rapport a z=0
  !point_sym(i) = numero du point symetrique du point i
  USE Gauss_points
  USE def_type_mesh

  IMPLICIT NONE
  TYPE(mesh_type), TARGET                                :: mesh
  INTEGER, DIMENSION(mesh%np) ,               INTENT(OUT):: point_sym
  INTEGER                                                :: i  

  !determination des points symetriques les uns par rapport aux autres
  DO i= 1, mesh%np
      point_sym(i) = find_point(mesh,mesh%rr(1,i), -mesh%rr(2,i))
  ENDDO
  write(*,*) 'symetrie du maillage en r : ', MAXVAL(ABS(mesh%rr(1,:)-mesh%rr(1,point_sym(:))))
  write(*,*) 'symetrie du maillage en z : ', MAXVAL(ABS(mesh%rr(2,:)+mesh%rr(2,point_sym(:))))

  END SUBROUTINE calcul_point_sym



  SUBROUTINE val_ener_sym_glob(mesh, list_mode, v, point_sym,  e_mode, e_mode_sym, e_mode_anti,type_sym)
  !type_sym = 1 pour un champ pair
  !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ  
    INTEGER, DIMENSION(mesh%np),                INTENT(IN) :: point_sym
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8), DIMENSION(3),                 INTENT(IN) :: type_sym
    !type_sym(r,theta,z)
    LOGICAL,                                    SAVE       :: once =  .true.
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti      !champ  

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    REAL(KIND=8) :: x, ray
    INTEGER      :: i, k, j, l, m, ni
    INTEGER      :: m_max_c


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me   
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0
    
    !calcul du champ anti-symetrique
    ! CN bug le 11/07/07
    champ_anti(:,1:2,:) = 0.5d0*( v(:,1:2,:) - type_sym(1)*v(point_sym(:),1:2,:))
    champ_anti(:,3:4,:) = 0.5d0*( v(:,3:4,:) - type_sym(2)*v(point_sym(:),3:4,:))
    champ_anti(:,5:6,:) = 0.5d0*( v(:,5:6,:) - type_sym(3)*v(point_sym(:),5:6,:))
    write(*,*) 'champ anti_sym calcule'
   
    !calcul des energies
    CALL val_ener(mesh, list_mode, v, e_mode)
    CALL val_ener(mesh, list_mode, champ_anti, e_mode_anti)
!   e_mode_anti = e_mode_anti
    e_mode_sym = e_mode - e_mode_anti


  END SUBROUTINE val_ener_sym_glob
  
  SUBROUTINE val_ener_sym_centrale(mesh, list_mode, v, e_mode, e_mode_sym, e_mode_anti)
  !type_sym = 1 pour un champ pair
  !type_sym = -1 pour un champ impair

    USE Gauss_points
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type), TARGET                                :: mesh
    INTEGER, DIMENSION(:),                      INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:)             ,INTENT(IN) :: v      !champ  
    REAL(KIND=8), DIMENSION(SIZE(list_mode))   ,INTENT(OUT):: e_mode, e_mode_sym, e_mode_anti !energie par mode
    REAL(KIND=8),DIMENSION(3)                              :: type_sym
    !type_sym(r,theta,z)
    INTEGER, ALLOCATABLE, DIMENSION(:),         SAVE       :: point_sym
    LOGICAL,                                    SAVE       :: once =  .true.
    REAL(KIND=8), DIMENSION(mesh%np,size(v,2),size(list_mode))  :: champ_anti      !champ  

    INTEGER,      DIMENSION(:,:), POINTER                  :: jj
    INTEGER,                      POINTER                  :: me
    REAL(KIND=8) :: x, ray
    INTEGER      :: i, k, j, l, m, ni
    INTEGER      :: m_max_c


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me   
    m_max_c = size(list_mode)
    e_mode      = 0.d0
    e_mode_sym  = 0.d0
    e_mode_anti = 0.d0
    
    IF (once) THEN 
       once = .false.
       ALLOCATE(point_sym(mesh%np))
       !determination des points symetriques les uns par rapport aux autres
       DO i= 1, mesh%np
          point_sym(i) = find_point(mesh,mesh%rr(1,i), -mesh%rr(2,i)) 
       ENDDO
       write(*,*) 'symetrie du maillage en r : ', MAXVAL(ABS(mesh%rr(1,:)-mesh%rr(1,point_sym(:))))
       write(*,*) 'symetrie du maillage en z : ', MAXVAL(ABS(mesh%rr(2,:)+mesh%rr(2,point_sym(:))))
    ENDIF
    

    !calcul du champ anti-symetrique
    DO i=1, size(list_mode)
       IF (mod(list_mode(i),2) == 0) THEN  !mode pair
          type_sym(1) = 1.d0
          type_sym(2) = 1.d0
          type_sym(3) = -1.d0
       ELSE
          type_sym(1) = -1.d0
          type_sym(2) = -1.d0
          type_sym(3) = 1.d0
       ENDIF      
       champ_anti(:,1:2,i) =  0.5d0*(v(:,1:2,i) - type_sym(1)*v(point_sym(:),1:2,i))
       champ_anti(:,3:4,i) =  0.5d0*(v(:,3:4,i) - type_sym(2)*v(point_sym(:),3:4,i))
       champ_anti(:,5:6,i) =  0.5d0*(v(:,5:6,i) - type_sym(3)*v(point_sym(:),5:6,i))
    ENDDO   
    write(*,*) 'champ anti_sym calcule'
   
    !calcul des energies
    CALL val_ener(mesh, list_mode, v, e_mode)
    CALL val_ener(mesh, list_mode, champ_anti, e_mode_anti)
    e_mode_sym = e_mode - e_mode_anti


  END SUBROUTINE val_ener_sym_centrale
  
  SUBROUTINE trace_profile(mesh, v, it, freq_plot, list_mode, nom_champ)

    USE chaine_caractere
    USE sub_plot
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type)                                 :: mesh
    REAL(KIND=8)   , DIMENSION(:,:,:), INTENT(INOUT):: v
    INTEGER,                           INTENT(IN)   :: it, freq_plot
    INTEGER, DIMENSION(:),             INTENT(IN)   :: list_mode
    CHARACTER(len=3),                  INTENT(IN)   :: nom_champ

    INTEGER                           :: l, lblank
    CHARACTER(len=3)                  :: tit, tmode 
    CHARACTER(len=1)                  :: tcomp
    INTEGER                           :: i, k


    WRITE(tit,'(i3)') it/freq_plot
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO

    DO i= 1, SIZE(v,3)

       WRITE(tmode,'(i3)') list_mode(i)
       lblank = eval_blank(3,tmode)
       DO l = 1, lblank - 1
          tmode(l:l) = '0'
       END DO

       DO k= 1, size(v,2)
          WRITE(tcomp,'(i1)') k
          CALL plot_scalar_field(mesh%jj, mesh%rr, &
               v(:,k,i) , nom_champ//tcomp//'_m='//tmode//'_'//tit//'.plt')          
       ENDDO

    END DO

  END SUBROUTINE trace_profile

FUNCTION norme_L2_champ_localise_par(mesh, list_mode, v, xmin, xmax, ymin, ymax) RESULT(norm)

!calcul la norme L2 du champ v dans le carre [xmin:xmax] x [ymin:ymax]

    USE def_type_mesh
    USE fem_tn_axi
    USE tn_parallele

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v
    REAL(KIND=8),                    INTENT(IN) :: xmin, xmax, ymin, ymax
    REAL(KIND=8), DIMENSION(size(v,1),size(v,2),size(v,3)) :: v_loc

    REAL(KIND=8)                        :: err1, s1, norm
    INTEGER                             :: k, nn

    !restriction du champ v au carre [xmin:xmax] x [ymin:ymax] dans v_loc
    DO nn = 1, mesh%np
       IF ((mesh%rr(1,nn) .le. xmax) .and. (mesh%rr(1,nn) .ge. xmin) &
           .and. (mesh%rr(2,nn) .le. ymax) .and. (mesh%rr(2,nn) .ge. ymin)) THEN
           !test sur le format
            IF (SIZE(v,2)==mesh%np) THEN   
               v_loc(:,nn,:) = v(:,nn,:)
            ELSE
               v_loc(nn,:,:) = v(nn,:,:)
            ENDIF   
       ELSE
            !test sur le format
            IF (SIZE(v,2)==mesh%np) THEN
               v_loc(:,nn,:) = 0.d0
            ELSE
               v_loc(nn,:,:) = 0.d0
            ENDIF
       ENDIF
    ENDDO

    err1 = 0.d0
    s1 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO k = 1, SIZE(list_mode)
          DO nn = 1,SIZE(v,1)
             s1 = 0.d0
             CALL ns_0(mesh , v_loc(nn,:,k), s1)
             err1 = err1 + s1**2
          END DO
       ENDDO
    ELSE
       DO k = 1, SIZE(list_mode)
          DO nn = 1,SIZE(v,2)
             s1 = 0.d0
             CALL ns_0(mesh , v_loc(:,nn,k), s1)
             err1 = err1 + s1**2
          END DO
       ENDDO
    END IF
    !norm = SQRT(err1)
    CALL integration_mode(SIZE(list_mode),err1, norm)
    norm = SQRT(norm)

  END FUNCTION norme_L2_champ_localise_par

  FUNCTION find_point(mesh,r,z)  RESULT(n)
  
  USE def_type_mesh
  USE fem_tn_axi
  USE tn_parallele

  IMPLICIT NONE
  TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
  REAL(KIND=8),                    INTENT(IN) :: r, z
  INTEGER                                     :: n, i, j
  REAL(KIND=8)                                :: dist, dist_min
  REAL(KIND=8), DIMENSION(1) :: jlg

  jlg = MINLOC((mesh%rr(1,:)-r)**2 + (mesh%rr(2,:)-z)**2)
  n = jlg(1)
  RETURN

  !June 27, 2008, JLG
  !dist_min = 1.d0
  !DO i =1, mesh%np
  !   dist = sqrt((mesh%rr(1,i)-r)**2 + (mesh%rr(2,i)-z)**2)
  !   IF (dist .LE. dist_min) THEN
  !      dist_min = dist
  !      n = i
  !   ENDIF
  !ENDDO 
  END FUNCTION find_point 

END MODULE post_processing
