!
!Authors Jean-Luc Guermond, Copyrights 1996
!
MODULE  periodic


  IMPLICIT NONE

  PUBLIC :: prep_periodic_bc, &
       prep_periodic, &
       prep_periodic_bloc, &
       prep_periodic_mhd_bc, &
       st_csr_periodic,  &
       st_csr_per,  &
       bc_periodic_M,    & 
       bc_per_M,    & 
       bc_periodic, &
       bc_per, &
       prep_periodic_H_p_phi_bc
  PRIVATE

CONTAINS

  SUBROUTINE prep_periodic_bc(dir, fil, mesh)
    !=========================================

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(len=64),       INTENT(IN) :: dir,fil
    TYPE(mesh_type)                     :: mesh

    INTEGER,      DIMENSION(:), POINTER :: list_loc, perlist_loc 
    INTEGER                             :: n, d_end, f_end, side1, side2
    REAL(KIND=8), DIMENSION(:), POINTER :: e

    ALLOCATE(e(SIZE(mesh%rr,1)))

    WRITE (*,*) 'Loading periodic-data file ...'

    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    OPEN  (30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')

    CALL read_until(30, 'data_periodic')
    READ  (30, *)  mesh%periodic%n_bord
    IF (mesh%periodic%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, mesh%periodic%n_bord
       READ  (30, *) side1, side2, e
       CALL list_periodic(mesh%nps, mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
            list_loc, perlist_loc)
       ALLOCATE (mesh%periodic%list(n)%DIL(SIZE(list_loc)), mesh%periodic%perlist(n)%DIL(SIZE(list_loc)))
       mesh%periodic%list(n)%DIL = list_loc
       mesh%periodic%perlist(n)%DIL = perlist_loc
    END DO

    CLOSE(30)
    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE prep_periodic_bc

  !jan 29 2007
  SUBROUTINE prep_periodic(dir, fil, mesh, periodic)
    !=========================================

    USE def_type_mesh
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(len=64),       INTENT(IN) :: dir,fil
    TYPE(mesh_type)                     :: mesh
    TYPE(periodic_type)                 :: periodic

    INTEGER,      DIMENSION(:), POINTER :: list_loc, perlist_loc 
    INTEGER                             :: n, d_end, f_end, side1, side2
    REAL(KIND=8), DIMENSION(:), POINTER :: e

    ALLOCATE(e(SIZE(mesh%rr,1)))
    WRITE (*,*) 'Loading periodic-data file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    OPEN  (30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    CALL read_until(30, 'data_periodic')
    READ  (30, *)  periodic%n_bord
    IF (periodic%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, periodic%n_bord
       READ  (30, *) side1, side2, e
       CALL list_periodic(mesh%nps, mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
            list_loc, perlist_loc)

       ALLOCATE (periodic%list(n)%DIL(SIZE(list_loc)), periodic%perlist(n)%DIL(SIZE(list_loc)))
       periodic%list(n)%DIL = list_loc
       periodic%perlist(n)%DIL = perlist_loc

       DEALLOCATE(list_loc,perlist_loc)
    END DO

    CLOSE(30)
    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE prep_periodic

  SUBROUTINE prep_periodic_bloc(dir, fil, mesh, periodic, nb_bloc)
    !=========================================
    USE chaine_caractere
    USE def_type_mesh
    IMPLICIT NONE
    CHARACTER(len=64),       INTENT(IN) :: dir,fil
    TYPE(mesh_type)                     :: mesh
    TYPE(periodic_type)                 :: periodic
    INTEGER,                 INTENT(IN) :: nb_bloc

    INTEGER,      DIMENSION(:), POINTER :: list_loc, perlist_loc
    INTEGER                             :: n, d_end, f_end, side1, side2, nsize, n_b
    INTEGER                             :: k, k_deb, k_fin 
    REAL(KIND=8), DIMENSION(2)          :: e

    WRITE (*,*) 'Loading periodic-data file ...'

    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    OPEN  (30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')

    CALL read_until(30, 'data_periodic')
    READ  (30, *)  periodic%n_bord
    IF (periodic%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, periodic%n_bord

       READ  (30, *) side1, side2, e

       CALL list_periodic(mesh%nps, mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
            list_loc, perlist_loc)

       !n_b = SIZE(list_loc)
       n_b = SIZE(perlist_loc)
       nsize = nb_bloc*n_b !SIZE(list_loc) !n_b


       ALLOCATE(periodic%list(n)%DIL(nsize), periodic%perlist(n)%DIL(nsize))

       DO k = 1, nb_bloc
          k_deb=(k-1)*n_b+1
          k_fin=k*n_b
          periodic%list(n)%DIL(k_deb:k_fin)    = list_loc(1:n_b)    + (k-1)*mesh%np ! First bloc
          periodic%perlist(n)%DIL(k_deb:k_fin) = perlist_loc(1:n_b) + (k-1)*mesh%np ! First bloc
       END DO

       DEALLOCATE(list_loc,perlist_loc)

    END DO

    CLOSE(30)
    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE prep_periodic_bloc
  !jan 29 2007

  !JLG+FL/Feb 2 2010
  SUBROUTINE prep_periodic_H_p_phi_bc(data_fichier, H_mesh, pmag_mesh, phi_mesh, H_p_phi_per)
    USE chaine_caractere
    USE def_type_mesh
    IMPLICIT NONE
    CHARACTER(len=64),       INTENT(IN) :: data_fichier
    TYPE(mesh_type)                     :: H_mesh, pmag_mesh, phi_mesh
    TYPE(periodic_type)                 :: H_p_phi_per

    INTEGER,      DIMENSION(:), POINTER :: b_list_loc, b_perlist_loc, &
         e_list_loc, e_perlist_loc, p_list_loc, p_perlist_loc
    INTEGER                             :: n, d_end, f_end, side1, side2, nsize, n_b, n_e, n_p
    REAL(KIND=8), DIMENSION(2)          :: e

    WRITE (*,*) 'Loading periodic-data file ...'

 
    OPEN  (30,FILE=data_fichier,FORM='formatted')

    CALL read_until(30, 'data_periodic')
    READ  (30, *)  H_p_phi_per%n_bord
    IF (H_p_phi_per%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, H_p_phi_per%n_bord

       READ  (30, *) side1, side2, e

       CALL list_periodic(H_mesh%nps, H_mesh%np, H_mesh%jjs, H_mesh%sides, H_mesh%rr, side1, side2, e, &
            b_list_loc, b_perlist_loc)

       CALL list_periodic(pmag_mesh%nps, pmag_mesh%np, pmag_mesh%jjs, pmag_mesh%sides, pmag_mesh%rr, side1, side2, e, &
            p_list_loc, p_perlist_loc)

       IF (phi_mesh%nps==0) THEN
          ALLOCATE(e_list_loc(1))
          n_e = 0
       ELSE
          CALL list_periodic(phi_mesh%nps, phi_mesh%np, phi_mesh%jjs, phi_mesh%sides, phi_mesh%rr, side1, side2, e, &
               e_list_loc, e_perlist_loc)
          n_e = SIZE(e_list_loc)
       END IF

       n_b = SIZE(b_list_loc)
       n_p = SIZE(p_list_loc)
       !n_e = SIZE(e_list_loc) 
       nsize = 3*n_b + n_p + n_e

       ALLOCATE(H_p_phi_per%list(n)%DIL(nsize), H_p_phi_per%perlist(n)%DIL(nsize))

       H_p_phi_per%list(n)%DIL(1:n_b)    = b_list_loc    ! First block
       H_p_phi_per%perlist(n)%DIL(1:n_b) = b_perlist_loc ! First block

       H_p_phi_per%list(n)%DIL(n_b+1:2*n_b)    = b_list_loc    + H_mesh%np ! Second block
       H_p_phi_per%perlist(n)%DIL(n_b+1:2*n_b) = b_perlist_loc + H_mesh%np ! Second block

       H_p_phi_per%list(n)%DIL(2*n_b+1:3*n_b)    = b_list_loc    + 2*H_mesh%np ! Third block
       H_p_phi_per%perlist(n)%DIL(2*n_b+1:3*n_b) = b_perlist_loc + 2*H_mesh%np ! Third block

       H_p_phi_per%list(n)%DIL(3*n_b+1:3*n_b+n_p)    = p_list_loc    + 3*H_mesh%np ! Forth block
       H_p_phi_per%perlist(n)%DIL(3*n_b+1:3*n_b+n_p) = p_perlist_loc + 3*H_mesh%np ! Forth block

       IF (n_e/=0) THEN
          H_p_phi_per%list(n)%DIL(3*n_b+n_p+1:)    = e_list_loc    + 3*H_mesh%np + pmag_mesh%np ! Fourth block
          H_p_phi_per%perlist(n)%DIL(3*n_b+n_p+1:) = e_perlist_loc + 3*H_mesh%np + pmag_mesh%np ! Fourth block
       END IF

       DEALLOCATE(b_list_loc, b_perlist_loc, e_list_loc, e_perlist_loc, p_list_loc, p_perlist_loc)
    END DO

    CLOSE(30)
    WRITE (*,*) 'Treatment of periodic-data done'

  END SUBROUTINE Prep_periodic_H_p_phi_bc

  SUBROUTINE prep_periodic_mhd_bc(dir, fil, H_mesh, phi_mesh, H_phi_mesh)
    !jan 29
    !=========================================
    USE chaine_caractere
    USE def_type_mesh
    IMPLICIT NONE
    CHARACTER(len=64),       INTENT(IN) :: dir,fil
    TYPE(mesh_type)                     :: H_mesh, phi_mesh
    TYPE(periodic_type)                 :: H_phi_mesh

    INTEGER,      DIMENSION(:), POINTER :: b_list_loc, b_perlist_loc, &
         e_list_loc, e_perlist_loc
    INTEGER                             :: n, d_end, f_end, side1, side2, nsize, n_b, n_e
    REAL(KIND=8), DIMENSION(2)          :: e

    WRITE (*,*) 'Loading periodic-data file ...'

    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    OPEN  (30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')

    CALL read_until(30, 'data_periodic')
    READ  (30, *)  H_phi_mesh%n_bord
    IF (H_phi_mesh%n_bord .GT. 20) THEN
       WRITE(*,*) 'PREP_MESH_PERIODIC: trop de bords periodiques'
       STOP
    END IF

    DO n= 1, H_phi_mesh%n_bord

       READ  (30, *) side1, side2, e

       CALL list_periodic(H_mesh%nps, H_mesh%np, H_mesh%jjs, H_mesh%sides, H_mesh%rr, side1, side2, e, &
            b_list_loc, b_perlist_loc)

       IF (phi_mesh%nps==0) THEN
          ALLOCATE(e_list_loc(1))
          n_e = 0
       ELSE
          CALL list_periodic(phi_mesh%nps, phi_mesh%np, phi_mesh%jjs, phi_mesh%sides, phi_mesh%rr, side1, side2, e, &
               e_list_loc, e_perlist_loc)
          n_e = SIZE(e_list_loc)
       END IF

       n_b = SIZE(b_list_loc)
       !n_e = SIZE(e_list_loc) 
       nsize = 3*n_b + n_e

       ALLOCATE(H_phi_mesh%list(n)%DIL(nsize), H_phi_mesh%perlist(n)%DIL(nsize))

       H_phi_mesh%list(n)%DIL(1:n_b)    = b_list_loc    ! First bloc
       H_phi_mesh%perlist(n)%DIL(1:n_b) = b_perlist_loc ! First bloc

       H_phi_mesh%list(n)%DIL(n_b+1:2*n_b)    = b_list_loc    + H_mesh%np ! Second bloc
       H_phi_mesh%perlist(n)%DIL(n_b+1:2*n_b) = b_perlist_loc + H_mesh%np ! Second bloc

       H_phi_mesh%list(n)%DIL(2*n_b+1:3*n_b)    = b_list_loc    + 2*H_mesh%np ! Third bloc
       H_phi_mesh%perlist(n)%DIL(2*n_b+1:3*n_b) = b_perlist_loc + 2*H_mesh%np ! Third bloc

       IF (n_e/=0) THEN
          H_phi_mesh%list(n)%DIL(3*n_b+1:)    = e_list_loc    + 3*H_mesh%np  ! Fourth bloc
          H_phi_mesh%perlist(n)%DIL(3*n_b+1:) = e_perlist_loc + 3*H_mesh%np  ! Fourth bloc
       END IF

       DEALLOCATE(b_list_loc, b_perlist_loc, e_list_loc, e_perlist_loc)
    END DO

    CLOSE(30)
    WRITE (*,*) 'Treatment of periodic-data done'

    !jan 29
  END SUBROUTINE prep_periodic_mhd_bc
  !END SUBROUTINE prep_periodic_bc
  !jan 29


  SUBROUTINE st_csr_periodic(mesh, ia, ja, ia_per, ja_per)
    !=================================================
    USE def_type_mesh
    IMPLICIT NONE

    INTEGER,      DIMENSION(:), POINTER     :: ia, ja, ia_per, ja_per
    TYPE(mesh_type)                         :: mesh

    CALL st_csr_prd(mesh%periodic%n_bord, mesh%periodic%list, mesh%periodic%perlist, &
         mesh%periodic%pnt, ia, ja, ia_per, ja_per) 

  END SUBROUTINE st_csr_periodic

  SUBROUTINE st_csr_per(periodic, ia, ja, ia_per, ja_per)
    !=================================================
    USE def_type_mesh
    IMPLICIT NONE

    INTEGER,      DIMENSION(:), POINTER     :: ia, ja, ia_per, ja_per
    TYPE(periodic_type)                     :: periodic 

    CALL st_csr_prd(periodic%n_bord, periodic%list, periodic%perlist, &
         periodic%pnt, ia, ja, ia_per, ja_per) 

  END SUBROUTINE st_csr_per

  SUBROUTINE bc_periodic_M(mesh, ia, ja, aa, ia_per, ja_per, aa_per) 
    !=============================================================
    USE def_type_mesh
    IMPLICIT NONE

    INTEGER,             DIMENSION(:), POINTER       :: ia_per, ja_per, ia, ja
    REAL(KIND=8),        DIMENSION(:), INTENT(IN)    :: aa 
    REAL(KIND=8),        DIMENSION(:), INTENT(OUT)   :: aa_per 
    TYPE(mesh_type)                         :: mesh

    CALL bc_prd_M(mesh%periodic%n_bord, mesh%periodic%list, mesh%periodic%perlist, &
         mesh%periodic%pnt, ia_per, ja_per, aa_per, ia, ja, aa)

  END SUBROUTINE bc_periodic_M

  !SUBROUTINE bc_per_M(periodic, ia, ja, aa, ia_per, ja_per, aa_per) 
  SUBROUTINE bc_per_M(periodic, mat, mat_per) 
    !=============================================================
    USE def_type_mesh
    USE matrix_type
    IMPLICIT NONE

    TYPE(matrice_bloc),                 INTENT(IN)     :: mat
    TYPE(matrice_bloc),                 INTENT(OUT)    :: mat_per
    TYPE(periodic_type)                              :: periodic

    CALL bc_prd_M(periodic%n_bord, periodic%list, periodic%perlist, &
         periodic%pnt, mat_per%ia, mat_per%ja, mat_per%aa, mat%ia, mat%ja, mat%aa)

  END SUBROUTINE bc_per_M

  SUBROUTINE bc_periodic(mesh, ff)
    !=========================
    USE def_type_mesh
    IMPLICIT NONE

    REAL(KIND=8),        DIMENSION(:), INTENT(INOUT)   :: ff
    TYPE(mesh_type)                         :: mesh

    CALL bc_prd(mesh%periodic%n_bord, mesh%periodic%list, mesh%periodic%perlist, ff)

  END SUBROUTINE bc_periodic

  SUBROUTINE bc_per(periodic, ff)
    !=========================
    USE def_type_mesh
    IMPLICIT NONE

    REAL(KIND=8),        DIMENSION(:), INTENT(INOUT)   :: ff
    TYPE(periodic_type)                                :: periodic

    CALL bc_prd(periodic%n_bord, periodic%list, periodic%perlist, ff)

  END SUBROUTINE bc_per

  SUBROUTINE list_periodic(nps,np, jjs, sides, rr, side1, side2, e, list_out, perlist_out)
    !============================================================================

    IMPLICIT NONE
    INTEGER,                      INTENT(IN)  :: nps, np 
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jjs 
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: sides 
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: rr
    INTEGER,                      INTENT(IN)  :: side1, side2 
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: e 
    INTEGER,      DIMENSION(:),   POINTER     :: list_out, perlist_out
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: list, perlist

    LOGICAL,      DIMENSION(np)               :: virgin 
    REAL(KIND=8), DIMENSION(SIZE(rr,1))       :: ri
    INTEGER :: ms, ns, i, j, long, inter
    REAL(KIND=8) :: r, epsilon = 1.d-9
    LOGICAL :: verif

    IF (ALLOCATED(list))    DEALLOCATE(list)
    IF (ALLOCATED(perlist)) DEALLOCATE(perlist)

    ALLOCATE (list(nps), perlist(nps))
    virgin = .TRUE.

    i = 0; j=0
    DO ms = 1, SIZE(sides)

       IF (sides(ms) .EQ. side1) THEN
          DO ns = 1, SIZE(jjs,1)
             IF (virgin(jjs(ns,ms))) THEN
                i = i + 1
                list(i) = jjs(ns,ms)
                virgin(jjs(ns,ms)) = .FALSE.
             END IF
          END DO
       ELSE IF (sides(ms) .EQ. side2) THEN
          DO ns = 1, SIZE(jjs,1)
             IF (virgin(jjs(ns,ms))) THEN
                j = j + 1
                perlist(j) = jjs(ns,ms)
                virgin(jjs(ns,ms)) = .FALSE.
             END IF
          END DO

       END IF

    END DO

    IF (i .NE. j) THEN
       WRITE(*,*) ' FEM_PERIODIC: side1 and side2 have', &
            ' different numbers of points' 
       STOP
    END IF
    long = i 

    DO i = 1, long 
       ri = rr(:,list(i))+e(:)
       verif = .FALSE.
       !if (i==2) stop
       DO j = i, long
          r = SUM(ABS(ri - rr(:,perlist(j))))
          !if (i==1) write(*,*) ' r',r,'j',  j
          IF (r .LE. epsilon ) THEN 
             inter = perlist(i)
             perlist(i) = perlist(j) 
             perlist(j) = inter 
             verif = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.verif) THEN
          WRITE(*,*) ' BUG dans  data_periodic ou le maillage:', &
               ' side1 + e /= side2'
          WRITE(*,*) ' i = ', i 
          !         STOP
       END IF
    END DO

    ALLOCATE (list_out(long))
    list_out(1:long) = list(1:long) 
    ALLOCATE (perlist_out(long))
    perlist_out(1:long) = perlist(1:long)


  END SUBROUTINE list_periodic

  SUBROUTINE st_csr_prd(n_bord, list, perlist, pnt, ia_in ,ja_in , ia, ja) 
    !=================================================================
    USE dyn_line
    IMPLICIT NONE

    INTEGER,                                     INTENT(IN) :: n_bord
    TYPE(dyn_int_line), DIMENSION(:),            INTENT(IN) :: list, perlist
    INTEGER,            DIMENSION(:),            POINTER    :: pnt
    INTEGER,            DIMENSION(:),            INTENT(IN) :: ia_in, ja_in
    INTEGER,            DIMENSION(:),            POINTER    :: ia 
    INTEGER,            DIMENSION(:),            POINTER    :: ja 

    LOGICAL,            DIMENSION(SIZE(ia_in)-1)            :: virgin, vierge, touch 
    INTEGER,            DIMENSION(SIZE(ia_in)-1)            :: tsil, tsilrep
    INTEGER,            DIMENSION(SIZE(ia_in))              :: ia_new 
    INTEGER,            DIMENSION(:),            POINTER    :: ja_new, pnt_new 

    INTEGER      :: n_syst, i, pi, l, n, nz, p, p2, ki, kpi, kim, kpim
    LOGICAL      :: NOTDONE

    IF (ASSOCIATED(pnt)) NULLIFY(pnt)
    IF (ASSOCIATED(ia))  NULLIFY(ia)
    IF (ASSOCIATED(ja))  NULLIFY(ja)

    ALLOCATE(ia(SIZE(ia_in)),ja(SIZE(ja_in)), pnt(SIZE(ja_in)))

    ia = ia_in
    ja = ja_in
    DO p = 1, SIZE(ja_in)
       pnt(p) = p
    END DO

    p = 0
    DO n = 1, n_bord
       p = p+SIZE(list(n)%DIL)
    END DO
    p = 50*p
    ALLOCATE(ja_new(SIZE(ja)+p), pnt_new(SIZE(ja)+p))

    n_syst = SIZE(ia) - 1
    ia_new(1) = ia(1)
    touch = .FALSE.

    DO n = 1, n_bord

       virgin = .TRUE.
       vierge = .TRUE.

       DO l = 1, SIZE(list(n)%DIL)
          virgin(list(n)%DIL(l)) = .FALSE.    ! This node is periodic 
          vierge(perlist(n)%DIL(l)) = .FALSE. ! This node is periodic 
          tsil(list(n)%DIL(l))   = l 
       END DO

       DO i = 1, n_syst 

          IF (virgin(i)) THEN ! This node is not on periodic boundary of type 1
             p2 = ia_new(i)
             IF (vierge(i)) THEN ! This node is not/periodic boundary of type 2
                DO p = ia(i), ia(i+1) - 1
                   ja_new(p2) = ja(p)
                   pnt_new(p) = p2
                   p2 = p2 + 1
                END DO
             ELSE                ! This node is on a periodic boundary of type 2
                p2 = ia_new(i)
                DO p = ia(i), ia(i+1) - 1
                   ja_new(p2) = ja(p)
                   p2 = p2 + 1
                END DO
             END IF
             ia_new(i+1) = p2
          ELSE IF (touch(i)) THEN           ! This node is periodic but touched
             pi = perlist(n)%DIL(tsil(i))   ! Don't do anything
             p2 = ia_new(i)
             DO p = ia(i), ia(i+1) - 1      ! It is a corner 
                ja_new(p2) = ja(p) 
                p2 = p2 + 1
             END DO
             ia_new(i+1) = p2
          ELSE                              ! This node is periodic not touched
             pi = perlist(n)%DIL(tsil(i))
             !  On ajoute les deux lignes en ordonnant les indices
             ki = ia(i);     kim = ia(i+1)
             kpi = ia(pi);   kpim= ia(pi+1)
             p = ia_new(i) 
             DO WHILE((ki.LT.kim) .OR. (kpi.LT.kpim))  
                IF (ki .EQ. kim) THEN
                   ja_new(p)    = ja(kpi)
                   pnt_new(kpi) = p
                   kpi = kpi + 1
                   p = p + 1
                ELSE IF (kpi.EQ.kpim) THEN
                   ja_new(p)    = ja(ki)
                   pnt_new(ki) = p
                   ki = ki + 1
                   p = p + 1
                ELSE IF (ja(kpi) .LT. ja(ki)) THEN
                   ja_new(p)    = ja(kpi)
                   pnt_new(kpi) = p
                   kpi = kpi + 1
                   p = p + 1
                ELSE IF(ja(ki) .LT. ja(kpi)) THEN    
                   ja_new(p)   = ja(ki)
                   pnt_new(ki) = p
                   ki = ki + 1
                   p = p + 1
                ELSE   ! Il y a egalite
                   ja_new(p)   = ja(ki)
                   pnt_new(ki) = p
                   pnt_new(kpi) = p 
                   ki = ki + 1
                   kpi = kpi + 1
                   p = p + 1    
                END IF
             END DO
             ia_new(i+1) = p 
          END IF
       END DO

       ia = ia_new
       nz = ia(n_syst+1)-ia(1)
       NULLIFY(ja); ALLOCATE(ja(nz))
       ja = ja_new(1:nz)
       DO p = 1, SIZE(ja_in)
          p2 = pnt_new(pnt(p))
          pnt(p) = p2
       END DO

       virgin = .TRUE.

       DO l = 1, SIZE(perlist(n)%DIL)
          virgin(perlist(n)%DIL(l)) = .FALSE.   ! This node is periodic 
          tsilrep(perlist(n)%DIL(l)) = l 
       END DO

       DO i = 1, n_syst

          IF (virgin(i)) THEN               ! This node is not periodic
             p2 = ia_new(i)
             DO p = ia(i), ia(i+1) - 1
                ja_new(p2) = ja(p)
                pnt_new(p) = p2
                p2 = p2 + 1
             END DO
             ia_new(i+1) = p2
          ELSE                              ! This node is periodic
             pi = list(n)%DIL(tsilrep(i))
             p2 = ia_new(i)
             IF (touch(i)) THEN             ! This line already touched
                NOTDONE = .TRUE.
                DO p = ia(i), ia(i+1) - 1   ! add column pi only
                   IF ((ja(pi) .LT. ja(p)) .AND. NOTDONE) THEN
                      ja_new(p2) = pi
                      ja_new(p2+1) = ja(p)
                      p2 = p2 + 2
                      NOTDONE = .FALSE.
                   ELSE
                      ja_new(p2) = ja(p)
                      p2 = p2 + 1
                   END IF
                END DO
                p = ia(i+1) - 1
                IF (ja(pi) .GT. ja(p) ) THEN ! Don't forget column pi
                   ja_new(p2) = pi
                   p2 = p2 + 1
                END IF
                ia_new(i+1) = p2
             ELSE                           ! This line not touched yet
                p = ia_new(i)               ! put columns i and pi only
                ia_new(i+1) = p + 2
                touch(i) = .TRUE.
                IF (i .LT. pi) THEN
                   ja_new(p) = i
                   ja_new(p+1) = pi 
                ELSE
                   ja_new(p) = pi 
                   ja_new(p+1) = i 
                END IF
             END IF
          END IF

       END DO

       ia = ia_new
       nz = ia(n_syst+1)-ia(1)
       NULLIFY(ja)
       ALLOCATE(ja(nz))
       ja = ja_new(1:nz)
       DO p = 1, SIZE(ja_in) 
          p2 = pnt_new(pnt(p))
          pnt(p) = p2
       END DO
    END DO

    NULLIFY(ja_new, pnt_new)

  END SUBROUTINE st_csr_prd


  SUBROUTINE bc_prd_M(n_bord, list, perlist, pnt, ia_per, ja_per, aa_per, ia, ja, aa) 
    !=====================================================================
    USE dyn_line
    IMPLICIT NONE
    INTEGER,                           INTENT(IN)    :: n_bord
    TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)    :: list, perlist
    INTEGER,             DIMENSION(:), INTENT(IN)    :: pnt, ia_per, ja_per, ia, ja
    REAL(KIND=8),        DIMENSION(:), INTENT(IN)    :: aa
    REAL(KIND=8),        DIMENSION(:), INTENT(OUT)   :: aa_per

    INTEGER      :: n_syst, n, i, pi, l, p

    n_syst = SIZE(ia) - 1
    aa_per = 0.d0

    DO p = ia(1), ia(n_syst + 1) - 1
       aa_per(pnt(p)) = aa_per(pnt(p)) + aa(p)
    END DO

    DO n = 1, n_bord

       DO l = 1, SIZE(list(n)%DIL)
          i  = list(n)%DIL(l)
          pi = perlist(n)%DIL(l)

          DO p = ia_per(pi), ia_per(pi+1) - 1 
             IF (ja_per(p) .EQ. pi) THEN 
                aa_per(p) =  aa_per(p) + 1.d0 
             ELSE IF (ja_per(p) .EQ. i) THEN 
                aa_per(p) =  - 1.d0 
             END IF
          END DO
       END DO

    END DO

  END SUBROUTINE bc_prd_M

  SUBROUTINE bc_prd(n_bord, list, perlist, ff)
    !===================================
    USE dyn_line
    IMPLICIT NONE

    INTEGER,                           INTENT(IN)    :: n_bord
    TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)    :: list, perlist
    REAL(KIND=8),        DIMENSION(:), INTENT(INOUT) :: ff

    INTEGER      :: i, pi, l, n

    DO n = 1, n_bord 
       DO l = 1, SIZE(list(n)%DIL)
          i = list(n)%DIL(l)
          pi = perlist(n)%DIL(l)

          ff(i) = ff(i) + ff(pi)
          ff(pi) = 0.d0

       END DO
    END DO

  END SUBROUTINE bc_prd


  FUNCTION last_c_leng (len_str, string) RESULT (leng)
    !===================================================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .EQ. ' ' ) THEN
          leng = i-1; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_c_leng

END MODULE periodic
