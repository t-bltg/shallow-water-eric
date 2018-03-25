!
!Authors Jean-Luc Guermond, Rapahel Laguerre, Copyrights 2005
!
MODULE tn_parallele

CONTAINS

  SUBROUTINE integration_mode(m_max_c,norm_loc,norm_tot)
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER,      INTENT(IN)  :: m_max_c
    REAL(KIND=8), INTENT(IN)  :: norm_loc
    REAL(KIND=8), INTENT(OUT) :: norm_tot
    INTEGER                   :: code, nb_procs
    REAL(KIND=8) :: pi
    pi = ACOS(-1.d0)
    CALL MPI_ALLREDUCE(norm_loc,norm_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, code)
! CN-AR Tue Jan 13 2009
    norm_tot = norm_tot*(2*pi)
  END SUBROUTINE integration_mode

  FUNCTION norme_max_champ_par(v) RESULT(norm)
    IMPLICIT NONE
    include 'mpif.h'
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm, norm_tot, pi
    INTEGER      :: code
    pi = ACOS(-1.d0)
    norm_loc = MAXVAL(ABS(v))
    CALL MPI_ALLREDUCE(norm_loc,norm_tot,1,MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, code)
    norm = norm_tot*(2*pi)
  END FUNCTION norme_max_champ_par

  FUNCTION norme_L2_champ_par(mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_NS_MHD
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm

    norm_loc = norme_L2_champ(mesh, list_mode, v)
    CALL integration_mode(SIZE(list_mode),norm_loc**2, norm)
    norm = SQRT(norm)

  END FUNCTION norme_L2_champ_par

  FUNCTION norme_H1_champ_par(mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_NS_MHD

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm

    norm_loc = norme_H1_champ(mesh, list_mode, v)
    CALL integration_mode(SIZE(list_mode),norm_loc**2, norm)
    norm = SQRT(norm)

  END FUNCTION norme_H1_champ_par

  FUNCTION norme_div_par(H_mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_NS_MHD

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm

    norm_loc = norme_div(H_mesh, list_mode, v)

    CALL integration_mode(SIZE(list_mode),norm_loc**2, norm)
    norm = SQRT(norm)

  END FUNCTION norme_div_par

  FUNCTION norme_curl_par(H_mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_NS_MHD

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8) :: norm_loc, norm

    norm_loc = norme_curl(H_mesh, list_mode, v)
    CALL integration_mode(SIZE(list_mode),norm_loc**2, norm)
    norm = SQRT(norm)

  END FUNCTION norme_curl_par

  SUBROUTINE moments(H_mesh, list_mode, v, dipole_out, quadripole_out)

    USE def_type_mesh

    IMPLICIT NONE
    include 'mpif.h'
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8), DIMENSION(3),      INTENT(OUT):: dipole_out
    REAL(KIND=8), DIMENSION(3,3),    INTENT(OUT):: quadripole_out
    REAL(KIND=8), DIMENSION(3)                  :: dipole
    REAL(KIND=8), DIMENSION(3,3)                :: quadripole
    REAL(KIND=8)         :: jr, ray, zed, pi
    INTEGER              :: m_max_c, mode, k, m, l, ni, i
    INTEGER              :: code
    REAL(KIND=8), DIMENSION(6) :: c

    pi = ACOS(-1.d0)
    m_max_c = SIZE(list_mode)
    dipole = 0
    quadripole = 0

    IF (SIZE(v,1)/=H_mesh%np .OR. SIZE(v,2)/=6 .OR. SIZE(v,3)/=m_max_c ) THEN
       WRITE(*,*) ' BUG in MOMENTS', SIZE(v,1), H_mesh%np
       STOP
    END IF
    DO m = 1, H_mesh%me     
       DO l = 1, H_mesh%gauss%l_G

          !--------On calcule le rayon et z du point gauss
          ray = 0
          zed = 0
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             zed = zed + H_mesh%rr(2,i)*H_mesh%gauss%ww(ni,l)
          END DO
          jr = ray * H_mesh%gauss%rj(l,m)

          DO k=1, m_max_c
             mode = list_mode(k)
             IF (mode /=0 .AND. mode /=1 .AND. mode /=2) CYCLE

             !--------Compute Curl------
             c = 0
             DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)
                !--------Composante r------
                c(1) = c(1) + ( mode/ray*v(i,6,k)*H_mesh%gauss%ww(ni,l) &
                     - v(i,3,k)*H_mesh%gauss%dw(2,ni,l,m)) 
                c(2) = c(2) + (-mode/ray*v(i,5,k)*H_mesh%gauss%ww(ni,l) &
                     - v(i,4,k)*H_mesh%gauss%dw(2,ni,l,m)) 
                !--------Composante theta------
                c(3) = c(3) + (v(i,1,k)*H_mesh%gauss%dw(2,ni,l,m) &
                     - v(i,5,k)*H_mesh%gauss%dw(1,ni,l,m)) 
                c(4) = c(4) + (v(i,2,k)*H_mesh%gauss%dw(2,ni,l,m) &
                     - v(i,6,k)*H_mesh%gauss%dw(1,ni,l,m))
                !--------Composante z------
                c(5) = c(5) + (v(i,3,k)*H_mesh%gauss%dw(1,ni,l,m) &
                     + v(i,3,k)*H_mesh%gauss%ww(ni,l)/ray &
                     - mode/ray*v(i,2,k)*H_mesh%gauss%ww(ni,l))
                c(6) = c(6) + (v(i,4,k)*H_mesh%gauss%dw(1,ni,l,m) &
                     + v(i,4,k)*H_mesh%gauss%ww(ni,l)/ray &
                     + mode/ray*v(i,1,k)*H_mesh%gauss%ww(ni,l))
             ENDDO

             !--------Compute dipole and quadripole------
             IF (mode == 0) THEN
                dipole(3) = dipole(3) + 2*pi*ray*c(3)*jr
                quadripole(1,1) = quadripole(1,1) + pi*(-zed*ray*c(3))*jr
                quadripole(1,2) = quadripole(1,2) + pi*(-zed*ray*c(1)+ray*ray*c(5))*jr
                quadripole(2,1) = quadripole(2,1) + pi*( zed*ray*c(1)-ray*ray*c(5))*jr
                quadripole(2,2) = quadripole(2,2) + pi*(-zed*ray*c(3))*jr
                quadripole(3,3) = quadripole(3,3)+2*pi*( zed*ray*c(3))*jr
             ELSE IF (mode == 1) THEN
                dipole(1) = dipole(1) + pi*(-zed*c(3) -zed*c(2) +ray*c(6))*jr
                dipole(2) = dipole(2) + pi*(-zed*c(4) +zed*c(1) -ray*c(5))*jr
                quadripole(1,3) = quadripole(1,3) + pi*(-zed*c(3)-zed*c(2)+ray*c(6))*zed*jr
                quadripole(2,3) = quadripole(2,3) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*zed*jr
                quadripole(3,1) = quadripole(3,1) + pi*(ray*ray*c(3))*jr
                quadripole(3,2) = quadripole(3,2) + pi*(ray*ray*c(4))*jr
             ELSE IF (mode == 2) THEN
                quadripole(1,1) = quadripole(1,1) + pi*(-zed*c(3)-zed*c(2)+ray*c(6))*ray*jr/2
                quadripole(1,2) = quadripole(1,2) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*ray*jr/2
                quadripole(2,1) = quadripole(2,1) + pi*(-zed*c(4)+zed*c(1)-ray*c(5))*ray*jr/2
                quadripole(2,2) = quadripole(2,2) + pi*( zed*c(3)+zed*c(2)-ray*c(6))*ray*jr/2
             END IF

          END DO
       END DO
    END DO

    !--------Collect from everybody------
    CALL MPI_ALLREDUCE(dipole,    dipole_out,    3,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, code)
    CALL MPI_ALLREDUCE(quadripole,quadripole_out,9,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, code)

    RETURN
  END SUBROUTINE moments

  END MODULE tn_parallele
