!  =================================
!  Entry points for VECTOR equations
!  =================================

!  ff  and  u0  are Scalar functions  (fs is defined only on the boundary)
!  gg  and  v0  are Vector fields     (gs is defined only on the boundary)

!  w   denotes either a scalar of a vector weigthing function
!  ws  denotes either a scalar of a vector weigthing function
!      defined only on the boundary

!  < _ , _ >     means integration over elements in the list  m0
! << _ , _ >>    means scalar product and integration over elements in  m0
!  < _ , _ >_s   means surface integration over boundary elements in  ms0
! << _ , _ >>_s  means scalar product and surface integration over b.e. ms0

!     D  is Nabla spatial derivative operator (grad, div, curl)

MODULE  fem_v

CONTAINS

SUBROUTINE qv_00 (mesh, gg,  v0)
!=================================

!  << w, g >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   INTEGER :: m, l, k, n

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) =  SUM(gg(k,jj(:,m)) * ww(:,l)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_00

!------------------------------------------------------------------------------

SUBROUTINE qvect_00 (mesh, gg,  v0)
!=================================

!  << w, g >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   INTEGER :: m, l, k, n


   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) =  SUM(gg(k,jj(:,m)) * ww(:,l)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qvect_00

SUBROUTINE qvect_00_bloc (mesh, gg,  v0)
!=================================

!  << w, g >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   INTEGER :: m, l, k, n

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) =  SUM(gg(k,jj(:,m)) * ww(:,l)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qvect_00_bloc
!------------------------------------------------------------------------------

SUBROUTINE qv_01 (m0, jj, ff,  v0)
!=================================

!  << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, n
   REAL(KIND=8), DIMENSION(k_d) :: gl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) =  SUM(ff(jj(:,m)) * dw(k,:,l,m)) * rj(l,m)
         ENDDO

         DO n = 1, n_w
            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_01

!------------------------------------------------------------------------------

SUBROUTINE qv_10 (m0, jj, ff,  v0)
!=================================

!  < D.w, f >   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l
   REAL(KIND=8) :: fl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         fl =  SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

         v0(:, jj(:,m)) = v0(:, jj(:,m)) + dw(:,:,l,m) * fl

      ENDDO
   ENDDO

END SUBROUTINE qv_10

!------------------------------------------------------------------------------

SUBROUTINE dv_11 (m0, jj, vv,  v0)
!=================================

!  < D.w, D.v >   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: vv
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l
   REAL(KIND=8) :: dvl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         dvl =  SUM(vv(:,jj(:,m)) * dw(:,:,l,m)) * rj(l,m)

         v0(:, jj(:,m)) = v0(:, jj(:,m)) + dw(:,:,l,m) * dvl

      ENDDO
   ENDDO

END SUBROUTINE dv_11

!------------------------------------------------------------------------------

SUBROUTINE qv_11_u (m0, jj, vv,  v0)
!===================================

!  SUM from k=1 to k_d  << (Dw_k), (Dv_k) >>   ===>   v0
!  Uncoupled version of the Laplacian of a VECTOR variable
!  -------------------------------------------------------

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: vv
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l
   REAL(KIND=8) :: dvl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         dvl =  SUM(vv(:,jj(:,m)) * dw(:,:,l,m)) * rj(l,m)

         v0(:, jj(:,m)) = v0(:, jj(:,m)) + dw(:,:,l,m) * dvl

      ENDDO
   ENDDO

END SUBROUTINE qv_11_u

!------------------------------------------------------------------------------

SUBROUTINE qv_001 (m0, jj, gg, vv,  v0)
!======================================

!  << w, (g.D)v >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, vv
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, i, n
   REAL(KIND=8) :: s
   REAL(KIND=8), DIMENSION(k_d) :: gl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d
            s = 0
            DO i = 1, k_d
              s = s + SUM(gg(i,jj(:,m)) * ww(:,l))   &
                    * SUM(vv(k,jj(:,m)) * dw(i,:,l,m))
            ENDDO
            gl(k) = s * rj(l,m)
         ENDDO

         DO n = 1, n_w
            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_001

!------------------------------------------------------------------------------

SUBROUTINE qv_00_s (ms0, js, is, gs,  v0)
!========================================

!  << ws, gs >>_s   ===>   v0    incremental accumulation of the boundary term

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: js, is
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, ms, ls, k, ns
   REAL(KIND=8), DIMENSION(k_d) :: gls

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) =  SUM(gs(k,is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         ENDDO

         DO ns = 1, n_ws
            v0(:, js(ns,ms)) = v0(:, js(ns,ms)) + ww(ns,ls) * gls
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_00_s

!------------------------------------------------------------------------------

SUBROUTINE qv_10_s (ms0, js, is, fs,  v0)
!========================================

!  < n.ws, fs >>_s   ===>   v0    incremental accumulation of the boundary term

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: js, is
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: fs
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, ms, ls, ns
   REAL(KIND=8) :: fls

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
      DO ls = 1, l_Gs

         fls =  SUM(fs(is(:,ms)) * wws(:,ls)) * rjs(ls,ms)

         DO ns = 1, n_ws
            v0(:, js(ns,ms)) = v0(:, js(ns,ms))   &
                                  + rnorms(:,ls,ms) * wws(ns,ls) * fls
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_10_s

SUBROUTINE qv_mass_grad (m0, jj, gg, ff,  v0)
!=================================

!  << w, g >>   
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, n
   REAL(KIND=8), DIMENSION(k_d) :: gl
   INTEGER, DIMENSION(n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(k_d,n_w) :: dwlm, ggkn 
   REAL(KIND=8), DIMENSION(n_w) :: ffn 

   v0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)
      j_loc(:) = jj(:,m)

!      ggkn(:,:) = gg(:,jj(:,m))
!      ffn(:)   = ff(jj(:,m))
      ggkn(:,:) = gg(:,j_loc(:))
      ffn(:)   = ff(j_loc(:))

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         gl = 0.
         DO k = 1, k_d
            DO n = 1, n_w
               gl(k) =  gl(k) + ggkn(k,n) * ww(n,l) &
                              + ffn(n) * dwlm(k,n)
            END DO
         ENDDO
         gl = gl * rj(l,m)

         DO n = 1, n_w
!            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
            DO k = 1, k_d
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * gl(k)
            ENDDO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_mass_grad

SUBROUTINE qv_mass_gradT (mesh, gg, ff,  v0)
!=================================

!  << w, g >>   
!  - << D.w, f >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
 
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: ggkn 
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8) :: fl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.

   DO m = 1, me
      j_loc(:) = jj(:,m)

!      ggkn(:,:) = gg(:,jj(:,m))
!      ffn(:)   = ff(jj(:,m))
      ggkn(:,:) = gg(:,j_loc(:))
      ffn(:)   = ff(j_loc(:))

      DO l = 1, l_G

         gl = 0.
         DO k = 1, k_d
            DO n = 1, n_w
               gl(k) =  gl(k) + ggkn(k,n) * ww(n,l) 
            END DO
         ENDDO
         gl = gl * rj(l,m)

         fl =  SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)


         DO n = 1, n_w
!            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
            DO k = 1, k_d
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * gl(k) - fl*dw(k,n,l,m)
            ENDDO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_mass_gradT

SUBROUTINE qv_gradT (mesh, ff,  v0)
!=================================
!  - << D.w, f >>   ===>   v0
   USE Gauss_points
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me
   INTEGER :: m, l, k, n
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8) :: fl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.d0

   DO m = 1, me
      j_loc(:) = jj(:,m)
      ffn(:)   = ff(j_loc(:))

      DO l = 1, l_G
         fl =  SUM(ffn* ww(:,l)) * rj(l,m)
         DO n = 1, n_w
            DO k = 1, k_d
               v0(k, j_loc(n)) = v0(k, j_loc(n)) - fl*dw(k,n,l,m)
            ENDDO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_gradT

!!$SUBROUTINE qv_diff_mass_gradT (m0, jj, alpha, uu, gg, ff,  v0)
!!$!=================================
!!$
!!$!  << w, g >>   + alpha*<<Dw,Du>>  - << D.w, f >>   ===>   v0
!!$
!!$   USE Gauss_points
!!$
!!$   IMPLICIT NONE
!!$   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
!!$   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
!!$   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu, gg
!!$   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
!!$   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
!!$   REAL(KIND=8), INTENT(IN)                  :: alpha
!!$  
!!$   INTEGER :: mm, m, l, k, kp, n
!!$   REAL(KIND=8), DIMENSION(k_d) :: gl
!!$   INTEGER, DIMENSION(n_w) :: j_loc 
!!$   REAL(KIND=8), DIMENSION(k_d,n_w) :: ggkn, uukn
!!$   REAL(KIND=8), DIMENSION(k_d,k_d) :: dul
!!$   REAL(KIND=8), DIMENSION(n_w) :: ffn 
!!$   REAL(KIND=8) :: fl
!!$
!!$   v0 = 0.
!!$
!!$   DO mm = 1, SIZE(m0);  m = m0(mm)
!!$      j_loc(:) = jj(:,m)
!!$
!!$!      ggkn(:,:) = gg(:,jj(:,m))
!!$!      ffn(:)   = ff(jj(:,m))
!!$      ggkn(:,:) = gg(:,j_loc(:))
!!$      uukn(:,:) = uu(:,j_loc(:))
!!$      ffn(:)   = ff(j_loc(:))
!!$
!!$      DO l = 1, l_G
!!$
!!$         gl = 0.
!!$         DO k = 1, k_d
!!$            DO n = 1, n_w
!!$               gl(k) =  gl(k) + ggkn(k,n) * ww(n,l) 
!!$            END DO
!!$         ENDDO
!!$         gl = gl * rj(l,m)
!!$         fl = SUM(ffn(:) * ww(:,l)) * rj(l,m)
!!$
!!$         DO k  = 1, k_d
!!$         DO kp = 1, k_d
!!$            dul(k,kp) = SUM(uukn(k,:) * dw(kp,:,l,m)) * rj(l,m) * alpha
!!$         ENDDO
!!$         ENDDO
!!$
!!$         DO n = 1, n_w
!!$!            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
!!$            DO k = 1, k_d
!!$               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * gl(k) &
!!$                               - fl*dw(k,n,l,m)  + SUM(dw(:,n,l,m)*dul(k,:))
!!$            ENDDO
!!$         ENDDO
!!$
!!$      ENDDO
!!$   ENDDO
!!$
!!$END SUBROUTINE qv_diff_mass_gradT

SUBROUTINE qv_diff_mass_adv_gradT (mesh, alpha, uu, gg, ff, vv, v0)
!=================================

!  << w, g >> + alpha*<<Dw,Du>> + << w, v.Dv >>   
!  - << D.w, f >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu, gg, vv
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   REAL(KIND=8), INTENT(IN)                  :: alpha
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, kp, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl, vkpl, vgdv
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, ggkn, uukn, vvkn, v0kn
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%k_d) :: dul, dvl
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8) :: fl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.

   DO m = 1, me

      j_loc(:) = jj(:,m)

      ggkn(:,:) = gg(:,j_loc(:))
      uukn(:,:) = uu(:,j_loc(:))
      vvkn(:,:) = vv(:,j_loc(:))
      ffn(:)   = ff(j_loc(:))
      
      v0kn = 0

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO k = 1, k_d
            gl(k)   = SUM(ggkn(k,:) * ww(:,l))
            vkpl(k) = SUM(vvkn(k,:) * ww(:,l))
         ENDDO
         fl = SUM(ffn(:) * ww(:,l))

         DO k  = 1, k_d
            DO kp = 1, k_d
               dul(k,kp) = SUM(uukn(k,:) * dwlm(kp,:))
               dvl(k,kp) = SUM(vvkn(k,:) * dwlm(kp,:))
            ENDDO
            vgdv(k) = SUM(vkpl(:)*dvl(k,:))
         ENDDO 

         gl = gl + vgdv ! source + Nl

         DO n = 1, n_w
            DO k = 1, k_d
               v0kn(k,n) = v0kn(k,n) & 
                           + (ww(n,l) * gl(k) &
                           - fl*dwlm(k,n)  &
                           + alpha*SUM(dwlm(:,n)*dul(k,:)))*rj(l,m)
            ENDDO
         ENDDO

      ENDDO
      v0(:, j_loc(:)) =  v0(:, j_loc(:)) + v0kn
   ENDDO

END SUBROUTINE qv_diff_mass_adv_gradT

SUBROUTINE qv_diff_mass_gradT (mesh, alpha, uu, gg, ff, v0)
!=================================

!  << w, g >> + alpha*<<Dw,Du>> 
!  - << D.w, f >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu, gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
   REAL(KIND=8), INTENT(IN)                  :: alpha
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, kp, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl, vkpl, vgdv
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, ggkn, uukn, v0kn
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%k_d) :: dul, dvl
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8) :: fl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.

   DO m = 1, me

      j_loc(:) = jj(:,m)

      ggkn(:,:) = gg(:,j_loc(:))
      uukn(:,:) = uu(:,j_loc(:))
      ffn(:)   = ff(j_loc(:))
      
      v0kn = 0

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO k = 1, k_d
            gl(k)   = SUM(ggkn(k,:) * ww(:,l))
         ENDDO
         fl = SUM(ffn(:) * ww(:,l))

         DO k  = 1, k_d
            DO kp = 1, k_d
               dul(k,kp) = SUM(uukn(k,:) * dwlm(kp,:))
            ENDDO
         ENDDO 

         !gl = gl + vgdv ! source + Nl

         DO n = 1, n_w
            DO k = 1, k_d
               v0kn(k,n) = v0kn(k,n) & 
                           + (ww(n,l) * gl(k) &
                           - fl*dwlm(k,n)  &
                           + alpha*SUM(dwlm(:,n)*dul(k,:)))*rj(l,m)
            ENDDO
         ENDDO

      ENDDO
      v0(:, j_loc(:)) =  v0(:, j_loc(:)) + v0kn
   ENDDO

END SUBROUTINE qv_diff_mass_gradT

SUBROUTINE qv_nl_ns (mesh, vv, v0)
!=================================

!   << w, v.Dv >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: vv
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, kp, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: vkpl, vgdv
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, vvkn, v0kn
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%k_d) :: dvl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.d0

   DO m = 1, me

      j_loc(:) = jj(:,m)
      vvkn(:,:) = vv(:,j_loc(:))
      v0kn = 0

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO k = 1, k_d
            vkpl(k) = SUM(vvkn(k,:) * ww(:,l))
         ENDDO

         DO k  = 1, k_d
            DO kp = 1, k_d
               dvl(k,kp) = SUM(vvkn(k,:) * dwlm(kp,:))
            ENDDO
            vgdv(k) = SUM(vkpl(:)*dvl(k,:))
         ENDDO 

         DO n = 1, n_w
            DO k = 1, k_d
               v0kn(k,n) = v0kn(k,n) + (ww(n,l) * vgdv(k))*rj(l,m)
            ENDDO
         ENDDO

      ENDDO
      v0(:, j_loc(:)) =  v0(:, j_loc(:)) + v0kn
   ENDDO

END SUBROUTINE qv_nl_ns

SUBROUTINE qv_mass_adv_gradT (mesh, gg, ff, vv, v0)
!=================================

!  << w, g >> + << w, v.Dv >> - << D.w, f >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, vv
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, kp, n
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl, vkpl, vgdv
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, ggkn, vvkn
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%k_d) :: dvl
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8) :: fl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.

   DO m = 1, me

      j_loc(:) = jj(:,m)

      ggkn(:,:) = gg(:,j_loc(:))
      vvkn(:,:) = vv(:,j_loc(:))
      ffn(:)   = ff(j_loc(:))

      !v0kn = 0

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO k = 1, k_d
            gl(k)   = SUM(ggkn(k,:) * ww(:,l))
            vkpl(k) = SUM(vvkn(k,:) * ww(:,l))
         ENDDO
         fl = SUM(ffn(:) * ww(:,l)) * rj(l,m)

         DO k  = 1, k_d
            DO kp = 1, k_d
               dvl(k,kp) = SUM(vvkn(k,:) * dwlm(kp,:))
            ENDDO
            vgdv(k) = SUM(vkpl(:)*dvl(k,:))
         ENDDO 

         gl = (gl + vgdv) * rj(l,m) ! source + Nl

         DO n = 1, n_w
            DO k = 1, k_d
               !v0kn(k,n) = v0kn(k,n) + ww(n,l) * gl(k) - fl*dwlm(k,n)
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * gl(k) - fl*dwlm(k,n)
            ENDDO
         ENDDO

      ENDDO
      !v0(:, j_loc(:)) =  v0(:, j_loc(:)) + v0kn
   ENDDO

END SUBROUTINE qv_mass_adv_gradT

SUBROUTINE qv_mass_adv_grad (mesh, gg, ff, uu,  v0)
!=================================

!  << w, g >>   
!  << w, u.Du >>   
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n, kp
   REAL(KIND=8) :: gl, ugduk, dkpuk, ukp, resul
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, uukn, ggkn, v0loc
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: ukpl 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.

   DO m = 1, me

      j_loc(:) = jj(:,m)

      uukn(:,:) = uu(:,jj(:,m))
      ggkn(:,:) = gg(:,jj(:,m))
      ffn(:)   = ff(jj(:,m))

      v0loc = 0.d0

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO kp = 1, k_d
            ukp  = 0.d0
            DO n = 1,n_w
               ukp   = ukp   + uukn(kp,n) * ww(n,l)
            END DO
            ukpl(kp)= ukp 
         ENDDO


         DO k = 1, k_d

            gl = 0.d0
            DO n = 1, n_w
               gl =  gl + ggkn(k,n) * ww(n,l) &
                        + ffn(n) * dwlm(k,n)
            END DO

            ugduk = 0.d0
            DO kp = 1, k_d
               dkpuk = 0.d0
               DO n = 1,n_w
                  dkpuk = dkpuk + uukn(k, n) * dwlm(kp,n)
               END DO
               ugduk = ugduk + ukpl(kp)*dkpuk 
            ENDDO

            resul = (gl + ugduk) * rj(l,m)
            DO n = 1, n_w
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * resul 
            ENDDO
         END DO
      END DO
   END DO
            !DO n = 1, n_w
            !   v0loc(k,n) = v0loc(k,n) + (gl + ugduk) * rj(l,m)* ww(n,l)
            !ENDDO

         !END DO
      !END DO

      !DO n = 1, n_w
      !   DO k = 1, k_d
      !      v0(k, j_loc(n)) = v0(k, j_loc(n)) + v0loc(k,n)
      !   ENDDO
      !ENDDO

   !ENDDO

END SUBROUTINE qv_mass_adv_grad

SUBROUTINE qv_massvar_adv_grad (mesh, source, rho, gg, ff, uu, v0)
!=================================
!  << w, rho *gg + source>>
!  << w, rho*u.D uu >>   
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: source, gg, uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: rho, ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n, kp
   REAL(KIND=8) :: gl, ugduk, dkpuk, ukp, resul, rhol
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, uukn, ggkn, srckn
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn, rhon
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: ukpl, gkpl, srckpl

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.d0

   DO m = 1, me

      j_loc(:) = jj(:,m)

      uukn(:,:) = uu(:,j_loc)
      ggkn(:,:) = gg(:,j_loc)
      ffn(:)    = ff(j_loc)
      rhon(:)   = rho(j_loc)
      srckn(:,:)   = source(:,j_loc)


      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO kp = 1, k_d
            ukpl(kp) = SUM(uukn(kp,:)*ww(:,l))
            gkpl(kp) = SUM(ggkn(kp,:)*ww(:,l))
            srckpl(kp) = SUM(srckn(kp,:)*ww(:,l))
         ENDDO
         rhol = SUM(rhon*ww(:,l))
         gkpl = rhol*gkpl + srckpl

         DO k = 1, k_d

            gl = SUM(ffn*dwlm(k,:))

            ugduk = 0.d0
            DO kp = 1, k_d
               ugduk = ugduk + ukpl(kp)*SUM(uukn(k,:)*dwlm(kp,:))
            ENDDO

            resul = (gkpl(k) + gl + rhol * ugduk) * rj(l,m)

            DO n = 1, n_w
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * resul 
            ENDDO

         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_massvar_adv_grad

SUBROUTINE qv_massvar_grad (mesh, gg, ff, rho, uu, v0)
!=================================

!  <<  w, g + rho*uu >>      
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, rho
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0
  
   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   INTEGER :: m, l, k, n, kp
   REAL(KIND=8) :: gl, ugduk, dkpuk, ukp, resul, rhol
   INTEGER, DIMENSION(mesh%gauss%n_w) :: j_loc 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dwlm, uukn, ggkn 
   REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: ffn 
   REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: ukpl 

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0.d0

   DO m = 1, me

      j_loc(:) = jj(:,m)

      uukn(:,:) = uu(:,jj(:,m))
      ggkn(:,:) = gg(:,jj(:,m))
      ffn(:)   = ff(jj(:,m))

      DO l = 1, l_G

         dwlm(:,:) = dw(:,:,l,m)

         DO kp = 1, k_d
            ukp  = 0.d0
            DO n = 1,n_w
               ukp   = ukp   + uukn(kp,n) * ww(n,l)
            END DO
            ukpl(kp)= ukp 
         ENDDO

         rhol = SUM( rho(jj(:,m))*ww(:,l) )
         ukpl = rhol*ukpl
         
         DO k = 1, k_d

            gl = 0.d0
            DO n = 1, n_w
               gl =  gl + ggkn(k,n) * ww(n,l) &
                        + ffn(n) * dwlm(k,n)
            END DO

            resul = (gl + ukpl(k)) * rj(l,m)

            DO n = 1, n_w
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * resul 
            ENDDO

         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_massvar_grad


SUBROUTINE qv_mass_grad_bis (m0, jj, gg, ff,  v0)
!=================================

!  << w, g >>   
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, n
   REAL(KIND=8), DIMENSION(k_d) :: gl

   v0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         gl = 0.
         DO k = 1, k_d
            DO n = 1, n_w
               gl(k) =  gl(k) + gg(k,jj(n,m)) * ww(n,l) &
                              + ff(jj(n,m)) * dw(k,n,l,m)
            END DO
         ENDDO
         gl = gl * rj(l,m)

         DO n = 1, n_w
!            v0(:, jj(n,m)) = v0(:, jj(n,m)) + ww(n,l) * gl
            DO k = 1, k_d
               v0(k, jj(n,m)) = v0(k, jj(n,m)) + ww(n,l) * gl(k)
            ENDDO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_mass_grad_bis

SUBROUTINE qvs_mass_grad (m0, jj, gg, ff,  v0)
!=================================

!  << w, g >>   
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, n
   REAL(KIND=8) :: gl

   v0 = 0.

   DO k = 1, k_d
   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         gl = 0.
         DO n = 1, n_w
            gl =  gl + gg(k,jj(n,m)) * ww(n,l) &
                     + ff(jj(n,m)) * dw(k,n,l,m)
         END DO
         gl = gl * rj(l,m)

         DO n = 1, n_w
               v0(k, jj(n,m)) = v0(k, jj(n,m)) + ww(n,l) * gl
         ENDDO

      ENDDO
   ENDDO
   ENDDO

END SUBROUTINE qvs_mass_grad

!=============================================================
SUBROUTINE cs_01 (mesh, gg,  u0)

!     < w, D x g >   ===>   u0
!     WORKS IN 2D ONLY

   USE Gauss_points

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: u0

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me
   REAL(KIND=8):: x, s
   INTEGER :: m, l, n

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   u0 = 0

   DO m = 1, me
      DO l = 1, l_G

         s = 0
         DO n = 1, n_w
            s = s + gg(2, jj(n,m)) * dw(1,n,l,m) &
                  - gg(1, jj(n,m)) * dw(2,n,l,m)
         ENDDO
	 s = s * rj(l,m)

         DO n = 1, n_w
              x = ww(n,l) * s
              u0(jj(n,m)) = u0(jj(n,m)) + x
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE cs_01

SUBROUTINE cv_01 (m0, jj, gg,  U0)
!=================================

!  << w, D x g >>   ===>   U0  ! Bloc vector 

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: U0

   INTEGER :: mm, m, n, l, i, i_b, k, k1, k2, bloc_size
   REAL(KIND=8) :: rgk

   bloc_size = SIZE(gg,2)

   U0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1     

            rgk = SUM(gg(k2, jj(:,m)) * dw(k1,:,l,m)   &
                    - gg(k1, jj(:,m)) * dw(k2,:,l,m)) * rj(l,m)

!!!!!            v0(k, jj(:,m)) = v0(k, jj(:,m)) + ww(:,l) * rgk

            DO n = 1, n_w;  i_b = jj(n, m)
                 i = i_b + (k-1)*bloc_size
                 U0(i) = U0(i) + ww(n,l) * rgk
            END DO

         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE cv_01

SUBROUTINE qv_01_s (ms0, js, is, gs,  v0)
!========================================

!  << ws x ns, gs >>_s   ===>   v0   incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: js, is
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs 
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, ms, ns, ls, k1, k2, k
   REAL(KIND=8), DIMENSION(3) :: gls, ngls


   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k,is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         END DO
         
         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1 
            ngls(k) = gls(k2)*rnorms(k1,ls,ms) - gls(k1)*rnorms(k2,ls,ms)
         END DO

         DO k = 1, k_d 
            DO ns = 1, n_ws
                v0(k,js(ns,ms)) = v0(k,js(ns,ms)) + wws(ns,ls) * ngls(k)
            END DO
         END DO

      ENDDO
   ENDDO

END SUBROUTINE qv_01_s

SUBROUTINE cv_00_s (ms0, js, is, gs,  U0)
!========================================

!  << ws x ns, gs >>_s   ===>   U0   incremental accumulation of boundary terms
!   U0  is a block vector

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: js, is
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs 
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: U0

   INTEGER :: mm, ms, ls, bloc_size, k1, k2, k
!   INTEGER, DIMENSION(n_ws) :: j
   INTEGER  :: j, ns
   REAL(KIND=8), DIMENSION(3) :: gls, ngls

   bloc_size = SIZE(U0)/3

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k,is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         END DO
         
         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1 
            ngls(k) = gls(k2)*rnorms(k1,ls,ms) - gls(k1)*rnorms(k2,ls,ms)
         END DO

         DO k = 1, k_d 
            DO ns = 1, n_ws
                j = js(ns,ms) + (k-1)*bloc_size
                U0(j) = U0(j) + wws(ns,ls) * ngls(k)
            END DO
         END DO

!         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1 
! 
!            j = js(:,ms) + (k-1)*bloc_size
! 
!            U0(j) = U0(j) + wws(:,ls) *   &
!                   (gls(k2)*rnorms(k1,ls,ms) - gls(k1)*rnorms(k2,ls,ms))
!         END DO

      ENDDO
   ENDDO

END SUBROUTINE cv_00_s

SUBROUTINE cc_00_s (ms0, js, is, gs,  U0)
!========================================

!  << ws x ns, gs x ns >>_s   ===>   U0   incremental accumulation of boundary terms
!   U0  is a block vector

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: js, is
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: U0

   INTEGER :: mm, ms, ls, bloc_size, k, k1, k2
   INTEGER, DIMENSION(n_ws) :: j
   REAL(KIND=8), DIMENSION(3) :: gls, hls

   bloc_size = SIZE(U0)/3

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k,is(:,ms)) * wws(:,ls)) * rjs(ls,ms)
         END DO

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            hls(k) = gls(k1)*rnorms(k2,ls,ms) - gls(k2)*rnorms(k1,ls,ms)
         END DO

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1

            j = js(:,ms) + (k-1)*bloc_size

            U0(j) = U0(j) + wws(:,ls) *   &
                   (hls(k2)*rnorms(k1,ls,ms) - hls(k1)*rnorms(k2,ls,ms))

         END DO

      ENDDO
   ENDDO

END SUBROUTINE cc_00_s

SUBROUTINE cv_10 (m0, jj, gg,  U0)                                              
!=================================                                              
                                                                                
!  << D x w,  g >>   ===>   U0                                                  
                                                                                
   USE Gauss_points
                                                                                
   IMPLICIT NONE                                                                
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0                              
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj                              
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg                              
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: U0                              
                                                                                
   INTEGER :: mm, m, l, k, k1, k2, bloc_size                                               
   REAL(KIND=8), DIMENSION(k_d) :: gl                                           
   REAL(KIND=8), DIMENSION(n_w) :: xx                                           
                                                                                
   bloc_size = SIZE(gg,2)
   U0 = 0                                                                       
                                                                                
   DO mm = 1, SIZE(m0);  m = m0(mm)                                             
      DO l = 1, l_G

         DO k = 1, k_d                                                          
            gl(k) = SUM(gg(k, jj(:,m)) * ww(:,l)) * rj(l,m)                     
         ENDDO                                                                  
                                                                                
         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1      
            xx = dw(k2,:,l,m)*gl(k1) - dw(k1,:,l,m) * gl(k2)                          
            U0((k-1)*bloc_size + jj(:,m)) = U0((k-1)*bloc_size + jj(:,m)) + xx
         ENDDO                                                                  
                                                                                
      ENDDO                                                                     
   ENDDO                                                                        
                                                                                
END SUBROUTINE cv_10                                                            

SUBROUTINE cc_10 (m0, jj, gg,  v0)                                              
!=================================                                              
                                                                                
!  << D x w,  g >>   ===>   v0                                                  
                                                                                
   USE Gauss_points
                                                                                
   IMPLICIT NONE                                                                
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0                              
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj                              
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg                              
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0                              
                                                                                
   INTEGER :: mm, m, l, k, k1, k2
   REAL(KIND=8), DIMENSION(k_d) :: gl                                           
   REAL(KIND=8), DIMENSION(n_w) :: xx                                           
                                                                                
   v0 = 0                                                                       
                                                                                
   DO mm = 1, SIZE(m0);  m = m0(mm)                                             
      DO l = 1, l_G

         DO k = 1, k_d                                                          
            gl(k) = SUM(gg(k, jj(:,m)) * ww(:,l)) * rj(l,m)                     
         ENDDO                                                                  
                                                                                
         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1      
            xx = dw(k2,:,l,m)*gl(k1) - dw(k1,:,l,m) * gl(k2)                          
            v0(k,jj(:,m)) = v0(k,jj(:,m)) + xx
         ENDDO                                                                  
                                                                                
      ENDDO                                                                     
   ENDDO                                                                        
                                                                                
END SUBROUTINE cc_10                                                            

SUBROUTINE cg_11 (m0, jj, ff,  u0)                                              
!=================================                                              
                                                                                
!  << D x w,  D f >>   ===>   u0                                                  
                                                                                
   USE Gauss_points
                                                                                
   IMPLICIT NONE                                                                
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0                              
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj                              
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff                          
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: u0                              
                                                                                
   INTEGER :: mm, l, m, n, k, k1, k2
                                             
   REAL(KIND=8), DIMENSION(k_d) :: gfl                                           
   REAL(KIND=8), DIMENSION(n_w) :: xx                                           
                                                                                
   u0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)                                             
      DO l = 1, l_G

         gfl = 0.
         DO k = 1, k_d                                                          
!            gfl(k) = SUM(ff(jj(:,m)) * dw(k,:,l,m)) * rj(l,m)
            DO n = 1, n_w
               gfl(k) = gfl(k) + ff(jj(n,m)) * dw(k,n,l,m)
            END DO
         ENDDO                                                                  
         gfl = gfl * rj(l,m)

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1      
            xx = dw(k2,:,l,m)*gfl(k1) - dw(k1,:,l,m) * gfl(k2)                          
            u0(k,jj(:,m)) = u0(k,jj(:,m)) + xx
         ENDDO                                                                  
                                                                                
      ENDDO                                                                     
   ENDDO                                                                        
                                                                                
END SUBROUTINE cg_11                                                            

SUBROUTINE gc_11 (m0, jj, gg,  u0)                                              
!=================================                                              
                                                                                
!  < D w,  D x g >   ===>   u0                                                  
                                                                                
   USE Gauss_points
                                                                                
   IMPLICIT NONE                                                                
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0                              
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj                              
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg                          
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0                              
                                                                                
   INTEGER :: mm, l, m, n, k, k1, k2
   REAL(KIND=8), DIMENSION(k_d) :: cgl
   REAL(KIND=8) :: xx                                           
                                                                                
   u0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)                                             
      DO l = 1, l_G

         cgl = 0.
         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1      
            DO n = 1, n_w
               cgl(k) = cgl(k) + dw(k1,n,l,m) * gg(k2,jj(n,m)) &
                               - dw(k2,n,l,m) * gg(k1,jj(n,m)) 
             END DO
         ENDDO                                                                  
         cgl = cgl * rj(l,m)

         DO n = 1, n_w 
            xx = 0.
            DO k = 1, k_d
               xx = xx + dw(k,n,l,m)*cgl(k)                          
            END DO
            u0(jj(n,m)) = u0(jj(n,m)) + xx
         ENDDO                                                                  
                                                                                
      ENDDO                                                                     
   ENDDO                                                                        
                                                                                
END SUBROUTINE gc_11                                                            

SUBROUTINE qv_11 (m0, jj, gg, v0)
!=======================================

!  << (D x w), (D x gg) >>  ===>   v0  

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg 
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: v0

   REAL(KIND=8), DIMENSION(3) :: rgl 
   INTEGER :: mm, m, l, n, k, k1, k2

   v0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            rgl(k) = SUM(gg(k2, jj(:,m)) * dw(k1,:,l,m)   &
                      -  gg(k1, jj(:,m)) * dw(k2,:,l,m)) * rj(l,m)
         ENDDO
          
         DO n = 1, n_w 
            DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
               v0(k,jj(n,m)) = v0(k,jj(n,m)) & 
                        + dw(k2,n,l,m)*rgl(k1) - dw(k1,n,l,m)*rgl(k2)
            END DO
         END DO 

      ENDDO
   ENDDO

END SUBROUTINE qv_11

SUBROUTINE cv_11 (m0, jj, gg, U0)
!=======================================

!  << (D x w), (D x gg) >>  ===>   U0    ! Bloc vector

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg 
   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: U0

   REAL(KIND=8), DIMENSION(3) :: rgl 
   INTEGER :: mm, m, l, n, k, k1, k2, i, i_b, bloc_size

   bloc_size = SIZE(gg,2)

   U0 = 0.

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            rgl(k) = SUM(gg(k2, jj(:,m)) * dw(k1,:,l,m)   &
                      -  gg(k1, jj(:,m)) * dw(k2,:,l,m)) * rj(l,m)
         ENDDO
          
         DO n = 1, n_w;  i_b = jj(n, m)
            DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
               i = i_b + (k-1)*bloc_size
               U0(i) = U0(i) + dw(k2,n,l,m)*rgl(k1) - dw(k1,n,l,m)*rgl(k2)
            END DO
         END DO 

      ENDDO
   ENDDO

END SUBROUTINE cv_11

SUBROUTINE cc_101 (m0, jj, gg, vv, U0)
!=======================================

!  << (D x w), gg x (D x vv) >>  ===>   U0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, vv
   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: U0

   REAL(KIND=8), DIMENSION(3) :: gl, rvl, grvl
   INTEGER :: mm, m, l, n, k, k1, k2, i, i_b, bloc_size

   bloc_size = SIZE(gg,2)

   DO mm = 1, SIZE(m0);  m = m0(mm)
      DO l = 1, l_G

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1

             gl(k) = SUM(gg(k, jj(:,m)) * ww(:,l)) * rj(l,m)                     

             rvl(k) = SUM(vv(k2, jj(:,m)) * dw(k1,:,l,m)   &
                        - vv(k1, jj(:,m)) * dw(k2,:,l,m))

         ENDDO

         DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
            grvl(k) = gl(k1)*rvl(k2) - gl(k2)*rvl(k1)
         ENDDO

         DO n = 1, n_w;  i_b = jj(n, m)
            DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
               i = i_b + (k-1)*bloc_size
               U0(i) = U0(i) + dw(k2,n,l,m)*grvl(k1) - dw(k1,n,l,m)*grvl(k2)
            END DO
         END DO

      ENDDO
   ENDDO

END SUBROUTINE cc_101

SUBROUTINE proj_v_p0 (m0, jj, gg,  v_e)
!=================================

!  << w, g >>   ===>   v_e

   USE Gauss_points

   IMPLICIT NONE
   INTEGER, DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v_e


   INTEGER :: mm, m, l, k, n
   REAL(KIND=8), DIMENSION(k_d) :: gl, int
   REAL(KIND=8) :: vol

   DO mm = 1, SIZE(m0);  m = m0(mm)
      vol = 0.
      int = 0.
      DO l = 1, l_G

         vol = vol + rj(l,m)

         gl = 0.
         DO k = 1, k_d
            DO n = 1, n_w
               gl(k) = gl(k) + gg(k,jj(n,m)) * ww(n,l)
            END DO
            int(k) = int(k) + gl(k)*rj(l,m) 
         ENDDO

      ENDDO

      v_e(:,m) = int / vol

   ENDDO

END SUBROUTINE proj_v_p0

SUBROUTINE qv_massvar_adv_grad_BDF (m0, jj, sig, sigo, sigoo, &
     pp, qpo, qpoo, uu, gg,  v0)
!====================================================

!  << w, g >>
!  << w, u.Du >>
!  + << w, Df >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: sig, sigo, sigoo, pp, qpo, qpoo
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu, gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   INTEGER :: mm, m, l, k, n, kp
   REAL(KIND=8) :: qol, qool, gl, rap, rapo, rapoo, ugduk, dkpuk, ukp, resul
   INTEGER,      DIMENSION(n_w) :: j_loc
   REAL(KIND=8), DIMENSION(k_d,n_w) :: dwlm, uukn, ggkn
   REAL(KIND=8), DIMENSION(n_w) :: ppn, qpon, qpoon, sign, sigon, sigoon

   REAL(KIND=8), DIMENSION(k_d) :: ukpl

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      j_loc(:) = jj(:,m)

      DO k = 1, k_d
         uukn(k,:) = sig(jj(:,m))*uu(k,jj(:,m))
      END DO
      ggkn(:,:) = gg(:,jj(:,m))
      ppn(:)    = pp(jj(:,m))
      qpon(:)   = qpo(jj(:,m))
      qpoon(:)  = qpoo(jj(:,m))
      sign(:) = sig(jj(:,m))
      sigon(:) = sigo(jj(:,m))
      sigoon(:) = sigoo(jj(:,m))

      DO l = 1, l_G


         dwlm(:,:) = dw(:,:,l,m)

         DO kp = 1, k_d
            ukp  = 0.
            DO n = 1,n_w
               ukp   = ukp   + uukn(kp,n) * ww(n,l)
            ENDDO
            ukpl(kp)= ukp
         ENDDO


         DO k = 1, k_d

            gl = 0
            qol = 0
            qool = 0
            DO n = 1, n_w
               gl   =  gl + ggkn(k,n) * ww(n,l) &
                        - ppn(n) * dwlm(k,n)
               qol  =  qol + qpon(n) * dwlm(k,n)
               qool =  qool + qpoon(n) * dwlm(k,n)
            ENDDO

            ugduk = 0.
            DO kp = 1, k_d
               dkpuk = 0.
               DO n = 1,n_w
                  dkpuk = dkpuk + uukn(k, n) * dwlm(kp,n)
               ENDDO
               ugduk = ugduk + ukpl(kp)*dkpuk

            ENDDO

            rap = 0
            rapo = 0
            rapoo = 0
            DO n = 1, n_w
               rap = rap + sign(n)*ww(n,l)
               rapo = rapo + sigon(n)*ww(n,l)
               rapoo = rapoo + sigoon(n)*ww(n,l)
            ENDDO

            resul = (gl - (4.d0/3)*(rap/rapo)*qol + (1.d0/3)*(rap/rapoo)*qool - ugduk) * rj(l,m)

            DO n = 1, n_w
               v0(k, j_loc(n)) = v0(k, j_loc(n)) + ww(n,l) * resul
            ENDDO

         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE qv_massvar_adv_grad_BDF



SUBROUTINE q_bloc_00 (mesh, gg, nb_bloc, bloc_size, v0)
!=================================

!  << w, g >>   ===>   v0   ! Bloc structure

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: gg
   INTEGER,                    INTENT(IN)  :: nb_bloc, bloc_size
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: v0

   INTEGER :: m, l, k, n
   REAL(KIND=8) :: gl

   TYPE(mesh_type), TARGET                     :: mesh
   INTEGER,      DIMENSION(:,:), POINTER       :: jj
   INTEGER,                      POINTER       :: me

   CALL gauss(mesh)
   jj => mesh%jj
   me => mesh%me

   v0 = 0

   DO m = 1, me
      DO l = 1, l_G

         DO k = 1, nb_bloc
            gl = 0
            DO n = 1, n_w
               gl = gl + gg((k-1)*bloc_size+jj(n,m)) * ww(n,l)
            ENDDO
            gl = gl * rj(l,m)

            DO n = 1, n_w
               v0((k-1)*bloc_size+jj(n,m)) = v0((k-1)*bloc_size+jj(n,m)) &
                                          + ww(n,l) * gl
            ENDDO
         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE q_bloc_00

END MODULE  fem_v
