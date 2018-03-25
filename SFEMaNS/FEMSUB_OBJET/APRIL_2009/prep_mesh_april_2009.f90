!
!Authors Jean-Luc Guermond, Luigi Quarapelle, Copyrights 1994, 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE prep_maill

  IMPLICIT NONE

  PUBLIC :: load_mesh, load_mesh_formatted, load_mesh_free_format, &
       load_dg_mesh_free_format, load_mesh_free_format_ordered
  PRIVATE

CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE load_dg_mesh_free_format(dir, fil, list_dom, list_inter, type_fe, mesh, mesh_formatted)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d_p1
    USE mod_gauss_points_2d_p2
    USE Dir_nodes
    !TEST june 11 2008
    !USE sub_plot
    !TEST june 11 2008


    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom, list_inter
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect, stat
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_ms
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    !TEST june 11 2008
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: const
    !INTEGER, DIMENSION(2) :: nface, ngauss
    !INTEGER :: n1, n2, face
    !TEST june 11 2008
    LOGICAL :: test, t1, t2
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, mop, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)
    IF (mesh_formatted) THEN
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    ELSE
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
          DO m = 1, me
             READ(30,*)
          END DO
          DO ms = 1, mes
             READ(30,*)
          END DO
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30) np, nw, me, nws, mes
          READ(30)
          READ(30)
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ(30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))
    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
       READ(30) jjs_lect, neighs_lect, sides_lect
       READ(30) rr_lect
    END IF


    !---Renumerotation------------------------------------------------------------
    ! Identify the status of faces
    ! stat = 1 (interface to be forgotten), stat = 2 (boundary), stat = 3 (real interface)
    ALLOCATE (stat(mes))
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       IF (neighs1==0) THEN
          WRITE(*,*) ' BUG in prep_mesh, neighs1=0 '
          STOP
       END IF
       IF (MINVAL(ABS(i_d_lect(neighs1) - list_dom))==0) THEN
          t1 = .TRUE.
       ELSE
          t1 = .FALSE.
       END IF
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
       END DO
       neighs2 = neigh_lect(n,neighs1)
       IF (neighs2==0) THEN
          IF (t1) THEN
             stat(ms) = 2  ! face on the boundary of the domain of interest 
          ELSE
             stat(ms) = 1  ! face does not touch the domain of interest 
          END IF
          CYCLE
       END IF
       ! neighs2 /=0 
       IF (MINVAL(ABS(i_d_lect(neighs2) - list_dom))==0) THEN 
          t2 = .TRUE.
       ELSE
          t2 = .FALSE.
       END IF

       IF (t1) THEN
          IF (t2) THEN
             IF (SIZE(list_inter)==0) THEN
                stat(ms) = 1 ! no inteface to treat
             ELSE IF (MINVAL(ABS(sides_lect(ms)-list_inter))==0) THEN
                stat(ms) = 3 ! real interface
             ELSE
                stat(ms) = 1 ! interface to be forgotten 
             END IF
          ELSE
             stat(ms) = 2 ! face at the boundary of the domain of interest
          END IF
       ELSE
          IF (t2) THEN
             stat(ms) = 2 ! face at the boundary of the domain of interest
          ELSE
             stat(ms) = 1 ! on an interface of no interest
          END IF
       END IF

    END DO

    ALLOCATE (nouv_nd(np),  nouv_els(mes), virgin_nd(np), virgin_ms(mes), nouv_el(me))
    nouv_nd = -1000
    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    mnouv = 0
    msnouv= 0
    nnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO m = 1, me ! Count new nodes from domain: i_d=dom
          IF (list_dom(dom) /= i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          mnouv = mnouv + 1  ! Nouvel element
          nouv_el(m) = mnouv
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
             END IF
          END DO
       END DO

       DO ms = 1, mes
          IF (stat(ms) /= 1 .AND. stat(ms) /=2 .AND. stat(ms) /=3) THEN
             WRITE(*,*) ' BUG in prep_mesh, stat out of bounds '
             STOP
          END IF

          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          IF (virgin_ms(ms)) THEN !New interface
             virgin_ms(ms) = .FALSE.
             msnouv = msnouv + 1
          END IF
          IF (stat(ms) ==3) THEN 
             ! Nodes and sides on the interface are virgin again
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
             virgin_ms(ms) = .TRUE.
          END IF
       END DO
    END DO
    mesh%me  = mnouv
    mesh%np  = nnouv
    mesh%mes = msnouv


    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))
    ALLOCATE(mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes),  mesh%sides(mesh%mes))
    ALLOCATE(mesh%rr(kd,mesh%np))

    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    nnouv = 0
    msnouv = 0
    DO dom = 1, SIZE(list_dom)
       !Loop on me and get nouv_el and nouv_nd right
       DO m = 1, me
          IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
                nouv_nd(i) = nnouv
             END IF
          END DO
       END DO

       !Loop again on me and update
       DO m = 1, me
          IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          DO n = 1, nw; i = jj_lect(n,m)
             IF (n .LE. nwneigh) THEN
                mop = neigh_lect(n,m) 
                IF (mop .LE. 0) THEN
                   mesh%neigh(n,nouv_el(m)) = 0
                ELSE IF (MINVAL(ABS(list_dom - i_d_lect(mop))) == 0) THEN
                   mesh%neigh(n,nouv_el(m)) = nouv_el(mop)
                ELSE
                   mesh%neigh(n,nouv_el(m)) = 0
                END IF
             END IF
             mesh%rr(:,nouv_nd(i)) = rr_lect(:,i)
          END DO
          mesh%i_d(nouv_el(m)) = i_d_lect(m)
          mesh%jj(:,nouv_el(m)) = nouv_nd(jj_lect(:,m))
       END DO

       !Loop on mes and get neighs_lect and nouv_els right
       DO ms = 1, mes
          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          IF (virgin_ms(ms)) THEN !New interface
             virgin_ms(ms) = .FALSE.
             msnouv = msnouv + 1
             nouv_els(ms) = msnouv
          END IF
          IF (stat(ms) ==3) THEN 
             ! Nodes and sides on the interface are virgin again
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
             virgin_ms(ms) = .TRUE.
          END IF

          ! Swapping problem
          neighs1 = neighs_lect(ms)
          IF ((ABS(i_d_lect(neighs1) - list_dom(dom)))==0) THEN
             t1 = .TRUE.
          ELSE
             t1 = .FALSE.
          END IF
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (neighs2==0) THEN
             CYCLE
          END IF
          ! neighs2 /=0 
          IF ((ABS(i_d_lect(neighs2) - list_dom(dom)))==0) THEN 
             t2 = .TRUE.
          ELSE
             t2 = .FALSE.
          END IF
          IF (.NOT.t1 .AND. t2) THEN
             neighs_lect(ms) = neighs2 !get things right (swap neighs) 
          END IF
       END DO

       !Loop again on mes and update
       DO ms = 1, mes

          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          mesh%jjs(:,nouv_els(ms))  = nouv_nd(jjs_lect(:,ms))
          mesh%neighs(nouv_els(ms)) = nouv_el(neighs_lect(ms))
          mesh%sides(nouv_els(ms))  = sides_lect(ms)
          IF (stat(ms)==3) THEN ! side is an interface to be kept
             neighs1 = neighs_lect(ms)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs1)) = 0
             neighs2 = neigh_lect(n,neighs1)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs2))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs2)) = 0
          END IF
       END DO

    END DO
    !---Fin renumerotation--------------------------------------------------------
    !TEST June 11, 2008
    !write(*,*) 'ok'
    !allocate(const(mesh%me))
    !const = 3.d0
    !do ms = 1, mesh%mes
    !   const(mesh%neighs(ms)) = mesh%i_d(mesh%neighs(ms))
    !end do
    !write(*,*) 'ok'
    !write(*,*) size(mesh%jj,2), size(const)
    !call plot_const_p1_label(mesh%jj, mesh%rr, const, 't.plt')

    !DO ms = 1, mesh%mes
    !   m = mesh%neighs(ms)
    !   !Il faut savoir quelle face est la bonne
    !   DO n = 1, 3
    !      IF (MINVAL(ABS(mesh%jjs(:,ms)-mesh%jj(n,m)))==0) CYCLE
    !      face = n 
    !   END DO
    !   n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
    !   IF (mesh%jjs(1,ms)==mesh%jj(n1,m) .AND. mesh%jjs(2,ms)==mesh%jj(n2,m)) THEN
    !      nface(1)=n1; nface(2)=n2
    !      ngauss(1)=1; ngauss(2)=2
    !   ELSE IF (mesh%jjs(1,ms)==mesh%jj(n2,m) .AND. mesh%jjs(2,ms)==mesh%jj(n1,m)) THEN
    !      nface(1)=n2; nface(2)=n1
    !      ngauss(1)=2; ngauss(2)=1
    !   ELSE
    !      WRITE(*,*) mesh%np
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(n1,m))
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(n2,m))  
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(face,m))
    !      WRITE(*,*)  mesh%rr(:,mesh%jjs(1,ms))
    !      WRITE(*,*)  mesh%rr(:,mesh%jjs(2,ms))
    !      STOP
    !   END IF
    !END DO
    !stop
    !write(100,*) mesh%neigh
    !write(100,*) mesh%jj
    !write(100,*) mesh%i_d
    !write(100,*) mesh%jjs
    !write(100,*) mesh%neighs
    !write(100,*) mesh%mes
    !write(100,*) mesh%sides
    !stop
    !TEST June 11, 2008


    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       CALL gauss_points_2d_p1(mesh)
    ELSE IF (nw==6 .AND. nws==3) THEN
       !CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
       !              mesh%jj, mesh%jjs, mesh%rr)  
       CALL gauss_points_2d_p2(mesh)
    ELSE IF (nw==4 .AND. nws==3) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    ELSE IF (nw==10 .AND. nws==6) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    END IF

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_dg_mesh_free_format

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_free_format_ordered(dir, fil, list_dom, type_fe, mesh, mesh_formatted)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d_p1
    USE mod_gauss_points_2d_p2
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)
    IF (mesh_formatted) THEN
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    ELSE
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !     READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN    ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) kd, nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF



    !---  Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO m = 1, me
          IF (ABS(list_dom(dom)-i_d_lect(m)) /= 0)  CYCLE ! i_d(m) dans la liste
          virgin_el(m) = .FALSE.
          mnouv = mnouv + 1   ! Nouvel element
          nouv_el(m) = mnouv
          ancien_el(mnouv) = m
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
                nouv_nd(i) = nnouv
                ancien_nd(nnouv) = i
             END IF
          END DO
       END DO
    END DO
    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---  Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    !---  Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO ms = 1, mes
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
          END DO
          neighs2 = neigh_lect(n,neighs1)
          test = .FALSE.
          IF (ABS(list_dom(dom)-i_d_lect(neighs1)) == 0) THEN
             test=.TRUE.
          ELSE IF (neighs2 /= 0) THEN
             IF (ABS(list_dom(dom)-i_d_lect(neighs2)) == 0) THEN
                test=.TRUE.
                neighs_lect(ms) = neighs2 ! On change de cote
             END IF
          END IF

          IF (.NOT.test) CYCLE
          !11 June 2007 
          IF (.NOT.virgin_el(ms)) CYCLE 
          !11 June 2007 
          virgin_el(ms) = .FALSE.
          msnouv = msnouv + 1
          nouv_els(ms) = msnouv
          ancien_els(msnouv) = ms
       END DO
    END DO
    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---  Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !     END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !     CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !     mesh%jj, mesh%jjs, mesh%rr)

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       CALL gauss_points_2d_p1(mesh)
    ELSE IF (nw==6 .AND. nws==3) THEN
       !CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
       !              mesh%jj, mesh%jjs, mesh%rr)  
       CALL gauss_points_2d_p2(mesh)
    ELSE IF (nw==4 .AND. nws==3) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE IF (nw==10 .AND. nws==6) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_free_format_ordered

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_free_format(dir, fil, list_dom, type_fe, mesh, mesh_formatted)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d_p1
    USE mod_gauss_points_2d_p2
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)
    IF (mesh_formatted) THEN
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    ELSE
       OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes

       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF

       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       CALL gauss_points_2d_p1(mesh)
    ELSE IF (nw==6 .AND. nws==3) THEN
       !CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
       !              mesh%jj, mesh%jjs, mesh%rr)  
       CALL gauss_points_2d_p2(mesh)
    ELSE IF (nw==4 .AND. nws==3) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE IF (nw==10 .AND. nws==6) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_free_format

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_formatted(dir, fil, list_dom, type_fe, mesh)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d_p1
    USE mod_gauss_points_2d_p2
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh

    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    WRITE (*,*) 'Loading mesh-file ...'

    OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    OPEN(UNIT=20,FILE='error_mesh', FORM='formatted',STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       READ(30,*) np, nw, me, nws, mes
       !READ(30) np, nw, me, nws, mes
       DO m = 1, me
          READ(30,*)
       END DO
       !READ(30)
       DO ms = 1, mes
          READ(30,*)
       END DO
       !READ(30)
       DO n = 1, np
          READ(30,*)
       END DO
       !READ(30)
    END IF

    READ  (30, *)  np,  nw,  me,  nws,  mes
    !READ  (30)  np,  nw,  me,  nws,  mes

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    DO m = 1, me
       READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
    END DO

    !READ(30) jj_lect, neigh_lect, i_d_lect

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    DO ms = 1, mes
       READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    END DO

    !READ(30) jjs_lect, neighs_lect, sides_lect

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF
       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    DO n = 1, np
       READ(30,*) rr_lect(:,n)
    END DO
    !READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)   
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       CALL gauss_points_2d_p1(mesh)
    ELSE IF (nw==6 .AND. nws==3) THEN
       !CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
       !              mesh%jj, mesh%jjs, mesh%rr)  
       CALL gauss_points_2d_p2(mesh)
    ELSE IF (nw==4 .AND. nws==3) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE IF (nw==10 .AND. nws==6) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_formatted

  SUBROUTINE load_mesh(dir, fil, list_dom, type_fe, mesh)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d_p1
    USE mod_gauss_points_2d_p2
    !USE gauss_points
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       !READ(30,*) np, nw, me, nws, mes
       READ(30) np, nw, me, nws, mes

       !DO m = 1, me
       !   READ(30,*)
       !END DO
       READ(30)

       !DO ms = 1, mes
       !   READ(30,*)
       !END DO

       ! DO m = 1, mes 
       !     READ(30)
       ! END DO
       READ(30)

       !DO n = 1, np
       !   READ(30,*)
       !END DO
       READ(30)
    END IF

    !READ  (30, *)  np,  nw,  me,  nws,  mes
    READ(30)  np,  nw,  me,  nws,  mes

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    !DO m = 1, me
    !   READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
    !END DO

    READ(30) jj_lect, neigh_lect, i_d_lect

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    !DO ms = 1, mes
    !   READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    !END DO


    !TEST
    !DO ms = 1, mes
    !READ(30) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    !END DO
    !TEST
    READ(30) jjs_lect, neighs_lect, sides_lect

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes

       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF

       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    !DO n = 1, np
    !   READ(30,*) rr_lect(:,n)
    !END DO

    READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       CALL gauss_points_2d_p1(mesh)
    ELSE IF (nw==6 .AND. nws==3) THEN
       !CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
       !              mesh%jj, mesh%jjs, mesh%rr)  
       CALL gauss_points_2d_p2(mesh)
    ELSE IF (nw==4 .AND. nws==3) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE IF (nw==10 .AND. nws==6) THEN
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh


  SUBROUTINE surf_nodes_i(jjs, j_s,  iis)

    !  generation of the surface element connectivity matrix  iis
    !  based on the surface node numbering, starting from the
    !  connectivity matrix  jjs  of the surface elements according
    !  to the volume node numbering, and from the array  j_s  of
    !  the boundary nodes according to the volume node numbering

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjs
    INTEGER, DIMENSION(:),   INTENT(IN)  :: j_s
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iis

    INTEGER :: ms, ls, j, i

    DO ms = 1, SIZE(jjs,2)
       DO ls = 1, SIZE(jjs,1)
          j = jjs(ls,ms)
          DO i = 1, SIZE(j_s)
             IF ( j == j_s(i) )  iis(ls,ms) = i
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE surf_nodes_i

END MODULE prep_maill




