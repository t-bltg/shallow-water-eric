MODULE boundary_conditions
  !in 2D
  !11=Paraboloid + friction (Eric T.)
  !12=Paraboloid + coriolis force (Eric T.)
  !13=Hyperbolic SGN model (Eric T. )
  USE input_data
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: bath
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: velocity
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: one_over_h
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: regul_h
  REAL(KIND=8)                            :: max_water_h, aspect_ratio
  !===Malpasset
  REAL(KIND=8), DIMENSION(2,12)   :: malpasset_rr
  INTEGER,      DIMENSION(12)     :: malpasset_m
  CHARACTER(LEN=20), DIMENSION(12):: malpasset_file, malpasset_title
  INTEGER, DIMENSION(:), ALLOCATABLE :: new_to_old
  INTEGER, DIMENSION(1) :: vmax
  !===
CONTAINS

  FUNCTION flux(un) RESULT(vv)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np):: vv
    REAL(KIND=8) :: hh, vel, reg, umax=3.d0
    INTEGER :: n
    SELECT CASE(inputs%type_test)
    CASE(1,2,3,4,5,6,7,8,9,10,11,12) ! added case 11,12 (Eric T.)
       IF (inputs%type_test==9) THEN !malpasset
          umax = 30*SQRT(inputs%gravity*MAXVAL(un(1,:)))
       END IF
       one_over_h = compute_one_over_h(un(1,:))
       DO n = 1, mesh%np
          velocity(1,n) = un(2,n)*one_over_h(n)
          velocity(2,n) = un(3,n)*one_over_h(n)
          !===Limit velocity
          !vel = SQRT(velocity(1,n)**2+velocity(2,n)**2)
          !reg = MIN(umax/vel,1.d0)
          !velocity(:,n) = reg*velocity(:,n)
          !un(2,n) = velocity(1,n)*un(1,n) !This is AK's fix.
          !un(3,n) = velocity(2,n)*un(1,n) !It makes Arnau's test long-time stable on fine grid.
          !===Limit velocity
          vv(1,1,n) = velocity(1,n)*un(1,n) ! u h
          vv(1,2,n) = velocity(2,n)*un(1,n) ! v h
          vv(2,1,n) = velocity(1,n)*un(2,n) ! u uh
          vv(2,2,n) = velocity(2,n)*un(2,n) ! v uh
          vv(3,1,n) = velocity(1,n)*un(3,n) ! u vh
          vv(3,2,n) = velocity(2,n)*un(3,n) ! v vh

       END DO
    CASE(13) ! case 13 for modifed SGN hyperbolic model

       one_over_h = compute_one_over_h(un(1,:))
       DO n = 1, mesh%np
          ! set up velocity vector here
          velocity(1,n) = un(2,n)*one_over_h(n) ! u
          velocity(2,n) = un(3,n)*one_over_h(n)*0.d0 ! v
          velocity(3,n) = un(4,n)*one_over_h(n) ! eta
          velocity(4,n) = un(5,n)*one_over_h(n) ! w
          ! set flux terms here without pressure
          vv(1,1,n) = velocity(1,n)*un(1,n)   ! u h
          vv(1,2,n) = velocity(2,n)*un(1,n)   ! v h
          vv(2,1,n) = velocity(1,n)*un(2,n)   ! u uh
          vv(2,2,n) = velocity(2,n)*un(2,n)   ! v uh
          vv(3,1,n) = velocity(1,n)*un(3,n)   ! u vh
          vv(3,2,n) = velocity(2,n)*un(3,n)   ! v vh
          vv(4,1,n) = velocity(1,n)*un(4,n)   ! u eta*h
          vv(4,2,n) = velocity(2,n)*un(4,n)   ! v eta*h
          vv(5,1,n) = velocity(1,n)*un(5,n)   ! u wh
          vv(5,2,n) = velocity(2,n)*un(5,n)   ! v wh
        END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux

  FUNCTION compute_one_over_h(un) RESULT(vv)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: vv
    REAL(KIND=8) :: vel_square
    INTEGER :: n
    regul_h = inputs%htiny
    DO n = 1, mesh%np
       vv(n) = 2*un(n)/(un(n)**2+max(un(n),regul_h(n))**2)
    END DO
  END FUNCTION compute_one_over_h

  SUBROUTINE init(un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: un
    INTEGER :: i, k
    REAL(KIND=8) :: ray1, ray2, ray3, scale, h0, a, hL, eta, omega, &
         hcone, rcone, scone, htop, radius, bx, q0, p, kappa, s, h1, h2, z, D
    !===For malpasset
    INTEGER, DIMENSION(3,mesh%me) :: jj_old
    REAL(KIND=8), DIMENSION(mesh%np) :: bath_old
    INTEGER :: m
    !===
    ALLOCATE(bath(mesh%np))
    ALLOCATE(one_over_h(mesh%np))
    ALLOCATE(regul_h(mesh%np))
    !ALLOCATE(velocity(k_dim,mesh%np)) !this is original
    ALLOCATE(velocity(4,mesh%np)) ! for modified hyperbolic SGN model, u v eta w

    !===
    SELECT CASE(inputs%type_test)
    CASE(1) !Arnau's case (revisited)
       inputs%gravity=1.d0
       max_water_h=0.2d0
       DO i = 1, mesh%np
          ray1 = 25.d0*(mesh%rr(1,i)-1.2d0)**2+50.d0*(mesh%rr(2,i)-0.7d0)**2
          ray2 = 50.d0*(mesh%rr(1,i)-0.45d0)**2+100.d0*(mesh%rr(2,i)-0.4d0)**2

          IF (ray1 .LE. 2.d0) THEN
             bath(i) = ((EXP(-ray1)-EXP(-1.d0))/(EXP(-2.d0)-EXP(-1.d0))-1.d0)*(-0.2d0)
             if (ABS(SQRT(ray1)-1.d0)<1.d-10) bath(i) = 0.2d0
          ELSE IF (ray2 .LE. 1.d0) THEN
             bath(i) = -0.5*(EXP(-ray2)-EXP(-1.d0))
          ELSE
             bath(i) = 0.d0
          END IF
       END DO
    CASE(2,3,6,7) !Bump 3.1 (3.1.3 sub, 3.1.4 trans, 3.1.5 critical)
       inputs%gravity=9.81d0
       IF (inputs%type_test==2) max_water_h=0.2d0
       IF (inputs%type_test==3) max_water_h=2.d0
       IF (inputs%type_test==6) max_water_h=2.d0
       IF (inputs%type_test==7) max_water_h=2.d0
       DO i = 1, mesh%np
          IF ((8.d0<mesh%rr(1,i)) .AND. (mesh%rr(1,i)<12.d0)) THEN
             bath(i) = 0.2d0*(mesh%rr(1,i)-8.d0)**3*(12.d0-mesh%rr(1,i))**3/64.d0
          ELSE
             bath(i) = 0.d0
          END IF
       END DO
    CASE(4) !Dam break wet/dry, 4.1.2 (Delestre)
       inputs%time = 1.d0
       inputs%gravity=9.81d0
       bath = 0.d0
       hL=0.005d0
       max_water_h=hL
    CASE(5) !Planar surface in a paraboloid
        inputs%gravity=9.81d0
        a=1.d0
        h0=0.1d0
        hL=4.d0
        eta=0.5d0
        max_water_h=hL
        omega=SQRT(2*inputs%gravity*h0)/a
        inputs%Tfinal = 6.d0*ACOS(-1.d0)/omega
        !inputs%Tfinal = 2*6.d0*ACOS(-1.d0)/omega
        bath = -h0*(1.d0-((mesh%rr(1,:)-hL/2.d0)**2+(mesh%rr(2,:)-hL/2.d0)**2)/a**2)
    CASE(8) !Solitary wave around an island
        inputs%gravity=9.81d0
        htop=0.625d0
        rcone=3.6d0
        scone=4.d0
        hcone=0.9d0
        h0=0.32
        max_water_h=h0
        DO i = 1, mesh%np
           radius = SQRT((mesh%rr(1,i)-15)**2 + (mesh%rr(2,i)-13)**2)
           IF (radius<rcone) THEN
              bath(i) = min(htop,hcone-radius/scone)
           ELSE
              bath(i) = 0.d0
           END IF
        END DO
    CASE(9) ! Malpasset
        inputs%gravity=9.81d0
        OPEN(unit=30,file='malpasset.wh',FORM='unformatted')
        READ(30) bath_old
        bath_old = MAX(bath_old,0.d0)
        CLOSE(30)
        max_water_h=1.d0 !!MAXVAL(bath_old)
        OPEN(unit=30,file='malpasset.jj',FORM='unformatted')
        READ(30) jj_old
        ALLOCATE(new_to_old(mesh%np))
        DO m = 1, mesh%me
           new_to_old(mesh%jj(:,m))=jj_old(:,m)
        END DO
        OPEN(unit=30,file='malpasset.bath',FORM='unformatted')
        READ(30) bath_old
        bath = bath_old(new_to_old)
        CLOSE(30)
        !===Find Gauges
        CALL malpasset_gauges
    CASE(11) ! planar surface in a  2d paraboloid with friction (Eric T.)
       inputs%gravity = 9.81d0
       a = 1.d0
       h0 = 0.1d0
       hL = 4.d0 ! L = 4 for this case
       max_water_h = hL
       kappa = 0.2d0
       p = SQRT((8*inputs%gravity*h0)/a**2)
       s = 1.d0/2.d0 * SQRT(p**2 - kappa**2)
       inputs%Tfinal = 2*2*ACOS(-1.d0)/s
       bath = -h0*( 1.d0-((mesh%rr(1,:)-hL/2.d0)**2+(mesh%rr(2,:)-hL/2.d0)**2)/a**2 )
    CASE(12) ! planar surface in a 2d paraboloid with coriolis force (Eric T.)
        inputs%gravity = 9.81d0
        a=3.d0
        h0=0.5d0
        hL = 0.d0 ! L = 0 for this case since disk is centered at the origin
        max_water_h = 6.d0
        bath = -h0*(1.d0-((mesh%rr(1,:)-hL/2.d0)**2+(mesh%rr(2,:)-hL/2.d0)**2)/a**2)
    CASE(13) ! modified hyperbolic SGN model (Eric T., 2/26/2018)
         inputs%gravity = 9.81d0
         h1 = 5.0d0 / 100.0d0
         h2 = 6.0d0 / 100.0d0
         max_water_h = h2
         bath = 0.d0
         D = SQRT(inputs%gravity * h2) ! constant wave velocity
         z = SQRT( (3.d0 * (h2 - h1) / (h2 * h1**2.d0)) )
         DO i = 1, mesh%np
           bath(i) = 0.d0
         END DO

    CASE DEFAULT
       WRITE(*,*) ' BUG in init'
       STOP
    END SELECT

    inputs%htiny=inputs%epsilon*max_water_h
    DO k = 1, inputs%syst_size
       un(k,:) = sol_anal(k,mesh%rr,inputs%time)
    END DO
    ray1=MAXVAL(mesh%rr(1,:))-MINVAL(mesh%rr(1,:))
    ray2=MAXVAL(mesh%rr(2,:))-MINVAL(mesh%rr(2,:))
    ray3=SQRT(SUM(mesh%gauss%rj))
    aspect_ratio=MAXVAL(un(1,:))/MAX(ray1,ray2,ray3)
  END SUBROUTINE init

  SUBROUTINE malpasset_gauges
    USE mesh_handling
    USE mesh_interpolation
    IMPLICIT NONE
    INTEGER :: m, n
    malpasset_rr(1,1) = 5500    !Trans A
    malpasset_rr(2,1) = 4400    !Trans A
    malpasset_rr(1,2) = 11900   !Trans B
    malpasset_rr(2,2) = 3250    !Trans B
    malpasset_rr(1,3) = 13000   !Trans C
    malpasset_rr(2,3) = 2700    !Trans C
    malpasset_rr(1,4) = 4947.46 !S6
    malpasset_rr(2,4) = 4289.71 !S6
    malpasset_rr(1,5) = 5717.30 !S7
    malpasset_rr(2,5) = 4407.61 !S7
    malpasset_rr(1,6) = 6775.14 !S8
    malpasset_rr(2,6) = 3869.23 !S8
    malpasset_rr(1,7) = 7128.20 !S9
    malpasset_rr(2,7) = 3162.00 !S9
    malpasset_rr(1,8) = 8585.30 !S10
    malpasset_rr(2,8) = 3443.08 !S10
    malpasset_rr(1,9) = 9674.97 !S11
    malpasset_rr(2,9) = 3085.89 !S11
    malpasset_rr(1,10)= 10939.15!S12
    malpasset_rr(2,10)= 3044.78 !S12
    malpasset_rr(1,11)= 11724.37!S13
    malpasset_rr(2,11)= 2810.41 !S13
    malpasset_rr(1,12)= 12723.70!S14
    malpasset_rr(2,12)= 2485.08 !S14
    DO n = 1, 12
       write(*,*) ' n', n
       DO m = 1, mesh%me
          IF (inside(mesh,m,malpasset_rr(:,n))) EXIT
       END DO
       malpasset_m(n)=m
       !write(*,*) malpasset_m(n)
       !write(*,*) mesh%rr(1,mesh%jj(:,malpasset_m(n))), malpasset_rr(1,n)
       !write(*,*) mesh%rr(2,mesh%jj(:,malpasset_m(n))), malpasset_rr(2,n)
       WRITE(*,*) bath(mesh%jj(:,malpasset_m(n))), SUM(bath(mesh%jj(:,malpasset_m(n))) &
            * FE_interpolation(mesh,malpasset_m(n),malpasset_rr(1:2,n)))
    END DO
    malpasset_file(1)  ='Elec_trans_A.plt'
    malpasset_file(2)  ='Elec_trans_B.plt'
    malpasset_file(3)  ='Elec_trans_C.plt'
    malpasset_file(4)  ='Gauge_S6.plt'
    malpasset_file(5)  ='Gauge_S7.plt'
    malpasset_file(6)  ='Gauge_S8.plt'
    malpasset_file(7)  ='Gauge_S9.plt'
    malpasset_file(8)  ='Gauge_S10.plt'
    malpasset_file(9)  ='Gauge_S11.plt'
    malpasset_file(10) ='Gauge_S12.plt'
    malpasset_file(11) ='Gauge_S13.plt'
    malpasset_file(12) ='Gauge_S14.plt'
    malpasset_title(1) ='%Elec trans A'
    malpasset_title(2) ='%Elec trans B'
    malpasset_title(3) ='%Elec trans C'
    malpasset_title(4) ='%Gauge S6'
    malpasset_title(5) ='%Gauge S7'
    malpasset_title(6) ='%Gauge S8'
    malpasset_title(7) ='%Gauge S9'
    malpasset_title(8) ='%Gauge S10'
    malpasset_title(9) ='%Gauge S11'
    malpasset_title(10)='%Gauge S12'
    malpasset_title(11)='%Gauge S13'
    malpasset_title(12)='%Gauge_S14'
  END SUBROUTINE malpasset_gauges

  FUNCTION sol_anal(k,rr,t) RESULT(vv) !CAREFUL HERE, rr comes from the calling subroutine
    USE mesh_handling
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN) :: rr
    REAL(KIND=8),                  INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: vv
    INTEGER :: i
    REAL(KIND=8) :: x, x0, speed, q0, hL, b, d, x1, x2, x3, &
         theta, Rcard, Scard, Tcard, Qcard, Dcard, tpio3, fpio3, a, omega, eta, h0, bernoulli, &
         xshock, h_pre_shock, h_post_shock, bath_shock, bathi, Ber_pre, Ber_post, &
         alpha, beta, chi, vel, xs, hcone, htop, radius, rcone, scone, bx, kappa, &
         s, htilde, p, f, delta, w_star, den, gamma, h1, h2, lambdaSGN, z
    !===Malpasset
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: h_old
    !===
    SELECT CASE(inputs%type_test)
    CASE(1)
       SELECT CASE(k)
       CASE(1)
          DO i = 1, SIZE(rr,2)
             vv(i)  = MAX(MAX(bath(i),0.2d0)-bath(i),0.d0)  !+ inputs%htiny
          END DO
       CASE DEFAULT
          vv(:) = 0.d0
       END SELECT
    CASE(2) ! 3.1.2
       SELECT CASE(k)
       CASE(1)
          DO i = 1, SIZE(rr,2)
             vv(i) = MAX(MAX(bath(i),0.1d0)-bath(i),0.d0) !+ inputs%htiny
          END DO
       CASE DEFAULT
          vv(:) = 0.d0
       END SELECT
    CASE(3) !Subcritical flow over a bump, Delestre 3.1.3
       q0 = 4.42d0
       hL = 2.d0
       tpio3=2*ACOS(-1.d0)/3
       fpio3=2*tpio3
       SELECT CASE(k)
       CASE(1)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                vv(i) = hL - bathi
             END DO
          ELSE
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                b = bathi-q0**2/(2*inputs%gravity*hL**2) -hL
                d = q0**2/(2*inputs%gravity)
                Rcard = (-27*d-2*b**3)/(54)
                Qcard = -b**2/9
                Dcard = Qcard**3+Rcard**2
                !write(*,*) 'Dcard', Dcard
                theta = ACOS(Rcard/SQRT(-Qcard**3))
                x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
                !write(*,*) x1**3+b*x1**2+d
                !x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
                !x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
                !write (*,*) x1, x2, x3
                vv(i) = x1
             END DO
          END IF
       CASE(2)
          IF (t.LE.1.d-14) THEN
             vv =0.d0
          ELSE
             vv = q0
          END IF
       CASE DEFAULT
          vv = 0.d0
       END SELECT
    CASE(4) !Dam break wet/dry, 4.1.2 (Delestre)
       x0=5.d0
       hL=0.005d0
       speed = SQRT(inputs%gravity*hL)
       SELECT CASE(k)
       CASE(1) ! water height
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                IF (rr(1,i) .LE. x0) THEN
                   vv(i) = hL
                ELSE
                   vv(i) = 0.d0 !+ inputs%htiny
                END IF
             END DO
          ELSE
             DO i = 1, SIZE(rr,2)
                x = rr(1,i)-x0
                IF (x .LE. -t*speed) THEN
                   vv(i) = hL
                ELSE IF (x.LE.2*t*speed) THEN
                   vv(i) = (4.d0/(9.d0*inputs%gravity))*(speed-x/(2*t))**2
                ELSE
                   vv(i) = 0.d0 + inputs%htiny
                END IF
             END DO
          END IF
       CASE(2)
          IF (t.LE.1.d-14) THEN
             vv = 0.d0
          ELSE
             DO i = 1, SIZE(rr,2)
                x = rr(1,i)-x0
                IF (x.LE.-t*speed) THEN
                   vv(i) = 0.d0
                ELSE IF (x.LE.2*t*speed) THEN
                   vv(i) = (8.d0/(27.d0*inputs%gravity))*(speed+x/t)*(speed-x/(2*t))**2
                ELSE
                   vv(i) = 0.d0
                END IF
             END DO
          END IF
       CASE DEFAULT
          vv(:) = 0.d0
       END SELECT
    CASE(5)
        a=1.d0
        h0=0.1d0
        hL=4.d0
        eta=0.5d0
        omega=SQRT(2*inputs%gravity*h0)/a
        SELECT CASE(k)
        CASE(1) ! water height
           DO i = 1, SIZE(rr,2)
              bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)
              vv(i) = MAX((eta*h0/a**2)*(2.d0*(rr(1,i)-hL/2.d0)*COS(omega*t) &
                   +2.d0*(rr(2,i)-hL/2.d0)*SIN(omega*t)) - bathi,0.d0) + inputs%htiny
           END DO
        CASE(2)
           DO i = 1, SIZE(rr,2)
              bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)
              vv(i) = MAX((eta*h0/a**2)*(2.d0*(rr(1,i)-hL/2.d0)*COS(omega*t) &
                   +2.d0*(rr(2,i)-hL/2.d0)*SIN(omega*t)) - bathi,0.d0)
              vv(i) = -eta*omega*SIN(omega*t)*vv(i)
           END DO
        CASE(3)
           DO i = 1, SIZE(rr,2)
              bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)
              vv(i) = MAX((eta*h0/a**2)*(2.d0*(rr(1,i)-hL/2.d0)*COS(omega*t) &
                   +2.d0*(rr(2,i)-hL/2.d0)*SIN(omega*t)) - bathi,0.d0)
              vv(i) = eta*omega*COS(omega*t)*vv(i)
           END DO
        END SELECT
    CASE(6) !Transcritical flow over a bump, Delestre ??
       x0 = 10.d0
       q0 = 1.53d0
       hL = 0.66d0
       tpio3=2*ACOS(-1.d0)/3
       fpio3=2*tpio3
       !TEST
       b = 0.2d0-.2d0-q0**2/(2*inputs%gravity*hL**2) -hL
       d = q0**2/(2*inputs%gravity)
       eta = -(27*d/4)**(1.d0/3)
       bernoulli = 0.2d0  - eta


       SELECT CASE(k)
       CASE(1)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                vv(i) = hL - bathi
             END DO
          ELSE
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                b = bathi - bernoulli
                d = q0**2/(2*inputs%gravity)
                Rcard = (-27*d-2*b**3)/(54)
                Qcard = -b**2/9
                Dcard = Qcard**3+Rcard**2
                IF (Dcard.GE.0.d0) THEN
                   !write(*,*) Rcard**2,Qcard**3,Qcard**3+Rcard**2
                   !write(*,*) Rcard + sqrt(Dcard)
                   Scard=-(-SQRT(Dcard)-Rcard)**(1.d0/3)
                   Tcard=-(SQRT(Dcard)-Rcard)**(1.d0/3)
                   x1 = Scard+Tcard-b/3.d0
                   x2 = -(Scard+Tcard)/2.d0 - b/3.d0
                   !write(*,*) x1**3+b*x1**2+d
                   vv(i) = x2
                   write(*,*) 'Dcard', Dcard, bath(i), x2
                ELSE
                   !write(*,*) 'ratio', Rcard/SQRT(-Qcard**3), 'Dcard', Dcard, '-Qcard', -Qcard
                   theta = ACOS(Rcard/SQRT(-Qcard**3))
                   x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
                   !write(*,*) x1**3+b*x1**2+d
                   x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
                   x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
                   !write (*,*) 'x1,x2,x3',x1, x2, x3
                   IF (rr(1,i)<x0) THEN
                      vv(i) = MAX(x1,x3)
                   ELSE
                      vv(i) = MIN(x1,x3)
                   END IF
                   !write(*,*) 'vv(i)', vv(i)
                END IF
             END DO
          END IF
          !stop
       CASE(2)
          IF (t.LE.1.d-14) THEN
             vv = 0.d0
          ELSE
             vv = q0
          END IF
       CASE DEFAULT
          vv = 0.d0
       END SELECT
    CASE(7) !critical flow over a bump, Delestre 3.1.5
       x0= 10.d0
       xshock = 11.7d0
       q0 = 0.18d0
       tpio3=2*ACOS(-1.d0)/3
       fpio3=2*tpio3
       d = q0**2/(2*inputs%gravity)
       eta = -(27*d/4)**(1.d0/3)
       Ber_Pre = 0.2d0  - eta
       IF (8.d0< xshock .AND. xshock<12.0d0) THEN
          bath_shock = 0.2d0*(xshock -8.d0)**3*(12.d0-xshock )**3/64.d0
          !bath_shock = 0.2d0 - 0.05*(xshock -10.d0)**2
       ELSE
          bath_shock = 0.d0
       END IF
       b = bath_shock - Ber_Pre
       Rcard = (-27*d-2*b**3)/(54)
       Qcard = -b**2/9
       Dcard = Qcard**3+Rcard**2
       IF (Dcard.GE.0.d0) THEN
          WRITE(*,*) 'BUG case(7)'
          !STOP
       ELSE
          theta = ACOS(Rcard/SQRT(-Qcard**3))
          x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
          x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
          x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
          h_pre_shock = x3
          h_post_shock = (-x3+SQRT(x3**2+8*q0**2/(inputs%gravity*x3)))/2
       END IF
       Ber_Post = (q0**2/(2*inputs%gravity))/h_post_shock**2 + h_post_shock + bath_shock
       !write(*,*) ' h_pre_shock,  h_post_shock', h_pre_shock,  h_post_shock
       !write(*,*) ' Ber_pre_shock,  Ber_post_shock', Ber_pre, Ber_post
       !write(*,*) ' H_L',  bath_shock + h_post_shock
       !stop
       SELECT CASE(k)
       CASE(1)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                !bathi =  MAX(0.2d0 - 0.05*(rr(1,i)-10.d0)**2,0d0)
                vv(i) = h_post_shock - bathi
             END DO
          ELSE
             DO i = 1, SIZE(rr,2)
                bathi = MAX(0.d0,0.2d0*(rr(1,i)-8.d0)**3*(12.d0-rr(1,i))**3/64.d0)
                !bathi = MAX(0.2d0 - 0.05*(rr(1,i)-10.d0)**2,0d0)
                b = bathi - ber_Pre
                d = q0**2/(2*inputs%gravity)
                Rcard = (-27*d-2*b**3)/(54)
                Qcard = -b**2/9
                Dcard = Qcard**3+Rcard**2
                IF (Dcard.GE.0.d0) THEN
                   Scard=-(-SQRT(Dcard)-Rcard)**(1.d0/3)
                   Tcard=-(SQRT(Dcard)-Rcard)**(1.d0/3)
                   x1 = Scard+Tcard-b/3.d0
                   x2 = -(Scard+Tcard)/2.d0 - b/3.d0
                   vv(i) = x2
                   !write(*,*) 'BUG case(7)',x2,  Dcard
                   !STOP
                ELSE
                   theta = ACOS(Rcard/SQRT(-Qcard**3))
                   x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
                   x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
                   x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
                   IF (rr(1,i)<x0) THEN
                      vv(i) = MAX(x1,x3)
                   ELSE IF (rr(1,i)<xshock) THEN
                      vv(i) = MIN(x1,x3)
                   ELSE
                      vv(i) = bath_shock + h_post_shock - bathi
                   END IF
                END IF
             END DO
          END IF
       CASE(2)
          IF (t.LE.1.d-14) THEN
             vv = q0 !0.d0
          ELSE
             vv = q0
          END IF
       CASE DEFAULT
          vv = 0.d0
       END SELECT
    CASE(8) !Solitary wave around an island
       h0 = 0.32
       alpha = 1*h0
       xs = 15.0d0-12.96d0
       vel = SQRT(3.d0*alpha/(4*h0**3))
       SELECT CASE(k)
       CASE(1)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                vv(i) = MAX(h0+alpha/(COSH(vel*(rr(1,i)-xs))**2)&
                     -bath(i),0.d0)
             END DO
          ELSE
             vv = h0
          END IF
       CASE(2)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                vv(i) = alpha/(COSH(vel*(rr(1,i)-xs))**2)
                vv(i) = vv(i)*SQRT(inputs%gravity/h0)
             END DO
          ELSE
             vv =0.d0
          END IF
       CASE(3)
          vv = 0.d0
       END SELECT
    CASE(9) ! Malpasset
        SELECT CASE(k)
        CASE(1)
           IF (t.LE.1.d-14) THEN
              OPEN(unit=30,file='malpasset.wh',FORM='unformatted')
              !DO i = 1, SIZE(rr,2)
              READ(30) h_old
              !END DO
              CLOSE(30)
              vv = MAX(h_old(new_to_old),0.d0)
           ELSE
              vv = 0.d0
           END IF
        CASE(2,3)
           vv = 0.d0
        END SELECT
    CASE(11) ! Added to accomodate for paraboid with friction (Eric T.)
     inputs%gravity = 9.81d0
     a = 1.d0
     h0 = 0.1d0
     hL = 4.d0 ! L = 4 for this case
     kappa = 0.2d0
     gamma = 0.5d0
     p = SQRT( (8 * inputs%gravity * h0)/a**2 )
     s = 1.d0/2.d0 * SQRT(p**2 - kappa**2)
     SELECT CASE(k)
     CASE(1) ! water height h
        DO i = 1, SIZE(rr,2)
           bathi = -h0*(1.d0- ((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2) /a**2 )

           htilde = -bathi - (gamma**2)/(2*inputs%gravity)*EXP(-kappa*t) &
           - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(1,i)-hL/2.d0) * &
           (2*s*COS(s*t) + kappa*SIN(s*t)) &
           - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(2,i)-hL/2.d0) * &
           (kappa*COS(s*t)-2.d0*s*SIN(s*t))

           vv(i) = MAX(htilde,0.d0)
        END DO
     CASE(2) ! u*h component of flow rate q
        DO i = 1, SIZE(rr,2)
          bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)

          htilde =  -bathi - (gamma**2)/(2*inputs%gravity)*EXP(-kappa*t) &
          - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(1,i)-hL/2.d0) * &
          (2*s*COS(s*t) + kappa*SIN(s*t)) &
          - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(2,i)-hL/2.d0) * &
          (kappa*COS(s*t)-2.d0*s*SIN(s*t))

           vv(i) = MAX(htilde,0.d0)

           vv(i) = vv(i) * gamma * EXP(-kappa*t/2.d0) * SIN(s*t)
        END DO
     CASE(3) ! v*h component of flow rate q
        DO i = 1, SIZE(rr,2)
          bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)

          htilde = -bathi - (gamma**2)/(2*inputs%gravity)*EXP(-kappa*t) &
          - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(1,i)-hL/2.d0) * &
          (2*s*COS(s*t) + kappa*SIN(s*t)) &
          - gamma/(2.d0*inputs%gravity)*EXP(-kappa*t/2.d0)*(rr(2,i)-hL/2.d0) * &
          (kappa*COS(s*t)-2.d0*s*SIN(s*t))

          vv(i) = MAX(htilde,0.d0)

          vv(i) = vv(i) * gamma * EXP(-kappa*t/2.d0) * COS(s*t)

        END DO
     END SELECT
    CASE(12) ! Added to accomodate for paraboid with coriolis force (Eric T.)
     inputs%gravity = 9.81d0
     a = 3.d0
     h0 = 0.5d0
     hL = 0.d0 ! L = 4 for this case
     f = 2.d0
     eta = 0.5d0
     den = f**2 * a**4 + SQRT(8.d0 * f**2 * inputs%gravity  * h0 * a**6 + f**4 * a**8)
     delta = 2.d0 * f * a**4 / (den)
     w_star = -1.d0/(SQRT(2.d0) * a**2) * SQRT(4.d0*inputs%gravity*h0* a**2 + den)
      SELECT CASE(k)
      CASE(1) ! water height
        DO i = 1, SIZE(rr,2)
            bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)

            htilde =  -bathi - eta/inputs%gravity * COS(w_star * t) * (1 - delta*f) * (rr(2,i)-hL/2.d0) &
            + eta/inputs%gravity * SIN(w_star * t) * (f/w_star - delta*w_star) * (rr(1,i)-hL/2.d0) &
            + eta**2 / (2.d0 * inputs%gravity) * (delta**2 *COS(w_star * t)*COS(w_star * t) + &
            SIN(w_star * t)*SIN(w_star * t)/(w_star**2) )

            vv(i) = MAX(htilde,0.d0)
          END DO
        CASE(2) ! u*h component of flow rate q
          DO i = 1, SIZE(rr,2)
            bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)

            htilde =  -bathi - eta/inputs%gravity * COS(w_star * t) * (1 - delta*f) * (rr(2,i)-hL/2.d0) &
            + eta/inputs%gravity * SIN(w_star * t) * (f/w_star - delta * w_star) * (rr(1,i)-hL/2.d0) &
            + eta**2 / (2.d0 * inputs%gravity) * (delta**2 *COS(w_star * t)*COS(w_star * t) + &
            SIN(w_star * t)*SIN(w_star * t)/(w_star**2) )

            vv(i) = MAX(htilde,0.d0)

            vv(i) = -vv(i) * delta * eta * COS(w_star * t)

          END DO
        CASE(3) ! v*h component of flow rate q
          DO i = 1, SIZE(rr,2)
            bathi = -h0*(1.d0-((rr(1,i)-hL/2.d0)**2+(rr(2,i)-hL/2.d0)**2)/a**2)

            htilde =  -bathi - eta/inputs%gravity * COS(w_star * t) * (1 - delta*f) * (rr(2,i)-hL/2.d0) &
            + eta/inputs%gravity * SIN(w_star * t) * (f/w_star - delta * w_star) * (rr(1,i)-hL/2.d0) &
            + eta**2 / (2.d0 * inputs%gravity) * (delta**2 *COS(w_star * t)*COS(w_star * t) + &
            SIN(w_star * t)*SIN(w_star * t)/(w_star**2) )

            vv(i) = MAX(htilde,0.d0)
            vv(i) = vv(i) * eta * SIN(w_star * t)/w_star

          END DO
        END SELECT
    CASE(13) ! Added to include hyperbolid SGN model (Eric T., 02/2018)
      ! here we are doing exact solitary wave solution from Favrie-Gavrilyuk paper
      ! initial constants go here
      inputs%gravity = 9.81d0
      h1 = 5.0d0 / 100.0d0
      h2 = 6.d0 / 100.0d0
      x0 = 1.75d0  ! we want largest solitary wave height starting here
      D = SQRT(inputs%gravity * h2) ! constant wave velocity
      z = SQRT( ( 3.0d0 * (h2 - h1) / (h2 * h1**2.0d0) ) )

      SELECT CASE(k)
      CASE(1) ! h water height
          ! for initial condition on water height
          !IF (t.LE.1.d-10) THEN
            DO i = 1, SIZE(rr,2)
              bathi = 0.d0
              htilde= h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
              vv(i) = max(htilde,0.d0)
          !  END DO
          ! ELSE ! exact solution
          !   DO i = 1, SIZE(rr,2)
          !     bathi = 0.d0
          !     htilde= h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
          !     vv(i) = max(htilde,0.d0)
            END DO
         !END IF

      CASE(2) ! u*h component of flow rate q
        ! for initial velocity u
        !IF (t.LE.1.d-10) THEN
          DO i = 1, SIZE(rr,2)
            bathi = 0.d0
            htilde = h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
            vv(i) = MAX(htilde,0.d0)
            vv(i) = vv(i) * D - h1
          END DO
        !ELSE
          ! DO i = 1, SIZE(rr,2)
          !   bathi = 0.d0
          !   htilde =  h1 + (h2 - h1)*1.0d0/COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0
          !   vv(i) = MAX(htilde,0.d0)
          !   vv(i) = vv(i) * D - h1
          ! END DO

       !END IF
      CASE(3) ! v*h component of flow rate q, just 0 for now
       ! for initial velocity u
         DO i = 1, SIZE(rr,2)
           bathi = 0.d0
           htilde =  h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
           vv(i) = MAX(htilde,0.d0)
           vv(i) = vv(i) * 0.d0
         END DO
      CASE(4) ! eta*h component of flow rate r
         ! initial condition from FAVRIE/GAVRILYUK paper
         IF (t.LE.1.d-10) THEN
           DO i = 1, SIZE(rr,2)
             bathi = 0.d0
             htilde = h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
             vv(i) = MAX(htilde,0.d0)
             vv(i) = vv(i)*vv(i)
           END DO
         END IF
      CASE(5) ! w*h component of flow rate r
        IF (t.LE.1.d-10) THEN
             DO i = 1, SIZE(rr,2)
               bathi = 0.d0
               htilde = h1 + (h2 - h1)*1.0d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.0d0)
               vv(i) = MAX(htilde,0.d0)
               vv(i) = -d * h1 * ( (h2 - h1) * z * 1.d0/(COSH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t))**2.d0) &
               * TANH(1.0d0/2.0d0*z*(rr(1,i)-x0-D*t)) )
             END DO
       END IF
      END SELECT

   CASE DEFAULT
      WRITE(*,*) ' BUG in sol_anal'
      STOP
   END SELECT

  END FUNCTION sol_anal

  SUBROUTINE compute_lambda(ul,ur,nij,lambda)
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: lambda
    REAL(KIND=8), INTENT(IN), DIMENSION(k_dim)  :: nij
    REAL(KIND=8), DIMENSION(inputs%syst_size)   :: ur, ul
    REAL(KIND=8) :: ht, vl, vr, hl, hr, lbdl, lbdr, lbd_dry, sql, sqr
    hl = ul(1)
    hr = ur(1)
    vl =  SUM(ul(2:k_dim+1)*nij)*(2*hl/(hl**2+max(hl,inputs%htiny)**2))
    vr =  SUM(ur(2:k_dim+1)*nij)*(2*hr/(hr**2+max(hr,inputs%htiny)**2))
    sql = SQRT(inputs%gravity*ABS(hl))
    sqr = SQRT(inputs%gravity*ABS(hr))

    SELECT CASE(inputs%type_test)
    CASE(1,2,3,4,5,6,7,8,11,12,13) ! added cases 11 and 12 (Eric T. )

       IF (hl.LE.inputs%htiny .AND. hr.LE.inputs%htiny) THEN
          !lambda = 0.d0
          lambda = MAX(ABS(vl - 2*sql),ABS(vr + 2*sqr))
       ELSE IF (hl.GE.inputs%htiny .AND. hr.GE.inputs%htiny) THEN
          !ht = (MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity)
          lbdl = vl - sql!*SQRT((1+max((ht-hl)/(2*hl),0.d0))*(1+max((ht-hl)/hl,0.d0)))
          lbdr = vr + sqr!*SQRT((1+max((ht-hr)/(2*hr),0.d0))*(1+max((ht-hr)/hr,0.d0)))
          lambda = MAX(ABS(lbdl),ABS(lbdr))
       ELSE IF  (hl.LE.inputs%htiny .AND. hr.GE.inputs%htiny) THEN
          !ht = (MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity)
          lbdl = vl - sql!*SQRT((1+max((ht-hl)/(2*hl),0.d0))*(1+max((ht-hl)/hl,0.d0)))
          lbdr = vr + sqr!*SQRT((1+max((ht-hr)/(2*hr),0.d0))*(1+max((ht-hr)/hr,0.d0)))
          lbd_dry =  vl - 2*sql
          lambda = MAX(ABS(lbdl),ABS(lbd_dry),ABS(lbdr))
       ELSE IF  (hl.GE.inputs%htiny .AND. hr.LE.inputs%htiny) THEN
          !ht = (MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity)
          lbdl = vl - sql!*SQRT((1+max((ht-hl)/(2*hl),0.d0))*(1+max((ht-hl)/hl,0.d0)))
          lbdr = vr + sqr!*SQRT((1+max((ht-hr)/(2*hr),0.d0))*(1+max((ht-hr)/hr,0.d0)))
          lbd_dry =  vr + 2*sqr
          lambda = MAX(ABS(lbdl),ABS(lbd_dry),ABS(lbdr))
       END IF
    CASE DEFAULT
       WRITE(*,*) 'BUG in compute_lambda'
       STOP
    END SELECT
    RETURN
  END SUBROUTINE compute_lambda


  SUBROUTINE compute_lambda_vacc(ul,ur,nij,lambda)
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: lambda
    REAL(KIND=8), INTENT(IN), DIMENSION(k_dim)  :: nij
    REAL(KIND=8), DIMENSION(inputs%syst_size)   :: ur, ul
    REAL(KIND=8) :: fl, fr, ht, vl, vr, hl, hr, lbdl, lbdr, lbd_dry, sql, sqr
    REAL(KIND=8) :: fh, a, c, Delta, x0, hmin, hmax, vmin, vmax, sqrmin, sqrmax
    REAL(KIND=8) :: ovhl, ovhr
    hl = ul(1)
    hr = ur(1)
    ovhl = 2*hl/(hl**2+max(hl,inputs%htiny)**2)
    ovhr = 2*hr/(hr**2+max(hr,inputs%htiny)**2)
    vl =  SUM(ul(2:k_dim+1)*nij)*ovhl
    vr =  SUM(ur(2:k_dim+1)*nij)*ovhr
    sql = SQRT(inputs%gravity*ABS(hl))
    sqr = SQRT(inputs%gravity*ABS(hr))

    SELECT CASE(inputs%type_test)
    CASE(1,2,3,4,5,6,7,8,9,10,11,12,13)
       !IF (hl.LE.inputs%htiny .AND. hr.LE.inputs%htiny) THEN
       !   lambda = 0.d0
       !   RETURN
       !END IF

       x0=(2.d0*SQRT(2.d0)-1.d0)**2
       IF(hl.LE.hr) THEN
          hmin = hl
          vmin = vl
          hmax = hr
          vmax = vr
       ELSE
          hmin = hr
          vmin = vr
          hmax = hl
          vmax = vl
       END IF

       !===Case 1
       fh = phi(x0*hmin,vl,hl,vr,hr)
       IF (0.d0.LE.fh) THEN
          ht = min(x0*hmin,(MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity))
          lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
          lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
          lambda = MAX(ABS(lbdl),ABS(lbdr))
          RETURN
       END IF
       !===Case 2
       fh = phi(x0*hmax,vl,hl,vr,hr)
       IF (0.d0.LE.fh) THEN
          sqrmin = SQRT(hmin)
          sqrmax = SQRT(hmax)
          a = 1/(2.d0*SQRT(2.d0))
          c = -hmin*a -sqrmin*sqrmax + sqrmin*(vr-vl)/(2*sqrt(inputs%gravity))
          Delta = hmin-4*a*c
          IF (Delta<0.d0) THEN
             WRITE(*,*) ' BUG in compute_lambda_vacc'
             STOP
          END IF
          ht = min(x0*hmax,((-sqrmin+sqrt(Delta))/(2*a))**2)
          lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
          lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
          lambda = MAX(ABS(lbdl),ABS(lbdr))
          RETURN
       END IF
       !===Case 3
       sqrmin = SQRT(hmin)
       sqrmax = SQRT(hmax)
       ht = sqrmin*sqrmax*(1+sqrt(2/inputs%gravity)*(vl-vr)/(sqrmin+sqrmax))
       lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
       lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
       lambda = MAX(ABS(lbdl),ABS(lbdr))

    CASE DEFAULT
       WRITE(*,*) 'BUG in compute_lambda'
       STOP
    END SELECT
    RETURN
  END SUBROUTINE compute_lambda_vacc


  FUNCTION phi(h,ul,hl,ur,hr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: h, ul, hl, ur, hr
    REAL(KIND=8)             :: vv, fl, fr, ovh, ovhl, ovhr
    IF (h>hl) THEN
       ovh = (2*h/(h**2+max(h,inputs%htiny)**2))
       ovhl = (2*hl/(hl**2+max(hl,inputs%htiny)**2))
       !fl = (h-hl)*SQRT((inputs%gravity/2)*(h+hl)/(h*hl))
       fl = (h-hl)*SQRT((inputs%gravity/2)*(ovh+ovhl))
    ELSE
       fl = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hl))
    END IF
    IF (h>hr) THEN
       ovh = (2*h/(h**2+max(h,inputs%htiny)**2))
       ovhr = (2*hr/(hr**2+max(hr,inputs%htiny)**2))
       !fr = (h-hr)*SQRT((inputs%gravity/2)*(h+hr)/(h*hr))
       fr = (h-hr)*SQRT((inputs%gravity/2)*(ovh+ovhr))
    ELSE
       fr = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hr))
    END IF
    vv = fl + fr + ur - ul
  END FUNCTION phi

END MODULE boundary_conditions
