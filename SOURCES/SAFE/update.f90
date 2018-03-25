MODULE update
  USE matrix_type
  USE space_dim
  USE input_data
  PUBLIC:: construct_matrices, euler, compute_dij
  TYPE(matrice_bloc), PUBLIC                  :: dij
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: lumped
  INTEGER, DIMENSION(:), POINTER, PUBLIC      :: diag
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: fix_roundoff
  TYPE(matrice_bloc), PUBLIC                  :: mass, pre_mass, mc_minus_ml
  PRIVATE
  TYPE(matrice_bloc), DIMENSION(k_dim):: cij
  TYPE(matrice_bloc)                  :: muij, resij, dijH
  INTEGER                             :: isolve, isolve_m
  LOGICAL                             :: if_roundoff_fix=.TRUE.
CONTAINS 

  SUBROUTINE construct_matrices
    USE st_matrix
    USE mesh_handling
    USE fem_s_M
    IMPLICIT NONE
    INTEGER :: m, p, ni, nj, i, j, d
    REAL(KIND=8), DIMENSION(k_dim) :: x
    !===mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)

    !===fix_roundoff
    ALLOCATE(fix_roundoff(mesh%np))

    !===lumped
    ALLOCATE(lumped(mesh%np))
    DO i = 1, mesh%np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO

    !===diag
    ALLOCATE(diag(mesh%np))
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          IF (i==mass%ja(p)) THEN
             diag(i) = p
             EXIT
          END IF
       END DO
    END DO

    !===Mass - lumped
    CALL st_csr(mesh%jj, mc_minus_ml%ia, mc_minus_ml%ja)
    ALLOCATE(mc_minus_ml%aa(SIZE(mc_minus_ml%ja)))
    mc_minus_ml%aa = mass%aa
    mc_minus_ml%aa(diag) = mc_minus_ml%aa(diag) - lumped

    !===Pre mass
    CALL st_csr(mesh%jj, pre_mass%ia, pre_mass%ja)
    ALLOCATE(pre_mass%aa(SIZE(pre_mass%ja)))
    DO i = 1, mesh%np
       pre_mass%aa(pre_mass%ia(i):pre_mass%ia(i+1)-1) = mass%aa(pre_mass%ia(i):pre_mass%ia(i+1)-1)/lumped(i)
    END DO

    !===dij
    CALL st_csr(mesh%jj, dij%ia, dij%ja)
    ALLOCATE(dij%aa(SIZE(dij%ja)))
    dij%aa = 0.d0

    !===muij
    CALL st_csr(mesh%jj, muij%ia, muij%ja)
    ALLOCATE(muij%aa(SIZE(muij%ja)))
    muij%aa = 0.d0

    !===dijH
    CALL st_csr(mesh%jj, dijH%ia, dijH%ja)
    ALLOCATE(dijH%aa(SIZE(dijH%ja)))
    dijH%aa = 0.d0

    !===cij = \int_K \GRAD(\phi_j) \phi_i \dif x
    DO d = 1, k_dim
       CALL st_csr(mesh%jj, cij(d)%ia, cij(d)%ja)
       ALLOCATE(cij(d)%aa(SIZE(cij(d)%ja)))
       cij(d)%aa = 0.d0
    END DO
    DO m = 1, mesh%me
       DO ni = 1, mesh%gauss%n_w  
          i = mesh%jj(ni, m)
          DO nj = 1, mesh%gauss%n_w  
             j = mesh%jj(nj, m)
             DO d = 1, k_dim
                x(d) = SUM(mesh%gauss%dw(d,nj,:,m) * mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
             END DO
             DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                IF (cij(1)%ja(p) == j) THEN  
                   DO d = 1, k_dim
                      cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                   END DO
                   EXIT
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !===entropy viscosity matrix
    CALL st_csr(mesh%jj, resij%ia, resij%ja)
    ALLOCATE(resij%aa(SIZE(resij%ja)))
    resij%aa = 0.d0

  END SUBROUTINE construct_matrices

  SUBROUTINE euler(un,unext)
    USE mesh_handling
    USE boundary_conditions
    !USE solve_sp
    USE pardiso_solve
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: rk
    REAL(KIND=8), DIMENSION(mesh%np)  :: ff
    REAL(KIND=8) :: xx, hmin
    INTEGER :: p, i, j, k, d
    LOGICAL, SAVE :: once=.TRUE.
    IF (once) THEN
       isolve=-1
       isolve_m=-1
       once=.FALSE.
    END IF

    !===Galerkin
    IF (inputs%viscosity_type=='galerkin') THEN
       dij%aa = 0.d0
       muij%aa = 0.d0
       CALL smb_2(un,rk)
       CALL divide_by_lumped(rk)
       IF (inputs%if_lumped) THEN
          unext = un+inputs%dt*rk
          !unext(1,:)  = max(un(1,:) + inputs%dt*rk(1,:),0.d0)
          !unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
       ELSE
          DO k = 1, inputs%syst_size
             CALL solve_pardiso(pre_mass%aa,pre_mass%ia,pre_mass%ja,rk(k,:),ff,isolve,2)
             isolve=ABS(isolve)
             !unext(k,:) = un(k,:)+inputs%dt*ff
             IF (k==1) THEN
                unext(k,:)  = max(un(k,:) + inputs%dt*ff,0.d0)
             ELSE
                unext(k,:) = un(k,:)+inputs%dt*ff
             END IF
          END DO
       END IF
       RETURN
    END IF

    !===Viscous solution
    CALL compute_dij(un)
    CALL compute_muij(un)
    IF (inputs%viscosity_type=='entropy_visc') THEN
       CALL entropy_residual(un)
       DO i = 1, mesh%np
          DO p = dij%ia(i), dij%ia(i+1) - 1
             j = dij%ja(p)
             dij%aa(p) = MIN(dij%aa(p),  1.d0*resij%aa(p))
             muij%aa(p)= MIN(muij%aa(p), 1.d0*resij%aa(p))
          END DO
       END DO
    END IF

    !===Compute right-hand side
    IF (inputs%viscous_type=='type1') THEN
       CALL smb_1(un,rk)
    ELSE IF (inputs%viscous_type=='type2') THEN
       dij%aa = MAX(dij%aa,muij%aa)
       IF (if_roundoff_fix) THEN
          CALL smb_2_roundoff(un,rk)
       ELSE
          CALL smb_2(un,rk)
       END IF
    ELSE
       WRITE(*,*) ' BUG in Euler, viscous_type'
       STOP
    END IF
    CALL divide_by_lumped(rk)

    !===Solve and update
    IF (inputs%if_lumped) THEN
       IF (inputs%viscous_type=='type1')  THEN
          unext = un+inputs%dt*rk
       ELSE IF (inputs%viscous_type=='type2')  THEN
          IF (if_roundoff_fix) THEN
             !===Trick with the roundoff works only with CALL smb_2_roundoff(un,rk)
             unext(1,:)  = un(1,:)*(1+inputs%dt*fix_roundoff/lumped) + inputs%dt*rk(1,:)
             unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
          ELSE
             !unext = un+inputs%dt*rk
             unext(1,:)  = max(un(1,:) + inputs%dt*rk(1,:),0.d0)
             unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
          END IF
       ELSE
          WRITE(*,*) ' BUG in euler, wrong inputs%viscous_type'
          STOP
       END IF
    ELSE
       DO k = 1, inputs%syst_size
          CALL solve_pardiso(pre_mass%aa,pre_mass%ia,pre_mass%ja,rk(k,:),ff,isolve,2)
          isolve=ABS(isolve)
          unext(k,:) = un(k,:)+inputs%dt*ff
       END DO
    END IF

    CALL check_Hmin(unext)
    !===Exit here if not entropy viscosity
    RETURN

  END SUBROUTINE euler

  SUBROUTINE compute_dij(un)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: norm_cij, lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8), DIMENSION(inputs%syst_size)     :: ur, ul
 
    !===Viscosity using compute_lambda
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p) !=== definition of cij is same as in the paper
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             ul=un(:,i)
             ur=un(:,j)
             !CALL compute_lambda(ul,ur,nij,lambda)
             CALL compute_lambda_vacc(ul,ur,nij,lambda)
             dij%aa(p) = norm_cij*lambda
          ELSE
             dij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO
 
    RETURN

  END SUBROUTINE compute_dij

  SUBROUTINE compute_muij(un)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: norm_cij, lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8), DIMENSION(k_dim)                :: ur, ul
 
    !===Viscosity using speed only
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p) !=== definition of cij is the same as in the paper
             END DO
             ul = un(2:,i)*(2*un(1,i)/(un(1,i)**2+max(un(1,i),inputs%htiny)**2))
             ur = un(2:,j)*(2*un(1,j)/(un(1,j)**2+max(un(1,j),inputs%htiny)**2))
             !ul=un(2:,i)/max(un(1,i),inputs%htiny)
             !ur=un(2:,j)/max(un(1,j),inputs%htiny)
             lambda=MAX(MAX(-SUM(nij*ul),0.d0),MAX(SUM(nij*ur),0.d0))
             !lambda=MAX(ABS(SUM(nij*ul)),ABS(SUM(nij*ur)))
             muij%aa(p) = lambda
          ELSE
             muij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    DO i = 1, mesh%np
       muij%aa(diag(i)) = -SUM(muij%aa(dij%ia(i):muij%ia(i+1)-1))
    END DO

    RETURN

  END SUBROUTINE compute_muij

  SUBROUTINE smb_1(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, Hstarij, Hstarji, ratij, ratji

    IF (inputs%if_alpha_limit) THEN
       CALL alpha_limit(un(1,:))
    END IF

    vv=flux(un)

    rk=0.d0
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij/un(1,i)
          ratji = Hstarji/un(1,j)
          IF (un(1,i).LE.inputs%htiny) THEN !Important for long-time WB
             ratij=0.d0
          ELSE 
             ratij = Hstarij/un(1,i)
          END IF
          IF (un(1,j).LE.inputs%htiny) THEN !Important for long-time WB
             ratji=0.d0
          ELSE
             ratji = Hstarji/un(1,j)
          END IF
          
          DO k = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j)*ratji + vv(k,d,i)*ratij)
             END DO
             rk(k,i) = rk(k,i) + xx + dij%aa(p)*(un(k,j)*ratji-un(k,i)*ratij)
          END DO
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - 0.5d0*inputs%gravity*(Hstarji**2 - Hstarij**2)*cij(k)%aa(p)
          END DO
       END DO
    END DO

    IF (inputs%type_test==9) THEN
       CALL friction(un,rk)
    END IF
  END SUBROUTINE smb_1

  SUBROUTINE friction(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    REAL(KIND=8), DIMENSION(mesh%np)  :: hloc_star, hstar, vel
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: fric
    REAL(KIND=8), PARAMETER :: alpha = 2.d0, eta=1.d0/3.d0, chi=1.d0

!!$    fric = inputs%gravity*inputs%mannings**2
!!$    hloc_star = (fric)**(1.d0-alpha)*(lumped)**(alpha*eta/2.d0)
!!$    hstar = un(1,:)**eta
!!$    vel = fric*SQRT(velocity(1,:)**2+velocity(2,:)**2)
!!$    DO i = 1, mesh%np
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          DO k = 2, inputs%syst_size
!!$             rk(k,i) = rk(k,i) &
!!$                  - mass%aa(p)*velocity(k-1,j)*vel(j)/(hstar(j) + hloc_star(j))
!!$          END DO
!!$       END DO
!!$    END DO

    fric = inputs%gravity*inputs%mannings**2
    hstar = MAX(un(1,:),inputs%htiny)**(1.d0+eta)
    !hstar = un(1,:)**(1.d0+eta)
    vel = fric*SQRT(velocity(1,:)**2+velocity(2,:)**2)
    hloc_star = chi*vel*inputs%dt !fric*SQRT(lumped) !t !
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 2, inputs%syst_size
             rk(k,i) = rk(k,i) &
                  - mass%aa(p)*un(k,j)*vel(j)/(hstar(j) + hloc_star(j))
          END DO
       END DO
    END DO


  END SUBROUTINE friction

  SUBROUTINE smb_2_roundoff(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, hmean, Hstarij, Hstarji, ratij, ratji

    IF (inputs%if_alpha_limit) THEN
       CALL alpha_limit(un(1,:))
    END IF

    vv=flux(un)
    rk=0.d0
    fix_roundoff = 0.d0
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          !IF (un(1,i).LE.0.d0) THEN 
          IF (un(1,i).LE.inputs%htiny) THEN !Important for long-time WB
             ratij=0.d0
          ELSE 
             ratij = Hstarij/un(1,i)
          END IF
          !IF (un(1,j).LE.0.d0) THEN
          IF (un(1,j).LE.inputs%htiny) THEN !Important for long-time WB
             ratji=0.d0
          ELSE
             ratji = Hstarji/un(1,j)
          END IF
          !===mass
          k=1
          xx = 0.d0
          DO d = 1, k_dim
             xx = xx - cij(d)%aa(p)*(velocity(d,j))
          END DO
          IF (i==j) THEN
             fix_roundoff(i) = fix_roundoff(i) + xx
          ELSE
             !===Fix roundoff error
             fix_roundoff(i) = fix_roundoff(i) - muij%aa(p) + (dij%aa(p)-muij%aa(p))*(-ratij)
             !fix_roundoff(i) = fix_roundoff(i) - (muij%aa(p) + (dij%aa(p)-muij%aa(p))*(ratij))
             rk(k,i) = rk(k,i) + (xx + muij%aa(p))*un(k,j) + max((dij%aa(p)-muij%aa(p))*Hstarji,0.d0)
          END IF
          !===rest
          DO k = 2, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j))
             END DO
             rk(k,i) = rk(k,i) + xx + dij%aa(p)*(un(k,j)*ratji-un(k,i)*ratij) &
                  + muij%aa(p)*(un(k,j)*(1.d0-ratji)-un(k,i)*(1.d0-ratij))
          END DO
          !===Topography contribution 
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - inputs%gravity*un(1,i)*(un(1,j)+bath(j)-un(1,i)-bath(i))*cij(k)%aa(p)
          END DO
       END DO
    END DO

    IF (inputs%type_test==9) THEN
       CALL friction(un,rk)
    END IF
  END SUBROUTINE smb_2_roundoff

  SUBROUTINE smb_2(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, hmean, Hstarij, Hstarji, ratij, ratji

    IF (inputs%if_alpha_limit) THEN
       CALL alpha_limit(un(1,:))
    END IF

    vv=flux(un)
    rk=0.d0
    fix_roundoff = 0.d0
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          IF (un(1,i).LE.0.d0) THEN
             ratij=0.d0
          ELSE 
             ratij = Hstarij/un(1,i) !Do not use one_over_h here
             !ratij = Hstarij*one_over_h(i) !Do not use one_over_h here
          END IF
          IF (un(1,j).LE.0.d0) THEN 
             ratji=0.d0
          ELSE
             ratji = Hstarji/un(1,j) !Do not use one_over_h here
             !ratji = Hstarji*one_over_h(j)
          END IF

!!$          k=1
!!$          xx = 0.d0
!!$          DO d = 1, k_dim
!!$             xx = xx - cij(d)%aa(p)*velocity(d,j)
!!$             !xx = xx - cij(d)%aa(p)*vv(k,d,j)
!!$          END DO
!!$          IF (i==j) THEN
!!$             fix_roundoff(i) = fix_roundoff(i) + xx 
!!$          ELSE
!!$             xx = xx*un(k,j)
!!$             fix_roundoff(i) = fix_roundoff(i) - muij%aa(p) - (dij%aa(p)-muij%aa(p))*ratij
!!$             !fix_roundoff(i) = fix_roundoff(i) - (muij%aa(p) + (dij%aa(p)-muij%aa(p))*ratij)*un(k,i)
!!$             rk(k,i) = rk(k,i) + xx + (muij%aa(p) + (dij%aa(p)-muij%aa(p))*ratji)*un(k,j)
!!$          END IF

          DO k = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j))
             END DO
             rk(k,i) = rk(k,i) + xx + dij%aa(p)*(un(k,j)*ratji-un(k,i)*ratij) &
                  + muij%aa(p)*(un(k,j)*(1.d0-ratji)-un(k,i)*(1.d0-ratij))
          END DO
          !===Topography contribution 
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - inputs%gravity*un(1,i)*(un(1,j)+bath(j)-un(1,i)-bath(i))*cij(k)%aa(p)
          END DO

       END DO
    END DO

    IF (inputs%type_test==9) THEN
       CALL friction(un,rk)
    END IF

  END SUBROUTINE smb_2

  SUBROUTINE alpha_limit(un)
    USE mesh_handling
    IMPLICIT NONE
    INTEGER :: i, j, p
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: alpha
    REAL(KIND=8) :: nump, numm, denom, num
    INTEGER :: exponent=2
    DO i = 1, mesh%np
       !===Keep full viscosity in dry regions
       IF (un(i).LE. inputs%htiny) THEN
          alpha(i) = 1.d0
          CYCLE
       END IF

       nump = 0.d0
       numm = 0.d0
       denom =0.d0
       num  = 0.d0
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          num = num + un(j) - un(i)
          denom = denom + ABS(un(j) - un(i))
       END DO
       IF (denom.LE.ABS(num)) THEN
          alpha(i) = 1.d0
          num = denom
       ELSE
          alpha(i) = ABS(num)/(denom)
       END IF
       alpha(i) = alpha(i)**exponent
    END DO
    !===Limit viscosity
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             dij%aa(p) = dij%aa(p)*MAX(alpha(i),alpha(j))
             muij%aa(p) = muij%aa(p)*MAX(alpha(i),alpha(j))
          ELSE
             dij%aa(p) = 0.d0
          END IF
       END DO
    END DO

    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
       muij%aa(diag(i)) = -SUM(muij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO

  END SUBROUTINE alpha_limit

  SUBROUTINE divide_by_lumped(rk)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rk
    INTEGER :: k
    DO k = 1, inputs%syst_size
       rk(k,:) = rk(k,:)/lumped
    END DO
  END SUBROUTINE divide_by_lumped

  SUBROUTINE check_hmin(h)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: h
    !WRITE(*,*) minval(h)
    IF (MINVAL(h(1,:))<0.d0) THEN
       WRITE(*,*) 'Min h<0, STOP', MINVAL(h)
       WRITE(*,*) 'MAXVAL(vel)', MAXVAL(ABS(velocity(1,:))), MAXVAL(ABS(velocity(2,:)))
       STOP
    END IF
  END SUBROUTINE check_hmin

  SUBROUTINE entropy_residual(un)
    USE mesh_handling
    USE pardiso_solve
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: rk, unext, Entprime
    REAL(KIND=8), DIMENSION(k_dim,mesh%np)  :: velnext
    REAL(KIND=8), DIMENSION(mesh%np)        :: res, Ent, Entnext, maxn, minn
    REAL(KIND=8), DIMENSION(mesh%np)  :: ff, rescale

    INTEGER :: k, i, j, p

    CALL smb_2(un,rk)
    CALL divide_by_lumped(rk)
    unext = un+inputs%dt*rk
    DO k = 1, k_dim
       velnext(k,:) = unext(k+1,:)*compute_one_over_h(unext(1,:))
    END DO
    Ent     = inputs%gravity*un(1,:)**2
    Entnext = inputs%gravity*unext(1,:)**2
    DO k = 1, k_dim
       Ent     = Ent     + velocity(k,:)*un(k+1,:)
       Entnext = Entnext + velnext(k,:)*unext(k+1,:)
    END DO
    Ent     = 0.5d0*Ent
    Entnext = 0.5d0*Entnext

    Entprime(1,:) = inputs%gravity*un(1,:)  ! -|u|^2/2+gh
    DO k = 1, k_dim 
       Entprime(1,:) = Entprime(1,:) -0.5d0*velocity(k,:)**2
       Entprime(k+1,:) = velocity(k,:)
    END DO

    res = lumped*(Entnext- Ent)/inputs%dt
    DO k = 1, inputs%syst_size
       res = res -lumped*rk(k,:)*Entprime(k,:)
    END DO

    CALL maxmin(ent,dij,maxn,minn)
    DO i = 1, mesh%np
       rescale = ABS(maxn-minn)/2 +1.d-10*inputs%gravity*max_water_h**2
       !IF (un(1,i).LE.100*inputs%htiny) then
       !   write(*,*) res(i), rescale(i), maxn(i), minn(i)
       !   end if
    END DO
    res = ABS(res)/rescale
    
    resij%aa = 0.d0
    DO i = 1, mesh%np
       DO p = resij%ia(i), resij%ia(i+1) - 1
          j = resij%ja(p)
          resij%aa(p) = max(res(i),res(j))
          !resij%aa(p) = max(res(i),res(j))/ABS(Ent(j)-Ent(i))
       END DO
    END DO
  END SUBROUTINE entropy_residual

  SUBROUTINE maxmin(un,mat,maxn,minn)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: maxn, minn
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
    INTEGER      :: i
    DO i = 1, SIZE(un)
       maxn(i) = MAXVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
       minn(i) = MINVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
    END DO
  END SUBROUTINE maxmin
END MODULE update
