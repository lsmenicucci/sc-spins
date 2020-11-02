MODULE spintronics
    IMPLICIT NONE

    ! Constants
    REAL, PARAMETER :: zpi = 8*ATAN(1.0)

    !
    ! Global state variables
    !

    ! Spin & position vectos
    REAL(8), ALLOCATABLE, DIMENSION(:)  :: sx, sy, sz
    REAL(8), ALLOCATABLE, DIMENSION(:)  :: rx, ry, rz
  
    ! Spin data
    REAL(8), ALLOCATABLE, DIMENSION(:)  :: spin 


    ! Exchange interaction
    REAL(8), ALLOCATABLE, DIMENSION(:, :, :)    :: V_exc
    INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: V_interacts_exc
    INTEGER, ALLOCATABLE, DIMENSION(:)          :: V_interacts_exc_count

    ! Dipolar interaction
    REAL(8), ALLOCATABLE, DIMENSION(:, :, :)    :: V_dip
    INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: V_interacts_dip
    INTEGER, ALLOCATABLE, DIMENSION(:)          :: V_interacts_dip_count

    ! Log
    REAL(8), ALLOCATABLE, DIMENSION(:)  :: energy_log

    CONTAINS

    !> Check the global parameters
    LOGICAL FUNCTION check_parameters(verbose)
        IMPLICIT NONE
            
        LOGICAL :: verbose

        INTEGER :: n
        INTEGER, DIMENSION(3) :: v_shape

        check_parameters = .FALSE.

        ! Check spin vectors
        IF(.NOT. ALLOCATED(sx)) RETURN
        IF(.NOT. ALLOCATED(sy)) RETURN
        IF(.NOT. ALLOCATED(sz)) RETURN
        IF(verbose) PRINT *, "Spin vectors allocated: OK"

        ! Check spin vectors
        IF(.NOT. ALLOCATED(rx)) RETURN
        IF(.NOT. ALLOCATED(ry)) RETURN
        IF(.NOT. ALLOCATED(rz)) RETURN
        IF(.NOT. ALLOCATED(spin)) RETURN
        IF(verbose) PRINT *, "Position vectors allocated: OK"
        
        ! Check potential related variables
        IF(.NOT. ALLOCATED(V_interacts_exc)) RETURN
        IF(.NOT. ALLOCATED(V_exc)) RETURN
        IF(verbose) PRINT *, "Exchange interaction vectors allocated: OK"

        IF(.NOT. ALLOCATED(V_interacts_dip_count)) RETURN
        IF(.NOT. ALLOCATED(V_interacts_dip)) RETURN
        IF(.NOT. ALLOCATED(V_dip)) RETURN
        IF(verbose) PRINT *, "Dipolar interaction vectors allocated: OK"

        ! Check size
        n = SIZE(sx)
        IF(SIZE(sy) /= n) RETURN   
        IF(SIZE(sz) /= n) RETURN
        IF(verbose) PRINT *, "Spin vectors have correct size: OK"

        IF(SIZE(V_exc) < 3*n) RETURN   
        IF(SIZE(V_interacts_exc)*3 /= SIZE(V_exc)) RETURN   
        IF(SIZE(V_interacts_exc_count) /= n) RETURN
        IF(verbose) PRINT *, "Exchange interaction vectors have correct size: OK"

        IF(SIZE(V_dip) < 3*n) RETURN   
        IF(SIZE(V_interacts_dip)*3 /= SIZE(V_dip)) RETURN   
        IF(SIZE(V_interacts_dip_count) /= n) RETURN
        IF(verbose) PRINT *, "Dipolar interaction vectors have correct size: OK"
        
        ! Check shape
        v_shape = SHAPE(V_exc)
        IF(v_shape(1) /= 3) RETURN

        v_shape = SHAPE(V_dip)
        IF(v_shape(1) /= 3) RETURN

        ! Check position vectos
        IF(SIZE(rx) /= n) RETURN
        IF(SIZE(ry) /= n) RETURN
        IF(SIZE(rz) /= n) RETURN
        IF(verbose) PRINT *, "Position vectors have correct size: OK"

        check_parameters = .TRUE.

    END FUNCTION

    !>  Expensivly calculates the total system's energy
    SUBROUTINE total_energy(E)
        IMPLICIT NONE 

        REAL(8), INTENT(OUT) :: E

        INTEGER :: n, i, j
        INTEGER :: neib 
        REAL(8) :: dx, dy, dz, dr
        REAL(8) :: E_exc, E_dip

        n = SIZE(sx)

        ! Initialize energy
        E_exc = 0
        E_dip = 0

        DO i = 1, n
            ! Calculate exchange energy
            DO j = 1, V_interacts_exc_count(i)
                neib = V_interacts_exc(j, i)

                ! Compute this term energy
                E_exc = E_exc + V_exc(1, j, i) * sx(i) * sx(neib)
                E_exc = E_exc + V_exc(2, j, i) * sy(i) * sy(neib)
                E_exc = E_exc + V_exc(3, j, i) * sz(i) * sz(neib)
            END DO

            ! Calculate dipolar energy
            DO j = 1, V_interacts_dip_count(i)
                neib = V_interacts_dip(j, i)

                dx = rx(i) - rx(neib)
                dy = ry(i) - ry(neib)
                dz = rz(i) - rz(neib)
                dr = SQRT(dx*dx + dy*dy + dz*dz)

                E_dip = E_dip + V_dip(1, j, i) * &
                    (3 * (sx(i) * dx) * (sx(neib) * dx)/dr**5 - sx(i) * sx(neib)/dr**3)
                E_dip = E_dip + V_dip(2, j, i) * &
                    (3 * (sy(i) * dy) * (sy(neib) * dy)/dr**5 - sy(i) * sy(neib)/dr**3)
                E_dip = E_dip + V_dip(3, j, i) * &
                    (3 * (sz(i) * dz) * (sz(neib) * dz)/dr**5 - sz(i) * sz(neib)/dr**3)
            END DO

        END DO

        E_exc = E_exc/2
        E_dip = E_dip/2

        ! Return the final value
        E = E_exc + E_dip

    END SUBROUTINE

    !> Generates a new spin on the site i
    SUBROUTINE new_spin(i)
        IMPLICIT NONE

        ! Spin site index
        INTEGER :: i

        REAL(8), DIMENSION(3) :: r
        REAL(8) :: theta, sxy, s 

        CALL RANDOM_NUMBER(r)

        ! Generate a random spin size
        s = (spin(i) - FLOOR(spin(i))) + FLOOR(spin(i)*r(1))

        ! Set sz
        sz(i) = s*(2*r(2) - 1)

        sxy = s * SQRT(1 - sz(i)**2)
        theta = r(3) * zpi

        ! Set sx and sy
        sx(i) = sxy * COS(theta)
        sy(i) = sxy * SIN(theta)

    END SUBROUTINE


    !> Perform n steps of metropolis simmulation on the current configuration
    SUBROUTINE metropolis(steps, beta, mean_energy, mean_mag_x, mean_mag_y, mean_mag_z)
        !f2py threadsafe
        IMPLICIT NONE   

        !> Number of steops
        INTEGER, INTENT(IN) :: steps
        REAL(8), INTENT(IN) :: beta
        REAL(8), INTENT(OUT) :: mean_energy
        REAL(8), INTENT(OUT) :: mean_mag_x
        REAL(8), INTENT(OUT) :: mean_mag_y
        REAL(8), INTENT(OUT) :: mean_mag_z


        ! Local variables
        INTEGER :: step, i, j, k, neib
        INTEGER :: n

        ! Iteration state variables
        REAL(8) :: current_energy, delta_energy, r
        REAL(8) :: sx_hold, sy_hold, sz_hold
        REAL(8) :: sx_delta, sy_delta, sz_delta
        REAL(8) :: dx, dy, dz, dr

        ! Check parameters
        IF(.NOT. check_parameters(.FALSE.)) THEN 
            PRINT *, "Current set of parameters is not valid"
            RETURN
        END IF

        ! Initialize array sizes
        n = SIZE(sx)

        ! Initialize energy
        CALL total_energy(current_energy)

        mean_energy = 0
        mean_mag_x = 0
        mean_mag_y = 0
        mean_mag_z = 0

        ! Perform simulation
        DO step = 1, steps
            ! Flip all the lattice
            DO i = 1, n
                ! Hold the current spin value
                sx_hold = sx(i)
                sy_hold = sy(i)
                sz_hold = sz(i)

                ! Generate a new spin
                CALL new_spin(i)

                ! Calculate the spin delta
                sx_delta = sx(i) - sx_hold
                sy_delta = sy(i) - sy_hold
                sz_delta = sz(i) - sz_hold

                ! Calculate the energy shift
                delta_energy = 0

                ! Calculate exchange energy shift
                DO j = 1,  V_interacts_exc_count(i)
                    neib = V_interacts_exc(j, i)

                    delta_energy = delta_energy + V_exc(1, j, i) * sx(neib) * sx_delta
                    delta_energy = delta_energy + V_exc(2, j, i) * sy(neib) * sy_delta
                    delta_energy = delta_energy + V_exc(3, j, i) * sz(neib) * sz_delta
                END DO
                
                ! Calculate dipolar energy shift
                DO j = 1, V_interacts_dip_count(i)
                    neib = V_interacts_dip(j, i)

                    dx = rx(i) - rx(neib)
                    dy = ry(i) - ry(neib)
                    dz = rz(i) - rz(neib)
                    dr = SQRT(dx*dx + dy*dy + dz*dz)

                    delta_energy = delta_energy + V_dip(1, j, i) * &
                        (3 * (sx_delta * dx) * (sx(neib) * dx)/dr**5 - sx_delta * sx(neib)/dr**3)
                    delta_energy = delta_energy + V_dip(2, j, i) * &
                        (3 * (sy_delta * dy) * (sy(neib) * dy)/dr**5 - sy_delta * sy(neib)/dr**3)
                    delta_energy = delta_energy + V_dip(3, j, i) * &
                        (3 * (sz_delta * dz) * (sz(neib) * dz)/dr**5 - sz_delta * sz(neib)/dr**3)
                END DO

                ! Update energy
                current_energy = current_energy + delta_energy

                ! Metropolis test
                IF (delta_energy > 0) THEN
                    CALL RANDOM_NUMBER(r)

                    ! Tests for rejection
                    IF (r > EXP(-beta * delta_energy)) THEN
                        sx(i) = sx_hold
                        sy(i) = sy_hold
                        sz(i) = sz_hold

                        ! Restore the energy
                        current_energy = current_energy - delta_energy
                    END IF
                END IF
            END DO

            mean_energy = mean_energy + current_energy/DBLE(steps)
            mean_mag_x = mean_mag_x + SUM(sx)/DBLE(steps)
            mean_mag_y = mean_mag_y + SUM(sy)/DBLE(steps)
            mean_mag_z = mean_mag_z + SUM(sz)/DBLE(steps)

        END DO
	END SUBROUTINE

    !> Calculates the spin commutator with H for a given configuration 
    SUBROUTINE spin_commutator(n, spin_x, spin_y, spin_z, commu_x, commu_y, commu_z)
        IMPLICIT NONE 

        INTEGER, INTENT(IN) :: n
        REAL(8), DIMENSION(n), INTENT(IN) :: spin_x
        REAL(8), DIMENSION(n), INTENT(IN) :: spin_y
        REAL(8), DIMENSION(n), INTENT(IN) :: spin_z
        REAL(8), DIMENSION(n), INTENT(OUT) :: commu_x
        REAL(8), DIMENSION(n), INTENT(OUT) :: commu_y
        REAL(8), DIMENSION(n), INTENT(OUT) :: commu_z
        
        ! Internal variables
        REAL(8) :: cx, cy, cz, rcx, rcy, rcz
        REAL(8) :: dx, dy, dz, dr
        INTEGER :: i, j, neib

        ! Initialize accumulator variable
        commu_x = 0
        commu_y = 0
        commu_z = 0

        DO i = 1, n
            ! Sum exchange contribution
            DO j = 1, V_interacts_exc_count(i)
                neib = V_interacts_exc(j, i)

                cx = spin_y(neib) * spin_z(i) - spin_y(i) * spin_z(neib)
                cy = spin_z(neib) * spin_x(i) - spin_z(i) * spin_x(neib)
                cz = spin_x(neib) * spin_y(i) - spin_x(i) * spin_y(neib)

                commu_x(i) = commu_x(i) + V_exc(1, j, i) * cx
                commu_y(i) = commu_y(i) + V_exc(2, j, i) * cy
                commu_z(i) = commu_z(i) + V_exc(3, j, i) * cz
            END DO

            ! Sum dipolar contribution
            DO j = 1, V_interacts_dip_count(i)
                neib = V_interacts_dip(j, i)

                dx = rx(i) - rx(neib)
                dy = ry(i) - ry(neib)
                dz = rz(i) - rz(neib)
                dr = SQRT(dx*dx + dy*dy + dz*dz)

                ! Compute terms
                cx = spin_y(neib) * spin_z(i) - spin_y(i) * spin_z(neib)
                cy = spin_z(neib) * spin_x(i) - spin_z(i) * spin_x(neib)
                cz = spin_x(neib) * spin_y(i) - spin_x(i) * spin_y(neib)

                rcx = dy * spin_z(i) - spin_y(i) * dz
                rcy = dz * spin_x(i) - spin_z(i) * dx
                rcz = dx * spin_y(i) - spin_x(i) * dy

                commu_x(i) = commu_x(i) + V_dip(1, j, i) * &
                        (3 * rcx * (spin_x(neib) * dx)/dr**5 - cx/dr**3)

                commu_y(i) = commu_y(i) + V_dip(2, j, i) * &
                        (3 * rcy * (spin_y(neib) * dy)/dr**5 - cy/dr**3)

                commu_z(i) = commu_z(i) + V_dip(3, j, i) * &
                        (3 * rcz * (spin_z(neib) * dz)/dr**5 - cz/dr**3)
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE integrate(time, dt)
        !f2py threadsafe
        IMPLICIT NONE

        REAL(8), INTENT(IN) :: time
        REAL(8), INTENT(IN) :: dt

        ! Internal variables
        REAL(8), DIMENSION(:, :), ALLOCATABLE :: com_sx, com_sy, com_sz
        REAL(8), DIMENSION(:, :), ALLOCATABLE :: kx, ky, kz
        REAL(8), DIMENSION(:), ALLOCATABLE :: pre_sx, pre_sy, pre_sz
        INTEGER :: steps, t, n, i

        ! Allocate buffer vectors
        n = SIZE(sx)

        ! Functions f(y_{i}), f(y_{i - 1}), f(y_{i - 2}), f(y_{i - 3})
        ALLOCATE(com_sx(4, n))
        ALLOCATE(com_sy(4, n))
        ALLOCATE(com_sz(4, n))

        ! Predictors P(y_{i}), P(y_{i - 1}), P(y_{i - 2}), P(y_{i - 3})
        ALLOCATE(pre_sx(n))
        ALLOCATE(pre_sy(n))
        ALLOCATE(pre_sz(n))

        ! Functions f(y_{i}), f(y_{i - 1}), f(y_{i - 2}), f(y_{i - 3})
        ALLOCATE(kx(4, n))
        ALLOCATE(ky(4, n))
        ALLOCATE(kz(4, n))

        ! Fill initial values with Runge-Kutta
        DO t = 1, 4
            ! Calculate k1
            CALL spin_commutator(n, sx, sy, sz, &
                                    kx(1, :), ky(1, :), kz(1, :))
            ! Calculate k2
            CALL spin_commutator(n, sx + dt*kx(1, :)/2.0, sy + dt*ky(1, :)/2.0, sz + dt*kz(1, :)/2.0, &
                                    kx(2, :), ky(2, :), kz(2, :))
            ! Calculate k3
            CALL spin_commutator(n, sx + dt*kx(2, :)/2.0, sy + dt*ky(2, :)/2.0, sz + dt*kz(2, :)/2.0, &
                                    kx(3, :), ky(3, :), kz(3, :))
            ! Calculate k4
            CALL spin_commutator(n, sx + dt*kx(3, :), sy + dt*ky(3, :), sz + dt*kz(3, :), &
                                    kx(4, :), ky(4, :), kz(4, :))

            sx = sx + (kx(1, :) + kx(2, :) + kx(3, :) + kx(4, :))*dt/6.0
            sy = sy + (ky(1, :) + ky(2, :) + ky(3, :) + ky(4, :))*dt/6.0
            sz = sz + (kz(1, :) + kz(2, :) + kz(3, :) + kz(4, :))*dt/6.0

            ! PRINT *, ""
            ! PRINT *, "step ", t
            ! PRINT *, sx
            ! PRINT *, sy
            ! PRINT *, sz
            ! PRINT *, (kx(1, :) + kx(2, :) + kx(3, :) + kx(4, :))*dt/6.0
            ! PRINT *, (ky(1, :) + ky(2, :) + ky(3, :) + ky(4, :))*dt/6.0
            ! PRINT *, (kz(1, :) + kz(2, :) + kz(3, :) + kz(4, :))*dt/6.0
            
            com_sx(5 - i,:) = kx(1, :)
            com_sy(5 - i,:) = ky(1, :)
            com_sz(5 - i,:) = kz(1, :)
        END DO

        ! Clear Runge-Kutta variables
        DEALLOCATE(kx)
        DEALLOCATE(ky)
        DEALLOCATE(kz)

        ! Start integration
        steps = INT(time/dt)

        ! Allocate log
        IF(ALLOCATED(energy_log)) DEALLOCATE(energy_log)
        ALLOCATE(energy_log(steps))

        DO t = 1, steps 
            ! Calculate the predicted value
            pre_sx = sx + & 
                (dt/24.0) * ( 55.0 * com_sx(1, :) - 59.0 * com_sx(2, :) + &
                              37.0 * com_sx(3, :) - 9.0 * com_sx(4, :) )
            pre_sy = sy + & 
                (dt/24.0) * ( 55.0 * com_sy(1, :) - 59.0 * com_sy(2, :) + &
                              37.0 * com_sy(3, :) - 9.0 * com_sy(4, :) )
            pre_sz = sz + & 
                (dt/24.0) * ( 55.0 * com_sz(1, :) - 59.0 * com_sz(2, :) + &
                              37.0 * com_sz(3, :) - 9.0 * com_sz(4, :) )

            ! Overwrite the last commutator value with the predictor
            CALL spin_commutator(n, pre_sx, pre_sy, pre_sz, com_sx(4, :), com_sy(4, :), com_sz(4, :))

            ! Compute the next spins
            sx = sx + & 
                (dt/24.0) * ( 9.0 * com_sx(4, :) + 19.0 * com_sx(1, :) - &
                              5.0 * com_sx(2, :) - com_sx(3, :) )
            sy = sy + & 
                (dt/24.0) * ( 9.0 * com_sy(4, :) + 19.0 * com_sy(1, :) - &
                              5.0 * com_sy(2, :) - com_sy(3, :) )
            sz = sz + & 
                (dt/24.0) * ( 9.0 * com_sz(4, :) + 19.0 * com_sz(1, :) - &
                              5.0 * com_sz(2, :) - com_sz(3, :) )

            ! Reorder the buffer values
            com_sx(4, :) = com_sx(3, :)
            com_sy(4, :) = com_sy(3, :)
            com_sz(4, :) = com_sz(3, :)

            com_sx(3, :) = com_sx(2, :)
            com_sy(3, :) = com_sy(2, :)
            com_sz(3, :) = com_sz(2, :)

            com_sx(2, :) = com_sx(1, :)
            com_sy(2, :) = com_sy(1, :)
            com_sz(2, :) = com_sz(1, :) 

            ! Calculate the commutator for the current configuration
            CALL spin_commutator(n, sx, sy, sz, com_sx(1, :), com_sy(1, :), com_sz(1, :))

            ! Calculate the energy
            energy_log(i) = SUM(sx)/SIZE(sx)
        END DO

    END SUBROUTINE

END MODULE spintronics
