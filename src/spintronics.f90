MODULE spintronics
    IMPLICIT NONE

    ! Constants
    REAL, PARAMETER :: zpi = 4*ATAN(1.0)

    !
    ! Global state variables
    !

    ! Spin & position vectos
    REAL, ALLOCATABLE, DIMENSION(:)  :: sx, sy, sz
    REAL, ALLOCATABLE, DIMENSION(:)  :: rx, ry, rz
    
    ! Exchange interaction
    REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: V_exc
    INTEGER, ALLOCATABLE, DIMENSION(:, :)   :: V_interacts_exc
    INTEGER, ALLOCATABLE, DIMENSION(:)      :: V_interacts_exc_count

    ! Dipolar interaction
    REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: V_dip
    INTEGER, ALLOCATABLE, DIMENSION(:, :)   :: V_interacts_dip
    INTEGER, ALLOCATABLE, DIMENSION(:)      :: V_interacts_dip_count

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
        IF(verbose) PRINT *, "Position vectors have correct size: OK"

        check_parameters = .TRUE.

    END FUNCTION

    !>  Expensivly calculates the total system's energy
    SUBROUTINE total_energy(E)
        IMPLICIT NONE 

        REAL, INTENT(OUT) :: E

        INTEGER :: n, i, j
        INTEGER :: neib 
        REAL :: dx, dy, dr
        REAL :: E_exc, E_dip

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
                dr = SQRT(dx*dx + dy*dy)

                E_dip = E_dip + V_dip(1, j, i) * &
                    (3 * (sx(i) * dx) * (sx(neib) * dx)/dr**5 - sx(i) * sx(neib)/dr**3)
                E_dip = E_dip + V_dip(2, j, i) * &
                    (3 * (sy(i) * dy) * (sy(neib) * dy)/dr**5 - sy(i) * sy(neib)/dr**3)
                E_dip = E_dip + V_dip(3, j, i) * &
                    (- sz(i) * sz(neib)/dr**3)
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

        REAL, DIMENSION(3) :: r
        REAL :: theta, sxy, s 

        CALL RANDOM_NUMBER(r)
        
        ! Set sz
        sz(i) = 2*r(1) - 1

        sxy = SQRT(1 - sz(i)**2)
        theta = r(2) * zpi

        ! Set sx and sy
        sx(i) = sxy * COS(theta)
        sy(i) = sxy * SIN(theta)

    END SUBROUTINE


    !> Perform n steps of metropolis simmulation on the current configuration
    SUBROUTINE metropolis(steps, beta, mean_energy)
        IMPLICIT NONE   

        !> Number of steops
        INTEGER, INTENT(IN) :: steps
        REAL, INTENT(IN) :: beta
        REAL, INTENT(OUT) :: mean_energy

        ! Local variables
        INTEGER :: step, i, j, k, neib
        INTEGER :: n

        ! Iteration state variables
        REAL :: current_energy, delta_energy, r
        REAL :: sx_hold, sy_hold, sz_hold
        REAL :: sx_delta, sy_delta, sz_delta
        REAL :: dx, dy, dr

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
                    dr = SQRT(dx*dx + dy*dy)

                    delta_energy = delta_energy + V_dip(1, j, i) * &
                        (3 * (sx_delta * dx) * (sx(neib) * dx)/dr**5 - sx_delta * sx(neib)/dr**3)
                    delta_energy = delta_energy + V_dip(2, j, i) * &
                        (3 * (sy_delta * dy) * (sy(neib) * dy)/dr**5 - sy_delta * sy(neib)/dr**3)
                    delta_energy = delta_energy + V_dip(3, j, i) * &
                        (- sz_delta * sz(neib)/dr**3)
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

            mean_energy = mean_energy + current_energy/REAL(steps)

        END DO
	END SUBROUTINE

END MODULE spintronics