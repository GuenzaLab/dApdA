!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       PROGRAM STRUCTURE_CALCULATOR
       IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       (N1)              a <----- b                        5'
!                                           |                     ^
!                                           |                     |
!                                           |                     |
!       (N2)              d <----- c                          3'
!
! This is a schematic representation of the dinucleotide system in question. Here
! N1 refers to the first or "top" nucleotide and N2 the "bottom" or second. The
! directionality is established by the 5' to 3' being the z axis, going positive
! towards the 5' end. The vectors for the nucleotide go in the direction from the
! backbone to the nucleic base.
!
! Sites "a" and "d" represents the geometric center between the C6 and N3 atoms
! for purines, and C4 and N1 atoms for pyrimmidines. Sites "b" and "c" represent
! the C4' carbon in the deoxyribose ring. The N1, N2 vectors along with the N3
! vector (the downward vector that connects the two nucleotides as the backbone)
! are used to find the stacking angle (tau).
!
! Experimental values from JACS 1992, Vol. 115, Num. 4,1205-14, which can be found
! condensed in Berova's Circular Dichroism: Principles and Applications, are used
! to establish the orientation (theta) of the electronic dipole transition moment
! (edtm) per band. The alpha angle, angle between points abc with reference to
! theta. In the case of no tilt/roll then alpha is 90; however, if the nucleotide
! is not perfectly linear then the theta must be added to the simple angle(abc).
! The case is similar to gamma which is the angle between points bcd.
!
! The following abreviations notation NxSy can be expanded as the xth nucleotide
! at the yth site. By arbitrary conventient, in purines S1=N3 and S2=C6 whereas
! in pyrimmidines S1=N1 and S2=C4.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       CHARACTER(128) :: line
       INTEGER(kind=16) :: limit
       INTEGER :: i, j, t, skip, residue, ios, nLines, nResidues, Time
       DOUBLE PRECISION :: pi, theta, a(4), b(4), c(4), d(4), R(4), alpha, beta
       DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:,:)

       pi=4.D0*DATAN(1.D0)
       theta = 90.d0*(2*pi/360.d0)

       nResidues=2

       ! File of all atomic sites
       OPEN(1,file="./DUMMYselected.gro")
       READ(1,*)
       READ(1,*)
! Allocate atoms array
       ALLOCATE(atoms(nResidues,8,3))

       ! Output parameter files
       OPEN(10,file="./DistanceDUMMY", status='REPLACE')
       OPEN(20,file="./TwistDUMMY", status='REPLACE')
       OPEN(30,file="./RollDUMMY", status='REPLACE')
       OPEN(40,file="./TiltDUMMY", status='REPLACE')



       Time=1
       ! Time Loop
       DO t=1,Time

              atoms = 0.d0
              DO residue=1,nResidues
                     ! Get the XYZ cartesian coordinates for site 1 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 1, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 1, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 1, 3)

                     ! Get the XYZ cartesian coordinates for site 2 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 2, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 2, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 2, 3)

                     ! Get the XYZ cartesian coordinates for site 3 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 3, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 3, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 3, 3)

                     ! Get the XYZ cartesian coordinates for site 4 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 4, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 4, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 4, 3)

                     ! Get the XYZ cartesian coordinates for site 5 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 5, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 5, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 5, 3)

                     ! Get the XYZ cartesian coordinates for site 6 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 6, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 6, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 6, 3)

                     ! Get the XYZ cartesian coordinates for site 4 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 7, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 7, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 7, 3)

                     ! Get the XYZ cartesian coordinates for site 2 atom
                     READ(1,'(A)')line
                     READ(line(24:30),'(F13.9)')atoms(residue, 8, 1)
                     READ(line(32:38),'(F13.9)')atoms(residue, 8, 2)
                     READ(line(40:46),'(F13.9)')atoms(residue, 8, 3)
              END DO

              DO residue=1,nResidues-1
                     ! Site 2    C
                     b(1) = 0.5*(atoms(residue, 1, 1)+atoms(residue, 2, 1))
                     b(2) = 0.5*(atoms(residue, 1, 2)+atoms(residue, 2, 2))
                     b(3) = 0.5*(atoms(residue, 1, 3)+atoms(residue, 2, 3))
                     b(4) = DSQRT(DOT_PRODUCT(b(1:3),b(1:3)))

                     ! Site 3
                     c(1) = 0.5*(atoms(residue+1, 1, 1)+atoms(residue+1, 2, 1))
                     c(2) = 0.5*(atoms(residue+1, 1, 2)+atoms(residue+1, 2, 2))
                     c(3) = 0.5*(atoms(residue+1, 1, 3)+atoms(residue+1, 2, 3))
                     c(4) = DSQRT(DOT_PRODUCT(c(1:3),c(1:3)))

                     ! Site 1    A
                     CALL ROTATE_IN_PLANE(a, atoms(residue, 1, :), atoms(residue, 2, :), &
                                                        & atoms(residue, 3, :), atoms(residue, 4, :), theta)
                     a(1) = b(1)+a(1)
                     a(2) = b(2)+a(2)
                     a(3) = b(3)+a(3)
                     a(4) = DSQRT(DOT_PRODUCT(a(1:3),a(1:3)))

                     ! Site 4
                     CALL ROTATE_IN_PLANE(d, atoms(residue+1, 1, :), atoms(residue+1, 2, :), &
                                                        & atoms(residue+1, 3, :), atoms(residue+1, 4, :), theta)
                     d(1) = c(1)+d(1)
                     d(2) = c(2)+d(2)
                     d(3) = c(3)+d(3)
                     d(4) = DSQRT(DOT_PRODUCT(d(1:3),d(1:3)))

                     ! Write distances to file
                     R = 0.d0
                     CALL DISTANCE(R, b(1:3), c(1:3))
                     !!WRITE(10, '(3(F16.8,A),F16.8)')R(1),',',R(2),',',R(3),',',R(4)
                     WRITE(10, '(F16.8)') R(4)

                     ! Write twist to file
                     alpha = 0.d0
                     CALL DIHEDRAL(alpha, a(1:3), b(1:3), c(1:3), d(1:3))
                     WRITE(20, '(F16.8)')alpha

                     ! Write roll to file
                     alpha = 0.d0
                     CALL DIHEDRAL(alpha, atoms(residue, 5, :), atoms(residue, 6, :), &
                                                        &atoms(residue, 7, :), atoms(residue, 8, :))

                     beta = 0.d0
                     CALL DIHEDRAL(beta, atoms(residue+1, 5, :), atoms(residue+1, 6, :), &
                                                        &atoms(residue+1, 7, :), atoms(residue+1, 8, :))
                     WRITE(30, '(F16.8, A, F16.8)')alpha,',',beta

                     ! Write tilt angles to file
                     alpha = 0.d0
                     CALL ANGLE(alpha, a, c, b)

                     beta = 0.40
                     CALL ANGLE(beta, d, b, c)
                     WRITE(40, '(F16.8, A, F16.8)')alpha,',',beta
              END DO
       END DO

       DEALLOCATE(atoms)

       CLOSE(1)
       CLOSE(10)
       CLOSE(20)
       CLOSE(30)
       CLOSE(40)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END PROGRAM STRUCTURE_CALCULATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ANGLE(Output, v1, Site1, Site2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(in) :: v1(4), Site1(3), Site2(3)
       DOUBLE PRECISION, INTENT(inout) :: Output
       DOUBLE PRECISION :: v2(4)

       v2(1) = Site1(1)-Site2(1)
       v2(2) = Site1(2)-Site2(2)
       v2(3) = Site1(3)-Site2(3)
       v2(4) = DSQRT(DOT_PRODUCT(v2(1:3), v2(1:3)))

       Output = DACOS(DOT_PRODUCT(v1(1:3), v2(1:3))/(v1(4)*v2(4)))
END SUBROUTINE ANGLE

SUBROUTINE DISTANCE(Output, Site1, Site2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(in) :: Site1(3), Site2(3)
       DOUBLE PRECISION, INTENT(inout) :: Output(4)

       Output(1) = Site1(1)-Site2(1)
       Output(2) = Site1(2)-Site2(2)
       Output(3) = Site1(3)-Site2(3)
       Output(4) = DSQRT(DOT_PRODUCT(Output(1:3), Output(1:3)))
END SUBROUTINE DISTANCE

SUBROUTINE DIHEDRAL(Output, Site1, Site2, Site3, Site4)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(in) :: Site1(3), Site2(3), Site3(3), Site4(3)
       DOUBLE PRECISION, INTENT(inout) :: Output
       DOUBLE PRECISION :: x, y, v1(4), v2(4), v3(4), N1(3), N2(3), N3(3)
       DOUBLE PRECISION :: pi
       pi = 4.d0*DATAN(1.d0)

! The results for atan2(y,x) are shown below. Notice however that these values
! are to produce angles from [-pi,pi]. These values are mapped by adding 2*pi to
! all negative values of tau to produce the angular range [0,2*pi]. In addition,
! to avoid the potential divergence producing error in the code, tau=0 when
! x=y=0.
!
!                   /   arctan(y/x)              x>0
!                   |        arctan(y/x)+pi/2       y>=0,x<0
!                   |        arctan(y/x)-pi/2       y<0,x<0
!        atan2(y,x) = |
!                   |        pi/2                     y<0,x=0
!                   |       -pi/2                     y>0,x=0
!                   \        undefine              y=0,x=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       TAU ANGLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reestablish vectors with proper orientation

       v1(1) = Site1(1)-Site2(1)
       v1(2) = Site1(2)-Site2(2)
       v1(3) = Site1(3)-Site2(3)
       v1(4) = DSQRT(DOT_PRODUCT(v1(1:3), v1(1:3)))

       v2(1) = Site2(1)-Site3(1)
       v2(2) = Site2(2)-Site3(2)
       v2(3) = Site2(3)-Site3(3)
       v2(4) = DSQRT(DOT_PRODUCT(v2(1:3), v2(1:3)))

       v3(1) = Site3(1)-Site4(1)
       v3(2) = Site3(2)-Site4(2)
       v3(3) = Site3(3)-Site4(3)
       v3(4) = DSQRT(DOT_PRODUCT(v3(1:3), v3(1:3)))

! The following creates an orthogonal frame to calculate tau.
! Where <v> = v/(vx**2+vy**2+vz**2)**0.5

! N1 = <v1 x v2>

       N1(1) =  (v1(2)*v2(3)-v1(3)*v2(2))
       N1(2) = -(v1(1)*v2(3)-v1(3)*v2(1))
       N1(3) =  (v1(1)*v2(2)-v1(2)*v2(1))

       N1(1) = N1(1)/(v1(4)*v2(4))
       N1(2) = N1(2)/(v1(4)*v2(4))
       N1(3) = N1(3)/(v1(4)*v2(4))

!N2 = <v2 x v3>

       N2(1) =  (v2(2)*v3(3)-v2(3)*v3(2))
       N2(2) = -(v2(1)*v3(3)-v2(3)*v3(1))
       N2(3) =  (v2(1)*v3(2)-v2(2)*v3(1))

       N2(1) = N2(1)/(v2(4)*v3(4))
       N2(2) = N2(2)/(v2(4)*v3(4))
       N2(3) = N2(3)/(v2(4)*v3(4))

!N3 = N1 x <v2>

       v2(1) = v2(1)/v2(4)
       v2(2) = v2(2)/v2(4)
       v2(3) = v2(3)/v2(4)

       N3(1) =  (N1(2)*v2(3)-N1(3)*v2(2))
       N3(2) = -(N1(1)*v2(3)-N1(3)*v2(1))
       N3(3) =  (N1(1)*v2(2)-N1(2)*v2(1))

! The vectors are projected onto the normal frame to create x and y
! components for the two argument atan2(y,x) function.
! x = N1*N2
! y = N3*N2

       x = DOT_PRODUCT(N1, N2)
       y = DOT_PRODUCT(N3, N2)
       Output = DATAN2(y,x)
END SUBROUTINE DIHEDRAL

SUBROUTINE ROTATE_IN_PLANE(Output, Site1, Site2, Site3, Site4, theta)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(in) :: Site1(3), Site2(3), Site3(3), Site4(3), theta
       DOUBLE PRECISION, INTENT(inout) :: Output(4)
       DOUBLE PRECISION :: N1(3), N2(3), RotAxis(4), R(3,3)

! Virtual atom is generated by rotating the end of vector N1 within the plane
! defined by the plane of N1, N2 along the rotation axis.
!
!                  ^ RotAxis
!                  |
!                  |
!  N1 <-----/
!                 /
!         v  N2

! First in plane vector
       N1(1) = Site1(1)-Site2(1)
       N1(2) = Site1(2)-Site2(2)
       N1(3) = Site1(3)-Site2(3)

! Second in plane vector
       N2(1) = Site3(1)-Site4(1)
       N2(2) = Site3(2)-Site4(2)
       N2(3) = Site3(3)-Site4(3)

! Normal Vector to define axis of rotation
       RotAxis(1) =  (N1(2)*N2(3)-N1(3)*N2(2))
       RotAxis(2) = -(N1(1)*N2(3)-N1(3)*N2(1))
       RotAxis(3) =  (N1(1)*N2(2)-N1(2)*N2(1))
       RotAxis(4) = DSQRT(DOT_PRODUCT(RotAxis(1:3), RotAxis(1:3)))

       RotAxis(1) = RotAxis(1)/RotAxis(4)
       RotAxis(2) = RotAxis(2)/RotAxis(4)
       RotAxis(3) = RotAxis(3)/RotAxis(4)

! Define rotation operator
       R(1,1)=(cos(theta)+(RotAxis(1)**2)*(1.0-cos(theta)))
       R(1,2)=((RotAxis(1)*RotAxis(2))*(1.0-cos(theta))-RotAxis(3)*sin(theta))
       R(1,3)=((RotAxis(1)*RotAxis(3))*(1.0-cos(theta))+RotAxis(2)*sin(theta))
       R(2,1)=((RotAxis(2)*RotAxis(1))*(1.0-cos(theta))+RotAxis(3)*sin(theta))
       R(2,2)=(cos(theta)+(RotAxis(2)**2)*(1.0-cos(theta)))
       R(2,3)=((RotAxis(2)*RotAxis(3))*(1.0-cos(theta))-RotAxis(1)*sin(theta))
       R(3,1)=((RotAxis(3)*RotAxis(1))*(1.0-cos(theta))-RotAxis(2)*sin(theta))
       R(3,2)=((RotAxis(3)*RotAxis(2))*(1.0-cos(theta))+RotAxis(1)*sin(theta))
       R(3,3)=(cos(theta)+(RotAxis(3)**2)*(1.0-cos(theta)))

! Rotate first in plane vector by theta
       Output(1) = R(1,1)*N1(1)+R(1,2)*N1(2)+R(1,3)*N1(3)
       Output(2) = R(2,1)*N1(1)+R(2,2)*N1(2)+R(2,3)*N1(3)
       Output(3) = R(3,1)*N1(1)+R(3,2)*N1(2)+R(3,3)*N1(3)
       Output(4) = DSQRT(DOT_PRODUCT(Output(1:3), Output(1:3)))
END SUBROUTINE ROTATE_IN_PLANE
