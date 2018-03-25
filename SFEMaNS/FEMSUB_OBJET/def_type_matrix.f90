!
!Authors: Jean-Luc Guermond, Copyright 2009
!
MODULE def_type_matrix
  TYPE csr_matrix
      REAL(KIND=8), POINTER, DIMENSION(:) :: aa
      INTEGER,      POINTER, DIMENSION(:) :: ia, ja
  END TYPE csr_matrix
END MODULE def_type_matrix
