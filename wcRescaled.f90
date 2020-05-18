!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!--------------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION :: v0, v1, v2, v3, y5, y6, y7, y8
      DOUBLE PRECISION :: A, B, G, aa, bb, gg, sigm, taua, taub, taug
      DOUBLE PRECISION :: m, C, C1, C2, C3, C4, C5, C6, C7


      ! VARIABLES

      v3 = U(1)
      v0 = U(2)
      v1 = U(3)
      v2 = U(4)

      y8 = U(5)
      y5 = U(6)
      y6 = U(7)
      y7 = U(8)


      ! PARAMETERS
      A = PAR(1)
      B = PAR(2)
      G = PAR(3)
      aa = PAR(4)
      bb = PAR(5)
      gg = PAR(6)
      m = PAR(7)
      C = PAR(8)
      C1 = PAR(9)
      C2 = PAR(10)
      C3 = PAR(12)
      C4 = PAR(13)
      C5 = PAR(14)
      C6 = PAR(15)
      C7 = PAR(16)

      taua=1/aa
      taub=1/bb
      taug=1/gg


      ! RIGHT-HAND SIDE

      F(1) = y8
      F(2) = y5*(taug/taua)
      F(3)= y6*(taug/taua)
      F(4)= y7*(taug/taub)

      F(5) = (G*sigm(c5*taua*y0-c6*taub*y2)- 2*y8- v3)
      F(6) = (A*sigm(m*A*taua+c2*taua*y1 - c4*taub*y2 -c7*taug*y3)- 2*y5- v0)*(taug/taua)
      F(7) = (A*sigm(c1*taua*y0)-2*y6- v1)*(taug/taua)
      F(8) = (B*sigm(c3*taua*y0) -2*y7- v2)*(taug/taub)



      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 

      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
     
      ! VARIABLES
      U(1) = 0.0 !v3
      U(2) = 0.150158 !v0
      U(3) = 0.9363391 !v1
      U(4) = 6.90745 !v2
      U(5) = 0.0d0 !y5
      U(6) = 0.0d0 !y6
      U(7) = 0.0d0 !y8
      U(8) = 0.0d0 !y9
      
      ! PARAMETERS
      PAR(1) = 5 !A
      PAR(2) = 40 !B
      PAR(3) = 35 !G
      PAR(4) = 100 !aa
      PAR(5) = 20 !bb
      PAR(6) = 350 !gg
      PAR(7) = 90 !m
      PAR(8) = 135 !C
      PAR(9) = 135 !C1
      PAR(10) = 108 !C2
      PAR(12) = 35 !C3
      PAR(13) = 25 !C4
      PAR(14) = 450 !C5
      PAR(15) = 121  !C6
      PAR(16) = 121 !C7

      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

!----------------------------------------------------------------------
      FUNCTION sigm(x)
!----------------------------------------
! Function for integration
!----------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: sigm, x
      sigm = 5/(1 + exp(0.56*(6-x)))
      RETURN
      END FUNCTION sigm
