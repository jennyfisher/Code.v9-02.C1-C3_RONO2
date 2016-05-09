! $Id: alpha_rno3.f,v 1.1 2016/05/09 14:06:27 jaf Exp $
      REAL*8 FUNCTION ALPHA_RNO3( XCARBN, PRESS_ATM, T_K )
!
!******************************************************************************
!  Function ALPHA_RNO3 returns organic nitrate yields YN = RKA/(RKA+RKB)
!  from RO2+NO reactions as a function of pressure and temperature for
!  species with this dependency (so far just ETNO3) (jaf, 5/9/16)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) XCARBN (REAL*8)    : Number of C atoms in RO2 - used to identify RO2
!  (2 ) PRESS_ATM (REAL*8) : Air pressure  [atm]
!  (3 ) T_K       (REAL*8) : Temperature   [K]
!    
!  NOTES: 
!  09 May 2016 - J. Fisher - Initial version
!
!  REFERENCES:
!  1. Butkovskaya, N., A. Kukui, and G. Le Bras. Pressure and Temperature
!     Dependence of Ethyl Nitrate Formation in the C2H5O2 + NO Reaction.
!     J. Phys. Chem. A, 114, 956-964, 2010.
!******************************************************************************
!
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP

      IMPLICIT NONE

      ! Arguments
      REAL*8, INTENT(IN) :: XCARBN, PRESS_ATM, T_K

      !=================================================================
      ! ALPHA_RNO3 begins here!
      !=================================================================

      ! Use carbon number to check that this species has a yield defined
      ! in this function.
      IF ( XCARBN /= 2 ) THEN
         WRITE(6,*) 'No P-T dependent yield found for species with '
         WRITE(6,*) 'C number = ',XCARBN
         WRITE(6,*) 'Check globchem.dat and alpha_rno3.F!'
         CALL GEOS_CHEM_STOP
      ENDIF

      ! EtNO3 -- turn into IF statement if additional species added
      ! Calculate yield based on Butkovskaya 2010
      ! 2.9488 is coefficient for P in atm from 3.88e-3 for P in Torr
      ! Divide by 100 to convert from % to yield
      ALPHA_RNO3 = 1e-2 * ( 2.9488 * PRESS_ATM + 
     &                  0.365 * (1e0 + 1.5e3 * (1e0/T_K - 1e0/298e0)) )

      ! Return to calling program
      END FUNCTION ALPHA_RNO3
