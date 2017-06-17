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
!  13 Dec 2016 - J. Fisher - Add MeNO3
!
!  REFERENCES:
!  1. Butkovskaya, N., A. Kukui, and G. Le Bras. Pressure and Temperature
!     Dependence of Ethyl Nitrate Formation in the C2H5O2 + NO Reaction.
!     J. Phys. Chem. A, 114, 956-964, 2010.
!  2. Williams, J. E., Le Bras, G., Kukui, A., Ziereis, H., &
!     Brenninkmeijer, C. A. M. (2014). The impact of the chemical
!     production of methyl nitrate from the NO +
!     CH3O2 reaction on the global distributions of
!     alkyl nitrates, nitrogen oxides and tropospheric ozone: a global
!     modelling study. Atmospheric Chemistry and Physics, 14(5),
!     2363-2382. http://doi.org/10.5194/acp-14-2363-2014
!  3. Flocke, F., Atlas, E., Madronich, S., Schauffler, S. M., Aikin,
!     K., Margitan, J. J., & Bui, T. P. (1998). Observations of methyl
!     nitrate in the lower stratosphere during STRAT: Implications for its
!     gas phase production mechanisms. Geophysical Research Letters,
!     25(11), 1891-1894. http://doi.org/10.1029/98GL01417 
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
      IF ( XCARBN /= 1 .AND. XCARBN /= 2 ) THEN
         WRITE(6,*) 'No P-T dependent yield found for species with '
         WRITE(6,*) 'C number = ',XCARBN
         WRITE(6,*) 'Check globchem.dat and alpha_rno3.F!'
         CALL GEOS_CHEM_STOP
      ENDIF

      IF ( XCARBN == 1 ) THEN
         ! MeNO3
         ! For now use midrange aircraft derived value from Flocke et
         ! al. (1998) as implemented in Williams et al. (2014). This is
         ! just a constant value, but it is too small to put into
         ! globchem.dat
         ! Update 16 Jun 2017
         ! MeNO3 too low everywhere. Using Carter & Atkinson formula
         ! made it way too high. Try smaller increase to 0.0003 which is
         ! upper limit of Flocke (= 1 order of magnitude smaller than
         ! Carter & Atkinson prediction extrapolated, and than lower
         ! limit of Butkovskaka)
         ALPHA_RNO3 = 0.0003
      ELSE IF ( XCARBN == 2 ) THEN
         ! EtNO3 
         ! Calculate yield based on Butkovskaya 2010
         ! 2.9488 is coefficient for P in atm from 3.88e-3 for P in Torr
         ! Divide by 100 to convert from % to yield
         ALPHA_RNO3 = 1e-2 * ( (2.9488e0 * PRESS_ATM )
     &                        *(1e0 + 1.5e3 * (1e0/T_K - 1e0/298e0)) )
      ENDIF

      ! Return to calling program
      END FUNCTION ALPHA_RNO3
