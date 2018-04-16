!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: rono2_mod
!
! !DESCRIPTION: Module RONO2\_MOD contains subroutines to calculate the ocean
!  source/sink of C1-C3(?) RONO2.
!\\
!\\
! !INTERFACE:
!
      MODULE RONO2_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_RONO2   
      PUBLIC :: READ_CHL
      PUBLIC :: OCEAN_SOURCE_RONO2
      PUBLIC :: OCEAN_SINK_RONO2
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Jacob, D.J., B.D. Field, E. Jin, I. Bey, Q. Li, J.A. Logan, and 
!        R.M. Yantosca, "Atmospheric budget of rono2", Geophys. Res. Lett., 
!        107(D11), 4100, 2002. 
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!
! !REVISION HISTORY:
!  14 Jun 2016 - J. Fisher - Initial version, based on acetone_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Avogadro's number
      REAL*8,  PARAMETER   :: AVO      = 6.022d23 

      ! Molecules RONO2 / kg C 
      REAL*8,  PARAMETER   :: XNUMOL_MeNO3 = AVO / 77d-3 
      REAL*8,  PARAMETER   :: XNUMOL_EtNO3 = AVO / 91d-3 

      ! Gas constant, J/mol/K
      REAL*8,  PARAMETER   :: R = 8.314d0

      ! Coefficients for Schmidt numbers, scaled to CO2
!      ! Vb: Schroeder additive method as described in Johnson et al. 2010
!      REAL*8,  PARAMETER     :: Vb_CO2   = 35d0 ! cm3/mol
!      REAL*8,  PARAMETER     :: Vb_MeNO3 = 63d0 ! cm3/mol
!      REAL*8,  PARAMETER     :: Vb_EtNO3 = 84d0 ! cm3/mol

      ! Vb: Le Bas additive  method as in Dahl et al. 2005
      ! Now use Kornilov & Klselev 2015 for MeNO3, EtNO3
      REAL*8,  PARAMETER     :: Vb_CO2   = 34d0 ! cm3/mol
      REAL*8,  PARAMETER     :: Vb_MeNO3 = 64d0 ! cm3/mol
      REAL*8,  PARAMETER     :: Vb_EtNO3 = 82.2d0 ! cm3/mol

      ! CO2 Schmidt number coefficients (see ocean_mercury_mod.F)
      REAL*8,  PARAMETER     :: ScCO2_A = 644.7d0
      REAL*8,  PARAMETER     :: ScCO2_B = -6.16d0
      REAL*8,  PARAMETER     :: ScCO2_C = 0.11d0

      ! Coefficients for Henry's Law constants, from Sander et al. 2015
      REAL*8, PARAMETER      :: H0_MeNO3 = 2.0d-2 ! mol/m3/Pa
      REAL*8, PARAMETER      :: H0_EtNO3 = 1.6d-2 ! mol/m3/Pa
      REAL*8, PARAMETER      :: H_T_MeNO3 = 4.7d3 ! K
      REAL*8, PARAMETER      :: H_T_EtNO3 = 5.4d3 ! K

!
! !MODULE VARIABLES:
!
      REAL*8,  ALLOCATABLE :: CHL(:,:)                              

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_source_rono2
!
! !DESCRIPTION: Subroutine OCEAN\_SOURCE\_RONO2 specifies the ocean source of 
!  rono2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SOURCE_RONO2( I, J, MeNO3, EtNO3, State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : CHECK_VALUE
      USE DIAG_MOD,           ONLY : AD11
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2, GET_YMID, GET_XMID
      USE TIME_MOD,           ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND11
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: I           ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN)    :: J           ! GEOS-Chem latitude index
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: MeNO3     ! MeNO3 emiss [atoms C/s]
      REAL*8,         INTENT(INOUT) :: EtNO3     ! EtNO3 emiss [atoms C/s]
!
! !REMARKS:

! 
! !REVISION HISTORY: 
!  14 Jun 2016 - J. Fisher - Initial version, based on OCEAN_SOURCE_ACET
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: U, USTAR, TC, CD, ScCO2
      REAL*8                 :: DTSRCE, AREA_CM2, FOCEAN, LAT, LON
      REAL*8                 :: SC_MeNO3, HSTAR_MeNO3, C_MeNO3
      REAL*8                 :: KG_MeNO3, KL_MeNO3, KKL_MeNO3
      REAL*8                 :: SC_EtNO3, HSTAR_EtNO3, C_EtNO3
      REAL*8                 :: KG_EtNO3, KL_EtNO3, KKL_EtNO3
      REAL*8                 :: OCEAN_MENO3, OCEAN_ETNO3
      REAL*8,  EXTERNAL      :: SFCWINDSQR
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! OCEAN_SOURCE_RONO2 begins here!
      !=================================================================

      ! Emission timestep in seconds
      DTSRCE   = GET_TS_EMIS() * 60d0 

      ! Fraction of (I,J) that is ocean
      FOCEAN   = 1d0 - State_Met%FRCLND(I,J)    

      ! Area of grid box (I,J) in cm^2
      AREA_CM2 = GET_AREA_CM2( I, J, 1 )    

      ! Latitude of grid box
      LAT = GET_YMID( I, J, 1)
      LON = GET_XMID( I, J, 1)

      !=================================================================
      ! Compute ocean source if 50% of the box is covered by water
      ! and we are in region where ocean fluxes defined (20S - 20N)
      !=================================================================
!      IF ( FOCEAN > 0.5d0 .and. ABS(LAT) <= 20d0 ) THEN 
      IF ( FOCEAN > 0.5d0 ) THEN

         ! Magnitude of resultant wind [m/s]
         ! SFCWINDSQR(I,J) is needed since this will compute the square
         ! of the surface wind correctly for all GEOS models (bmy, 9/5/01)
         U        = SQRT( SFCWINDSQR( I, J,
     &                                State_Met%U10M,
     &                                State_Met%V10M ) )

         ! Updated KG to the Johnson [2010] parameterization (evf, 5/13/2011)
         ! USTAR = friction velocity (U* in Johnson [2010]
         USTAR = SQRT(6.1d-4 + U*6.3d-5)*U
         ! CD = drag coefficient
         CD = (USTAR/U)**2

         ! TC is temperature in Celsius
         ! Also make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01)
         TC       = State_Met%TS(I,J) - 273.15d0                           
!!         TC       = MIN( MAX( TC, -5d0 ), 30d0 )

         ! Schmidt number of CO2 in seawater, with temperature
         ! correction based on Poissant et al. (2000) - see Hg code
         ScCO2 = ScCO2_A + TC * ( ScCO2_B + TC * ( ScCO2_C ) )

         ! Henry's Law constants from Sander et al. 2015 
         ! (Hcp, in mol/m3/Pa)
         HSTAR_MeNO3 = H0_MeNO3 * 
     &           EXP(-H_T_MeNO3 * (1./298.15d0 - 1./State_Met%TS(I,J)))
         HSTAR_EtNO3 = H0_EtNO3 * 
     &           EXP(-H_T_EtNO3 * (1./298.15d0 - 1./State_Met%TS(I,J)))

         ! Need the dimensionless version, Hcc - Hcp * RT 
         ! this is [conc aq]/[conc gas phase]
         HSTAR_MeNO3  = HSTAR_MeNO3 * R * State_Met%TS(I,J)
         HSTAR_EtNO3  = HSTAR_EtNO3 * R * State_Met%TS(I,J)

         ! Want Henry constant expressed as [conc gas-phase]/[conc aq]
         ! Take reciprocal.
         HSTAR_MeNO3 = 1d0/HSTAR_MeNO3
         HSTAR_EtNO3 = 1d0/HSTAR_EtNO3

         ! SC is Schmidt # for RONO2 [unitless], scaled to CO2
         SC_MeNO3 =  ScCO2 * (Vb_MeNO3 / Vb_CO2)**0.6
         SC_EtNO3 =  ScCO2 * (Vb_EtNO3 / Vb_CO2)**0.6

         ! KL is conductance for mass transfer in liquid phase 
         ! (Nightingale et al 2000b), which has units of [cm/hr]
         KL_MeNO3 = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/SC_MeNO3 )
         KL_EtNO3 = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/SC_EtNO3 )

         !KG is airside transfer velocity (Johnson 2010)
         ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
         KG_MeNO3 = 1d-3 + (USTAR/(13.3*SC_MeNO3**(0.5) + 
     &                      CD**(-0.5)-5+LOG(SC_MeNO3)/0.8))
         KG_EtNO3 = 1d-3 + (USTAR/(13.3*SC_EtNO3**(0.5) + 
     &                      CD**(-0.5)-5+LOG(SC_EtNO3)/0.8))
         KG_MeNO3 = KG_MeNO3 * 360000d0
         KG_EtNO3 = KG_EtNO3 * 360000d0

         ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
         ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
         KKL_MeNO3  = 1d0 / ( 1d0/KL_MeNO3 + 
     &                        1d0/( HSTAR_MeNO3 * KG_MeNO3 ) ) 
         KKL_EtNO3  = 1d0 / ( 1d0/KL_EtNO3 + 
     &                        1d0/( HSTAR_EtNO3 * KG_EtNO3 ) ) 
         KKL_MeNO3  = KKL_MeNO3 / 3600d0
         KKL_EtNO3  = KKL_EtNO3 / 3600d0

         ! Turn off ocean uptake over snow/ice (ALBEDO > 0.4)
         IF ( State_Met%ALBD(I,J) > 0.4d0 ) KKL_MeNO3 = 0d0
         IF ( State_Met%ALBD(I,J) > 0.4d0 ) KKL_EtNO3 = 0d0
      
         !===============================================================
         ! Assume constant seawater concentration by latitude (C), pM
         ! 
!         ! For now, use values from Dahl et al. 2007 distributions.
!         ! Separate into low- and high- chlorophyll (0.1 mg/m3).
!!         IF ( ABS(LAT) <= 10.0 ) THEN
!!
!!            ! Roughly "Warm Pool"
!!            IF ( CHL(I,J) < 0.1 ) THEN
!!                C_MeNO3 = 25d0
!!                C_EtNO3 = 5d0
!!            ! Roughly "Equatorial" & "Sub-equatorial" / "HC"
!!            ELSE
!!                C_MeNO3 = 300d0
!!                C_EtNO3 = 60d0
!!            ENDIF
!!
!!         ! Roughly "Gyre"
!!         ELSE IF ( ABS(LAT) <= 20.0 ) THEN
!!            C_MeNO3 = 5d0
!!            C_EtNO3 = 2d0
!!         ELSE
!!            C_MeNO3 = 0d0
!!            C_EtNO3 = 0d0
!!         ENDIF
!!
!!
!         ! Test applying chlorophyll-based cutoffs everywhere globally
!         ! Maintain higher values in low-Chl sub-equatorial/warm pool
!         IF ( CHL(I,J) >= 0.1 ) THEN
!            C_MeNO3 = 300d0
!            C_EtNO3 = 60d0
!         ELSE
!            IF ( ABS(LAT) <= 10.0 ) THEN
!                C_MeNO3 = 25d0
!                C_EtNO3 = 5d0
!            ELSE
!                C_MeNO3 = 5d0
!                C_EtNO3 = 2d0
!            ENDIF
!         ENDIF
!! end test
!@!         ! New test - tie to ocean nitrite (NO2-). Based on Dahl et al.
!@!         ! 2008. May explain relationship with chlorophyll as co-variate
!@!         ! Nitrite measurements very limited, so for now using ODV
!@!         ! datasets GEOSECS and Reid-Mantyla to identify regions with
!@!         ! elevated nitrite by lat/lon. This only allows two options.
!@!         ! Might eventually tie tropical amount to chlorophyll as better
!@!         ! proxy...
!@!         IF (( LAT < -40.0 ) .OR. 
!@!     &      ( LAT > 50.0  ) .OR.
!@!     &      ((LAT > 40.0) .AND. (LON < -130.0)) .OR.
!@!     &      ((LAT > 40.0) .AND. (LON > 130.0)) .OR.
!@!     &      ((LAT > -10.0) .AND. (LAT < 5.0) .AND. (LON < -80.0)) .OR.
!@!     &      ((LAT > -10.0) .AND. (LAT < 5.0) .AND. (LON > 130.0)) .OR.
!@!     &      ((LAT > -15.0) .AND. (LAT < -5.0) .AND. (LON > -25.0) .AND.
!@!     &       (LON < 30.0)) ) THEN
!@!            C_MeNO3 = 350d0
!@!            C_EtNO3 = 60d0
!@!         ELSE
!@!            C_MeNO3 = 0d0
!@!            C_EtNO3 = 0d0
!@!         ENDIF

         ! Combine the two... use nitrite based cut-off for polar waters
         ! but chlorophyll as more spatially resolved proxy for
         ! equatorial waters...
         ! Note relationships in Dahl et al. 2007 - MeNO3 ~ 6xEtNO3
         ! There is no observational rationale for NPac << SPac, EqPac
         ! except to fit the aircraft obs. Same for SPac < EqPac
         IF ( LAT < -40.0 ) THEN
            C_MeNO3 = 200d0
            C_EtNO3 = 33d0
         ELSE IF (( LAT > 50.0  ) .OR.
     &       ((LAT > 40.0) .AND. (LON < -130.0)) .OR.
     &       ((LAT > 40.0) .AND. (LON >  130.0)) ) THEN
            C_MeNO3 =120d0 
            C_EtNO3 = 20d0
         ELSE IF ( (LAT > -15.0) .AND. (LAT < 10.0) ) THEN
            IF ( CHL(I,J) >= 0.1 ) THEN
                C_MeNO3 = 400d0
                C_EtNO3 = 67d0
            ELSE
                C_MeNO3 = 25d0
                C_EtNO3 = 5d0
            ENDIF
         ELSE
            C_MeNO3 = 0d0
            C_EtNO3 = 0d0
         ENDIF


         ! convert from pM to kg/cm3
         C_MeNO3 = C_MeNO3 * 77d-18
         C_EtNO3 = C_EtNO3 * 91d-18

         !correct for the fraction of the grid cell that is ocean
         !and compute the flux ( kg/cm2/s)(evf, 5/11/11)
         OCEAN_MeNO3 = C_MeNO3 * KKL_MeNO3 * FOCEAN 
         OCEAN_EtNO3 = C_EtNO3 * KKL_EtNO3 * FOCEAN 
          
         ! Convert to kg / box / step (evf, 5/11/11)
         OCEAN_MeNO3  = OCEAN_MeNO3 * DTSRCE * AREA_CM2
         OCEAN_EtNO3  = OCEAN_EtNO3 * DTSRCE * AREA_CM2

      ELSE

         ! If there is less than 50% water in the grid box, zero 
         ! the ocean source from RONO2 
         OCEAN_MeNO3 = 0d0
         OCEAN_EtNO3 = 0d0
         
      ENDIF

      ! Add ocean source to total biogenic source in [molec/box/s]
      MeNO3 = MeNO3 + ( OCEAN_MeNO3 * XNUMOL_MeNO3 / DTSRCE )
      EtNO3 = EtNO3 + ( OCEAN_EtNO3 * XNUMOL_EtNO3 / DTSRCE )
            
      !=================================================================
      ! ND11 diag -- save RONO2 from the ocean in [molec/cm2/s]
      !=================================================================
      IF ( ND11 > 0 ) THEN
         AD11(I,J,6) = AD11(I,J,6) + ( OCEAN_MeNO3 * XNUMOL_MeNO3 ) / 
     &                               ( AREA_CM2   * DTSRCE   )
         AD11(I,J,8) = AD11(I,J,8) + ( OCEAN_EtNO3 * XNUMOL_EtNO3 ) / 
     &                               ( AREA_CM2   * DTSRCE   )
      ENDIF

      END SUBROUTINE OCEAN_SOURCE_RONO2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_sink_rono2
!
! !DESCRIPTION: Subroutine OCEAN\_SINK\_RONO2 applies the ocean sink to global
!  RONO2 concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SINK_RONO2( MeNO3, EtNO3, State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD11
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: MeNO3(IIPAR,JJPAR) !  mass [kg]
      REAL*8,         INTENT(INOUT) :: EtNO3(IIPAR,JJPAR) !  mass [kg]
! 
! !REVISION HISTORY: 
!  14 Jun 2016 - J. Fisher - Initial version, based on OCEAN_SINK_ACET
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J
      REAL*8                 :: AREA_CM2, DTCHEM, FOCEAN
      REAL*8                 :: U, USTAR, TC, CD, ScCO2
      REAL*8                 :: SC_MeNO3, HSTAR_MeNO3, C_MeNO3
      REAL*8                 :: KG_MeNO3, KL_MeNO3, KKL_MeNO3
      REAL*8                 :: FLUX_MeNO3, CG_MeNO3, OCEAN_MeNO3
      REAL*8                 :: SC_EtNO3, HSTAR_EtNO3, C_EtNO3
      REAL*8                 :: KG_EtNO3, KL_EtNO3, KKL_EtNO3
      REAL*8                 :: FLUX_EtNO3, CG_EtNO3, OCEAN_EtNO3
      REAL*8                 :: MeNO3_LOSS, EtNO3_LOSS
!
! !EXTERNAL FUNCTIONS:
!
      REAL*8, EXTERNAL  :: BOXVL, SFCWINDSQR
!
! !DEFINED PARAMETERS:
!

      !=================================================================
      ! OCEAN_SINK_RONO2 begins here! 
      !
      ! Compute RONO2 lost to ocean sink and subtract from STT
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, FOCEAN, U, USTAR, CD, TC )
!$OMP+PRIVATE( HSTAR_MeNO3, Sc_MeNO3, KL_MeNO3, KG_MeNO3, KKL_MeNO3)
!$OMP+PRIVATE( CG_MeNO3, FLUX_MeNO3, OCEAN_MeNO3, MeNO3_LOSS)
!$OMP+PRIVATE( HSTAR_EtNO3, Sc_EtNO3, KL_EtNO3, KG_EtNO3, KKL_EtNO3)
!$OMP+PRIVATE( CG_EtNO3, FLUX_EtNO3, OCEAN_EtNO3, EtNO3_LOSS)

      ! Loop over grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box area in cm2
         AREA_CM2 = GET_AREA_CM2( I, J, 1 )

         ! Fraction of grid box that is ocean
         FOCEAN = 1d0 - State_Met%FRCLND(I,J)

         !===========================================================
         ! Only compute ocean sink if there is more than 50% ocean
         ! in the grid box, and if it is not ice (albedo > 0.4)
         ! (mje, rvm, bmy, 11/26/01)
         !===========================================================
         IF ( FOCEAN > 0.5d0 .and.
     &        State_Met%ALBD(I,J) <= 0.4d0 ) THEN

            ! Magnitude of surface wind [m/s]
            ! SFCWINDSQR(I,J) is needed since this will compute the 
            ! square of the surface wind correctly for all GEOS models 
            U = SQRT( SFCWINDSQR( I, J,
     &                            State_Met%U10M,
     &                            State_Met%V10M ) )

            ! USTAR = friction velocity (U* in Johnson [2010]
            USTAR  = SQRT(6.1d-4 + U*6.3d-5)*U
            ! CD = drag coefficient
            CD = (USTAR/U)**2
         
            ! TC is temperature in Celsius
            ! Bug fix: make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01) 
            TC    = State_Met%TS(I,J) - 273.15d0  
!!            TC    = MIN( MAX( TC, -5d0 ), 30d0 )

            ! Henry's Law constants from Sander et al. 2015 
            ! (Hcp, in mol/m3/Pa)
            HSTAR_MeNO3 = H0_MeNO3 * 
     &            EXP(-H_T_MeNO3 * (1./298.15d0 - 1./State_Met%TS(I,J)))
            HSTAR_EtNO3 = H0_EtNO3 * 
     &            EXP(-H_T_EtNO3 * (1./298.15d0 - 1./State_Met%TS(I,J)))

            ! Need the dimensionless version, Hcc - Hcp * RT 
            ! this is [conc aq]/[conc gas phase]
            HSTAR_MeNO3  = HSTAR_MeNO3 * R * State_Met%TS(I,J)
            HSTAR_EtNO3  = HSTAR_EtNO3 * R * State_Met%TS(I,J)

            ! Want Henry constant expressed as [conc gas-phase]/[conc aq]
            ! Take reciprocal.
            HSTAR_MeNO3 = 1d0/HSTAR_MeNO3
            HSTAR_EtNO3 = 1d0/HSTAR_EtNO3

            ! SC is Schmidt # for RONO2 [unitless], scaled to CO2
            SC_MeNO3 =  ScCO2 * (Vb_MeNO3 / Vb_CO2)**0.6
            SC_EtNO3 =  ScCO2 * (Vb_EtNO3 / Vb_CO2)**0.6

            ! KL is conductance for mass transfer in liquid phase 
            ! (Nightingale et al 2000b), which has units of [cm/hr]
            KL_MeNO3 = ( 0.24d0*U*U + 0.061d0*U ) * SQRT(600d0/SC_MeNO3)
            KL_EtNO3 = ( 0.24d0*U*U + 0.061d0*U ) * SQRT(600d0/SC_EtNO3)

            ! KG is airside transfer velocity (Johnson 2010)
            ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
            KG_MeNO3 = 1d-3 + (USTAR/(13.3*SC_MeNO3**(0.5) + 
     &                         CD**(-0.5)-5+LOG(SC_MeNO3)/0.8))
            KG_EtNO3 = 1d-3 + (USTAR/(13.3*SC_EtNO3**(0.5) + 
     &                         CD**(-0.5)-5+LOG(SC_EtNO3)/0.8))
            KG_MeNO3 = KG_MeNO3 * 360000d0
            KG_EtNO3 = KG_EtNO3 * 360000d0

            ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
            ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
            KKL_MeNO3  = 1d0 / ( 1d0/KL_MeNO3 + 
     &                           1d0/( HSTAR_MeNO3 * KG_MeNO3 ) ) 
            KKL_EtNO3  = 1d0 / ( 1d0/KL_EtNO3 + 
     &                           1d0/( HSTAR_EtNO3 * KG_EtNO3 ) ) 
            KKL_MeNO3  = KKL_MeNO3 / 3600d0
            KKL_EtNO3  = KKL_EtNO3 / 3600d0

            ! CG is the gas concentration of RONO2 [kg/cm3]
            CG_MeNO3    = MeNO3(I,J) / BOXVL(I,J,1,State_Met)
            CG_EtNO3    = EtNO3(I,J) / BOXVL(I,J,1,State_Met)

            ! FLUX is the air-to-sea flux of RONO2 in [kg/cm2/s].
            FLUX_MeNO3  = KKL_MeNO3 * CG_MeNO3 / HSTAR_MeNO3    
            FLUX_EtNO3  = KKL_EtNO3 * CG_EtNO3 / HSTAR_EtNO3    

            !========================================================
            ! Ocean loss of RONO2 consists of the following terms:
            !
            ! (1) FLUX, the air-to-sea flux of RONO2 in [kg/cm2/s]
            ! (2) AREA_CM2, the grid box surface area
            ! (3) DTCHEM, the number of seconds per chemistry timestep
            ! (4) FOCEAN, the fraction of the grid box that is ocean.
            !
            ! The units of the resultant ocean loss (in OCEAN_RONO2) are:
            !
            !     kg    | AREA_CM2 cm2 |  DTCHEM s           kg 
            !  ---------+--------------+--------------- = ------------
            !   cm2 * s |   grid box   | chem timestep     box * step
            !========================================================
            OCEAN_MeNO3 = ( FLUX_MeNO3 * AREA_CM2 * DTCHEM * FOCEAN )
            OCEAN_EtNO3 = ( FLUX_EtNO3 * AREA_CM2 * DTCHEM * FOCEAN )
               
         ELSE

            ! If there is less than 50% water in the grid box, or  
            ! if there is ice on the ocean, then zero the ocean sink
            OCEAN_MeNO3 = 0d0
            OCEAN_EtNO3 = 0d0

         ENDIF

         ! Set to zero for safety
         MeNO3_LOSS = 0d0
         EtNO3_LOSS = 0d0

         ! Apply exponential loss to RONO2 mass
         IF ( ABS( MeNO3(I,J) ) .gt. 1d-99 )
     &      MeNO3_LOSS = MeNO3(I,J) *
     &                   ( 1d0 - EXP(-OCEAN_MeNO3/MeNO3(I,J)) )
         IF ( ABS( EtNO3(I,J) ) .gt. 1d-99 )
     &      EtNO3_LOSS = EtNO3(I,J) *
     &                   ( 1d0 - EXP(-OCEAN_EtNO3/EtNO3(I,J)) )

         ! Cap the loss -- cannot remove more than in the atmosphere
         MeNO3_LOSS = MIN( MeNO3_LOSS, MeNO3(I,J) )
         EtNO3_LOSS = MIN( EtNO3_LOSS, EtNO3(I,J) )

         ! Apply loss to atmospheric RONO2 mass
         MeNO3(I,J) = MeNO3(I,J) - MeNO3_LOSS
         EtNO3(I,J) = EtNO3(I,J) - EtNO3_LOSS

        !===========================================================
        ! Diagnostics: save RONO2 lost to ocean in [molec/cm2/s]
        !===========================================================
        IF ( ND11 > 0 ) THEN
            AD11(I,J,7) = AD11(I,J,7) + 
     &                    ( MeNO3_LOSS * XNUMOL_MeNO3 ) / 
     &                    ( AREA_CM2 * DTCHEM )
            AD11(I,J,9) = AD11(I,J,9) + 
     &                    ( EtNO3_LOSS * XNUMOL_EtNO3 ) / 
     &                    ( AREA_CM2 * DTCHEM )
        ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE OCEAN_SINK_RONO2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_chl
!
! !DESCRIPTION: Subroutine READ\_CHL reads MODIS chlorophyll distributions (mg/m3)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CHL( THISMONTH )

      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT, GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local Variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: AS
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
! 
! !REVISION HISTORY: 
!  14 Jun 2016 - J. Fisher - Initial version based on ocean_mercury_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC

      ! First time, allocate array
      IF ( FIRST ) THEN
         ALLOCATE( CHL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHL' )
         CHL = 0d0
         FIRST = .FALSE.
      ENDIF

      ! Chl file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/Chl_2003.'        // 
     &           GET_NAME_EXT()     // '.'         //
     &           GET_RES_EXT()
#else
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201007/Chl_2003.geos.'   //
     &           GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - Reading Chlorophyll from:  ', a )  

      ! TAU0 values (uses year 2003)
      TAU = GET_TAU0( THISMONTH, 1, 2003 )      

#if defined (GRID2x25)
      ! Read data
      CALL READ_BPCH2( FILENAME, 'CHLO-A-$',    0,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.FALSE. )
#else 
      CALL READ_BPCH2( FILENAME, 'CHLO-A-$',    1,
     &                 TAU,       IIPAR,        JJPAR,
     &                 1,         ARRAY(:,:,1), QUIET=.FALSE. )
#endif

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), CHL )


      END SUBROUTINE READ_CHL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_rono2
!
! !DESCRIPTION: Subroutine CLEANUP\_RONO2 deallocates module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_RONO2
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - R. Yantosca - Initial version
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  19 Mar 2012 - M. Payer    - Removed deallocation of JO1D and XRESP. They
!                              are obsolete after removal of READ_JO1D and
!                              READ_RESP routines.
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( ALLOCATED( CHL       ) ) DEALLOCATE( CHL       )

      END SUBROUTINE CLEANUP_RONO2   
!EOC
      END MODULE RONO2_MOD
