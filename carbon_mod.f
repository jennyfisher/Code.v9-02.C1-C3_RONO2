 ! $Id: carbon_mod.f,v 1.9 2005/02/10 19:53:23 bmy Exp $
      MODULE CARBON_MOD
!
!******************************************************************************
!  Module CARBON_MOD contains arrays and routines for performing a 
!  carbonaceous aerosol simulation.  Original code taken from Mian Chin's 
!  GOCART model and modified accordingly. (rjp, bmy, 4/2/04, 1/18/05)
!
!  4 Aerosol species : Organic and Black carbon 
!                    : hydrophilic (soluble) and hydrophobic of each
!
!  For secondary organic aerosol (SOA) simulation orginal code developed
!  by Chung and Seinfeld [2002] and Hong Liao from John Seinfeld's group 
!  at Caltech was taken and further modified accordingly (rjp, bmy, 7/15/04)
!  This simulation introduces additional following species:
!     ALPH, LIMO, ALCO, SOA1, SOA2, SOA3, SOG1, SOG2, SOG3
!
!  Module Variables:
!  ============================================================================
!  (1 ) ANTH_BLKC        (REAL*8 ) : BC anthropogenic emissions       [kg C ]
!  (2 ) ANTH_ORGC        (REAL*8 ) : OC anthropogenic emissions       [kg C ]
!  (3 ) APROD            (REAL*8 ) : Aerosol mass ratio               
!  (4 ) BCCONV           (REAL*8 ) : Hydrophilic BC from Hydrophobic  [kg C ]
!  (5 ) BIOB_BLKC        (REAL*8 ) : BC biomass emissions             [kg C ]
!  (6 ) BIOB_ORGC        (REAL*8 ) : OC biomass emissions             [kg C ]
!  (7 ) BIOF_BLKC        (REAL*8 ) : BC biofuel emissions             [kg C ]
!  (8 ) BIOF_ORGC        (REAL*8 ) : OC biofuel emissions             [kg C ]
!  (9 ) BIOG_ALPH        (REAL*8 ) : A-PINENE biogenic emissions      [kg]
!  (10) BIOG_LIMO        (REAL*8 ) : LIMONENE biogenic emissions      [kg]
!  (11) BIOG_ALCO        (REAL*8 ) : ALCOHOL biogenic emissions       [kg]
!  (12) BIOG_TERP        (REAL*8 ) : TERPENE biogenic emissions       [kg]
!  (13) BIOG_SESQ        (REAL*8 ) : SESQTERPENE biogenic emissions   [kg]
!  (14) DIUR_ORVC        (REAL*8 ) : Diurnal variation of NVOC        [kg C]
!  (15) DRYBCPI          (INTEGER) : Index for BCPI in drydep array 
!  (16) DRYOCPI          (INTEGER) : Index for OCPI in drydep array 
!  (17) DRYBCPO          (INTEGER) : Index for BCPO in drydep array 
!  (18) DRYOCPO          (INTEGER) : Index for OCPO in drydep array 
!  (19) DRYALPH          (INTEGER) : Index for ALPH in drydep array 
!  (20) DRYLIMO          (INTEGER) : Index for LIMO in drydep array 
!  (21) DRYALCO          (INTEGER) : Index for ALCO in drydep array 
!  (22) DRYSOG1          (INTEGER) : Index for SOG1 in drydep array 
!  (23) DRYSOG2          (INTEGER) : Index for SOG2 in drydep array 
!  (24) DRYSOG3          (INTEGER) : Index for SOG3 in drydep array 
!  (25) DRYSOA1          (INTEGER) : Index for SOA1 in drydep array 
!  (26) DRYSOA2          (INTEGER) : Index for SOA2 in drydep array 
!  (27) DRYSOA3          (INTEGER) : Index for SOA3 in drydep array 
!  (28) EF_BLKC          (REAL*8 ) : Emission factors for BC          [kg/kg]
!  (29) EF_ORGC          (REAL*8 ) : Emission factors for OC          [kg/kg]
!  (30) GEIA_ORVC        (REAL*8 ) : NVOC emissions from GEIA         [kg C ]
!  (31) I1_NA            (REAL*8 ) : Starting lon index for N. America region
!  (32) I2_NA            (REAL*8 ) : Ending   lon index for N. America region
!  (33) J1_NA            (REAL*8 ) : Starting lat index for N. America region
!  (34) J2_NA            (REAL*8 ) : Ending   lat index for N. America region
!  (35) GPROD            (REAL*8 ) : Gas mass ratio                   
!  (36) MHC              (INTEGER) : Number of VOC clasees
!  (37) NDAYS            (INTEGER) : Array w/ number of days per month
!  (38) NPROD            (INTEGER) : Number of products by oxdation
!  (39) OCCONV           (REAL*8 ) : Hydrophilic OC from Hydrophobic  [kg C ]
!  (40) ORVC_SESQ        (REAL*8 ) : SESQTERPENE concentration        [kg]
!  (41) ORVC_TERP        (REAL*8 ) : MONOTERPENES concentration       [kg]
!  (42) TCOSZ            (REAL*8 ) : Summing array for SUNCOS
!  (43) TERP_ORGC        (REAL*8 ) : Lumped terpene emissions         [kg C ]
!  (44) SMALLNUM         (REAL*8 ) : A small positive number
!  (45) USE_MONTHLY_ANTH (LOGICAL) : Toggles monthly or annual anthro emissions
!  (46) USE_MONTHLY_BIOB (LOGICAL) : Toggles monthly or annual biomass emiss.
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMCARBON         : Driver program for carbon aerosol chemistry
!  (2 ) CHEM_BCPO          : Chemistry routine for hydrophobic BC (aka EC)
!  (3 ) CHEM_BCPI          : Chemistry routine for hydrophilic BC (aka EC)
!  (4 ) CHEM_OCPO          : Chemistry routine for hydrophobic OC
!  (5 ) CHEM_OCPI          : Chemistry routine for hydrophilic OC
!  (6 ) SOA_CHEMISTRY      : Driver routine for SOA chemistry
!  (7 ) SOA_EQUIL          : Determines mass of SOA
!  (8 ) ZEROIN             : Finds root of an equation by bisection method
!  (9 ) SOA_PARA           : Gets SOA yield parameters
!  (10) CHEM_NVOC          : Computes oxidation of HC by O3, OH, NO3
!  (11) SOA_PARTITION      : Partitions mass of SOA gas & aerosol tracers
!  (12) SOA_LUMP           : Returns organic gas & aerosol back to STT
!  (13) SOA_DEPO           : Performs dry deposition of SOA tracers & species
!  (14) EMISSCARBON        : Driver routine for carbon aerosol emissions
!  (15) BIOGENIC_OC        : Computes biogenic OC [each time step]
!  (16) ANTHRO_CARB_TBOND  : Computes anthropogenic OC/EC [annual data]
!  (17) ANTHRO_CARB_COOKE  : Computes anthropogenic OC/EC [monthly data]
!  (18) BIOMASS_CARB_TBOND : Computes biomass burning OC/EC [annual data]
!  (19) BIOMASS_CARB_GEOS  : Computes biomass burning OC/EC [monthly data]
!  (20) EMITHIGH           : Computes complete mixing of emission within PBL
!  (21) OHNO3TIME          : Computes the sum of the cosine of SZA
!  (22) GET_OH             : Returns monthly-mean OH conc.  at grid box (I,J,L)
!  (23) GET_NO3            : Returns monthly-mean O3 conc.  at grid box (I,J,L)
!  (24) GET_O3             : Returns monthly-mean NO3 conc. at grid box (I,J,L)
!  (25) INIT_CARBON        : Allocates and zeroes all module arrays
!  (26) CLEANUP_CARBON     : Deallocates all module arrays
!
!  NOTE: Choose either (16) or (17) for ANTHROPOGENIC emission
!        Choose either (18) or (19) for BIOMASS BURNING emission.
!
!  GEOS-CHEM modules referenced by carbon_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f      : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f        : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f       : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f  : Module containing GEOS-CHEM data & met field dirs
!  (5 ) drydep_mod.f     : Module containing routines for dry deposition
!  (6 ) error_mod.f      : Module containing I/O error and NaN check routines
!  (7 ) global_no3_mod.f : Module containing routines to read 3-D NO3 field
!  (8 ) global_oh_mod.f  : Module containing routines to read 3-D OH  field
!  (9 ) global_o3_mod.f  : Module containing routines to read 3-D O3  field
!  (10) grid_mod.f       : Module containing horizontal grid information
!  (11) logical_mod.f    : Module containing GEOS-CHEM logical switches
!  (12) pressure_mod.f   : Module containing routines to compute P(I,J,L)
!  (13) time_mod.f       : Module containing routines for computing time & date
!  (14) tracer_mod.f     : Module containing GEOS-CHEM tracer array STT etc. 
!  (15) tracerid_mod.f   : Module containing pointers to tracers & emissions
!  (16) transfer_mod.f   : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Added code from the Caltech group for SOA chemistry (rjp, bmy, 7/15/04)
!  (2 ) Now references "directory_mod.f", "logical_mod.f", "tracer_mod.f".
!        (bmy, 7/20/04)
!  (3 ) Now read data from carbon_200411/ subdir of DATA_DIR.  Also added
!        some extra debug output.  Now read T. Bond yearly emissions as 
!        default, but overwrite N. America with the monthly Cooke/RJP 
!        emissions.  Added module variables I1_NA, I2_NA, J1_NA, J2_NA.
!        (rjp, bmy, 12/1/04)
!  (4 ) Now can read seasonal or interannual BCPO, OCPO biomass emissions.
!        Also parallelize loop in OHNO3TIME. (rjp, bmy, 1/18/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "carbon_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMCARBON
      PUBLIC :: EMISSCARBON
      PUBLIC :: CLEANUP_CARBON

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      LOGICAL             :: USE_MONTHLY_BIOB = .TRUE.

      INTEGER, PARAMETER  :: MHC              = 5  
      INTEGER, PARAMETER  :: NPROD            = 3  
      INTEGER             :: DRYBCPI, DRYOCPI, DRYBCPO, DRYOCPO
      INTEGER             :: DRYALPH, DRYLIMO, DRYALCO
      INTEGER             :: DRYSOG1, DRYSOG2, DRYSOG3
      INTEGER             :: DRYSOA1, DRYSOA2, DRYSOA3
      INTEGER             :: I1_NA,   J1_NA
      INTEGER             :: I2_NA,   J2_NA
      REAL*8,  PARAMETER  :: SMALLNUM         = 1d-20

      ! Arrays
      REAL*8, ALLOCATABLE :: ANTH_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: ANTH_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOB_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOB_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOF_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOF_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: EF_BLKC(:,:)
      REAL*8, ALLOCATABLE :: EF_ORGC(:,:)
      REAL*8, ALLOCATABLE :: TERP_ORGC(:,:)
      REAL*8, ALLOCATABLE :: BCCONV(:,:,:)
      REAL*8, ALLOCATABLE :: OCCONV(:,:,:)
      REAL*8, ALLOCATABLE :: BIOG_ALPH(:,:)
      REAL*8, ALLOCATABLE :: BIOG_LIMO(:,:)
      REAL*8, ALLOCATABLE :: BIOG_ALCO(:,:)
      REAL*8, ALLOCATABLE :: BIOG_TERP(:,:)
      REAL*8, ALLOCATABLE :: BIOG_SESQ(:,:)    
      REAL*8, ALLOCATABLE :: DIUR_ORVC(:,:)
      REAL*8, ALLOCATABLE :: GEIA_ORVC(:,:)
      REAL*8, ALLOCATABLE :: TCOSZ(:,:)
      REAL*8, ALLOCATABLE :: ORVC_SESQ(:,:,:)
      REAL*8, ALLOCATABLE :: ORVC_TERP(:,:,:)
      REAL*8, ALLOCATABLE :: GPROD(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: APROD(:,:,:,:,:)

      ! Days per month (based on 1998)
      INTEGER             :: NDAYS(12) = (/ 31, 28, 31, 30, 31, 30, 
     &                                      31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCARBON
!
!******************************************************************************
!  Subroutine CHEMCARBON is the interface between the GEOS-CHEM main 
!  program and the carbon aerosol chemistry routines that calculates
!  dry deposition, chemical conversion between hydrophilic and 
!  hydrophobic, and SOA production. (rjp, bmy, 4/1/04, 7/20/04)
!
!  NOTES:
!  (1 ) Added code from the Caltech group for SOA chemistry.  Also now 
!        reference "global_oh_mod.f", "global_o3_mod.f", "global_no3_mod.f".
!        (rjp, bmy, 7/8/04)
!  (2 ) Now reference LSOA and LEMIS from CMN_SETUP.  Now only call OHNO3TIME
!        if it hasn't been done before w/in EMISSCARBON. (rjp, bmy, 7/15/04)
!  (3 ) Now reference LSOA, LEMIS, LPRT from "logical_mod.f".  Now reference
!        STT and ITS_AN_AEROSOL_SIM from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,     ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD, ONLY : GET_GLOBAL_NO3
      USE GLOBAL_O3_MOD,  ONLY : GET_GLOBAL_O3
      USE LOGICAL_MOD,    ONLY : LSOA, LEMIS, LPRT
      USE TIME_MOD,       ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: FIRSTCHEM = .TRUE.
      INTEGER       :: N, THISMONTH

      !=================================================================
      ! CHEMCARBON begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_CARBON

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'BCPI' )
                  DRYBCPI = N
               CASE ( 'OCPI' )
                  DRYOCPI = N
               CASE ( 'BCPO' )
                  DRYBCPO = N
               CASE ( 'OCPO' )
                  DRYOCPO = N
               CASE ( 'ALPH' )
                  DRYALPH = N
               CASE ( 'LIMO' )
                  DRYLIMO = N
               CASE ( 'ALCO' )
                  DRYALCO = N
               CASE ( 'SOG1' )
                  DRYSOG1 = N
               CASE ( 'SOG2' )
                  DRYSOG2 = N
               CASE ( 'SOG3' )
                  DRYSOG3 = N
               CASE ( 'SOA1' )
                  DRYSOA1 = N
               CASE ( 'SOA2' )
                  DRYSOA2 = N
               CASE ( 'SOA3' )
                  DRYSOA3 = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Do chemistry for carbon aerosol tracers 
      !=================================================================

      ! Chemistry for hydrophobic BC
      IF ( IDTBCPO > 0 ) THEN
         CALL CHEM_BCPO( STT(:,:,:,IDTBCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPO' )
      ENDIF

      ! Chemistry for hydrophilic BC
      IF ( IDTBCPI > 0 ) THEN
         CALL CHEM_BCPI( STT(:,:,:,IDTBCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPI' )
      ENDIF

      ! Chemistry for hydrophobic OC
      IF ( IDTOCPO > 0 ) THEN
         CALL CHEM_OCPO( STT(:,:,:,IDTOCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPO' )
      ENDIF

      ! Chemistry for hydrophilic OC
      IF ( IDTOCPI > 0 ) THEN 
         CALL CHEM_OCPI( STT(:,:,:,IDTOCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPI' )
      ENDIF

      !=================================================================
      ! Do chemistry for secondary organic aerosols 
      !=================================================================
      IF ( LSOA ) THEN

         ! Read offline OH, NO3, O3 fields from disk
         IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Current month
            THISMONTH = GET_MONTH()

            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_OH(  THISMONTH )
               CALL GET_GLOBAL_NO3( THISMONTH )
               CALL GET_GLOBAL_O3(  THISMONTH )
            ENDIF

            ! Compute time scaling arrays for offline OH, NO3
            ! but only if it hasn't been done in EMISSCARBON
            IF ( LSOA .and. ( .not. LEMIS ) ) THEN
               CALL OHNO3TIME
               IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARB: a OHNO3TIME' )
            ENDIF
         ENDIF

         ! Compute SOA chemistry
         CALL SOA_CHEMISTRY
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_SOA' )

      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMCARBON

!-----------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPO converts hydrophobic BC to hydrophilic BC and
!  calculates the dry deposition of hydrophobic BC. (rjp, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophobic BC tracer 
!
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_BC 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTBCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,       J,   L
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)     
      REAL*8                 :: DTCHEM, FLUX, KBC, FREQ
      REAL*8                 :: TC0,    CNEW, RKT, AREA_CM2, BL_FRAC
      REAL*8,  PARAMETER     :: BC_LIFE = 1.15D0


      !=================================================================
      ! CHEM_BCPO begins here!
      !=================================================================

      ! Return if BCPO isn't defined
      IF ( IDTBCPO == 0 .or. DRYBCPO == 0 ) RETURN

      ! Initialize
      KBC    = 1.D0 / ( 86400d0 * BC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !
      ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5          
      ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial BC mass [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq
         FREQ = 0d0

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic at all levels
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! BC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYBCPO) * BL_FRAC

         ENDIF

         ! Amount of BCPO left after chemistry and drydep [kg]
         RKT  = ( KBC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of BCPO converted to BCPI [kg/timestep]
         BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]  
             ! XNUMOL is the ratio [molec tracer/kg tracer]   
             FLUX     = TC0 - CNEW - BCCONV(I,J,L) 
             FLUX     = FLUX * XNUMOL(IDTBCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-philic BC from H_phobic BC [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
             AD07_BC(I,J,L) = AD07_BC(I,J,L) + BCCONV(I,J,L)
         ENDIF

         ! Store new concentration back into tracer array
         TC(I,J,L) = CNEW
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPO,1) = AD44(I,J,DRYBCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPO

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic BC.
!  (rjp, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTBCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC, AREA_CM2
      REAL*8                 :: TC0,    CNEW, CCV,     FREQ
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_BCPI begins here!
      !=================================================================

      ! Return if BCPI isn't defined
      IF ( IDTBCPI == 0 .or. DRYBCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for ND44 diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, FREQ, BL_FRAC, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic BC [kg]
         TC0 = TC(I,J,L)

         ! H-philic BC that used to be H-phobic BC [kg]
         CCV = BCCONV(I,J,L)
         
         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency
            FREQ = DEPSAV(I,J,DRYBCPI) * BL_FRAC
            
            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT)) 
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! Comment out for now
            !CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !     + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep flux [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 .and. FREQ > 0d0 ) THEN
  
               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [molec/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTBCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, omit the exponential to save on clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Save new concentration of H-philic IC in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero out the BCCONV array for the next iteration
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0.d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPI,1) = AD44(I,J,DRYBCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPI

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_OCPO converts hydrophobic OC to hydrophilic OC and
!  calculates the dry deposition of hydrophobic OC. (rjp, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophobic OC tracer [kg]
!
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_OC 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTOCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: DTCHEM, FLUX, KOC,  BL_FRAC
      REAL*8                 :: TC0,    FREQ, CNEW, RKT, AREA_CM2
      REAL*8,  PARAMETER     :: OC_LIFE = 1.15D0

      !=================================================================
      ! CHEM_OCPO begins here!
      !=================================================================

      ! Return if OCPO isn't defined
      IF ( IDTOCPO == 0 .or. DRYOCPO == 0 ) RETURN

      ! Initialize
      KOC    = 1.D0 / ( 86400d0 * OC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !    Hydrophobic --> Hydrophilic,  k  = 1.0e-5          
      !    Aerosols are dry-deposited,   kd = DEPSAV (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial OC [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq 
         FREQ = 0d0

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! OC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYOCPO) * BL_FRAC

         ENDIF

         ! Amount of OCPO left after chemistry and drydep [kg]
         RKT  = ( KOC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of OCPO converted to OCPI [kg/timestep]
         OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
             ! XNUMOL is the ratio [molec tracer/kg tracer]     
             FLUX     = TC0 - CNEW - OCCONV(I,J,L)
             FLUX     = FLUX * XNUMOL(IDTOCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-Philic OC from H-phobic [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
            AD07_OC(I,J,L) = AD07_OC(I,J,L) + OCCONV(I,J,L) 
         ENDIF

         ! Store modified OC concentration back in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPO,1) = AD44(I,J,DRYOCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF   

      ! Return to calling program
      END SUBROUTINE CHEM_OCPO

!-----------------------------------------------------------------------

      SUBROUTINE CHEM_OCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic OC.
!  (rjp, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!  (1 ) Remove reference to "CMN", it's obsolete (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTOCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC
      REAL*8                 :: TC0, CNEW, CCV, FREQ, AREA_CM2
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_OCPI begins here!
      !=================================================================
      IF ( IDTOCPI == 0 .or. DRYOCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic OC [kg]
         TC0 = TC(I,J,L)

         ! H-philic OC that used to be H-phobic OC [kg]
         CCV = OCCONV(I,J,L)

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DRYOCPI) * BL_FRAC

            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT))
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !       + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTOCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero OCCONV array for next timestep
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPI,1) = AD44(I,J,DRYOCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF    

      ! Return to calling program
      END SUBROUTINE CHEM_OCPI

!------------------------------------------------------------------------------

      SUBROUTINE SOA_CHEMISTRY
!
!******************************************************************************
!  Subroutine SOA_CHEMISTRY performs SOA formation.  This code is from the
!  Caltech group (Hong Liao, Serena Chung, et al) and was modified for 
!  GEOS-CHEM. (rjp, bmy, 7/8/04, 7/20/04)
!
!  Procedure:
!  ============================================================================
!  (1 ) Read in NO3, O3, OH in CHEM_SOA
!  (2 ) Scales these fields using OHNO3TIME in sulfate_mod.f (see GET_OH)
!  (3 ) Calculate reaction rates (Serena's OCHEMPARAETER)
!  (4 ) CALCULATE DELHC
!  (5 ) get T0M gas products
!  (6 ) equilibrium calculation
!
!  There are total of 37 tracers considered in this routine:
!
!  4 classes of primary carbonaceous aerosols:  
!     BCPI = Hydrophilic black carbon
!     OCPI = Hydrophilic organic carbon
!     BCPO = Hydrophobic black carbon
!     OCPO = Hydrophobic organic carbon
!
!  5 reactive biogenic hydrocarbon groups (NVOC):
!     ALPH = a-pinene, b-pinene, sabinene, carene, terpenoid ketones
!     LIMO = limonene
!     TERP = a-terpinene, r-terpinene, terpinolene
!     ALCO = myrcene, terpenoid alcohols, ocimene
!     SESQ = sesquiterpenes
!
!  NOTE: TERP and SESQ are not tracers because of their high reactivity
!
!  28 organic oxidation products by O3+OH and NO3: 
!     6 ( 3 gases + 3 aerosols ) from each of first four NVOC  = 24 
!     4 ( 2 gases + 2 aerosols ) from sesquiterpenes oxidation = 4
!
!  NOTE: We aggregate these into 6 tracers according to HC classes
!     SOG1 = lump of gas products of first three (ALPH+LIMO+TERP) HC oxidation.
!     SOG2 = gas product of ALCO oxidation
!     SOG3 = gas product of SESQ oxidation 
!     SOA1 = lump of aerosol products of first 3 (ALPH+LIMO+TERP) HC oxidation.
!     SOA2 = aerosol product of ALCO oxidation
!     SOA3 = aerosol product of SESQ oxidation 
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE DAO_MOD,      ONLY : T, AD, AIRVOL, SUNCOS
      USE DIAG_MOD,     ONLY : AD07_HC 
      USE TRACER_MOD,   ONLY : STT 
      USE TRACERID_MOD
      USE TIME_MOD,     ONLY : GET_TS_CHEM, GET_MONTH

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER :: I,        J,    L,     N,     JHC,  IPR
      REAL*8  :: RTEMP,    VOL,  FAC,   MPOC,  MNEW, MSOA_OLD
      REAL*8  :: MPRODUCT, CAIR, LOWER, UPPER, TOL,  VALUE
      REAL*8  :: KO3(MHC), KOH(MHC), KNO3(MHC)
      REAL*8  :: ALPHA(NPROD,MHC),   KOM(NPROD,MHC)
      REAL*8  :: GM0(NPROD,MHC),     AM0(NPROD,MHC)
      REAL*8  :: ORG_AER(NPROD,MHC), ORG_GAS(NPROD,MHC)

      !=================================================================
      ! SOA_CHEMISTRY begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,        J,        L,     JHC,   IPR,   GM0,  AM0  )
!$OMP+PRIVATE( VOL,      FAC,      RTEMP, KO3,   KOH,   KNO3, CAIR )
!$OMP+PRIVATE( MPRODUCT, MSOA_OLD, VALUE, UPPER, LOWER, MNEW, TOL  )
!$OMP+PRIVATE( ORG_AER,  ORG_GAS,  ALPHA, KOM,   MPOC              )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Volume of grid box [m3]
         VOL   = AIRVOL(I,J,L)    

         ! conversion factor from kg to ug/m3
         FAC   = 1.D9 / VOL       

         ! air conc. in kg/m3
         CAIR  = AD(I,J,L) / VOL  

         ! Temperature [K]
         RTEMP = T(I,J,L)         

         ! Get SOA yield parameters
         CALL SOA_PARA( RTEMP, KO3, KOH, KNO3, ALPHA, KOM )

         ! Partition mass of gas & aerosol tracers 
         ! according to 5 VOC classes & 3 oxidants
         CALL SOA_PARTITION( I, J, L, GM0, AM0 )

         ! Compute oxidation of hydrocarbons by O3, OH, NO3
         CALL CHEM_NVOC( I, J, L, KO3, KOH, KNO3, ALPHA, GM0 )

         !==============================================================
         ! Equilibrium calculation between GAS (SOG) and Aerosol (SOA)
         !==============================================================

         ! Total VOC oxidation products (gas and aerosol) [kg]
         MPRODUCT = 0.D0   

         ! Initialize aerosol-only total [kg]
         MSOA_OLD = 0.D0  
     
         ! Compute total VOC and aerosol-only total
         DO JHC = 1, MHC
         DO IPR = 1, NPROD
            MPRODUCT = MPRODUCT + GM0(IPR,JHC) + AM0(IPR,JHC)
            MSOA_OLD = MSOA_OLD + AM0(IPR,JHC)
         ENDDO
         ENDDO

         ! No need to proceed because there is no
         ! products that need to be re-equilibrated.
         IF ( ( MPRODUCT / AD(I,J,L) ) <= 29.D-18 ) CYCLE

         ! Individual SOA's; units: [ug/m3]
         DO JHC = 1, MHC
         DO IPR = 1, NPROD
            ORG_GAS(IPR,JHC) = GM0(IPR,JHC) * FAC
            ORG_AER(IPR,JHC) = AM0(IPR,JHC) * FAC
         ENDDO
         ENDDO
                
         ! Primary organic aerosol concentrations [ug/m3]
         ! We carry carbon mass only in STT arry and here multiply 1.4
         ! to account for the mass of other chemical component attached.
         MPOC = (STT(I,J,L,IDTOCPI) + STT(I,J,L,IDTOCPO)) * FAC
         MPOC = MPOC * 1.4D0

         !==============================================================
         ! Solve for MNEW by solving for SOA=0
         !==============================================================
         IF ( ( MPOC / ( CAIR*1.D9 ) ) <= 2.1D-18 ) THEN
            VALUE = 0.D0
            UPPER = 0.D0

            DO JHC = 1, MHC
            DO IPR = 1, NPROD
               VALUE = VALUE + KOM(IPR,JHC) *
     &                 (ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC))

               UPPER = UPPER + ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC)
            ENDDO
            ENDDO

            IF ( VALUE <= 1.D0 ) THEN
               MNEW  = 0.D0
            ELSE
               LOWER = 1.D-18 * ( CAIR * 1.D9 )
               TOL   = 1.D-18 
               MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)
            ENDIF

         ELSE

            UPPER = MPOC
            DO JHC = 1, MHC
            DO IPR = 1, NPROD
               UPPER = UPPER + ORG_GAS(IPR,JHC) + ORG_AER(IPR,JHC)
            ENDDO
            ENDDO

            LOWER = MPOC
            TOL   = 1.D-9*MPOC
            MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)

         ENDIF

         !==============================================================
         ! Equilibrium partitioning into new gas and aerosol 
         ! concentrations for individual contributions of SOA
         !==============================================================
         IF ( MNEW > 0.D0 ) THEN
             DO JHC = 1, MHC
             DO IPR = 1, NPROD

                ORG_AER(IPR,JHC) = KOM(IPR,JHC)*MNEW /
     &                            (1.D0 + KOM(IPR,JHC) * MNEW ) *
     &                            (ORG_AER(IPR,JHC) + ORG_GAS(IPR,JHC))

                IF ( KOM(IPR,JHC).NE.0D0 ) THEN
                    ORG_GAS(IPR,JHC) = ORG_AER(IPR,JHC) * 1.D8 /
     &                                 ( KOM(IPR,JHC) * MNEW * 1.D8 )
                ELSE
                   ORG_GAS(IPR,JHC) = 0.D0
                ENDIF

             ENDDO
             ENDDO

             ! STORE PRODUCT INTO T0M 
             DO JHC = 1, MHC
             DO IPR = 1, NPROD
                GM0(IPR,JHC) = ORG_GAS(IPR,JHC) / FAC
                AM0(IPR,JHC) = ORG_AER(IPR,JHC) / FAC
             ENDDO
             ENDDO

         !==============================================================
         ! Mnew=0.D0, all SOA evaporates to the gas-phase
         !==============================================================
         ELSE

             DO JHC = 1, MHC
             DO IPR = 1, NPROD
                GM0(IPR,JHC) = GM0(IPR,JHC) + AM0(IPR,JHC)
                AM0(IPR,JHC) = 1.D-18 * AD(I,J,L)
             ENDDO
             ENDDO

         ENDIF

         ! Lump SOA
         CALL SOA_LUMP( I, J, L, GM0, AM0 )
       
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================       
      ! Calculate dry-deposition
      !=================================================================
      CALL SOA_DEPO( STT(:,:,:,IDTALPH), DRYALPH, IDTALPH )
      CALL SOA_DEPO( STT(:,:,:,IDTLIMO), DRYLIMO, IDTLIMO )
      CALL SOA_DEPO( STT(:,:,:,IDTALCO), DRYALCO, IDTALCO )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG1), DRYSOG1, IDTSOG1 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG2), DRYSOG2, IDTSOG2 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOG3), DRYSOG3, IDTSOG3 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA1), DRYSOA1, IDTSOA1 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA2), DRYSOA2, IDTSOA2 )
      CALL SOA_DEPO( STT(:,:,:,IDTSOA3), DRYSOA3, IDTSOA3 )

      ! Return to calling program
      END SUBROUTINE SOA_CHEMISTRY

!------------------------------------------------------------------------------

      FUNCTION SOA_EQUIL( MASS, MPOC, AEROSOL, GAS, KOM ) 
     &         RESULT( SOA_MASS )
!
!******************************************************************************
!  Subroutine SOA_EQUIL solves SOAeqn=0 to determine Mnew (= mass)
!  See Eqn (27) on page 70 of notes.  Originally written by Serena Chung at
!  Caltech, and modified for inclusion into GEOS-CHEM. (rjp, bmy, 7/8/04)
!
!  This version does NOT assume that the gas and aerosol phases are in 
!  equilibrium before chemistry; therefore, gas phase concentrations are 
!  needed explicitly.  The gas and aerosol phases are assumed to be in 
!  equilibrium after chemistry.
! 
!  Note: Unlike FUNCTION SOA, this function assumes no reactions.  It only 
!  considers the partitioning of existing products of VOC oxidation.
!
!  HC_JHC + OXID_IOXID - > 
!    alpha(1,IOXID,JHC) [SOAprod_gas(1,IOXID,JHC)+SOAprod(1,IOXID,JHC)]+
!    alpha(2,IOXID,JHC) [SOAprod_gas(2,IOXID,JHC)+SOAprod(2,IOXID,JHC)]
!
!  SOAprod_gas(IPR,IOXID,JHC) <--> SOAprod(IPR,IOXID,JHC)   
!                                           (aerosol phase)
!
!  w/ equilibrium partitioning:
!
!                                   SOAprod(IPR,IOXID,JHC)
!    SOAprod_gas(IPR,IOXID,JHC) = ------------------------
!                                     Kom(IPR,IOXID,JHC)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MASS    (REAL*8) : Pre-existing aerosol mass                [ug/m3]
!  (2 ) MPOC    (REAL*8) : POA Mass                                 [ug/m3]
!  (3 ) AEROSOL (REAL*8) : Aerosol concentration                    [ug/m3]
!  (4 ) GAS     (REAL*8) : Gas-phase concentration                  [ug/m3]
!  (5 ) KOM     (REAL*8) : Equilibrium gas-aerosol partition coeff. [m3/ug]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: MASS, MPOC         
      REAL*8, INTENT(IN) :: AEROSOL(NPROD,MHC)
      REAL*8, INTENT(IN) :: GAS(NPROD,MHC)
      REAL*8, INTENT(IN) :: KOM(NPROD,MHC)

      ! Local variables
      INTEGER            :: JHC,   IPR
      REAL*8             :: VALUE, SOA_MASS

      !=================================================================
      ! SOA_EQUIL begins here!
      !=================================================================

      ! Equation (39) on page 139 of notes:
      VALUE = 0.D0

      ! 
      DO JHC = 1, MHC
      DO IPR = 1, NPROD
         VALUE = VALUE + KOM(IPR,JHC)                        / 
     &                   ( 1.D0 + KOM(IPR,JHC) * MASS      ) * 
     &                   ( GAS(IPR,JHC) + AEROSOL(IPR,JHC) )
      ENDDO                
      ENDDO   

      ! Compute SOA mass
      SOA_MASS = VALUE + ( 1.D5 * MPOC ) / ( 1.D5 * MASS ) - 1.0D0

      ! Return to calling program
      END FUNCTION SOA_EQUIL

!------------------------------------------------------------------------------

      FUNCTION ZEROIN(AX,BX,TOL,MPOC,AEROSOL,GAS,KOM) RESULT( MNEW )
!
!******************************************************************************
! NOTE: This function may be problematic -- it uses GOTO's, which are not
! good for parallelization. (bmy, 7/8/04)
!
! shc I got this code from http://www.netlib.org
!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0d0)
!
!
!  output..
!
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!  NOTES:
!  (1 ) Change dabs to ABS and dsign to SIGN, in order to avoid conflicts
!        with intrinsic function names on the PGI compiler. (bmy, 12/2/04)
!******************************************************************************
!
      real*8, intent(in) :: ax,bx,tol
      REAL*8, INTENT(IN) :: Mpoc      
      REAL*8, INTENT(IN) :: Aerosol(NPROD,MHC), Gas(NPROD,MHC)
      REAL*8, INTENT(IN) :: Kom(NPROD,MHC)

      !local variables
      real*8             :: MNEW
      real*8             :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
c
c  compute eps, the relative machine precision
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d0) go to 10
c
c initialization
c
      a  = ax
      b  = bx
      fa = SOA_equil( A, MPOC, Aerosol, GAS, Kom )
      fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom ) 
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d

   30 if (ABS(fc) .ge. ABS(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*ABS(b) + 0.5d0*tol
      xm = 0.5D0*(c - b)
      if (ABS(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (ABS(e) .lt. tol1) go to 70
      if (ABS(fa) .le. ABS(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = ABS(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - ABS(tol1*q))) go to 70
      if (p .ge. ABS(0.5d0*e*q)) go to 70

      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (ABS(d) .gt. tol1) b = b + d
      if (ABS(d) .le. tol1) b = b + SIGN(tol1, xm)

      fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom ) 
      if ((fb*(fc/ABS(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 MNEW = b

      ! Return to calling program
      END FUNCTION ZEROIN

!------------------------------------------------------------------------------

      FUNCTION RTBIS( X1,   X2,      XACC, 
     &                MPOC, AEROSOL, GAS, KOM ) RESULT( ROOT )
!
!******************************************************************************
!  Function RTBIS finds the root of the function SOA_EQUIL via the bisection
!  method.  Original algorithm from "Numerical Recipes" by Press et al, 
!  Cambridge UP, 1986.  Modified for inclusion into GEOS-CHEM. (bmy, 7/8/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X1      (REAL*8) : Endpoint #1 
!  (2 ) X2      (REAL*8) : Endpoint #2
!  (3 ) XACC    (REAL*8) : Desired accuracy of solution
!  (4 ) MPOC    (REAL*8) : POA Mass                                 [ug/m3]
!  (5 ) AEROSOL (REAL*8) : Aerosol concentration                    [ug/m3]
!  (6 ) GAS     (REAL*8) : Gas-phase concentration                  [ug/m3]
!  (7 ) KOM     (REAL*8) : Equilibrium gas-aerosol partition coeff. [m3/ug]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      REAL*8, INTENT(IN) :: X1, X2, XACC, MPOC
      REAL*8, INTENT(IN) :: AEROSOL(NPROD,MHC)
      REAL*8, INTENT(IN) :: GAS(NPROD,MHC)
      REAL*8, INTENT(IN) :: KOM(NPROD,MHC)

      ! Local variables
      INTEGER, PARAMETER :: JMAX = 100
      INTEGER            :: J
      REAL*8             :: ROOT, DX, F, FMID, XMID

      !=================================================================
      ! RTBIS begins here!
      !=================================================================

      ! Compute value of function SOA_EQUIL at endpoints
      FMID = SOA_EQUIL( X2, MPOC, AEROSOL, GAS, KOM )
      F    = SOA_EQUIL( X1, MPOC, AEROSOL, GAS, KOM )

      ! Test if we are bracketing a root
      IF ( F * FMID >= 0d0 ) THEN
         CALL ERROR_STOP( 'Root must be bracketed!', 
     &                    'RTBIS ("carbon_mod.f")' )
      ENDIF

      ! Set initial root and interval
      IF ( F < 0d0 ) THEN
         ROOT = X1
         DX   = X2 - X1
      ELSE
         ROOT = X2
         DX   = X1 - X2
      ENDIF

      ! Loop until max iteration count
      DO J = 1, JMAX

         ! Halve the existing interval
         DX   = DX * 0.5D0

         ! Compute midpoint of new interval
         XMID = ROOT + DX

         ! Compute value of function SOA_EQUIL at new midpoint
         FMID = SOA_EQUIL( XMID, MPOC, AEROSOL, GAS, KOM )

         ! We have found the root!
         IF ( FMID <= 0D0 ) ROOT = XMID

         ! We have reached the tolerance, so return
         IF ( ABS( DX ) < XACC .OR. FMID == 0.D0 ) RETURN
      ENDDO

      ! Stop with error condition
      CALL ERROR_STOP( 'Too many bisections!', 
     &                 'RTBIS ("carbon_mod.f")' )

      ! Return to calling program
      END FUNCTION RTBIS

!------------------------------------------------------------------------------

      SUBROUTINE SOA_PARA( TEMP, KO3, KOH, KNO3, RALPHA, KOM )
!
!******************************************************************************
!  Subroutine SOA_PARA gves mass-based stoichiometric coefficients for semi-
!  volatile products from the oxidation of hydrocarbons.  It calculates 
!  secondary organic aerosol yield parameters.  Temperature effects are
!  included.  Original code from the CALTECH group and modified for inclusion
!  to GEOS-CHEM. (rjp, bmy, 7/8/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TEMP   (REAL*8) : Temperature [k]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) KO3    (REAL*8) : Rxn rate for HC oxidation by O3         [cm3/molec/s]
!  (3 ) KOH    (REAL*8) : Rxn rate for HC oxidation by OH         [cm3/molec/s]
!  (4 ) KNO3   (REAL*8) : Rxn rate for HC oxidation by NO3        [cm3/molec/s]
!  (5 ) RALPHA (REAL*8) : Mass-based stoichiometric coefficients  [unitless]
!  (6 ) KOM    (REAL*8) : Equilibrium gas-aerosol partition coeff [m3/ug]
!
!  References:
!  ============================================================================
!  PHOTO-OXIDATION RATE CONSTANTS OF ORGANICS come from:
!   (1) Atkinson, el al., Int. J. Chem.Kinet., 27: 941-955 (1995)
!   (2) Shu and Atkinson, JGR 100: 7275-7281 (1995)
!   (3) Atkinson, J. Phys. Chem. Ref. Data 26: 215-290 (1997)   
!   (4) Some are reproduced in Table 1 of Griffin, et al., JGR 104: 3555-3567
!   (5) Chung and Seinfeld (2002)
!
!  ACTIVATION ENERGIES come from:
!   (6) Atkinson, R. (1994) Gas-Phase Tropospheric Chemistry of Organic 
!        Compounds.  J. Phys. Chem. Ref. Data, Monograph No.2, 1-216. 
!   (7) They are also reproduced in Tables B.9 and B.10 of Seinfeld and 
!        Pandis (1988).
!  
!  NOTES:
!  (1 ) Now use temporary variables TMP1, TMP2, TMP3 to pre-store the values
!        of exponential terms outside of DO-loops (bmy, 7/8/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: TEMP
      REAL*8, INTENT(OUT) :: KO3(MHC), KOH(MHC), KNO3(MHC)
      REAL*8, INTENT(OUT) :: RALPHA(NPROD,MHC), KOM(NPROD,MHC)

      ! Local variables
      INTEGER             :: IPR,  JHC,  J
      REAL*8              :: TMP1, TMP2, TMP3, OVER

      ! Activation Energy/R [K] for O3, OH, NO3 (see Refs #6-7)
      REAL*8, PARAMETER   :: ACT_O3     =  732.0d0   
      REAL*8, PARAMETER   :: ACT_OH     = -400.0d0   
      REAL*8, PARAMETER   :: ACT_NO3    = -490.0d0   

      ! Heat of vaporization (from CRC Handbook of Chemistry & Physics)
      REAL*8, PARAMETER   :: HEAT_VAPOR = 5.d3     

      ! Reciprocal reference temperatures at 298K and 310K
      REAL*8, PARAMETER   :: REF298     = 1d0 / 298d0
      REAL*8, PARAMETER   :: REF310     = 1d0 / 310d0

      !=================================================================
      ! SOA_PARA begins here!
      !=================================================================

      ! Photo-oxidation rates of O3 [cm3/molec/s] (See Refs #1-4)
      KO3(1) = 56.15d-18
      KO3(2) = 200.d-18
      KO3(3) = 7707.d-18
      KO3(4) = 422.5d-18
      KO3(5) = ( 11600.D0 + 11700.D0 ) / 2.D0 * 1.D-18

      ! Photo-oxidation rates of OH [cm3/molec/s] (See Refs #1-4)
      KOH(1) = 84.4d-12
      KOH(2) = 171.d-12
      KOH(3) = 255.d-12
      KOH(4) = 199.d-12
      KOH(5) = ( 197.d0 + 293.d0 ) / 2.d0 * 1.d-12

      ! Photo-oxidation rate of NO3 [cm3/molec/s] (See Refs #1-4)
      KNO3(1) = 6.95d-12
      KNO3(2) = 12.2d-12
      KNO3(3) = 88.7d-12
      KNO3(4) = 14.7d-12
      KNO3(5) = ( 19.d0 + 35.d0 ) / 2.d0 * 1.d-12

      !=================================================================
      ! Temperature Adjustments of KO3, KOH, KNO3
      !=================================================================

      ! Reciprocal temperature [1/K]
      OVER = 1.0d0 / TEMP

      ! Compute the exponentials once outside the DO loop
      TMP1 = EXP( ACT_O3  * ( REF298 - OVER ) )
      TMP2 = EXP( ACT_OH  * ( REF298 - OVER ) )
      TMP3 = EXP( ACT_NO3 * ( REF298 - OVER ) )

      ! Multiply photo-oxidation rates by exponential of temperature
      DO JHC = 1, MHC 
         KO3(JHC)  = KO3(JHC)  * TMP1
         KOH(JHC)  = KOH(JHC)  * TMP2
         KNO3(JHC) = KNO3(JHC) * TMP3
      ENDDO

      !=================================================================
      ! SOA YIELD PARAMETERS
      ! 
      ! Aerosol yield parameters for photooxidation of biogenic organics
      ! The data (except for C7-C10 n-carbonyls, aromatics, and higher 
      ! ketones are from: 
      !
      ! (7) Tables 1 and 2 of Griffin, et al., Geophys. Res. Lett. 
      !      26: (17)2721-2724 (1999)
      !
      ! These parameters neglect contributions of the photooxidation 
      ! by NO3. 
      !
      ! For the aromatics, the data are from
      ! (8) Odum, et al., Science 276: 96-99 (1997).
      !=================================================================

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
      RALPHA(1,1) = 0.067d0            
      RALPHA(1,2) = 0.239d0

      ! Average of TERPINENES and TERPINOLENE
      RALPHA(1,3) = 0.0685d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE      
      RALPHA(1,4) = 0.06675d0

      ! Average of BETA-CARYOPHYLLENE and ALPHA-HUMULENE      
      RALPHA(1,5) = 1.0d0

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE      
      RALPHA(2,1) = 0.35425d0
      RALPHA(2,2) = 0.363d0

      ! Average of TERPINENES and TERPINOLENE
      RALPHA(2,3) = 0.2005d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
      RALPHA(2,4) = 0.135d0

      ! Not applicable
      RALPHA(2,5) = 0.0d0

      ! Using BETA-PINENE for all species for NO3 oxidation
      ! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
      RALPHA(3,:) = 1.d0           

      !=================================================================
      ! Equilibrium gas-particle partition coefficients of 
      ! semi-volatile compounds [ug-1 m**3]
      !=================================================================
      KOM(1,1) = 0.1835d0

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
      KOM(1,2) = 0.055d0
      KOM(1,3) = 0.133d0

      ! Average of TERPINENES and TERPINOLENE
      KOM(1,4) = 0.22375d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
      KOM(1,5) = ( 0.04160d0 + 0.0501d0 ) / 2.d0

      ! Average of BETA-CARYOPHYLLENE and and ALPHA-HUMULENE
      KOM(2,1) = 0.004275d0

      ! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
      KOM(2,2) = 0.0053d0
      KOM(2,3) = 0.0035d0

      ! Average of TERPINENES and TERPINOLENE
      KOM(2,4) = 0.0082d0

      ! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
      KOM(2,5) = 0.0d0

      ! NOT APPLICABLE -- using BETA-PINENE for all species
      ! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
      KOM(3,:) = 0.0163d0
                           
      !=================================================================
      ! Temperature Adjustments of KOM
      !=================================================================

      ! Reciprocal temperature [1/K]
      OVER = 1.0D0 / TEMP

      ! Divide TEMP by 310K outside the DO loop
      TMP1 = ( TEMP / 310.D0 )

      ! Compute the heat-of-vaporization exponential term outside the DO loop
      TMP2 = EXP( HEAT_VAPOR * ( OVER - REF310 ) )

      ! Multiply KOM by the temperature and heat-of-vaporization terms
      DO JHC = 1, 5
      DO IPR = 1, 3
         KOM(IPR,JHC) = KOM(IPR,JHC) * TMP1 * TMP2
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE SOA_PARA

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NVOC( I, J, L, KO3, KOH, KNO3, ALPHA, GM0 )
!
!******************************************************************************
!  Subroutine CHEM_NVOC computes the oxidation of Hydrocarbon by O3, OH, and 
!  NO3.  This comes from the Caltech group (Hong Liao, Serena Chung, et al)
!  and was incorporated into GEOS-CHEM. (rjp, bmy, 7/6/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I     (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J     (INTEGER) : GEOS-CHEM latitude index
!  (3 ) L     (INTEGER) : GEOS-CHEM altitude index
!  (4 ) KO3   (REAL*8 ) : Rxn rate for HC oxidation by O3        [cm3/molec/s]
!  (5 ) KOH   (REAL*8 ) : Rxn rate for HC oxidation by OH        [cm3/molec/s]
!  (6 ) KNO3  (REAL*8 ) : Rxn rate for HC oxidation by NO3       [cm3/molec/s]
!  (7 ) ALPHA (REAL*8 ) : Mass-based stoichiometric coefficients [unitless]  
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) GM0   (REAL*8 ) : Gas mass for each HC and its oxidation product [kg]
!  
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 kmodules
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD
      USE TIME_MOD,     ONLY : GET_TS_CHEM, GET_MONTH

#     include "CMN_SIZE"     ! Size parameters

      INTEGER, INTENT(IN)    :: I, J, L
      REAL*8,  INTENT(IN)    :: KO3(MHC), KOH(MHC), KNO3(MHC)
      REAL*8,  INTENT(IN)    :: ALPHA(NPROD,MHC)
      REAL*8,  INTENT(INOUT) :: GM0(NPROD,MHC)

      ! LOCAL VARIABLES
      INTEGER                :: JHC, IPR
      REAL*8                 :: DELHC(NPROD), CHANGE(MHC), NMVOC(MHC)
      REAL*8                 :: OHMC, TTNO3, TTO3, DTCHEM, RK
      REAL*8                 :: OVER, DO3, DOH, DNO3

      !=================================================================
      ! CHEM_NVOC begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Get offline OH, NO3, O3 concentrations [molec/cm3]
      OHMC    = GET_OH(  I, J, L ) 
      TTNO3   = GET_NO3( I, J, L )
      TTO3    = GET_O3(  I, J, L ) 

      ! 5 classes NVOC concentrations followed by 4 primary species
      NMVOC(1) = STT(I,J,L,IDTALPH)
      NMVOC(2) = STT(I,J,L,IDTLIMO)
      NMVOC(3) = ORVC_TERP(I,J,L)
      NMVOC(4) = STT(I,J,L,IDTALCO)
      NMVOC(5) = ORVC_SESQ(I,J,L)

      ! Initialize DELHC so that the values from the previous
      ! time step are not carried over.
      DELHC(:) = 0.D0

      !=================================================================
      ! Change in NVOC concentration due to photooxidation [kg]
      !=================================================================
      DO JHC = 1, MHC
         RK          = KO3(JHC)*TTO3 + KOH(JHC)*OHMC + KNO3(JHC)*TTNO3
         CHANGE(JHC) = NMVOC(JHC) * ( 1.D0 - DEXP( -RK * DTCHEM ) )

         ! In case that the biogenic hydrocarbon is the limiting reactant
         IF ( CHANGE(JHC) >= NMVOC(JHC) ) CHANGE(JHC) = NMVOC(JHC)
      
         ! NMVOC concentration after oxidation reactions
         NMVOC(JHC) = NMVOC(JHC) - CHANGE(JHC)

         IF( CHANGE(JHC) > 1.D-16 ) THEN
            OVER  = 1.D0 / RK
            DO3   = CHANGE(JHC) * KO3(JHC)  * TTO3  * OVER ![=]KG
            DOH   = CHANGE(JHC) * KOH(JHC)  * OHMC  * OVER ![=]KG
            DNO3  = CHANGE(JHC) * KNO3(JHC) * TTNO3 * OVER ![=]KG
         ELSE
            DO3   = 0.D0
            DOH   = 0.D0
            DNO3  = 0.D0
         ENDIF
                  
         ! VOC change by photooxidation of O3 and OH [kg]
         DELHC(1) =  DO3 + DOH 
         DELHC(2) =  DO3 + DOH 
    
         ! VOC change by photooxidation of NO3 [kg]
         DELHC(3) = DNO3
                          
         !Individual Gas-Phase Products:
         DO IPR = 1, NPROD
            GM0(IPR,JHC) = GM0(IPR,JHC) + ALPHA(IPR,JHC) * DELHC(IPR) 
         ENDDO
      ENDDO                     

      !=================================================================
      ! Store Hydrocarbon remaining after oxidation rxn back into STT
      !=================================================================
      STT(I,J,L,IDTALPH) = MAX( NMVOC(1), 1.D-32 )
      STT(I,J,L,IDTLIMO) = MAX( NMVOC(2), 1.D-32 )
      ORVC_TERP(I,J,L)   = MAX( NMVOC(3), 1.D-32 )
      STT(I,J,L,IDTALCO) = MAX( NMVOC(4), 1.D-32 )
      ORVC_SESQ(I,J,L)   = MAX( NMVOC(5), 1.D-32 )

      ! Return to calling program
      END SUBROUTINE CHEM_NVOC

!------------------------------------------------------------------------------

      SUBROUTINE SOA_PARTITION( I, J, L, GM0, AM0 )
!
!******************************************************************************
!  Subroutine SOA_PARTITION partitions the mass of gas and aerosol 
!  tracers according to five Hydrocarbon species and three oxidants.
!  (rjp, bmy, 7/7/04, 7/20/04)
!
!  NOTE: GPROD and APROD are mass ratios of individual oxidation 
!        products of gas/aerosol to the sum of all. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I   (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J   (INTEGER) : GEOS-CHEM latitude index
!  (3 ) L   (INTEGER) : GEOS-CHEM altitude index
! 
!  Arguments as Output:
!  ============================================================================
!  (4 ) GM0 (REAL*8 ) : Gas mass for each HC and its oxidation product     [kg]
!  (5 ) AM0 (REAL*8 ) : Aerosol mass for each HC and its oxidation product [kg]
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! Refrences to F90 modules
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSOG1, IDTSOG2, IDTSOG3, 
     &                         IDTSOA1, IDTSOA2, IDTSOA3

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L
      REAL*8,  INTENT(OUT) :: GM0(NPROD,MHC), AM0(NPROD,MHC)

      ! Local variables
      INTEGER              :: JHC, IPR

      !=================================================================
      ! SOA_PARTITION begins here!
      !=================================================================

      ! Initialize
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = 0D0
         AM0(IPR,JHC) = 0D0
      ENDDO
      ENDDO

      ! Partition the lump of first three HC (ALPH + LIMO + TERP)
      ! oxidation products.  These are grouped together because they
      ! have the same molecular weight.
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG1) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA1) * APROD(I,J,L,IPR,JHC)
      ENDDO                       
      ENDDO

      ! Alcohol
      JHC = 4
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG2) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA2) * APROD(I,J,L,IPR,JHC)
      ENDDO

      ! Sesqterpene
      JHC = 5
      DO IPR = 1, NPROD
         GM0(IPR,JHC) = STT(I,J,L,IDTSOG3) * GPROD(I,J,L,IPR,JHC)
         AM0(IPR,JHC) = STT(I,J,L,IDTSOA3) * APROD(I,J,L,IPR,JHC)
      ENDDO
      
      ! Return to calling program
      END SUBROUTINE SOA_PARTITION

!------------------------------------------------------------------------------

      SUBROUTINE SOA_LUMP( I, J, L, GM0, AM0 )
!
!******************************************************************************
!  Subroutine SOA_LUMP returns the organic gas and aerosol back to the
!  STT array.  (rjp, bmy, 7/7/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I   (INTEGER) : Longitude index
!  (2 ) J   (INTEGER) : Latitude index
!  (3 ) L   (INTEGER) : Altitude index
!  (4 ) GM0 (REAL*8 ) : Gas mass for each HC and its oxidation product     [kg]
!  (5 ) AM0 (REAL*8 ) : Aerosol mass for each HC and its oxidation product [kg]
! 
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD07_HC 
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44
      
      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L
      REAL*8,  INTENT(IN)  :: GM0(NPROD,MHC), AM0(NPROD,MHC)

      ! Local variables
      INTEGER              :: JHC, IPR
      REAL*8               :: GASMASS, AERMASS

      !=================================================================
      ! SOA_LUMP begins here!
      !=================================================================

      ! Initialize
      GASMASS = 0D0
      AERMASS = 0D0

      ! Compute total gas & aerosol mass
      DO JHC = 1, 3
      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO                       
      ENDDO

      !----------------------------
      ! SOA1 net Production (kg)
      !----------------------------
      IF ( ND07 > 0 ) THEN
         AD07_HC(I,J,L,1) = AD07_HC(I,J,L,1)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA1) )
      ENDIF

      STT(I,J,L,IDTSOG1) = MAX( GASMASS, 1.D-32 )
      STT(I,J,L,IDTSOA1) = MAX( AERMASS, 1.D-32 )

      IF ( STT(I,J,L,IDTSOG1) > 1.0E-6 ) THEN
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG1)
         ENDDO                       
         ENDDO
      ELSE
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC)= 0.111111111d0
         ENDDO
         ENDDO   
      ENDIF

      IF ( STT(I,J,L,IDTSOA1) > 1.0E-6 ) THEN
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA1)
         ENDDO                       
         ENDDO
      ELSE
         DO JHC = 1, 3
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.111111111d0
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! Lump of products of fourth Hydrocarbon class (ALCOHOL)
      !=================================================================

      JHC     = 4
      GASMASS = 0D0
      AERMASS = 0D0
      
      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO

      !---------------------------
      ! SOA2 net Production (kg)
      !---------------------------
      IF ( ND07 > 0 ) THEN
         AD07_HC(I,J,L,2) = AD07_HC(I,J,L,2)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA2) )
      ENDIF

      STT(I,J,L,IDTSOG2) = MAX(GASMASS, 1.D-32)
      STT(I,J,L,IDTSOA2) = MAX(AERMASS, 1.D-32)

      IF ( STT(I,J,L,IDTSOG2) > 1.0E-6 ) THEN
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG2)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.333333333d0
         ENDDO
      ENDIF
      
      IF ( STT(I,J,L,IDTSOA2) > 1.0E-6 ) THEN
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA2)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = 0.333333333d0
         ENDDO
      ENDIF

      !=================================================================
      ! Lump of products of fifth Hydrocarbon class (SESQTERPINE)
      !=================================================================
      JHC     = 5
      GASMASS = 0D0
      AERMASS = 0D0

      DO IPR = 1, NPROD
         GASMASS = GASMASS + GM0(IPR,JHC)
         AERMASS = AERMASS + AM0(IPR,JHC)               
      ENDDO

      !---------------------------
      ! SOA3 net Production (kg)
      !---------------------------
      IF ( ND07 > 0 ) THEN
         AD07_HC(I,J,L,3) = AD07_HC(I,J,L,3)
     &                    + ( AERMASS - STT(I,J,L,IDTSOA3) )
      ENDIF

      STT(I,J,L,IDTSOG3) = MAX(GASMASS, 1.D-32)
      STT(I,J,L,IDTSOA3) = MAX(AERMASS, 1.D-32)

      IF ( STT(I,J,L,IDTSOG3) > 1.0E-6 ) THEN
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = GM0(IPR,JHC) / STT(I,J,L,IDTSOG3)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.5D0
         ENDDO
      ENDIF

      IF ( STT(I,J,L,IDTSOA3) > 1.0E-6 ) THEN
         DO IPR = 1, NPROD
            APROD(I,J,L,IPR,JHC) = AM0(IPR,JHC) / STT(I,J,L,IDTSOA3)
         ENDDO                       
      ELSE
         DO IPR = 1, NPROD
            GPROD(I,J,L,IPR,JHC) = 0.5D0
         ENDDO
      ENDIF

      ! make sure there is no second oxidation product 
      ! for SESQTERPENE by OH + O3
      GPROD(I,J,L,2,JHC) = 0.D0
      APROD(I,J,L,2,JHC) = 0.D0

      ! Return to calling program
      END SUBROUTINE SOA_LUMP

!------------------------------------------------------------------------------

      SUBROUTINE SOA_DEPO( TC, DEPID, TRID )
!
!******************************************************************************
!  Subroutine SOA_DEPO computes dry-deposition of a particular SOA species.
!  (rjp, bmy, 7/8/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC    (REAL*8 ) : Array of SOA tracer 
!  (2 ) DEPID (INTEGER) : Dry deposition ID # (from DEPVEL) 
!  (3 ) TRID  (INTEGER) : GEOS-CHEM tracer number 
! 
!  NOTES:
!  (1 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)
      INTEGER, INTENT(IN)    :: DEPID, TRID

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC
      REAL*8                 :: TC0, CNEW, FREQ, AREA_CM2
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! SOA_DEPO begins here!
      !=================================================================

      ! Return if tracer ID or tracer ID is undefined
      IF ( TRID == 0 .OR. DEPID == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial SOA [kg]
         TC0 = TC(I,J,L)

         ! PBLFRAC is only defined up to the tropopause
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DEPID) * BL_FRAC

            ! Amount of SOA[G] left after drydep [kg]
            CNEW = TC0 * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
               FLUX = ( TC0 - CNEW ) 
               FLUX = FLUX * XNUMOL(TRID) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DEPID,1) = AD44(I,J,DEPID,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF    

      ! Return to calling program
      END SUBROUTINE SOA_DEPO

!------------------------------------------------------------------------------

      SUBROUTINE EMISSCARBON
!
!******************************************************************************
!  Subroutine EMISSCARBON is the interface between the GEOS-CHEM model
!  and the CARBONACEOUS AEROSOL emissions (rjp, bmy, 1/24/02, 12/1/04)
!
!  NOTES:
!  (1 ) Now references LSOA from "CMN_SETUP".  Also now call OHNO3TIME since
!        biogenic emissions also have a diurnal variation. (rjp, bmy, 7/15/04)
!  (2 ) Now references LSOA and LPRT from "logical_mod.f".  Now references
!        STT from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Bug fix: removed "," from FORMAT 111.  Also added extra DEBUG_MSG
!        output after calling emissions routines. (bmy, 11/19/04)
!  (4 ) Now always call ANTHRO_CARB_TBOND and ANTHRO_CARB_COOKE.  This will
!        read the T. Bond et al [2004] emissions but overwrite the North
!        America region with monthly-mean emissions from Cooke et al [1999] 
!        with imposed seasonality from R. Park [2003].  (bmy, 12/1/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,    ONLY : AD07
      USE DAO_MOD,     ONLY : PBL
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LSOA,      LPRT
      USE TIME_MOD,    ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,  ONLY : STT
      USE TRACERID_MOD

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND07

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, MONTH, N
      REAL*8               :: BCSRC(IIPAR,JJPAR,2)
      REAL*8               :: OCSRC(IIPAR,JJPAR,2)

      !=================================================================
      ! EMISSCARBON begins here!
      !
      ! Read carbonaceous aerosols from disk and compute hydrophilic 
      ! and hydrophobic fractions. NOTE, CARBON AEROSOLS HAVE TO BE 
      ! ORDERED AS Hydrophilic(BC[1], OC[2]) Hydrophobic(BC[3], OC[4]).
      !=================================================================      

      !--------------------------
      ! Read time-invariant data
      !--------------------------
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100 )

         ! Echo info about ANTHRO emissionls
         WRITE( 6, 110 )
         WRITE( 6, 111 )
         WRITE( 6, 112 )
         WRITE( 6, 113 )

         ! Monthly or annual BIOMASS emissions?
         IF ( USE_MONTHLY_BIOB ) THEN
            WRITE( 6, 120 )
         ELSE
            WRITE( 6, 130 )
         ENDIF
         
         ! Write spacer
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'C A R B O N   A E R O S O L   E M I S S I O N S'    )
 110     FORMAT( 'w/ ANTHROPOGENIC emissions from Bond et al [2004]'  )
 111     FORMAT( '   overwritten with North American emissions from'  )
 112     FORMAT( '   Cooke et al [1999] having imposed seasonality'   )
 113     FORMAT( '   following Park et al [2003]'                     )
 120     FORMAT( 'w/ BIOMASS emissions from GEOS-CHEM inventory'      )
 130     FORMAT( 'w/ BIOMASS emissions from Bond et al [2004]'        )


         ! Initialize arrays
         CALL INIT_CARBON       
         
         ! Read annual mean anthro emissions from T. Bond [2004]
         CALL ANTHRO_CARB_TBOND
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a A_CRB_TBOND' )

         ! Read annual mean biomass emissions if necessary
         IF ( .not. USE_MONTHLY_BIOB ) THEN
            CALL BIOMASS_CARB_TBOND
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a B_CRB_TBOND' )
         ENDIF

         ! Reset flag
         FIRST = .FALSE.
      ENDIF

      ! Compute time scaling arrays which are used both for
      ! biogenic emission and offline OH (rjp, bmy, 7/15/04)
      IF ( LSOA ) THEN
         CALL OHNO3TIME
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: after OHNO3TIME' )
      ENDIF
      
      !--------------------------
      ! Read monthly-mean data
      !--------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN
      
         ! Current month
         MONTH = GET_MONTH()

         ! Overwrite the T. Bond [2004] emissions over North America
         ! with monthly mean anthro emissions from Cooke et al [1999] 
         ! having imposed seasonality by R. Park [2003]
         CALL ANTHRO_CARB_COOKE( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a A_CRB_COOKE' )

         ! Read monthly mean biomass emissions
         IF ( USE_MONTHLY_BIOB ) THEN
            CALL BIOMASS_CARB_GEOS( MONTH )
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: a B_CRB_COOKE' )
         ENDIF

      ENDIF

      !--------------------------
      ! Compute biogenic OC
      !--------------------------
      CALL BIOGENIC_OC
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSCARB: after BIOGENIC_OC' )

      !=================================================================
      ! Sum up BC and OC sources. 
      ! N=1 is HYDROPHILIC; N=2 is HYDROPHOBIC.
      !
      ! COMMENT: Maybe someday we'll want to play with the different 
      ! emission height for different source type.  For example the
      ! carbon from biomass burning could be emitted to the higher 
      ! altitude due to the thermal bouyancy and shallow convection.
      ! The current setting to use EMITHIGH seems rather inefficient 
      ! but robust for sensitivity studies for emission height 
      ! variation on carbon concentrations, so please keep using the 
      ! current setup until we decide otherwise. (rjp, 4/2/02)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Total HYDROPHILIC BC source [kg]
         BCSRC(I,J,1) = ANTH_BLKC(I,J,1) + 
     &                  BIOF_BLKC(I,J,1) + 
     &                  BIOB_BLKC(I,J,1)   

         ! Total HYDROPHOBIC BC source [kg]
         BCSRC(I,J,2) = ANTH_BLKC(I,J,2) +
     &                  BIOF_BLKC(I,J,2) +
     &                  BIOB_BLKC(I,J,2)  

         IF ( LSOA ) THEN

            ! Total HYDROPHILIC OC source [kg]
            ! (Don't use archived TERP_ORGC if LSOA=T)
            OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &                     BIOF_ORGC(I,J,1) + 
     &                     BIOB_ORGC(I,J,1)

         ELSE

            ! Total HYDROPHILIC OC source [kg]
            ! (Use archived TERP_ORGC for if LSOA=F)
            OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &                     BIOF_ORGC(I,J,1) + 
     &                     BIOB_ORGC(I,J,1) + 
     &                     TERP_ORGC(I,J)

         ENDIF

         ! Total HYDROPHOBIC OC source [kg]
         OCSRC(I,J,2) = ANTH_ORGC(I,J,2) + 
     &                  BIOF_ORGC(I,J,2) + 
     &                  BIOB_ORGC(I,J,2) 
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Sum up all carbon tracers throughout the boundary layer
      CALL EMITHIGH( BCSRC, OCSRC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARB: after EMITHIGH' )

      !=================================================================
      ! ND07 diagnostic: Carbon aerosol emissions [kg/timestep]
      !=================================================================
      IF ( ND07 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Anthropogenic BC source
            AD07(I,J,1) = AD07(I,J,1)        +
     &                    ( ANTH_BLKC(I,J,1) + ANTH_BLKC(I,J,2) )
            
            ! Biogenic BC source
            AD07(I,J,2) = AD07(I,J,2)        +
     &                    ( BIOB_BLKC(I,J,1) + BIOB_BLKC(I,J,2) )

            ! Biofuel BC source
            AD07(I,J,3) = AD07(I,J,3)        +
     &                    ( BIOF_BLKC(I,J,1) + BIOF_BLKC(I,J,2) )

            ! Anthropogenic OC source
            AD07(I,J,4) = AD07(I,J,4)        +
     &                    ( ANTH_ORGC(I,J,1) + ANTH_ORGC(I,J,2) )

            ! Biomass OC source
            AD07(I,J,5) = AD07(I,J,5)        +
     &                    ( BIOB_ORGC(I,J,1) + BIOB_ORGC(I,J,2) )

            ! Biofuel OC source
            AD07(I,J,6) = AD07(I,J,6)        + 
     &                    ( BIOF_ORGC(I,J,1) + BIOF_ORGC(I,J,2) )

            ! Terpene source
            AD07(I,J,7) = AD07(I,J,7)        + TERP_ORGC(I,J)

            IF ( LSOA ) THEN

               ! ALPHA-PINENE
               AD07(I,J,8)  = AD07(I,J,8)    + BIOG_ALPH(I,J)

               ! LIMONENE
               AD07(I,J,9)  = AD07(I,J,9)    + BIOG_LIMO(I,J)

               ! TERPENE
               AD07(I,J,10) = AD07(I,J,10)   + BIOG_TERP(I,J)

               ! ALCOHOL
               AD07(I,J,11) = AD07(I,J,11)   + BIOG_ALCO(I,J)

               ! SESQTERPENE
               AD07(I,J,12) = AD07(I,J,12)   + BIOG_SESQ(I,J)

            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARB: after ND07' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE EMISSCARBON

!------------------------------------------------------------------------------

      SUBROUTINE BIOGENIC_OC
!
!******************************************************************************
!  Subroutine BIOGENIC_OC emits secondary organic carbon aerosols.
!  Also modified for SOA tracers. (rjp, bmy, 4/1/04, 11/15/04)
!
!  Terpene emissions as a source of OC:  TERP.GEIA90.a1.2x2.5.*
!  Assuming 10% yield of OC(hydrophilic) from terpene emission.
!
!  NOTES:
!  (1 ) Now separate computation for FULLCHEM and OFFLINE runs (bmy, 7/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Now references LSOA
!        from "logical_mod.f". (bmy, 7/20/04)
!  (3 ) Now reads data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,       ONLY : SUNCOS
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LSOA
      USE TIME_MOD,      ONLY : GET_MONTH,   GET_TS_CHEM, 
     &                          GET_TS_EMIS, ITS_A_NEW_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I, J, IJLOOP, THISMONTH
      REAL*4                :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                :: CONVERT(NVEGTYPE)
      REAL*8                :: GMONOT(NVEGTYPE)
      REAL*8                :: FD2D(IIPAR,JJPAR)
      REAL*8                :: TMMP, EMMO, VALUE
      REAL*8                :: XTAU, STEPS_PER_MON
      REAL*8, PARAMETER     :: FC1 = 136.2364D0 / 120.11D0
      REAL*8, PARAMETER     :: FC2 = 154.2516D0 / 120.11D0
      REAL*8, PARAMETER     :: FC3 = 204.3546D0 / 180.165D0
      REAL*8, PARAMETER     :: FC4 = 152.D0     / 120.11D0
      CHARACTER(LEN=255)    :: FILENAME

      ! Fraction of yield of OC (hydrophilic) from terpene emission
      REAL*8, PARAMETER     :: FBIOG = 1.0d-1

      ! External functions
      REAL*8,  EXTERNAL     :: XLTMMP
      REAL*8,  EXTERNAL     :: EMMONOT

      !=================================================================
      ! BIOGENIC_OC begins here!
      !=================================================================

      ! Get ISOPRENE baseline emissions (first-time only)
      IF ( FIRST ) THEN
         CALL RDISOPT( CONVERT )
         CALL RDMONOT( GMONOT  )
         CALL SETBASE( CONVERT, GMONOT )

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! If secondary organic aerosols are turned off ...
      ! Compute biogenic organic carbon as 0.1 * MONOTERPENES
      !=================================================================
      IF ( .not. LSOA ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, TMMP, EMMO )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! 1-D loop index
            IJLOOP         = ( (J-1) * IIPAR ) + I

             ! Surface temperature [K]
            TMMP           = XLTMMP(I,J,IJLOOP)

            ! EMMO = [kg C/box/time-step] from monoterpenes
            EMMO           = EMMONOT( IJLOOP, TMMP, 1.d0 )

            ! Fraction of EMMO that converts into OC [kg/box/timestep]
            TERP_ORGC(I,J) = EMMO * FBIOG
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! If secondary organic aerosols are turned on ...
      ! Then use CALTECH algorithm 
      !=================================================================
      ELSE 

         ! Current month
         THISMONTH     = GET_MONTH()

         ! Number of emission timesteps per month
         STEPS_PER_MON = ( ( 1440 * NDAYS(THISMONTH) ) / GET_TS_EMIS() )

         !-----------------------------------------
         ! Read data from disk if it's a new month
         !-----------------------------------------
         IF ( ITS_A_NEW_MONTH() ) THEN
         
            ! Get TAU0 value to index the punch file
            XTAU  = GET_TAU0( THISMONTH, 1, 1990 )

            ! Filename for carbon aerosol from fossil fuel use
            FILENAME = TRIM( DATA_DIR ) //
     &                 'carbon_200411/NVOC.geos.' // GET_RES_EXT()

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - BIOGENIC_OC: Reading ', a )

            ! Read NVOC emission in kg/month
            CALL READ_BPCH2( FILENAME, 'NVOCSRCE', 35, 
     &                       XTAU,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast to REAL*8 and resize
            CALL TRANSFER_2D( ARRAY(:,:,1), GEIA_ORVC )

            ! from kgC/month to kgC/timestep
            GEIA_ORVC(:,:) = GEIA_ORVC(:,:) / STEPS_PER_MON
         ENDIF

         !------------------------------
         ! Get TERP_ORGC and DIUR_ORVC
         !------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, TMMP, EMMO )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! 1-D loop index
            IJLOOP         = ( (J-1) * IIPAR ) + I

            ! Surface temperature [K]
            TMMP           = XLTMMP(I,J,IJLOOP)

            ! Monoterpene emission [kg C/box/timestep]
            TERP_ORGC(I,J) = EMMONOT( IJLOOP, TMMP, 1.d0 )

            !---------------------------------------------------
            ! Impose a diurnal variation on NVOC during the day
            !---------------------------------------------------
            IF ( SUNCOS(IJLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN
               DIUR_ORVC(I,J) = GEIA_ORVC(I,J)*
     &                          ( SUNCOS(IJLOOP) / TCOSZ(I,J) ) *
     &                          ( 1440d0        / GET_TS_CHEM() )

               ! Make sure ORVC is not negative
               DIUR_ORVC(I,J) = MAX( DIUR_ORVC(I,J), 0d0 )

            ELSE

               ! At night, ORVC goes to zero
               DIUR_ORVC(I,J) = 0d0

            ENDIF

            !===========================================================
            ! For SOA (Hong Liao, 02/13/04)
            ! 
            ! The input emission data files have units of 
            ! [kg C/box/timestep]
            !
            ! The variable scale is used to convert to the relevant 
            ! units as follows:
            ! (1) Convert from kg C/step to kg compound/step
            ! (2) Multiply by the fraction of monoterpenes that 
            !      contributes to the particular species of interest.
            !
            ! The fraction of monoterpenes is from Table 4 of Griffin
            !  et al., Geophys. Res. Lett. 26 (17): 2721-2724 (1999)
            !===========================================================

            ! ALPHA-PINENE (0.35)
            ! BETA-PINENE lumped with ALPHA-PINENE (0.23)
            ! SABINENE lumped with ALPHA-PINENE (0.05)
            ! D3-CARENE lumped with ALPHA-PINENE (0.04)
            BIOG_ALPH(I,J) = TERP_ORGC(I,J) * FC1 * 0.67D0

            ! TERPENOID KETONE is lumped with SABINENE
            ! Then SABINENE is lumped with ALPHA-PINENE
            BIOG_ALPH(I,J) = BIOG_ALPH(I,J)
     &                     + ( DIUR_ORVC(I,J) * FC4 * 0.04D0 ) !using campher

            ! LIMONENE
            BIOG_LIMO(I,J) = TERP_ORGC(I,J) * FC1 * 0.23D0

            ! TERPINENE is lumped with TERPINOLENE
            BIOG_TERP(I,J) = TERP_ORGC(I,J) * FC1 * 0.03D0

            ! MYRCENE is lumped with TERPENOID ALCOHOL (0.05)
            ! OCIMENE lumped with TERPENOID ALCOHOL (0.02)
            BIOG_ALCO(I,J) = TERP_ORGC(I,J) * FC1 * 0.07D0

            ! Other reactive volatile organic carbon emissions
            BIOG_ALCO(I,J) = BIOG_ALCO(I,J)
     &                     + ( DIUR_ORVC(I,J) * FC2 * 0.09D0 ) !using LINALOOL

            ! We do not transport SESQ (C15H24) 
            ! because its chemical lifetime is short (reactive)
            BIOG_SESQ(I,J) = DIUR_ORVC(I,J) * FC3 * 0.05D0

         ENDDO
         ENDDO

      ENDIF

      ! Return to calling program
      END SUBROUTINE BIOGENIC_OC

!------------------------------------------------------------------------------

      SUBROUTINE ANTHRO_CARB_TBOND
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_TBOND computes annual mean anthropogenic and 
!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!  fractions. (rjp, bmy, 4/2/04, 11/15/04)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_YR
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON (aka ELEMENTAL CARBON)
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON 
      REAL*8, PARAMETER   :: FHO = 0.5d0 

      !=================================================================
      ! ANTHRO_CARB_TBOND begins here!
      !=================================================================

      ! Number of emission timesteps per year
      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )

      ! Get TAU0 value to index the punch file
      XTAU         = GET_TAU0( 1, 1, 2001 )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/year].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !=================================================================

      ! Filename for carbon aerosol from fossil fuel use
      FILENAME = TRIM( DATA_DIR )                        // 
     &           'carbon_200411/BCOC_TBond_fossil.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_TBOND: Reading ', a )

      ! Read BLCK emission
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR
         
         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
        
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/year].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,1) =          FHO *   FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kgC/timestep]
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion as tracer #34 in [kg C/year].  Then convert to 
      ! [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'carbon_200411/BCOC_TBond_biofuel.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic BLACK CARBON from biofuels [kg C /timestep]
         BIOF_BLKC(I,J,1) =          FHB *   FD2D(I,J) / STEPS_PER_YR
         
         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biofuel combustion as tracer #35 in
      ! [kg C/year].  Convert to [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================
      
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_TBOND

!------------------------------------------------------------------------------

      SUBROUTINE ANTHRO_CARB_COOKE( THISMONTH )
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_COOKE computes monthly mean anthropogenic and 
!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!  fractions. (rjp, bmy, 4/2/04, 12/1/04)
!
!  Emissions data comes from the Cooke et al. [1999] inventory and 
!  seasonality imposed by Park et al. [2003].  The data has units of 
!  [kg C/month].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR.  Now only apply
!        Cooke/RJP emissions over the North American region (i.e. the region
!        bounded by indices I1_NA, J1_NA, I2_NA, J2_NA).  (rjp, bmy, 12/1/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_MON
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON aerosol
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON aerosol
      REAL*8, PARAMETER   :: FHO = 0.5d0

      !=================================================================
      ! ANTHRO_CARB_COOKE begins here!
      !=================================================================

      ! Return if we are running on a nested grid (e.g. China) which
      ! does not cover the North America region (rjp, bmy, 12/1/04)
      IF ( I1_NA + J1_NA + I2_NA + J2_NA == 0 ) RETURN

      ! Number of emission timesteps per month
      STEPS_PER_MON = ( ( 1440 * NDAYS( THISMONTH ) ) / GET_TS_EMIS() )
      
      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( THISMONTH, 1, 1998 )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/month].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !
      ! The ANTH_BLKC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of ANTH_BLKC over North America below.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                    //
     &           'carbon_200411/BCOC_anthsrce.geos.' // 
     &            GET_RES_EXT()
       
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_COOKE: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA

         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON
      ENDDO
      ENDDO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/month].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      ! 
      ! The ANTH_ORGC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of ANTH_ORGC over North America below.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA
         
         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,1) = FHO * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON
      ENDDO
      ENDDO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion over Canada and the US as tracer #34 in [kg C/year].  
      ! Then convert to [kg C/timestep] and store in BIOF_BLKC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !
      ! The BIOF_BLKC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of BIOF_BLKC over North America below.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                   //
     &           'carbon_200411/BCOC_biofuel.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data7
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA
         
         ! Hydrophilic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON

      ENDDO
      ENDDO

      !=================================================================
      ! Read ORGANIC CARBON emission from biofuel combustion over 
      ! Canada and the US as tracer #35 in [kg C/year].  Then convert 
      ! to [kg C/timestep] and store in BIOF_ORGC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !
      ! The BIOF_ORGC array is initialized with the Bond et al [2004]
      ! emissions in READ_ANTHRO_TBOND on the very first timestep.
      ! Overwrite the contents of BIOF_ORGC over North America below.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

      DO J = J1_NA, J2_NA
      DO I = I1_NA, I2_NA

         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_COOKE

!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_TBOND
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_TBOND computes annual mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04, 11/15/04)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"   ! Size parameters 

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_YR     
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of carbonaceous aerosols
      REAL*8, PARAMETER   :: FHB = 0.2d0
      REAL*8, PARAMETER   :: FHO = 0.5d0

      !=================================================================
      ! BIOMASS_CARB_TBOND begins here!
      !=================================================================

      ! Number of emission timesteps per year
      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )

      ! Filename containing biomass emissions
      FILENAME = TRIM( DATA_DIR )                         //
     &           'carbon_200411/BCOC_TBond_biomass.geos.' // 
     &            GET_RES_EXT()

      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( 1, 1, 2001 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - BIOMASS_CARB_TBOND: Reading ', a )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from  
      ! biomass burning as tracer #34 in [kg C/year].  Then 
      ! convert to [kg C/timestep] and store in BIOB_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biomass burning as tracer #35 in 
      ! [kg C/year].  Then convert to [kg C/timestep] and store in 
      ! BIOF_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE BIOMASS_CARB_TBOND

!------------------------------------------------------------------------------
! Prior to 1/11/05:
! Now read biomass burning directly from files w/o applying emission factors
! (rjp, bmy, 1/11/05)
!      SUBROUTINE BIOMASS_CARB_GEOS( THISMONTH )
!!
!!******************************************************************************
!!  Subroutine BIOMASS_CARB_TBOND computes annual mean biomass burning 
!!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!!  (rjp, bmy, 4/2/04, 7/20/04)
!!
!!  Emissions data comes from the Bond et al [2004] inventory and has units
!!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!!
!!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!!
!!  NOTES:
!!  (1 ) Now references DATA_DIR from "directory_mod.f".  Also removed CMN,
!!        it's obsolete. (bmy, 7/20/04)
!!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE DIRECTORY_MOD, ONLY : DATA_DIR
!      USE GRID_MOD,      ONLY : GET_AREA_M2
!      USE TIME_MOD,      ONLY : GET_TS_EMIS
!      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!
!#     include "CMN_SIZE"     ! Size parameters
!
!      !-------------------
!      ! Arguments
!      !-------------------
!      INTEGER, INTENT(IN) :: THISMONTH
!
!      !-------------------
!      ! Local variables
!      !-------------------
!
!      ! Hydrophilic fraction of BLACK CARBON
!      REAL*8, PARAMETER   :: FHB = 0.2d0
!
!      ! Hydrophilic fraction of ORGANIC CARBON
!      REAL*8, PARAMETER   :: FHO = 0.5d0 
!
!      ! Black Carbon aerosols emission factor 
!      ! from biomass burning (0.002 kgC/kg)
!      REAL*8, PARAMETER   :: FEC = 2.d-3
!
!      ! Organic Carbon aerosols emission factor 
!      ! from biomass burning (0.014 kgC/kg)
!      REAL*8, PARAMETER   :: FOC = 1.4d-2
!
!      INTEGER             :: I, J
!      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
!      REAL*8              :: FD2D(IIPAR,JJPAR)
!      REAL*8              :: XTAU, BIOCARB
!      REAL*8              :: STEPS_PER_MON, AREA_M2
!      CHARACTER(LEN=255)  :: FILENAME
!      LOGICAL, SAVE       :: FIRSTEMIS = .TRUE.
!      LOGICAL             :: EMFAC_VEG = .TRUE.
!
!      !=================================================================
!      ! BIOMASS_CARB_GEOS begins here!
!      !=================================================================
!
!      ! Number of emission timesteps per month
!      STEPS_PER_MON = ( 1440 * NDAYS( THISMONTH ) ) / GET_TS_EMIS()
!
!      ! Only do the following on the first timestep
!      IF ( FIRSTEMIS ) THEN
!
!         ! Read vegetation factors?
!         IF ( EMFAC_VEG ) THEN
!
!            !===========================================================
!            ! If EMFAC_VEG=T, then read carbon aerosol emission factor 
!            ! [kg/kg] from biomass burning sources depedning on surface 
!            ! vegetation type compiled by [rjp, 2003].  Emission factors 
!            ! are from Andreae and Merlet [2001].
!            !===========================================================
!            FILENAME = TRIM( DATA_DIR )               //
!!-----------------------------------------------------------------------------
!! Prior to 11/15/04:
!!     &                'carbon_200404/emis_fac.EC-OC.' // GET_RES_EXT()
!!-----------------------------------------------------------------------------
!     &                'carbon_200411/emis_fac.EC-OC.' // GET_RES_EXT()
!
!            ! TAU value for reading from the bpch files
!            XTAU = GET_TAU0( 1, 1, 1985 )
!
!            ! Echo info
!            WRITE( 6, 100 ) TRIM( FILENAME )
! 100        FORMAT( '     - BIOMASS_CARB_GEOS: Reading ', a )
!
!            !------------------
!            ! BLACK CARBON
!            !------------------
!            CALL READ_BPCH2( FILENAME, 'EMISFAC', 81,  
!     &                       XTAU,      IGLOB,    JGLOB,
!     &                       1,         ARRAY,    QUIET=.TRUE. ) 
!
!            ! Cast to REAL*8 and resize
!            CALL TRANSFER_2D( ARRAY(:,:,1), EF_BLKC )
!
!            !------------------
!            ! ORGANIC CARBON
!            !------------------
!            CALL READ_BPCH2( FILENAME, 'EMISFAC', 82,  
!     &                       XTAU,      IGLOB,    JGLOB,     
!     &                       1,         ARRAY,    QUIET=.TRUE. ) 
!            
!            ! Cast to REAL*8 and resize
!            CALL TRANSFER_2D( ARRAY(:,:,1), EF_ORGC )
!
!         ELSE
!
!            !===========================================================
!            ! If EMFAC_VEG=F, then use these emission factors:
!            !   BLACK CARBON   : 0.002 [kg C/kg]
!            !   ORGANIC CARBON : 0.014 [kg C/kg] 
!            ! Thus OC/BC = 7. [Chin et al., 2000]
!            !============================================================
!            EF_BLKC(:,:) = FEC
!            EF_ORGC(:,:) = FOC
!           
!         ENDIF
!        
!         ! Reset first-time flag
!         FIRSTEMIS = .FALSE.          
!      ENDIF
!
!      !=================================================================
!      ! Read TOTAL biomass burning [g/cm2/month] as tracer #33
!      ! Convert to [kg C/box/timestep] and store in BIOCARB.  
!      ! 
!      ! Then compute HYDROPHILIC and HYDROPHOBIC fractions of
!      ! BLACK CARBON and ORGANIC CARBON.
!      !=================================================================
!
!      ! Filename
!      FILENAME = TRIM( DATA_DIR )                        //
!     &           'biomass_200110/bioburn.seasonal.geos.' //
!     &           GET_RES_EXT()
!
!      ! Get TAU value for reading the punch file
!      XTAU = GET_TAU0( THISMONTH, 1, 1985 )
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
!
!      ! Read data
!      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 33, 
!     &                 XTAU,      IGLOB,     JGLOB,     
!     &                 1,         ARRAY,     QUIET=.TRUE. )
!
!      ! Cast to REAL*8 and resize
!      CALL TRANSFER_2D ( ARRAY(:,:,1), FD2D )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, AREA_M2, BIOCARB )
!        DO J = 1, JJPAR
!
!           ! Surface area [m2]
!           AREA_M2 = GET_AREA_M2( J )
!
!           DO I = 1, IIPAR
!
!              ! Convert total biomass burned from [g/cm2/month] to
!              ! [kg C/box/timestep].  The factor 10 comes from 
!              ! 1e4 cm2/m2 divided by 1000 g/kg.  
!              BIOCARB = FD2D(I,J) * AREA_M2 * 10.d0 / STEPS_PER_MON
!
!              ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
!              BIOB_BLKC(I,J,1) =          FHB   * EF_BLKC(I,J) * BIOCARB
!
!              ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
!              BIOB_BLKC(I,J,2) = ( 1.D0 - FHB ) * EF_BLKC(I,J) * BIOCARB
!
!              ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
!              BIOB_ORGC(I,J,1) =          FHO   * EF_ORGC(I,J) * BIOCARB
!
!              ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
!              BIOB_ORGC(I,J,2) = ( 1.D0 - FHO ) * EF_ORGC(I,J) * BIOCARB
!
!           ENDDO
!        ENDDO
!!$OMP END PARALLEL DO  
!
!        ! Return to calling program
!        END SUBROUTINE BIOMASS_CARB_GEOS
!
!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_GEOS( THISMONTH )
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_GEOS computes monthly mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04, 7/20/04, 1/11/05)
!
!  Emissions data comes from the Duncan et al. [2001] inventory and has units
!  of [kg C/month].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f".  Also removed CMN,
!        it's obsolete. (bmy, 7/20/04)
!  (2 ) Now read data from "carbon_200411" subdir of DATA_DIR (bmy, 11/15/04)
!  (3 ) Now read BCPO, OCPO biomass burning data directly from files instead
!        of computing from emission factors. (rjp, bmy, 1/11/05)
!******************************************************************************
!
!      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE LOGICAL_MOD,   ONLY : LBBSEA
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR, GET_MONTH, 
     &                          GET_TAU,        GET_YEAR, 
     &                          GET_TS_EMIS
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters

      !-------------------
      ! Arguments
      !-------------------
      INTEGER, INTENT(IN) :: THISMONTH

      !-------------------
      ! Local variables
      !-------------------

      ! Hydrophilic fraction of BLACK CARBON
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON
      REAL*8, PARAMETER   :: FHO = 0.5d0 

      INTEGER             :: I, J, THISYEAR
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: FD2D_BC(IIPAR,JJPAR)
      REAL*8              :: FD2D_OC(IIPAR,JJPAR)
      REAL*8              :: XTAU, BIOBC, BIOOC, TAU
      REAL*8              :: STEPS_PER_MON, AREA_M2
      CHARACTER(LEN=4  )  :: CYEAR
      CHARACTER(LEN=255)  :: BC_FILE, OC_FILE

      !=================================================================
      ! BIOMASS_CARB_GEOS begins here!
      !=================================================================

      ! Number of emission timesteps per month
      STEPS_PER_MON = ( 1440 * NDAYS(THISMONTH) ) / GET_TS_EMIS()

      ! Create a string for the 4-digit year
      THISYEAR = GET_YEAR()
      WRITE( CYEAR, '(i4)' ) THISYEAR

      !=================================================================
      ! Read BC/OC biomass burning [kg C/box/month] as tracer #34, 35
      ! Convert to [kg C/box/timestep] and store in BIO[BC]OC
      ! 
      ! Then compute HYDROPHILIC and HYDROPHOBIC fractions of
      ! BLACK CARBON and ORGANIC CARBON.
      !=================================================================

      ! Use seasonal or interannual emisisons?
      IF ( LBBSEA ) THEN

         !-----------------------------------------
         ! Use seasonal biomass emissions
         !-----------------------------------------

         ! File name for seasonal BCPO biomass emissions
         BC_FILE = TRIM( DATA_DIR )                             //
     &             'biomass_200110/BCPO.bioburn.seasonal.geos.' //
     &             GET_RES_EXT()

         ! File name for seasonal OCPO biomass emissions
         OC_FILE = TRIM( DATA_DIR )                             //
     &             'biomass_200110/OCPO.bioburn.seasonal.geos.' //
     &             GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )      

      ELSE

         !-----------------------------------------
         ! Use interannual biomass emissions for
         ! years between 1996 and 2002
         !-----------------------------------------

         ! File name for interannual BCPO biomass burning emissions
         BC_FILE = TRIM( DATA_DIR )                                //
     &             'biomass_200110/BCPO.bioburn.interannual.geos.' // 
     &             GET_RES_EXT() // '.' // CYEAR

         ! File name for interannual BCPO biomass burning emissions
         OC_FILE = TRIM( DATA_DIR )                                //
     &             'biomass_200110/OCPO.bioburn.interannual.geos.' // 
     &             GET_RES_EXT() // '.' // CYEAR

         ! Use TAU0 value on the 1st of this month to index bpch file
         XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )
  
      ENDIF

!-----------------------------------------------------------------------
! Prior to 1/11/05:
!      ! Echo info
!      !WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - BIOMASS_CARB_GEOS: Reading ', a )
!-----------------------------------------------------------------------

      !--------------
      ! BCPO biomass
      !--------------

      ! Echo info
      WRITE( 6, 100 ) TRIM( BC_FILE )
 100  FORMAT( '     - BIOMASS_CARB_GEOS: Reading ', a )

      ! Read BC emission data [kg/mon]
      CALL READ_BPCH2( BC_FILE, 'BIOBSRCE', 34, 
     &                 XTAU,     IGLOB,     JGLOB,     
     &                 1,        ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), FD2D_BC )

      !--------------
      ! OCPO biomass
      !--------------

      ! Echo info
      WRITE( 6, 100 ) TRIM( OC_FILE )

      ! Read OC emission data [kg/mon]
      CALL READ_BPCH2( OC_FILE, 'BIOBSRCE', 35, 
     &                 XTAU,     IGLOB,     JGLOB,     
     &                 1,        ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), FD2D_OC )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, BIOBC, BIOOC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! kg/mon -> kg/timestep
         BIOBC = FD2D_BC(I,J) / STEPS_PER_MON  
         BIOOC = FD2D_OC(I,J) / STEPS_PER_MON

         ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,1) =          FHB   * BIOBC

         ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,2) = ( 1.D0 - FHB ) * BIOBC

         ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,1) =          FHO   * BIOOC

         ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,2) = ( 1.D0 - FHO ) * BIOOC

      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

        ! Return to calling program
        END SUBROUTINE BIOMASS_CARB_GEOS

!------------------------------------------------------------------------------

      SUBROUTINE EMITHIGH( BCSRC, OCSRC )
!
!******************************************************************************
!  Subroutine EMITHIGH mixes tracer completely from the surface to the PBL
!  top. (rjp, bmy, 4/2/04, 7/8/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) BCSRC (REAL*8) : Array which holds Total BC (H-phobic & H-philic)
!  (2 ) OCSRC (REAL*8) : Array which holds Total OC (H-phobic & H-philic)
!
!  NOTES:
!  (1 ) Now also mix ALPH, LIMO, ALCO tracers (rjp, bmy, 7/8/04)
!  (2 ) Now reference STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : PBL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTBCPO, IDTOCPI, IDTOCPO, 
     &                         IDTALPH, IDTLIMO, IDTALCO
      USE PRESSURE_MOD, ONLY : GET_PEDGE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! SCALE_HEIGHT

      ! Arguments
      REAL*8, INTENT(IN) :: BCSRC(IIPAR,JJPAR,2)
      REAL*8, INTENT(IN) :: OCSRC(IIPAR,JJPAR,2)

      ! Local variables
      LOGICAL            :: IS_BCPO, IS_OCPO, IS_BCPI, IS_OCPI
      LOGICAL            :: IS_ALPH, IS_LIMO, IS_ALCO
      INTEGER            :: I,       J,       L
      REAL*8             :: BLTOP,   BLTHIK,  FTOT,    PB
      REAL*8             :: PT,      DELP,    FEMIS

      !=================================================================
      ! EMITHIGH begins here!
      !=================================================================

      ! Define logical flags for expediency
      IS_BCPI = ( IDTBCPI > 0 )
      IS_OCPI = ( IDTOCPI > 0 ) 
      IS_BCPO = ( IDTBCPO > 0 )
      IS_OCPO = ( IDTOCPO > 0 )
      IS_ALPH = ( IDTALPH > 0 )
      IS_LIMO = ( IDTLIMO > 0 )
      IS_ALCO = ( IDTALCO > 0 )

      !=================================================================
      ! Compute FEMIS -- fraction of box (I,J,L) w/in the PBL
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FTOT, BLTOP, BLTHIK, L, PB, PT, DELP, FEMIS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize
         FTOT  = 0d0
         FEMIS = 0d0

#if   defined( GEOS_4 )

         ! BLTOP = pressure at PBL top [hPa]
         ! Use barometric law since PBL is in [m]
         BLTOP  = GET_PEDGE(I,J,1) * EXP( -PBL(I,J) / SCALE_HEIGHT )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = GET_PEDGE(I,J,1) - BLTOP

#else

         ! BLTOP = pressure of PBL top [hPa]
         BLTOP  = GET_PEDGE(I,J,1) - MAX( PBL(I,J), 1.D0 )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = MAX( PBL(I,J), 1.D0 )

#endif

         !==============================================================
         ! Loop thru tropospheric levels
         !==============================================================
         DO L = 1, LLTROP

            ! Pressure at edges of grid box(I,J,L) [hPa]
            PB   = GET_PEDGE(I,J,L)
            PT   = GET_PEDGE(I,J,L+1)

            ! Thickness of grid box (I,J,L) [hPa]
            DELP = PB - PT

            ! FEMIS is the fraction of the PBL 
            ! which is occupied by this level L
            IF ( BLTOP <= PT )  THEN
               FEMIS = DELP / BLTHIK

            ELSEIF ( BLTOP > PT .AND. BLTOP <= PB ) THEN
               FEMIS = ( PB - BLTOP ) / BLTHIK

            ELSEIF ( BLTOP > PB ) THEN
               CYCLE      
 
            ENDIF
            
            ! Fraction of data partitioned into 
            ! each level should sum to 1.0
            FTOT = FTOT + FEMIS

            !===========================================================
            ! Partition organic tracers equally throughout the 
            ! boundary layer -- store back into STT array
            !===========================================================

            ! Hydrophilic BLACK CARBON
            IF ( IS_BCPI ) THEN
               STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) + 
     &                              ( FEMIS * BCSRC(I,J,1) )
            ENDIF

            ! Hydrophilic ORGANIC CARBON
            IF ( IS_OCPI ) THEN
               STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) + 
     &                              ( FEMIS * OCSRC(I,J,1) )
            ENDIF
            
            ! Hydrophobic BLACK CARBON
            IF ( IS_BCPO ) THEN
               STT(I,J,L,IDTBCPO) = STT(I,J,L,IDTBCPO) + 
     &                              ( FEMIS * BCSRC(I,J,2) )
            ENDIF

            ! Hydrophobic ORGANIC CARBON
            IF ( IS_OCPO ) THEN
               STT(I,J,L,IDTOCPO) = STT(I,J,L,IDTOCPO) + 
     &                              ( FEMIS * OCSRC(I,J,2) )
            ENDIF

            ! ALPHA-PINENE
            IF ( IS_ALPH ) THEN
               STT(I,J,L,IDTALPH) = STT(I,J,L,IDTALPH) + 
     &                              ( FEMIS * BIOG_ALPH(I,J) )
            ENDIF

            ! LIMONENE
            IF ( IS_LIMO ) THEN
               STT(I,J,L,IDTLIMO) = STT(I,J,L,IDTLIMO) + 
     &                              ( FEMIS * BIOG_LIMO(I,J) )

               ORVC_TERP(I,J,L)   = ORVC_TERP(I,J,L) + 
     &                              ( FEMIS * BIOG_TERP(I,J) )
            ENDIF

            ! ALCOHOL and SESQTERPENE (not a tracer)
            IF ( IS_ALCO ) THEN
               STT(I,J,L,IDTALCO) = STT(I,J,L,IDTALCO) + 
     &                              ( FEMIS * BIOG_ALCO(I,J) )
               
               ORVC_SESQ(I,J,L)   = ORVC_SESQ(I,J,L) + 
     &                              ( FEMIS * BIOG_SESQ(I,J) )
            ENDIF

         ENDDO

         ! Error check
         IF ( ABS( FTOT - 1.d0 ) > 1.d-3 ) THEN
            CALL ERROR_STOP( 'Check vertical. distribution!',
     &                       'EMITHIGH ("carbon_mod.f")' )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMITHIGH

!------------------------------------------------------------------------------

      SUBROUTINE OHNO3TIME
!
!******************************************************************************
!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 1/18/05)
!
!  NOTES:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (bmy, 1/18/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC, 
     &                     GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TCOSZ(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            !-------------------------------------------
            ! Prior to 1/18/05:
            !! IJLOOP is the 1-D loop index for SUNCOS
            !IJLOOP = 0
            !-------------------------------------------

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               !---------------------------------------
               ! Prior to 1/18/05:
               ! Now use analytic function for IJLOOP
               !IJLOOP = IJLOOP + 1
               !---------------------------------------
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!------------------------------------------------------------------------------

      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_OH returns OH from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  Imposes a diurnal variation on
!  OH for offline simulations. (bmy, 7/9/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DAO_MOD,       ONLY : SUNCOS
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,  ONLY : IDOH

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !---------------------
         ! Coupled simulation
         !---------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Take OH from the SMVGEAR array CSPEC
         ! OH is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
         ELSE
            OH_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                     ( 1440d0        / GET_TS_CHEM() )

            ! Make sure OH is not negative
            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
         ELSE

            ! At night, OH goes to zero
            OH_MOLEC_CM3 = 0d0

         ENDIF

      ELSE

         !---------------------
         ! Invalid sim type!
         !---------------------        
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_OH

!------------------------------------------------------------------------------

      FUNCTION GET_NO3( I, J, L ) RESULT( NO3_MOLEC_CM3 ) 
!
!******************************************************************************
!  Function GET_NO3 returns NO3 from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  For offline runs, the concentration
!  of NO3 is set to zero during the day. (rjp, bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!        zero during the day. (rjp, bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : CSPEC, JLOP
      USE DAO_MOD,        ONLY : AD,    SUNCOS
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_NO3_MOD, ONLY : NO3
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDNO3

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: NO3_MOLEC_CM3
      REAL*8,  PARAMETER  :: XNUMOL_NO3 = 6.022d23 / 62d-3
 
      ! External functions
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_NO3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !----------------------
         ! Fullchem simulation
         !----------------------

         ! 1-D SMVGEAR grid box index
         JLOOP = JLOP(I,J,L)

         ! Take NO3 from the SMVGEAR array CSPEC
         ! NO3 is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
         ELSE
            NO3_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !==============================================================  
         ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
         ! in [v/v].  Convert these to [molec/cm3] as follows:
         !
         !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
         !  ------- = --------- * -------- * ---------------- =  kg NO3 
         !  vol air   moles air      1        kg air/mole air
         !
         ! And then we convert [kg NO3] to [molec NO3/cm3] by:
         !  
         !  kg NO3   molec NO3   mole NO3     1     molec NO3
         !  ------ * --------- * -------- * ----- = --------- 
         !     1     mole NO3     kg NO3     cm3       cm3
         !          ^                    ^
         !          |____________________|  
         !            this is XNUMOL_NO3
         !
         ! If at nighttime, use the monthly mean NO3 concentration from
         ! the NO3 array of "global_no3_mod.f".  If during the daytime,
         ! set the NO3 concentration to zero.  We don't have to relax to 
         ! the monthly mean concentration every 3 hours (as for HNO3) 
         ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02) 
         !==============================================================

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test if daylight
         IF ( SUNCOS(JLOOP) > 0d0 ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0d0
              
         ELSE

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L) * AD(I,J,L) * ( 62d0/28.97d0 ) 
               
            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL(I,J,L)
                  
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )

      ELSE

         !----------------------
         ! Invalid sim type!
         !----------------------       
         CALL ERROR_STOP( 'Invalid NSRCX!','GET_NO3 (sulfate_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_NO3

!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,       ONLY : SUNCOS, AD
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_O3_MOD, ONLY : O3
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,  ONLY : IDO3

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! XNUMOLAIR

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: O3_MOLEC_CM3
      REAL*8,  PARAMETER  :: XNUMOL_O3 = 6.022d23 / 48d-3
 
      ! External functions
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_O3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Get O3 from CSPEC [molec/cm3]
         ! O3 data will only be defined below the tropopause
         IF ( JLOOP  > 0 ) THEN
            O3_MOLEC_CM3 = CSPEC(JLOOP,IDO3)
         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF
         
      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !--------------------
         ! Offline simulation
         !--------------------

         ! Get O3 [v/v] for this gridbox & month
         ! O3 data will only be defined below the tropopause
         IF ( L <= LLTROP ) THEN

            ! Get O3 [v/v] and convert it to [kg]
            O3_MOLEC_CM3 = O3(I,J,L) * AD(I,J,L) * ( 48d0/28.97d0 )
               
            ! Convert O3 from [kg] to [molec/cm3]
            O3_MOLEC_CM3 = O3_MOLEC_CM3 * XNUMOL_O3 / BOXVL(I,J,L)
         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            O3_MOLEC_CM3 = O3_MOLEC_CM3                     *        
     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) )   *
     &                     ( 1440d0        / GET_TS_CHEM() )

            ! Make sure OH is not negative
            O3_MOLEC_CM3 = MAX( O3_MOLEC_CM3, 0d0 )

         ELSE
            O3_MOLEC_CM3 = 0d0
         ENDIF

      ELSE

         !--------------------
         ! Invalid sim type!
         !--------------------
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_O3 (sulfate_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------

      SUBROUTINE INIT_CARBON
!
!******************************************************************************
!  Subroutine INIT_CARBON initializes all module arrays. 
!  (rjp, bmy, 4/1/04, 12/1/04)
!
!  NOTES:
!  (1 ) Also added arrays for secondary organic aerosols (rjp, bmy, 7/8/04)
!  (2 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (3 ) Now reference LSOA from "logical_mod.f" not CMN_SETUP.  Now call
!        GET_BOUNDING_BOX from "grid_mod.f" to compute the indices I1_NA,
!        I2_NA, J1_NA, J2_NA which define the N. America region. (bmy, 12/1/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR, ERROR_STOP
      USE GRID_MOD,    ONLY : GET_BOUNDING_BOX
      USE LOGICAL_MOD, ONLY : LSOA

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS, INDICES(4)
      REAL*8        :: COORDS(4)

      !=================================================================
      ! INIT_CARBON begins here!
      !=================================================================
      
      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN

      ALLOCATE( ANTH_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_BLKC' )
      ANTH_BLKC = 0d0

      ALLOCATE( ANTH_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_ORGC' )
      ANTH_ORGC = 0d0

      ALLOCATE( BIOB_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_BLKC' )
      BIOB_BLKC = 0d0

      ALLOCATE( BIOB_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_ORGC' )
      BIOB_ORGC = 0d0

      ALLOCATE( BIOF_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_BLKC' )
      BIOF_BLKC = 0d0

      ALLOCATE( BIOF_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_ORGC' )
      BIOF_ORGC = 0d0

      ALLOCATE( TERP_ORGC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TERP_ORGC' )
      TERP_ORGC = 0d0

      ALLOCATE( BCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV' )
      BCCONV = 0d0

      ALLOCATE( OCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
      OCCONV = 0d0

      !=================================================================
      ! These only have to be allocated if we are
      ! reading in monthly mean biomass burning
      !=================================================================
      IF ( USE_MONTHLY_BIOB ) THEN

         ALLOCATE( EF_BLKC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_BLKC' )
         EF_BLKC = 0d0

         ALLOCATE( EF_ORGC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_ORGC' )
         EF_ORGC = 0d0

      ENDIF

      !=================================================================
      ! SOA arrays only have to be allocated if LSOA = T
      !=================================================================
      IF ( LSOA ) THEN

         ALLOCATE( BIOG_ALPH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_ALPH' )
         BIOG_ALPH = 0d0

         ALLOCATE( BIOG_LIMO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_LIMO' )
         BIOG_LIMO = 0d0

         ALLOCATE( BIOG_ALCO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_ALCO' )
         BIOG_ALCO = 0d0

         ALLOCATE( BIOG_TERP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_TERP' )
         BIOG_TERP = 0d0

         ALLOCATE( BIOG_SESQ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOG_SESQ' )
         BIOG_SESQ = 0d0

         ALLOCATE( DIUR_ORVC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DIUR_ORVC' )
         DIUR_ORVC = 0d0

         ALLOCATE( GEIA_ORVC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEIA_ORVC' )
         GEIA_ORVC = 0d0

         ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
         TCOSZ = 0d0

         ALLOCATE( ORVC_TERP( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORVC_TERP' )
         ORVC_TERP = 0d0

         ALLOCATE( ORVC_SESQ( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORVC_SESQ' )
         ORVC_SESQ = 0d0

         ALLOCATE( GPROD( IIPAR, JJPAR, LLPAR, NPROD, MHC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GPROD' )
         GPROD = 0D0

         ALLOCATE( APROD( IIPAR, JJPAR, LLPAR, NPROD, MHC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'APROD' )
         APROD = 0D0

      ENDIF

      !=================================================================
      ! Compute indices which define the N. America region so that we 
      ! can overwrite T. Bond emissions w/ Cooke/RJP emissions 
      !=================================================================

#if   defined( GRID1x1 ) && defined( NESTED_NA )

      ! For 1x1 N. America nested grid: set indices to grid extent
      I1_NA = 1
      J1_NA = 1
      I2_NA = IIPAR
      J2_NA = JJPAR

#elif defined( GRID1x1 ) && defined( NESTED_CH )

      ! For 1x1 China nested grid: we don't cover N. America region
      ! Setting these to zero will turn off Cooke/RJP emissions
      I1_NA = 0
      J1_NA = 0
      I2_NA = 0
      J2_NA = 0

#else

      ! Definition of the N. American bounding box
      ! with LL corner (10N,165W) and UR corner (90N,40W)
      !            Lon_LL  Lat_LL  Lon_UR  Lat_UR
      COORDS = (/ -165d0,  10d0,  -40d0,   90d0  /)
      
      ! Get the indices corresponding to the lon/lat values in COORDS
      CALL GET_BOUNDING_BOX( COORDS, INDICES )

      ! Copy values from INDEX array to scalars
      I1_NA = INDICES(1)
      J1_NA = INDICES(2)
      I2_NA = INDICES(3)
      J2_NA = INDICES(4)

#endif
              
      ! Reset IS_INIT before exiting
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_CARBON

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_CARBON
!
!******************************************************************************
!  Subroutine CLEANUP_CARBON deallocates all module arrays 
!  (rjp, bmy, 4/1/04, 7/8/04)
!
!  NOTES:
!  (1 ) Now deallocate arrays for secondary organic aerosols (rjp, bmy, 7/8/04)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_CARBON begins here!
      !=================================================================
      IF ( ALLOCATED( ANTH_BLKC ) ) DEALLOCATE( ANTH_BLKC )
      IF ( ALLOCATED( ANTH_ORGC ) ) DEALLOCATE( ANTH_ORGC )
      IF ( ALLOCATED( BIOB_BLKC ) ) DEALLOCATE( BIOB_BLKC )
      IF ( ALLOCATED( BIOB_ORGC ) ) DEALLOCATE( BIOB_ORGC )
      IF ( ALLOCATED( BIOF_BLKC ) ) DEALLOCATE( BIOF_BLKC )
      IF ( ALLOCATED( BIOF_ORGC ) ) DEALLOCATE( BIOF_ORGC )
      IF ( ALLOCATED( TERP_ORGC ) ) DEALLOCATE( TERP_ORGC )
      IF ( ALLOCATED( BCCONV    ) ) DEALLOCATE( BCCONV    )
      IF ( ALLOCATED( OCCONV    ) ) DEALLOCATE( OCCONV    )
      IF ( ALLOCATED( EF_BLKC   ) ) DEALLOCATE( EF_BLKC   )
      IF ( ALLOCATED( EF_ORGC   ) ) DEALLOCATE( EF_ORGC   )
      IF ( ALLOCATED( BIOG_ALPH ) ) DEALLOCATE( BIOG_ALPH )
      IF ( ALLOCATED( BIOG_LIMO ) ) DEALLOCATE( BIOG_LIMO )
      IF ( ALLOCATED( BIOG_ALCO ) ) DEALLOCATE( BIOG_ALCO )
      IF ( ALLOCATED( BIOG_TERP ) ) DEALLOCATE( BIOG_TERP )
      IF ( ALLOCATED( BIOG_SESQ ) ) DEALLOCATE( BIOG_SESQ )
      IF ( ALLOCATED( DIUR_ORVC ) ) DEALLOCATE( DIUR_ORVC )
      IF ( ALLOCATED( GEIA_ORVC ) ) DEALLOCATE( GEIA_ORVC )
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ     )
      IF ( ALLOCATED( GPROD     ) ) DEALLOCATE( GPROD     )
      IF ( ALLOCATED( APROD     ) ) DEALLOCATE( APROD     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_CARBON

!------------------------------------------------------------------------------

      END MODULE CARBON_MOD
