!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: turner_flux_mod
!
! !DESCRIPTION: Module TURNER\FLUX\_MOD contains variables and routines to 
!  read the Alex Turner's methane fluxes and convert to equivalent
!  emissions of C2H6 and C3H8 if required
!\\
!\\
! !INTERFACE: 
!
      MODULE TURNER_FLUX_MOD
! 
! !USES:
!
      IMPLICIT NONE
#     include "netcdf.inc"
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_TURNER_FLUX
      PUBLIC  :: EMISS_TURNER_FLUX
      PUBLIC  :: GET_TURNER_ANTHRO
      PUBLIC  :: GET_TURNER_BIOFUEL
!
! !PRIVATE MEMBER FUNCTIONS:
      PRIVATE :: INIT_TURNER_FLUX
      PRIVATE :: TOTAL_TURNER_TG
!
! !REMARKS:
!     
! !REVISION HISTORY:
!  13 May 2014 - J. Fisher  - Initial version, based on nei2011_anthro_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Arrays for emissions (lon/lat) - anthro & biofuel
      REAL*8,  ALLOCATABLE        :: C2H6an(:,:), C2H6bf(:,:) 
      REAL*8,  ALLOCATABLE        :: C3H8an(:,:), C3H8bf(:,:)  
      REAL*8                      :: T_C2H6an_MON, T_C2H6bf_MON
      REAL*8                      :: T_C3H8an_MON, T_C3H8bf_MON

      ! Shadow logical variables from Input_Opt
      LOGICAL                      :: LTZOMPA

      ! Variables for regridding
      INTEGER                       :: IIIPAR0
      INTEGER                       :: JJJPAR0
      REAL*8,  ALLOCATABLE          :: XEDGE_TURNER(:)
      REAL*8,  ALLOCATABLE          :: YEDGE_TURNER(:)
      REAL*8,  ALLOCATABLE          :: XEDGE_MODELG(:)
      REAL*8,  ALLOCATABLE          :: YEDGE_MODELG(:)
!
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_turner_anthro
!
! !DESCRIPTION: Function GET\_TURNER\_ANTHRO returns the anthropogenic
!  emission for GEOS-Chem grid box (I,J) and tracer N (C2H6 or C3H8),
!  scaled to Turner CH4 fluxes
!  Emissions can be returned in units of [kg/s] or [molec/cm2/s].
!\\
! !INTERFACE:
!
      FUNCTION GET_TURNER_ANTHRO( I, J, N ) RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TRACERID_MOD, ONLY : IDTC2H6, IDTC3H8
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N
 
      ! ALWAYS returns emissions in [molec/cm2/s]

! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE  

! !REMARKS:
!
! !REVISION HISTORY: 
!  13 May 2016 - J. Fisher - initial version  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_TURNER_ANTHRO begins here!
      !=================================================================

      IF ( N == IDTC2H6 ) THEN
         ! [molec/cm2/s]
         VALUE = C2H6an(I,J)
      ELSE IF ( N == IDTC3H8 ) THEN
         ! [molec/cm2/s]
         VALUE = C3H8an(I,J)
      ELSE
         ! Otherwise return a negative value to indicate
         ! that there are no Turner-derived emissions for tracer N
         VALUE = -1d0
         RETURN
      ENDIF

      ! Return to calling program
      END FUNCTION GET_TURNER_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_turner_biofuel
!
! !DESCRIPTION: Function GET\_TURNER\_BIOFUEL returns the biofuel
!  emission for GEOS-Chem grid box (I,J) and tracer N (C2H6 or C3H8),
!  scaled to Turner CH4 fluxes
!  Emissions can be returned in units of [kg/s] or [molec/cm2/s].
!\\
! !INTERFACE:
!
      FUNCTION GET_TURNER_BIOFUEL( I, J, N ) RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TRACERID_MOD, ONLY : IDTC2H6, IDTC3H8
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N
 
      ! ALWAYS returns emissions in [molec/cm2/s]

! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE  

! !REMARKS:
!
! !REVISION HISTORY: 
!  13 May 2016 - J. Fisher - initial version  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_TURNER_BIOFUEL begins here!
      !=================================================================

      IF ( N == IDTC2H6 ) THEN
         ! [molec/cm2/s]
         VALUE = C2H6bf(I,J)
      ELSE IF ( N == IDTC3H8 ) THEN
         ! [molec/cm2/s]
         VALUE = C3H8bf(I,J)
      ELSE
         ! Otherwise return a negative value to indicate
         ! that there are no Turner-derived emissions for tracer N
         VALUE = -1d0
         RETURN
      ENDIF

      ! Return to calling program
      END FUNCTION GET_TURNER_BIOFUEL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_turner_flux
!
! !DESCRIPTION: Subroutine EMISS\_TURNER\_FLUX reads the Turner CH4
!  flux fields at .5 x 0.67 and 4x5 resolution and converts to
!  appropriate C2H6 and C3H8 emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_TURNER_FLUX( am_I_Root, Input_Opt, &
                                             State_Chm, RC         )
!
! !USES:
!
      USE DIRECTORY_MOD,     ONLY : DATA_DIR
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE GRID_MOD,          ONLY : GET_XOFFSET
      USE GRID_MOD,          ONLY : GET_YOFFSET
      USE GLOBAL_GRID_MOD,   ONLY : GET_XEDGE_G, GET_YEDGE_G
      USE GRID_MOD,          ONLY : GET_XEDGE, GET_YEDGE
      USE REGRID_A2A_MOD,    ONLY : DO_REGRID_A2A, MAP_A2A
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH, GET_DAY
      USE TIME_MOD,          ONLY : GET_HOUR, GET_DAY_OF_WEEK
      USE TRACER_MOD,        ONLY : XNUMOL

      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE NCDF_MOD,          ONLY : NC_READ
      USE TRACERID_MOD,      ONLY : IDTC2H6, IDTC3H8

      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
      
      USE m_netcdf_io_open     
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      USE m_netcdf_io_get_dimlen
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

!
! !REVISION HISTORY:
!  13 May 2016 - J. Fisher   - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, A, B, I0, J0, SNo
      INTEGER                    :: KLM, SPECIES_ID(18)
      INTEGER                    :: st3d(3), ct3d(3)
      INTEGER                    :: fId
      REAL*8                     :: ARRAYtmp(72,46,1)
      REAL*8                     :: ARRAYan(72,46)
      REAL*8                     :: ARRAYbf(72,46)
      REAL*8, TARGET             :: GEOS_NATIVE(I01x01,J01x01)
      CHARACTER(LEN=255)         :: DATA_DIR_TURNER
      CHARACTER(LEN=255)         :: FILENAME, LLFILENAME
      CHARACTER(LEN=24)          :: SPCLIST(18)
      CHARACTER(LEN=8)           :: SId
      CHARACTER(LEN=5)           :: SNAME
      CHARACTER(LEN=3)           :: TTMON
      CHARACTER(LEN=2)           :: FDAY, FMON
      CHARACTER(LEN=4)           :: FYR
      REAL*8, POINTER            :: INGRID(:,:) => NULL()
      REAL*8, POINTER            :: OUTGRID(:,:) => NULL()

      ! Scaling relative to CH4 (see Tzompa et al. 2016)
      ! This is [kg SPECIES]/[kg CH4]
      REAL*4, PARAMETER          :: ScC2H6an=8.0d-2
      REAL*4, PARAMETER          :: ScC2H6bf=2.86d-1
      REAL*4, PARAMETER          :: ScC3H8an=9.32d-2
      REAL*4, PARAMETER          :: ScC3H8bf=6.63d-2

      ! For grid
      REAL*4                     :: DEG2RAD
      !=================================================================
      ! EMISS_TURNER_FLUX begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_TURNER_FLUX( am_I_Root, Input_Opt, RC )

         DEG2RAD = (4. * ATAN(1.) ) /180.

         ! Define TURNER grid box lat and lon edges
         XEDGE_TURNER( 1 ) = -180.d0
         DO I = 2,I01x01 +1
            XEDGE_TURNER( I ) = XEDGE_TURNER( I-1 ) + 1.d-1
         END DO

         YEDGE_TURNER( 1 ) = -89.975d0
         DO J = 2, J01x01+1
            YEDGE_TURNER( J ) = YEDGE_TURNER( J-1 ) + 1.d-1
         END DO

         DO J = 1,J01x01+1
            YEDGE_TURNER( J ) = SIN( YEDGE_TURNER( J ) * DEG2RAD)
         END DO

         ! Define global grid box lat and lon edges at model resolution
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE_G ( I )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE_G ( J )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#else

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE( I, 1, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE( 1, J, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#endif
         FIRST = .FALSE.
      ENDIF

      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) // &
                  'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic01x01.nc'

      ! Base data directory - for now 4x5 files only!
      DATA_DIR_TURNER = '/short/m19/geos-chem/data/' // &
                        'ExtData/HEMCO/TURNER/v2016-05/4x5/'

      ! Filename - biofuel & oil/gas fluxes not time-varying
      FILENAME = TRIM( DATA_DIR_TURNER ) // &
                 'CH4_flux_OilGas_Biofuel.nc'

      ! Allocate start and count arrays
      st3d = (/1, 1, 1/)          !Start lon/lat/time
      ct3d = (/IIPAR, JJPAR, 1/)  !Count lon/lat/time

      ! Open netCDF file for reading
      CALL Ncop_Rd(fId,  TRIM(FILENAME))
      
 100  FORMAT( '     - EMISS_TURNER_FLUX:  Reading ', a, ' -> ', a )

      ! Read oil&gas and biofuels from netCDF files
      ! Units are in kg/m2/s
      WRITE( 6, 100 )  TRIM( FILENAME ), '2-OilGas'
      Call NcRd(ARRAYtmp, fId, 'emissions_cat02', st3d, ct3d )

      ! Cast to 2-D array
      ARRAYan = ARRAYtmp(:,:,1)

      WRITE( 6, 100 )  TRIM( FILENAME ), '6-Biofuel'
      Call NcRd(ARRAYtmp, fId, 'emissions_cat06', st3d, ct3d )

      ! Cast to 2-D array
      ARRAYbf = ARRAYtmp(:,:,1)

      ! Close netCDF file
      CALL NcCl( fId )

!! Put regridding code here

      ! Convert from [kg CH4/m2/s] to [kg SPECIES/cm2/s] using scale
      ! factors, then convert to [atom C/cm2/s]
      ! [kg C2H6/m2/s] * [1 m2/1d4 cm2]* [1 mole/30d-3 kg] * [Na molec/mole] 
      !                * [2 atom C / molec C2H6]
      C2H6an = ( ARRAYan * 6.0225d23 * 2d0 / 30d1 ) * ScC2H6an
      C2H6bf = ( ARRAYbf * 6.0225d23 * 2d0 / 30d1 ) * ScC2H6bf
      ! [kg C3H8/m2/s] * [1 m2/1d4 cm2]* [1 mole/44d-3 kg] * [Na molec/mole] 
      !                * [3 atom C / molec C3H8]
      C3H8an = ( ARRAYan * 6.0225d23 * 3d0 / 44d1 ) * ScC3H8an
      C3H8bf = ( ARRAYbf * 6.0225d23 * 3d0 / 44d1 ) * ScC3H8bf

      !--------------------------
      ! Print emission totals
      !--------------------------

      CALL TOTAL_TURNER_Tg

      ! Return to calling program
      END SUBROUTINE EMISS_TURNER_FLUX

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_turner_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_TURNER\_TG prints the totals for the 
!  emissions of C2H6 and C3H8 scaled from Turner CH4 fluxes.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_TURNER_TG
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE TIME_MOD,     ONLY : ITS_A_NEW_MONTH
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTC2H6, IDTC3H8

!
! !REVISION HISTORY: 
!   13 May 2016 - J. Fisher - initial version
!   24 May 2016 - J. Fisher - updated since files are time invariant
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: II, JJ
      REAL*8              :: T_C2H6an, T_C3H8an
      REAL*8              :: T_C2H6bf, T_C3H8bf
      REAL*8              :: tmpArea(IIPAR, JJPAR)
  
      !=================================================================
      ! TOTAL_TURNER_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'T U R N E R - D E R I V E D   E M I S S I O N S', / )
      
      tmpArea = 0d0
      DO II = 1, IIPAR
      DO JJ = 1, JJPAR
         tmpArea(II,JJ) = GET_AREA_CM2(II,JJ,1)
      ENDDO
      ENDDO
         
      ! Units are in [atom C/cm2/s]
      ! Total C2H6  [Tg C]
      IF ( IDTC2H6 /= 0 ) THEN
          T_C2H6an = SUM( C2H6an*tmpArea ) * &
                     ( 365d0 * 86400d-9 ) / XNUMOL(IDTC2H6)
          T_C2H6bf = SUM( C2H6bf*tmpArea ) * &
                     ( 365d0 * 86400d-9 ) / XNUMOL(IDTC2H6)
      ENDIF
      ! Total C3H8  [Tg C]
      IF ( IDTC3H8 /= 0 ) THEN
          T_C3H8an = SUM( C3H8an*tmpArea ) * &
                     ( 365d0 * 86400d-9 ) / XNUMOL(IDTC3H8)
          T_C3H8bf = SUM( C3H8bf*tmpArea ) * &
                     ( 365d0 * 86400d-9 ) / XNUMOL(IDTC3H8)
      ENDIF

      ! Print anthro totals in [Tg C]
      WRITE( 6, 111 ) 'C2H6 ', T_C2H6an, '[ Tg C/yr ]'
      WRITE( 6, 111 ) 'C3H8 ', T_C3H8an, '[ Tg C/yr ]'

      ! Print biofuel totals in [Tg C]
      WRITE( 6, 121 ) 'C2H6 ', T_C2H6bf, '[ Tg C/yr ]'
      WRITE( 6, 121 ) 'C3H8 ', T_C3H8bf, '[ Tg C/yr ]'

      ! Format statement
 111  FORMAT( 'Turner-derived anthro ', a5, '(base year 2009-2011): ', &
              f11.5, 1x, a11 )
 121  FORMAT( 'Turner-derived biofuel ', a5, '(base year 2009-2011): ', &
              f11.5, 1x, a11 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_TURNER_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_turner_flux
!
! !DESCRIPTION: Subroutine INIT\_TURNER\_FLUX allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TURNER_FLUX( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
      USE CMN_SIZE_MOD    ! Size parameters

      USE GLOBAL_GRID_MOD, ONLY : GET_IIIPAR
      USE GLOBAL_GRID_MOD, ONLY : GET_JJJPAR
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  13 May 2016 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_TURNER_FLUX begins here!
      !=================================================================
      ! Assume success
      RC        =  GIGC_SUCCESS
      
      ! Return if LTZOMPA is false
      IF ( .not. Input_Opt%LTZOMPA ) RETURN

      ! Get global longitude/latitude extent [# of boxes]
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

      ! Nested grids utilize global longitude and latidude extent
      ! parameters IIIPAR and JJJPAR (from global_grid_mod)
      IIIPAR0 = GET_IIIPAR()
      JJJPAR0 = GET_JJJPAR()

#else

      ! Global grids utilize window longitude and latidude extent
      ! parameters IIPAR and JJPAR (from CMN_SIZE_mod)
      IIIPAR0 = IIPAR
      JJJPAR0 = JJPAR

#endif
      
     ! Allocate array to hold TURNER grid box lon edges
      ALLOCATE( XEDGE_TURNER( I01x01+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE_TURNER' )
      XEDGE_TURNER = 0.d0

      ! Allocate array to hold TURNER grid box lat edges
      ALLOCATE( YEDGE_TURNER( J01x01+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_TURNER' )
      YEDGE_TURNER = 0.d0

      ! Allocate array to hold GEOS-Chem grid box lon edges
      ALLOCATE( XEDGE_MODELG( IIIPAR0+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE_MODELG' )
      XEDGE_MODELG = 0.d0

      ! Allocate array to hold GEOS-Chem grid box lat edges
      ALLOCATE( YEDGE_MODELG( JJJPAR0+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_MODELG' )
      YEDGE_MODELG = 0.d0

      ! Emission arrays
      ALLOCATE( C2H6an( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H6an' )
      C2H6an = 0d0

      ALLOCATE( C2H6bf( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H6bf' )
      C2H6bf = 0d0

      ALLOCATE( C3H8an( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C3H8an' )
      C3H8an = 0d0

      ALLOCATE( C3H8bf( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C3H8bf' )
      C3H8bf = 0d0

      END SUBROUTINE INIT_TURNER_FLUX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_turner_flux
!
! !DESCRIPTION: Subroutine CLEANUP\_TURNER\_FLUX deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TURNER_FLUX
!
! !REVISION HISTORY: 
!  13 May 2016 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TURNER_FLUX begins here!
      !=================================================================
      IF ( ALLOCATED( C2H6an   ) ) DEALLOCATE( C2H6an  )
      IF ( ALLOCATED( C2H6bf   ) ) DEALLOCATE( C2H6bf  )
      IF ( ALLOCATED( C3H8an   ) ) DEALLOCATE( C3H8an  )
      IF ( ALLOCATED( C3H8bf   ) ) DEALLOCATE( C3H8bf  )

      END SUBROUTINE CLEANUP_TURNER_FLUX
!EOC
      END MODULE TURNER_FLUX_MOD
