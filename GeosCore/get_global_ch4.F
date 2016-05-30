! $Id: get_global_ch4.f,v 1.1 2009/09/16 14:06:26 bmy Exp $
      SUBROUTINE GET_GLOBAL_CH4( am_I_Root, Input_Opt )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_CH4 retrieves monthly surface concentrations of CH4 
!  from observations
!
!  Arguments as Input:
!  ===========================================================================
!
!  Arguments as Output:
!  ===========================================================================
!
!  NOTES:
!  (1 ) GET_GLOBAL_CH4 only has to be called at the start of the new year,
!        as long as A3090S, A0030S, A0030N, A3090N are saved in the
!        calling program (bmy, 1/3/01)
!  (2 ) Also need to compute yearly gradients for CH4 beyond 1997 --
!        will do this later (bmy, 1/3/01)
!  (3 ) Bug fix: add missing comma to FORMAT statement (bmy, 3/23/03)
!  (4 ) Place WRITE statments w/in an !$OMP CRITICAL block, so as to make
!        sure that only one processor at a time writes them.  Also now use
!        F90 REPEAT intrinsic function.  Also replaced old CH4 gradient values
!        with updated values for 1983-2001.  Use data for 2001 as a proxy for
!        years past 2001, since data for those years has not been reported
!        yet. (mje, bmy, 7/7/03)
!  (5 ) Split off from module "global_ch4_mod.f".  Updated for IPCC future
!        emissions scenarios. (swu, bmy, 5/30/06)     
!  (6 ) Add the preindustrial CH4 scenarios.  Also set 2001 as the default
!        in case we are running 2030 or 2050 met but present-day emissions.
!        (swu, havala, bmy, 1/25/08)
!  (7 ) Updated CH4 vales with version 2008-07-02 for 1983-2007. Also use
!        2007 for years past 2007 (jaf, 4/15/09) 
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LFUTURE instead of LFUTURE
!                              from logical_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  13 Jul 2013 - L. Murray - Replace zonal bands with new 2-D interpolated
!        maps, kriged from NOAA flask data.
!******************************************************************************
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCENARIO
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput

      USE REGRID_A2A_MOD,       ONLY : DO_REGRID_A2A
      USE DIRECTORY_MOD,        ONLY : DATA_DIR_1x1
      USE GRID_MOD,             ONLY : GET_YMID
      USE TIME_MOD,             ONLY : GET_YEAR, GET_MONTH

      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_close

      USE COMODE_MOD
      USE CMN_SIZE_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE                 :: PREVMONTH = -1
      INTEGER, SAVE                 :: PREVYEAR  = -1
      CHARACTER(LEN=2)              :: FUTURE_SCENARIO
      CHARACTER(LEN=256)            :: FILENAME, LLFILENAME
      INTEGER                       :: fileID, J
      REAL*4                        :: array_sp_1x1(I1x1,J1x1-1,1,1)
      REAL*8                        :: array_dp_1x1(I1x1,J1x1-1)
      INTEGER                       :: THISYEAR, THISMONTH

      !=================================================================
      ! GET_GLOBAL_CH4 begins here!
      !=================================================================

      ! What year do we want to use for methane?
      IF ( CH4_YEAR .eq. -1 ) THEN 
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = CH4_YEAR 
      ENDIF
      THISMONTH = GET_MONTH()
 
      ! Is this year and/or month different than what has been read?
      IF ( PREVYEAR .ne. THISYEAR .or. PREVMONTH .ne. THISMONTH ) THEN

         ! Get methane surface concentration at 1x1
         FILENAME = '../noaa.gmd.methane.nc'
100      FORMAT( '         => Reading from file: ', a )
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) TRIM( filename )
         ENDIF
       
         call NcOp_Rd( fileID, TRIM( FILENAME ) )
         call NcRd( array_sp_1x1, fileID, 'SFC_CH4', 
     &        (/   1,       1, THISMONTH, THISYEAR - 1978 /), ! Start
     &        (/ I1x1, J1x1-1,         1,               1 /)  ) ! Count
         array_dp_1x1(:,:) = array_sp_1x1(:,:,1,1)
         call NcCl( fileID )
       
         ! Regrid to model resolution
         LLFILENAME = TRIM( DATA_DIR_1x1) //
     &        'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'
         CALL DO_REGRID_A2A( LLFILENAME, I1x1, J1x1-1, 
     &        array_dp_1x1, SURF_CH4, 0, netCDF=.TRUE. )

      !=================================================================
      ! Print the latitudinal CH4 gradient for this month to stdout
      !=================================================================
!$OMP CRITICAL
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            
105         FORMAT( 'GET_GLOBAL_CH4: YEAR = ', i4, ' MONTH = ', i2 )
106         FORMAT( f7.2, ':  ', f7.1, ' ppbv' )     
            DO J=JJPAR,1,-1 
               WRITE( 6, 106   ) GET_YMID(1,J,1), 
     &              SUM(SURF_CH4(:,J)) / REAL(IIPAR)
            ENDDO
            
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         ENDIF
!$OMP END CRITICAL
         
         PREVYEAR = THISYEAR
         PREVMONTH = THISMONTH
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_CH4
