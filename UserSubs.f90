   ! NOTE: This source file contains dummy placeholders for ALL of the
   !       user-specified routines available in FAST.  These routines
   !       are as follows:
   !          Routine       Description
   !          ------------  ---------------------------------------------------
   !          PitchCntrl()  User-specified blade pitch control (either
   !                        independent or rotor-collective) model.
   !          UserGen()     User-specified generator torque and power model.
   !          UserHSSBr()   User-specified high-speed shaft brake model.
   !          UserPtfmLd()  User-specified platform loading model.
   !          UserRFrl()    User-specified rotor-furl spring/damper model.
   !          UserTeet()    User-specified rotor-teeter spring/damper model.
   !          UserTFin()    User-specified tail fin aerodynamics model.
   !          UserTFrl()    User-specified tail-furl spring/damper model.
   !          UserVSCont()  User-specified variable-speed torque and power
   !                        control model.
   !          UserYawCont() User-specified nacelle-yaw control model.
   !       In order to interface FAST with your own user-specified routines,
   !       you can develop your own logic within these dummy placeholders and
   !       recompile FAST; OR comment out the appropriate dummy placeholders,
   !       create your own routines in their own source files, and recompile
   !       FAST while linking in these additional source files.  For example,
   !       the executable version of FAST that is distributed with the FAST
   !       archive is linked with the example PitchCntrl() routine contained in
   !       source file PitchCntrl_ACH.f90 and the example UserGen() and
   !       UserVSCont() routines contained in source file UserVSCont_KP.f90;
   !       thus, the dummy placeholders for routines PitchCntrl(), UserGen(),
   !       and UserVSCont() are commented out within this source file.  The
   !       example pitch controller was written by Craig Hansen (ACH) and the
   !       example generator and variable speed controllers were written by
   !       Kirk Pierce (KP).  Please see the aforementioned source files for
   !       additional information on these example user-specified routines.

   ! NOTE: If you (the user) wants to access the current value of ANY of the
   !       output parameters available as outputs from FAST from your
   !       user-defined routines, then do the following:
   !          (1) USE MODULE Output() in your routine.
   !          (2) Access the output parameter by typing "AllOuts(OutName)",
   !              where OutName is the PRIMARY name of the output parameter.
   !              For example, to access the current value of the in-plane
   !              bending moment at the root of blade 1 (in kN·m), type in
   !              "AllOuts(RootMxc1)", since RootMxc1 is the primary name of
   !              this output parameter--RootMIP1 will not work in place of
   !              RootMxc1, since it is a SECONDARY name.  Also, you CANNOT use
   !              the prefixes ("-", "_", "m", or "M") in front of OutName to
   !              reverse the sign of the selected output channel.
   !       Note that OutName DOES NOT have to be one of the output parameters
   !       you listed in OutList from your primary input file.  Also note that
   !       this technique WILL also work for user-defined routines written for
   !       ADAMS datasets extracted using the FAST-to-ADAMS preprocessor.

!=======================================================================
SUBROUTINE PitchCntrl ( BlPitch, ElecPwr, HSS_Spd, GBRatio, TwrAccel, NumBl, ZTime, DT, DirRoot, BlPitchCom_out )

   ! This pitch control routine was written by Eric Anderson of UC Davis. It is based on the baseline NREL
   ! 5MW pitch control scheme outlined in the NREL 5MW specifications with the addition the derating scheme
   ! outlined in Eric's dissertaion. Derating is achieved by using subroutine ControlParameters() to scale
   ! relevant control parameters.
USE								precision
USE								EAControl 	! contains variables: TimeDRStart, TimeDREnd, DerateFactor, TEmShutdown, maxOverspeed, EmergencyShutdown, GenSpeedF, PC_RefSpd, PC_MinPit, VS_Rgn2_K, and VS_RtPwr. See EAControl module in FAST_Mods.f90 for variable descriptions.


IMPLICIT                        NONE


   ! Passed variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

REAL(ReKi), INTENT(IN )      :: BlPitch   (NumBl)                               ! Current values of the blade pitch angles, rad.
REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(IN )      :: ElecPwr                                         ! Electrical power, watts.
REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
REAL(ReKi), INTENT(OUT)      :: BlPitchCom_out(NumBl)                               ! Commanded blade pitch angles (demand pitch angles), rad.
REAL(ReKi), INTENT(IN )      :: TwrAccel                                        ! Tower Acceleration, m/s^2.
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.

 ! Local Variables:

REAL(ReKi)                      :: ElapTime                                        ! Elapsed time since the last call to the controller, sec.
REAL(ReKi)                      :: GK                                              ! Current value of the gain correction factor, used in the gain scheduling law of the pitch controller, (-).
REAL(ReKi), SAVE                :: IntSpdErr                                       ! Current integral of speed error w.r.t. time, rad.
REAL(ReKi), SAVE                :: LastTimePC                                      ! Last time the pitch  controller was called, sec.
REAL(ReKi), PARAMETER           :: OnePlusEps    = 1.0 + EPSILON(OnePlusEps)       ! The number slighty greater than unity in single precision.
REAL(ReKi), PARAMETER           :: PC_DT         = 0.00125  !JASON:THIS CHANGED FOR ITI BARGE:      0.0001                    ! Communication interval for pitch  controller, sec.
REAL(ReKi), PARAMETER           :: PC_KI         =       0.008068634               ! Integral gain for pitch controller at rated pitch (zero), (-).
REAL(ReKi), PARAMETER           :: PC_KK         =       0.1099965                 ! Pitch angle were the the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero), rad.
REAL(ReKi), PARAMETER           :: PC_KP         =       0.01882681                ! Proportional gain for pitch controller at rated pitch (zero), sec.
REAL(ReKi), PARAMETER           :: PC_MaxPit     =       1.570796                  ! Maximum pitch setting in pitch controller, rad.
REAL(ReKi), PARAMETER           :: PC_MaxRat     =       0.1396263                 ! Maximum pitch  rate (in absolute value) in pitch  controller, rad/s.
REAL(ReKi), SAVE                :: PitCom    (3)                                   ! Commanded pitch of each blade the last time the controller was called, rad.
REAL(ReKi)                      :: PitComI                                         ! Integral term of command pitch, rad.
REAL(ReKi)                      :: PitComP                                         ! Proportional term of command pitch, rad.
REAL(ReKi)                      :: PitComT                                         ! Total command pitch based on the sum of the proportional and integral terms, rad.
REAL(ReKi)                      :: PitRate   (3)                                   ! Pitch rates of each blade based on the current pitch angles and current pitch command, rad/s.
REAL(ReKi)                      :: SpdErr                                          ! Current speed error, rad/s.
LOGICAL                   	 :: Initialize        = .TRUE.                                 ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
INTEGER(ReKi)                   :: K                                               ! Loops through blades.

!=======================================================================
   !Initialize variables:

   
IF ( Initialize ) THEN  ! .TRUE. if we're on the first call to the subroutine
   WRITE(*,*)  'Running with pitch control programmed by Eric Anderson '// &
              'in subroutine PitchCntrl(), which can be found in UserSubs.f90 '
	Initialize = .FALSE.
   ! Initialize the SAVEd variables:
   PitCom     = BlPitch                         ! This will ensure that the variable speed controller picks the correct control region and the pitch controller pickes the correct gain on the first call. If pitchCtrl()is called before UsrVSCtrl() initializing it here will work, if not I need to do something else.
   GK         = 1.0/( 1.0 + PitCom(1)/PC_KK )   ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   IntSpdErr  = PitCom(1)/( GK*PC_KI )          ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   LastTimePC = ZTime - PC_DT                    	! This will ensure that the pitch  controller is called on the first pass 
ENDIF

!=======================================================================
   ! Pitch control:

   ! Compute the elapsed time since the last call to the controller:
ElapTime = ZTime - LastTimePC

   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the pitch controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when PC_DT = DT, even in the presence of
   !       numerical precision errors.

IF ( ( ZTime*OnePlusEps - LastTimePC ) >= PC_DT )  THEN

	CALL updateControlParameters( HSS_Spd, ZTime )
	
	IF ( EmergencyShutdown ) THEN
		PitComT = 3.1415926535/2
	ELSE
	   ! Compute the gain scheduling correction factor based on the previously
	   !   commanded pitch angle for blade 1:
		  GK = 1.0/( 1.0 + PitCom(1)/PC_KK )

	   ! Compute the current speed error and its integral w.r.t. time; 
		  SpdErr    = GenSpeedF - PC_RefSpd                                 ! Current speed error
		  IntSpdErr = IntSpdErr + SpdErr*ElapTime                           ! Current integral of speed error w.r.t. time
	  
		! saturate the integral term using the pitch angle limits:
		  IntSpdErr = MIN( MAX( IntSpdErr, PC_MinPit/( GK*PC_KI ) ) , PC_MaxPit/( GK*PC_KI ))    ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits

		! Compute the pitch commands associated with the proportional and integral gains: 
		  PitComP   = GK*PC_KP*   SpdErr                                    ! Proportional term
		  PitComI   = GK*PC_KI*IntSpdErr                                    ! Integral term (saturated)


		! Superimpose the individual commands to get the total pitch command;
		  PitComT   = PitComP + PitComI                                     ! Overall command (unsaturated)
	  
	  
		!   saturate the overall command using the pitch angle limits:
		  PitComT   = MIN( MAX( PitComT, PC_MinPit ), PC_MaxPit )           ! Saturate the overall command using the pitch angle limits
	ENDIF


   ! Saturate the overall commanded pitch using the pitch rate limit:
   ! NOTE: Since the current pitch angle may be different for each blade
   !       (depending on the type of actuator implemented in the structural
   !       dynamics model), this pitch rate limit calculation and the
   !       resulting overall pitch angle command may be different for each
   !       blade.
      DO K = 1,NumBl ! Loop through all blades
         PitRate(K) = ( PitComT - BlPitch(K) )/ElapTime                 ! Pitch rate of blade K (unsaturated)
         PitRate(K) = MIN( MAX( PitRate(K), -PC_MaxRat ), PC_MaxRat )   ! Saturate the pitch rate of blade K using its maximum absolute value
         PitCom (K) = BlPitch(K) + PitRate(K)*ElapTime                  ! Saturate the overall command of blade K using the pitch rate limit
      ENDDO          ! K - all blades


   ! Reset the value of LastTimePC to the current value:
      LastTimePC = ZTime
ENDIF

BlPitchCom_out = PitCom                   ! Pass the most recent blade pitch command out of the subroutine
 

RETURN
END SUBROUTINE PitchCntrl
!=======================================================================
SUBROUTINE UserGen ( HSS_Spd, GBRatio, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! generator torque and power model.  Modify this code to create your
!   ! own model.
!
!   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
!   !       runtime from this user-defined routine, then do the following:
!   !          (1) USE MODULE DOFs().
!   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
!   !              .FALSE. depending on whether you want to turn-on or turn-off
!   !              the DOF, respectively.  Turning off the DOF forces the
!   !              current RATE to remain fixed.  If the rate is currently zero,
!   !              the current POSITION will remain fixed as well.
!   !       Note that this technique WILL NOT work for user-defined routines
!   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
!   !       preprocessor.
!
!
USE                             Precision
USE                            NWTC_Library

!
!
IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
REAL(ReKi), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
CALL ProgAbort ( ' Error: There User defined Generator [ subroutine UserGen() in UserSubs.f90]. Please select a different option for GenModel' )
!
GenTrq  = 0.0 + DelGenTrq  ! Make sure to add the pertubation on generator torque, DelGenTrq.  This is used only for FAST linearization (it is zero otherwise).
!
!
!   ! The generator efficiency is either additive for motoring,
!   !   or subtractive for generating power.
!
IF ( GenTrq > 0.0 )  THEN
   ElecPwr = GenTrq*HSS_Spd*GenEff
ELSE
   ElecPwr = GenTrq*HSS_Spd/GenEff
ENDIF
!
!
!
RETURN
END SUBROUTINE UserGen
!=======================================================================
SUBROUTINE UserHSSBr ( GenTrq, ElecPwr, HSS_Spd, GBRatio, NumBl, ZTime, DT, DirRoot, HSSBrFrac )

	! wrote a very simple High speed brake torque control routine. If the high speed shaft 
	! comes to a complete stop the high speed shaft brake will be engaged at full torque. 
	! Otherwise the brake will not be used. The intention is to keep the rotor stationary 
	! after an emergency shutdown, but not use the brake durring the emergency shutdown or 
	! in normal operation. (Eric Anderson)
	
	
   ! This is a dummy routine for holding the place of a user-specified
   ! HSS brake model.  This routine must specify the fraction
   ! (HSSBrFrac) of full torque to be applied to the HSS by the HSS
   ! brake.  The magnitude of the full torque (HSSBrFrac = 1.0) equals
   ! HSSBrTqF from the primary input file.  Modify this code to create
   ! your own model.

   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision
USE								EAControl 	! contains LOGICAL variable EmergencyShutdown.
! USE								DriveTrain


IMPLICIT                        NONE


   ! Passed Variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(IN )      :: ElecPwr                                         ! Electrical power (account for losses), watts.
REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
REAL(ReKi), INTENT(IN )      :: GenTrq                                          ! Electrical generator torque, N-m.
REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
REAL(ReKi), INTENT(OUT)      :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-).
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.

REAL(ReKi), SAVE			:: brakeStartTime									! Time when HHS Brake is initiated.
LOGICAL, SAVE				:: brakeOff = .TRUE.
REAL(ReKi), PARAMETER		:: HSSBrDT = 0.6									! Time it takes for HSS brake to reach full deployment once deployed.
	
IF  ( ( EmergencyShutdown ) .AND. ( GenSpeedF < 1 ) ) THEN ! If emergency shutdown has been initiated and generator speed is almost zero.

	IF ( brakeOff ) THEN
		brakeStartTime = ZTime
		brakeOff = .FALSE.
		WRITE(*,*)  'HSS Brake initiated at T =',ZTime,' GenSpeedF =',GenSpeedF
		HSSBrFrac = 0.0 
	ELSEIF ( (ZTime-brakeStartTime) < HSSBrDT ) THEN
		HSSBrFrac = (ZTime-brakeStartTime)/HSSBrDT
	ELSE
		HSSBrFrac = 1.0 !Engage brake
	ENDIF
ELSE
	HSSBrFrac = 0.0   
ENDIF


RETURN
END SUBROUTINE UserHSSBr
!=======================================================================
SUBROUTINE UserPtfmLd ( X, XD, ZTime, DirRoot, PtfmAM, PtfmFt )


   ! This is a dummy routine for holding the place of a user-specified
   ! platform loading model.  Modify this code to create your own model.
   ! The local variables and associated calculations below provide a
   ! template for making this user-specified platform loading model
   ! include linear 6x6 damping and stiffness matrices.  These are
   ! provided as an example only and can be modified or deleted as
   ! desired by the user without detriment to the interface (i.e., they
   ! are not necessary for the interface).

   ! The platform loads returned by this routine should contain contributions
   !   from any external load acting on the platform other than loads
   !   transmitted from the wind turbine.  For example, these loads should
   !   contain contributions from foundation stiffness and damping [not
   !   floating] or mooring line restoring and damping [floating], as well as
   !   hydrostatic and hydrodynamic contributions [offshore].  The platform
   !   loads will be applied on the platform at the instantaneous platform
   !   reference position within FAST and ADAMS.

   ! This routine assumes that the platform loads are transmitted through a
   !   medium like soil [foundation] and/or water [offshore], so that added
   !   mass effects are important.  Consequently, the routine assumes that the
   !   total platform load can be written as:
   !
   ! PtfmF(i) = SUM( -PtfmAM(i,j)*XDD(j), j=1,2,..,6) + PtfmFt(i) for i=1,2,...,6
   !
   ! where,
   !   PtfmF(i)    = the i'th component of the total load applied on the
   !                 platform; positive in the direction of positive motion of
   !                 the i'th DOF of the platform
   !   PtfmAM(i,j) = the (i,j) component of the platform added mass matrix
   !                 (output by this routine)
   !   XDD(j)      = the j'th component of the platform acceleration vector
   !   PtfmFt(i)   = the i'th component of the portion of the platform load
   !                 associated with everything but the added mass effects;
   !                 positive in the direction of positive motion of the i'th
   !                 DOF of the platform (output by this routine)

   ! The order of indices in all arrays passed to and from this routine is as
   !   follows:
   !      1 = Platform surge / xi-component of platform translation (internal DOF index = DOF_Sg)
   !      3 = Platform sway  / yi-component of platform translation (internal DOF index = DOF_Sw)
   !      3 = Platform heave / zi-component of platform translation (internal DOF index = DOF_Hv)
   !      4 = Platform roll  / xi-component of platform rotation    (internal DOF index = DOF_R )
   !      5 = Platform pitch / yi-component of platform rotation    (internal DOF index = DOF_P )
   !      6 = Platform yaw   / zi-component of platform rotation    (internal DOF index = DOF_Y )

   ! NOTE: The added mass matrix returned by this routine, PtfmAM, must be
   !       symmetric.  FAST and ADAMS will abort otherwise.
   !
   !       Please also note that the hydrostatic restoring contribution to the
   !       hydrodynamic force returned by this routine should not contain the
   !       effects of body weight, as is often done in classical marine
   !       hydrodynamics.  The effects of body weight are included within FAST
   !       and ADAMS.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: PtfmAM (6,6)                                    ! Platform added mass matrix, kg, kg-m, kg-m^2.
REAL(ReKi), INTENT(OUT)      :: PtfmFt   (6)                                    ! The 3 components of the portion of the platform force (in N  ) acting at the platform reference and the 3 components of the portion of the platform moment (in N-m  ) acting at the platform reference associated with everything but the added-mass effects; positive forces are in the direction of motion.
REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame.
REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame.
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.


   ! Local Variables:

REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: J                                               ! Generic index.



Damp  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

Stff  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

PtfmAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
PtfmAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

PtfmFt(1)   = 0.0
PtfmFt(2)   = 0.0
PtfmFt(3)   = 0.0
PtfmFt(4)   = 0.0
PtfmFt(5)   = 0.0
PtfmFt(6)   = 0.0

DO J = 1,6
   DO I = 1,6
      PtfmFt(I) = PtfmFt(I) - Damp(I,J)*XD(J) - Stff(I,J)*X(J)
   ENDDO
ENDDO



RETURN
END SUBROUTINE UserPtfmLd
!=======================================================================
SUBROUTINE UserRFrl ( RFrlDef, RFrlRate, ZTime, DirRoot, RFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! rotor-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the rotor-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_RFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the rotor-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: RFrlDef                                         ! Rotor-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: RFrlMom                                         ! Rotor-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: RFrlRate                                        ! Rotor-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



RFrlMom = 0.0



RETURN
END SUBROUTINE UserRFrl
!=======================================================================
SUBROUTINE UserTeet ( TeetDef, TeetRate, ZTime, DirRoot, TeetMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! teeter spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the teeter DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_Teet) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the teeter hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TeetDef                                         ! Rotor-teeter angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TeetMom                                         ! Rotor-teeter restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TeetRate                                        ! Rotor-teeter angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TeetMom = 0.0



RETURN
END SUBROUTINE UserTeet
!=======================================================================
SUBROUTINE UserTFin ( TFrlDef , TFrlRate, ZTime   , DirRoot, &
                      TFinCPxi, TFinCPyi, TFinCPzi,          &
                      TFinCPVx, TFinCPVy, TFinCPVz,          &
                      TFinAOA , TFinQ   ,                    &
                      TFinCL  , TFinCD  ,                    &
                      TFinKFx , TFinKFy                        )


   ! This is a dummy routine for holding the place of user-specified
   ! computations for tail fin aerodynamic loads.  Modify this code to
   ! create your own logic.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: TFinAOA                                         ! Angle-of-attack between the relative wind velocity and tail fin chordline, rad.
REAL(ReKi), INTENT(OUT)      :: TFinCD                                          ! Tail fin drag            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(OUT)      :: TFinCL                                          ! Tail fin lift            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(IN )      :: TFinCPVx                                        ! Absolute Velocity of the tail center-of-pressure along tail fin chordline pointing toward tail fin trailing edge, m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVy                                        ! Absolute Velocity of the tail center-of-pressure normal to plane of tail fin pointing towards suction surface   , m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVz                                        ! Absolute Velocity of the tail center-of-pressure in plane of tail fin normal to chordline and nominally upward  , m/s.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Improve the description of input arguments TFinCPxi, TFinCPyi, and
!jmj   TFinCPzi:
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi), INTENT(OUT)      :: TFinKFx                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge, N.
REAL(ReKi), INTENT(OUT)      :: TFinKFy                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface   , N.
REAL(ReKi), INTENT(OUT)      :: TFinQ                                           ! Dynamic pressure of the relative wind velocity, Pa.
REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFinAOA = 0.0
TFinCL  = 0.0
TFinCD  = 0.0
TFinQ   = 0.0
TFinKFx = 0.0
TFinKFy = 0.0



RETURN
END SUBROUTINE UserTFin
!=======================================================================
SUBROUTINE UserTFrl ( TFrlDef, TFrlRate, ZTime, DirRoot, TFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! tail-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the tail-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_TFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the tail-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TFrlMom                                         ! Tail-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFrlMom = 0.0



RETURN
END SUBROUTINE UserTFrl
!=======================================================================
SUBROUTINE UserVSCont ( HSS_Spd, GBRatio, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )

   ! This torque control routine was written by Eric Anderson of UC Davis. The torque control scheme 
   ! is based on the controller described in the NREL 5MW specifications, with the addition of a 
   ! derating scheme described in Eric's dissertation. Derating is achieved by using subroutine 
   ! ControlParameters() to scale relevant control parameters. WARNING: This routine does not include the 
   ! logic required to run linearizations.

USE								precision
USE  							Output
USE								EAControl 	! contains variables: TimeDRStart, TimeDREnd, DerateFactor, TEmShutdown, maxOverspeed, EmergencyShutdown, GenSpeedF, PC_RefSpd, PC_MinPit, VS_Rgn2_K, and VS_RtPwr. See EAControl module in FAST_Mods.f90 for variable descriptions.

IMPLICIT                        NONE

   ! Passed Variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
REAL(ReKi), INTENT(OUT)      :: GenTrq                                   		! Electrical generator torque, N-m.
REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.

! Local Variables:
REAL(ReKi)                      :: ElapTime                                        ! Elapsed time since the last call to the controller, sec.
REAL(ReKi), SAVE                :: LastGenTrq                                      ! Commanded electrical generator torque the last time the controller was called, N-m.
REAL(ReKi), SAVE                :: LastTimeVS                                  	! Last time the torque controller was called, sec.
REAL(ReKi), PARAMETER           :: OnePlusEps    = 1.0 + EPSILON(OnePlusEps)       ! The number slighty greater than unity in single precision.
REAL(ReKi)                      :: TrqRate                                         ! Torque rate based on the current and last torque commands, N-m/s.
REAL(ReKi), PARAMETER           :: VS_Rgn3MP     =       0.01745329                ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed, rad. -- chosen to be 1.0 degree above PC_MinPit
REAL(ReKi), PARAMETER           :: VS_CtInSp     =      70.16224                   ! Transitional generator speed (HSS side) between regions 1 and 1 1/2, rad/s.
REAL(ReKi), PARAMETER           :: VS_DT         = 0.00125  !JASON:THIS CHANGED FOR ITI BARGE:      0.0001                    ! Communication interval for torque controller, sec.
REAL(ReKi), PARAMETER           :: VS_MaxRat     =   15000.0                       ! Maximum torque rate (in absolute value) in torque controller, N-m/s.
REAL(ReKi), PARAMETER           :: VS_MaxTq      =   47402.91                      ! Maximum generator torque in Region 3 (HSS side), N-m. -- chosen to be 10% above VS_RtTq = 43.09355kNm
REAL(ReKi), PARAMETER           :: VS_Rgn2Sp     =      91.21091                   ! Transitional generator speed (HSS side) between regions 1 1/2 and 2, rad/s.
REAL(ReKi)           			 :: VS_RtGnSp                         				! Rated generator speed (HSS side), rad/s. -- cequal to 99% of PC_RefSpd
REAL(ReKi)                		 :: VS_Slope15                                      ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
REAL(ReKi)                		 :: VS_Slope25                                      ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
REAL(ReKi), PARAMETER           :: VS_SlPc       =      10.0                       ! Rated generator slip percentage in Region 2 1/2, %.
REAL(ReKi)                		 :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator, rad/s.
REAL(ReKi)                		 :: VS_TrGnSp                                       ! Transitional generator speed (HSS side) between regions 2 and 2 1/2, rad/s.
REAL(ReKi)    					 :: BlPitchCom                  					! Commanded blade pitch angle for blade 1 (demand pitch angle), rad.
! REAL(ReKi), PARAMETER           :: R2D           =      57.295780                  ! Factor to convert radians to degrees.

LOGICAL, SAVE					:: Initialize1 = .TRUE.					!Flag used to initialize some saved variables on the first call to this subroutine
LOGICAL, SAVE					:: Initialize2 = .TRUE.					!Flag used to initialize some saved variables on the first call to this subroutine


!=======================================================================
   ! Initialize saved variables on first call to subroutine
	IF ( Initialize1 ) THEN
	   WRITE(*,*)  'Running with torque control programmed by Eric Anderson '// &
              'in subroutine UserVSCont(), which can be found in UserSubs.f90 '
	Initialize1 = .FALSE.
		! NOTE: LastGenTrq, though SAVEd, is initialized below for simplicity, not here.
		LastTimeVS = ZTime - VS_DT                    ! This will ensure that the torque controller is called on the first pass 
	ENDIF
!=======================================================================
	! Variable-speed torque control:

   ! Compute the elapsed time since the last call to the controller:
   ElapTime = ZTime - LastTimeVS


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the torque controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when VS_DT = DT, even in the presence of
   !       numerical precision errors.

IF ( ( ZTime*OnePlusEps - LastTimeVS ) >= VS_DT )  THEN

	! Get up to date control parameters
	CALL updateControlParameters( HSS_Spd, ZTime )
	
	IF ( EmergencyShutdown ) THEN
		GenTrq = 0
	ELSE

	   ! Determine some torque control parameters not specified directly:
		VS_RtGnSp = 0.99*PC_RefSpd
		VS_SySp    = VS_RtGnSp/( 1.0 +  0.01*VS_SlPc )
		VS_Slope15 = ( VS_Rgn2_K*VS_Rgn2Sp*VS_Rgn2Sp )/( VS_Rgn2Sp - VS_CtInSp )
		VS_Slope25 = ( VS_RtPwr/VS_RtGnSp           )/( VS_RtGnSp - VS_SySp   )
		VS_TrGnSp = ( VS_Slope25 - SQRT( VS_Slope25*( VS_Slope25 - 4.0*VS_Rgn2_K*VS_SySp ) ) )/( 2.0*VS_Rgn2_K ) !Transition speed from region 2 to region 2.5

		BlPitchCom =  AllOuts(PtchPMzc1)/R2D
	   ! Compute the generator torque, which depends on which region we are in:

		  IF ( (   GenSpeedF >= VS_RtGnSp ) .OR. (  BlPitchCom >= (VS_Rgn3MP + PC_MinPit ) ) ) THEN ! We are in region 3 - power is constant
			 GenTrq = VS_RtPwr/GenSpeedF
		  ELSEIF ( GenSpeedF <= VS_CtInSp )  THEN                                    ! We are in region 1 - torque is zero
			 GenTrq = 0.0
		  ELSEIF ( GenSpeedF <  VS_Rgn2Sp )  THEN                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
			 GenTrq = VS_Slope15*( GenSpeedF - VS_CtInSp )
		  ELSEIF ( GenSpeedF <  VS_TrGnSp )  THEN                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
			 GenTrq = VS_Rgn2_K*GenSpeedF*GenSpeedF
		  ELSE                                                                       ! We are in region 2 1/2 - simple induction generator transition region
			 GenTrq = VS_Slope25*( GenSpeedF - VS_SySp   )
		  ENDIF

	   ! Saturate the commanded torque using the maximum torque limit:
		  GenTrq  = MIN( GenTrq                    , VS_MaxTq  )   ! Saturate the command using the maximum torque limit

		!Initialize saved variables on first call to subroutine
		  IF ( Initialize2 ) THEN 
			Initialize2 = .FALSE.
			LastGenTrq = GenTrq                 ! Initialize the value of LastGenTrq on the first pass only
		  ENDIF
	  ENDIF
	  
	! Saturate the commanded torque using the torque rate limit:  
      TrqRate = ( GenTrq - LastGenTrq )/ElapTime               ! Torque rate (unsaturated)
      TrqRate = MIN( MAX( TrqRate, -VS_MaxRat ), VS_MaxRat )   ! Saturate the torque rate using its maximum absolute value
      GenTrq  = LastGenTrq + TrqRate*ElapTime                  ! Saturate the command using the torque rate limit

   ! Reset the values of LastTimeVS and LastGenTrq to the current values:

      LastTimeVS = ZTime
      LastGenTrq = GenTrq

ELSE
	GenTrq = LastGenTrq
ENDIF
   
IF ( GenTrq > 0.0 )  THEN
   ElecPwr = GenTrq*HSS_Spd*GenEff
ELSE
   ElecPwr = GenTrq*HSS_Spd/GenEff
ENDIF
   
RETURN
END SUBROUTINE UserVSCont
!=======================================================================
SUBROUTINE UserYawCont ( YawPos, YawRate, WindDir, YawError, NumBl, ZTime, DT, DirRoot, YawPosCom, YawRateCom )


   ! This is a dummy routine for holding the place of a user-specified
   ! nacelle-yaw controller.  Modify this code to create your own device.


   ! As indicated, the yaw controller must always specify a command (demand)
   !   yaw angle, YawPosCom, AND command (demand) yaw rate, YawRateCom.
   !   Normally, you should correlate these commands so that the commanded yaw
   !   angle is the integral of the commanded yaw rate, or likewise, the
   !   commanded yaw rate is the derivative of the commanded yaw angle.  FAST
   !   WILL NOT compute these correlations for you and DOES NOT check to
   !   ensure that they are correlated.  In some situations, it is desirable to
   !   set one of the commands (either yaw angle OR yaw rate) to ZERO depending
   !   on the desired transfer function of FAST's built-in actuator model (see
   !   below for a discussion of FAST's built-in actuator model).  In general,
   !   the commanded yaw angle and rate SHOULD NEVER be defined independent of
   !   each other with BOTH commands NONZERO.


   ! The yaw controller's effect on the FAST model depends on whether or not
   !   the yaw DOF is enabled as follows:
   !
   ! YawDOF = False - If the yaw DOF is disabled, then the commanded yaw angle
   !                  and rate will be the ACTUAL yaw angle and yaw rate used
   !                  internally by FAST (in general, you should ensure these
   !                  are correlated).  In this case, any desired actuator
   !                  effects should be built within this controller.  Also in
   !                  this case, FAST WILL NOT compute the correlated yaw
   !                  acceleration, but assume that it is ZERO.  If the
   !                  commanded yaw rate is zero while the commanded yaw angle
   !                  is changing in time, then the yaw controller's effect
   !                  on yaw angle is the identical to routine PitchCntrl()'s
   !                  effect on pitch angle (i.e., routine PitchCntrl()
   !                  commands changes in pitch angle with no associated
   !                  changes in pitch rate or pitch acceleration).  For yaw
   !                  control, this situation should be avoided however, since
   !                  yaw-induced gyroscopic pitching loads on the turbine
   !                  brought about by the yaw rate may be significant.
   !
   ! YawDOF = True  - If the yaw DOF is enabled, then the commanded yaw angle
   !                  and rate, YawPosCom and YawRateCom, become the neutral
   !                  yaw angle, YawNeut, and neutral yaw rate, YawRateNeut, in
   !                  FAST's built-in second-order actuator model defined by
   !                  inputs YawSpr and YawDamp.


   ! Description of FAST's built-in actuator model:
   !
   ! In the time-domain, FAST's built-in actuator model is defined as follows:
   !
   ! YawIner*YawAccel + YawDamp*YawRate + YawSpr*YawPos
   !                             = YawDamp*YawRateNeut + YawSpr*YawNeut + YawTq
   !
   ! so that the transmitted torque is:
   !
   ! YawMom = YawSpr*( YawPos - YawNeut ) + YawDamp*( YawRate - YawRateNeut )
   !
   ! where,
   !   YawSpr      = nacelle-yaw spring constant (defined in FAST's primary
   !                 input file)
   !   YawDamp     = nacelle-yaw damping constant (defined in FAST's primary
   !                 input file)
   !   YawIner     = instantaneous inertia of the nacelle and rotor about the
   !                 yaw axis
   !   YawNeut     = the commanded (neutral) yaw angle = YawPosCom
   !   YawRateNeut = the commanded (neutral) yaw rate  = YawRateCom
   !   YawPos      = yaw angle (position)
   !   YawRate     = yaw rate
   !   YawAccel    = yaw acceleration
   !   YawTq       = torque about the yaw axis applied by external forces above
   !                 the yaw bearing, such as wind loading
   !   YawMom      = torque transmitted through the yaw bearing
   !
   ! If the commanded yaw angle and rate are correlated (so that the commanded
   !   yaw angle is the integral of the commanded yaw rate, or likewise, the
   !   commanded yaw rate is the derivative of the commanded yaw angle), then
   !   FAST's built-in second-order actuator model will have the following
   !   characteristic transfer function:
   !
   !               YawDamp*s + YawSpr             2*Zeta*OmegaN*s + OmegaN^2
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2
   !
   ! where,
   !   T(s)    = the transfer function of FAST's built-in 2nd order actuator
   !             model
   !   OmegaN  = SQRT(YawSpr/YawIner) = yaw actuator natural frequency
   !   Zeta    = YawDamp/(2*SQRT(YawSpr*YawIner)) = yaw actuator damping ratio
   !             in fraction of critical
   !
   ! If only the yaw angle is commanded, and YawRateCom is zeroed, then the
   !   charecteristic transfer function of FAST's built-in second-order
   !   actuator model simplifies to:
   !
   !                     YawSpr                            OmegaN^2
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2
   !
   ! If only the yaw rate is commanded, and YawPosCom is zeroed, then the
   !   charecteristic transfer function of FAST's built-in second-order
   !   actuator model simplifies to:
   !
   !                    YawDamp                         2*Zeta*OmegaN
   ! T(s) = -------------------------------- = --------------------------------
   !        YawIner*s^2 + YawDamp*s + YawSpr   s^2 + 2*Zeta*OmegaN*s + OmegaN^2


   ! NOTE: If you (the user) wants to switch on-or-off the yaw DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_Yaw) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF acts is like
   !              setting YawDOF to False.
   !       This technique is useful, for example, if the yaw bearing has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).

REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
REAL(ReKi), INTENT(IN )      :: WindDir                                         ! Current horizontal hub-height wind direction (positive about the zi-axis), rad.
REAL(ReKi), INTENT(IN )      :: YawError                                        ! Current nacelle-yaw error estimate (positve about the zi-axis), rad.
REAL(ReKi), INTENT(IN )      :: YawPos                                          ! Current nacelle-yaw angular position, rad.
REAL(ReKi), INTENT(OUT)      :: YawPosCom                                       ! Commanded nacelle-yaw angular position (demand yaw angle), rad.
REAL(ReKi), INTENT(IN )      :: YawRate                                         ! Current nacelle-yaw angular rate, rad/s.
REAL(ReKi), INTENT(OUT)      :: YawRateCom                                      ! Commanded nacelle-yaw angular rate (demand yaw rate), rad/s.
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



YawPosCom  = 0.0
YawRateCom = 0.0



RETURN
END SUBROUTINE UserYawCont
!=======================================================================
SUBROUTINE updateControlParameters( HSS_Spd, ZTime )

USE								precision
USE								EAControl 	! contains variables: TimeDRStart, TimeDREnd, DerateFactor, TEmShutdown, maxOverspeed, EmergencyShutdown, GenSpeedF, PC_RefSpd, PC_MinPit, VS_Rgn2_K, and VS_RtPwr. See EAControl module in FAST_Mods.f90 for variable descriptions.

IMPLICIT                        NONE

   ! Passed Variables:
REAL(ReKi), INTENT(IN)       		:: HSS_Spd                     ! Current  HSS (generator) speed, rad/s.
REAL(ReKi), INTENT(IN )      		:: ZTime                        ! Current simulation time, sec.

	! Local variables storing baseline control parameters (FROM the controller described in the NREL 5MW specifications)
REAL(ReKi), PARAMETER				:: PC_MinPit_baseline     = 0.0        		! Minimum pitch setting in NREL 5MW baseline pitch controller, rad.
REAL(ReKi), PARAMETER           	:: PC_RefSpd_baseline     = 122.9096        ! Desired (reference) HSS speed for NREL 5MW baseline pitch controller, rad/s.
REAL(ReKi), PARAMETER           	:: VS_Rgn2K_baseline      = 2.332287        ! Generator torque constant in Region 2 (HSS side) for NREL 5MW baseline controller, N-m/(rad/s)^2.
REAL(ReKi), PARAMETER           	:: VS_RtGnSp_baseline     = 121.6805        ! Rated generator speed (HSS side) for NREL 5MW baseline controller, rad/s. -- chosen to be 99% of PC_RefSpd
REAL(ReKi), PARAMETER			    :: VS_RtPwr_baseline      = 5296610.0       ! Rated generator generator power in Region 3 for NREL 5MW baseline controller, Watts. -- chosen to be 5MW divided by tthe electrical generator efficiency of 94.ReKi%
REAL(ReKi), PARAMETER			    :: VS_Rgn2Sp_baseline	  = 91.21091        ! Transitional generator speed (HSS side) between regions 1 1/2 and 2 for NREL 5MW baseline controller, rad/s.
REAL(ReKi), PARAMETER			    :: VS_CtInSp_baseline	  = 70.16224        ! Transitional generator speed (HSS side) between regions 1 and 1 1/2 for NREL 5MW baseline controller, rad/s.

	!Local variables used to filter generator speed
REAL(ReKi)                      :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
REAL(ReKi), SAVE                :: LastTime                                        ! Last time this subroutine was called, sec.
REAL(ReKi), PARAMETER           :: CornerFreq    =       1.570796                  ! Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1Hz = 0.25Hz = 1.570796rad/s)

	!Local variables used for derate calculations
REAL(ReKi), PARAMETER				:: pDR = 0.2							!- poles of the second order derate input filter.
REAL(ReKi), SAVE                   :: FF_pwrFactor = 1.0 					! The derate factor. A fraction of 1, where 1 is not derated.

LOGICAL, SAVE					:: Initialize = .TRUE.					!Flag used to initialize some saved variables on the first call to this subroutine


!=======================================================================
   ! Initialize saved variables on first call to subroutine
	IF ( Initialize ) THEN
		WRITE(*,*)  'First call to subroutine updateControlParameters(), programmed by '// &
					'Eric Anderson. Subroutine can be found in UserSubs.f90 '
		Initialize = .FALSE.
		GenSpeedF  = HSS_Spd            ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
		LastTime   = ZTime               ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
	ENDIF
   
!=======================================================================

   ! Filter the HSS (generator) speed measurement:
   ! NOTE: This is a very simple recursive, single-pole, low-pass filter with
   !       exponential smoothing.
		! Update the coefficient in the recursive formula based on the elapsed time
		!   since the last call to the controller:
		Alpha     = EXP( ( LastTime - ZTime )*CornerFreq )
		! Apply the filter:
		GenSpeedF = ( 1.0 - Alpha )*HSS_Spd + Alpha*GenSpeedF
!=======================================================================
	!Derate Calculations
	IF (ZTime >= TimeDREnd) THEN
		!return turbine to normal operation
		FF_pwrFactor = 1.0 - DerateFactor + DerateFactor*(1.0 - pDR*(ZTime - TimeDREnd)*EXP(-pDR*(ZTime - TimeDREnd)) - EXP(-pDR*(ZTime - TimeDREnd)))
	ELSEIF (ZTime >= TimeDRStart) THEN
		!Derate turbine
		FF_pwrFactor = 1.0 - DerateFactor*(1.0 - pDR*(ZTime - TimeDRStart)*EXP(-pDR*(ZTime - TimeDRStart)) - EXP(-pDR*(ZTime - TimeDRStart)))
	ENDIF
!=======================================================================
	! Set pitch control parameters
	PC_RefSpd = PC_RefSpd_baseline*FF_pwrFactor
	PC_MinPit = PC_MinPit_baseline !correct this later
!=======================================================================
	! Set torque control parameters	
	VS_Rgn2_K = VS_Rgn2K_baseline/(FF_pwrFactor**2) 	! Region 2 torque constant
	VS_RtPwr = VS_RtPwr_baseline*FF_pwrFactor		! Rated power
!=======================================================================
	! Check to see if emergency shutdown should be initiated
	IF ((ZTime > 30.0) .AND. (EmergencyShutdown == .FALSE.)) THEN ! If simulation has run long enough to pass the initial transient behavior and an emergency shutdown hasn't been requested yet. 
		IF ( ( ZTime > TEmShutdown ) .OR. ( (100*(HSS_Spd-PC_RefSpd_baseline)/PC_RefSpd_baseline) .GE. maxOverspeed)) THEN ! Should an emergency shutdown be requested now?
			EmergencyShutdown = .TRUE.
			WRITE(*,*)  'Emergency shutdown requested at T =',ZTime, ' Overspeed =',(100*(HSS_Spd-PC_RefSpd_baseline)/PC_RefSpd_baseline) 
		ENDIF
	ENDIF
!=======================================================================
! Reset the value of LastTime to the current value:
   LastTime = ZTime
	
RETURN
END SUBROUTINE updateControlParameters
