/**
 * Emotiv Development Kit (EDK) API
 * Copyright (c) 2009 Emotiv Systems, Inc.
 *
 * The main interface that allows interactions between external programs and the Emotiv detection engine.
 *
 * None of the API functions are thread-safe.
 *
 * This header file is designed to be includable under C and C++ environment.
 *
 */


#ifndef EDK_H
#define EDK_H

#include <string>
#include <sstream>
#include "edkErrorCode.h"
#include "EmoStateDLL.h"


#ifndef EDK_STATIC_LIB

#ifdef EDK_EXPORTS
#ifdef _WIN32
#define EDK_API __declspec(dllexport)
#else
#define EDK_API
#endif
#else
#ifdef _WIN32
#define EDK_API __declspec(dllimport)
#else
#define EDK_API
#endif
#endif

#else

#define EDK_API extern

#endif


extern "C"
{
	//! Expressiv Suite threshold type enumerator
	typedef enum EE_ExpressivThreshold_enum {
		EXP_SENSITIVITY
	} EE_ExpressivThreshold_t;

	//! Expressiv Suite training control enumerator
	typedef enum EE_ExpressivTrainingControl_enum {
		EXP_NONE = 0, EXP_START, EXP_ACCEPT, EXP_REJECT, EXP_ERASE, EXP_RESET
	} EE_ExpressivTrainingControl_t;

	//! Expressiv Suite signature type enumerator
	typedef enum EE_ExpressivSignature_enum {
		EXP_SIG_UNIVERSAL = 0, EXP_SIG_TRAINED
	} EE_ExpressivSignature_t;

	 //! Cognitiv Suite training control enumerator
	typedef enum EE_CognitivTrainingControl_enum {
		COG_NONE = 0, COG_START, COG_ACCEPT, COG_REJECT, COG_ERASE, COG_RESET
	} EE_CognitivTrainingControl_t;

//DEPLOYMENT::STABLE_RELEASE::REMOVE_START

	//! Cognitiv Suite level enumerator
	//@@ This constant has been obsoleted
	typedef enum EE_CognitivLevel_enum {
		COG_LEVEL1 = 1, COG_LEVEL2, COG_LEVEL3, COG_LEVEL4
	} EE_CognitivLevel_t;

//DEPLOYMENT::STABLE_RELEASE::REMOVE_END

	//! Handle to internal EmoState structure allocated by EE_EmoStateCreate()
	typedef void*         EmoStateHandle;

	//! Handle to internal event structure allocated by EE_EmoEngineEventCreate()
	typedef void*         EmoEngineEventHandle;

	//! Handle to internal event structure allocated by EE_OptimizationParamCreate()
	typedef void*         OptimizationParamHandle;

	typedef void*		  DataHandle;

	//! EmoEngine event types
	typedef enum EE_Event_enum {
		EE_UnknownEvent		= 0x0000,
		EE_EmulatorError	= 0x0001,
		EE_ReservedEvent	= 0x0002,
		EE_UserAdded		= 0x0010,
		EE_UserRemoved		= 0x0020,
		EE_EmoStateUpdated	= 0x0040,
		EE_ProfileEvent		= 0x0080,
		EE_CognitivEvent	= 0x0100,
		EE_ExpressivEvent	= 0x0200,
		EE_InternalStateChanged = 0x0400,
		EE_AllEvent			= EE_UserAdded | EE_UserRemoved | EE_EmoStateUpdated | EE_ProfileEvent |
							  EE_CognitivEvent | EE_ExpressivEvent | EE_InternalStateChanged
	} EE_Event_t;

	//! Expressiv-specific event types
	typedef enum EE_ExpressivEvent_enum {
		EE_ExpressivNoEvent = 0, EE_ExpressivTrainingStarted, EE_ExpressivTrainingSucceeded,
		EE_ExpressivTrainingFailed, EE_ExpressivTrainingCompleted, EE_ExpressivTrainingDataErased,
		EE_ExpressivTrainingRejected, EE_ExpressivTrainingReset
	} EE_ExpressivEvent_t;
	
	//! Cognitiv-specific event types
	typedef enum EE_CognitivEvent_enum {
		EE_CognitivNoEvent = 0, EE_CognitivTrainingStarted, EE_CognitivTrainingSucceeded,
		EE_CognitivTrainingFailed, EE_CognitivTrainingCompleted, EE_CognitivTrainingDataErased,
		EE_CognitivTrainingRejected, EE_CognitivTrainingReset,
		EE_CognitivAutoSamplingNeutralCompleted, EE_CognitivSignatureUpdated
	} EE_CognitivEvent_t;

//DEPLOYMENT::NON_PREMIUM_RELEASE::REMOVE_START
	typedef enum EE_DataChannels_enum {
		ED_COUNTER = 0, ED_INTERPOLATED, ED_RAW_CQ,
		ED_AF3, ED_F7, ED_F3, ED_FC5, ED_T7, 
		ED_P7, ED_O1, ED_O2, ED_P8, ED_T8, 
		ED_FC6, ED_F4, ED_F8, ED_AF4, ED_GYROX, 
		ED_GYROY, ED_TIMESTAMP, ED_ES_TIMESTAMP, ED_FUNC_ID, ED_FUNC_VALUE, ED_MARKER, 
		ED_SYNC_SIGNAL
	} EE_DataChannel_t;
//DEPLOYMENT::NON_PREMIUM_RELEASE::REMOVE_END

	//! Input sensor description
	typedef struct InputSensorDescriptor_struct {
		EE_InputChannels_t channelId;  // logical channel id
		int                fExists;    // does this sensor exist on this headset model
		const char*        pszLabel;   // text label identifying this sensor
		double             xLoc;       // x coordinate from center of head towards nose
		double             yLoc;       // y coordinate from center of head towards ears
		double             zLoc;       // z coordinate from center of head toward top of skull
	} InputSensorDescriptor_t;


	//! Initializes the connection to EmoEngine. This function should be called at the beginning of programs that make use of EmoEngine, most probably in initialization routine or constructor.
	/*!	
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if a connection is established

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EngineConnect(const char* strDevID = "Emotiv Systems-5");

	
	//! Initializes the connection to a remote instance of EmoEngine.
	/*!
		Blocking call

		\param szHost - A null-terminated string identifying the hostname or IP address of the remote EmoEngine server
		\param port - The port number of the remote EmoEngine server
					- If connecting to the Emotiv Control Panel, use port 3008
					- If connecting to the EmoComposer, use port 1726
	
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if a connection is established

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EngineRemoteConnect(const char* szHost, unsigned short port, const char* strDevID = "Emotiv Systems-5");

	
	//! Terminates the connection to EmoEngine. This function should be called at the end of programs which make use of EmoEngine, most probably in clean up routine or destructor.
	/*!
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if disconnection is achieved

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EngineDisconnect();


	//! Controls the output of logging information from EmoEngine (which is off by default). This should only be enabled if instructed to do so by Emotiv developer support for the purposes of collecting diagnostic information.
	/*!
	    \param szFilename - The path of the logfile
		\param fEnable - If non-zero, then diagnostic information will be written to logfile.
		               - If zero, then all diagnostic information is suppressed (default).
	    \param nReserved - Reserved for future use.

	    \return EDK_ERROR_CODE
		        - EDK_ERROR_CODE = EDK_OK if the command succeeded
	*/
	EDK_API int
		EE_EnableDiagnostics(const char* szFilename, int fEnable, int nReserved);

	
	//! Returns a handle to memory that can hold an EmoEngine event. This handle can be reused by the caller to retrieve subsequent events.
	/*!
		\return EmoEngineEventHandle
	*/
	EDK_API EmoEngineEventHandle
		EE_EmoEngineEventCreate();


	//! Returns a handle to memory that can hold a profile byte stream. This handle can be reused by the caller to retrieve subsequent profile bytes.
	/*!
		\return EmoEngineEventHandle
	*/
	EDK_API EmoEngineEventHandle
		EE_ProfileEventCreate();

	
	//! Frees memory referenced by an event handle.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate() or EE_ProfileEventCreate()
	*/
	EDK_API void
		EE_EmoEngineEventFree(EmoEngineEventHandle hEvent);

	
	//! Returns a handle to memory that can store an EmoState. This handle can be reused by the caller to retrieve subsequent EmoStates.
	/*!
		\return EmoStateHandle
	*/
	EDK_API EmoStateHandle
		EE_EmoStateCreate();

	
	//! Frees memory referenced by an EmoState handle.
	/*!
		\param hState - a handle returned by EE_EmoStateCreate()
	*/
	EDK_API void
		EE_EmoStateFree(EmoStateHandle hState);


	//! Returns the event type for an event already retrieved using EE_EngineGetNextEvent.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate()
	
		\return EE_Event_t
	*/
	EDK_API EE_Event_t
		EE_EmoEngineEventGetType(EmoEngineEventHandle hEvent);

	
	//! Returns the Cognitiv-specific event type for an EE_CognitivEvent event already retrieved using EE_EngineGetNextEvent.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate()
	
		\return EE_CognitivEvent_t
	*/
	EDK_API EE_CognitivEvent_t
		EE_CognitivEventGetType(EmoEngineEventHandle hEvent);


	//! Returns the Expressiv-specific event type for an EE_ExpressivEvent event already retrieved using EE_EngineGetNextEvent.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate()
	
		\return EE_ExpressivEvent_t
	*/
	EDK_API EE_ExpressivEvent_t
		EE_ExpressivEventGetType(EmoEngineEventHandle hEvent);
	

	//! Retrieves the user ID for EE_UserAdded and EE_UserRemoved events.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate()
		\param pUserIdOut - receives the user ID associated with the current event
	
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful.

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EmoEngineEventGetUserId(EmoEngineEventHandle hEvent, unsigned int *pUserIdOut);

	
	//! Copies an EmoState returned with a EE_EmoStateUpdate event to memory referenced by an EmoStateHandle.
	/*!
		\param hEvent - a handle returned by EE_EmoEngineEventCreate() and populated with EE_EmoEngineGetNextEvent()
		\param hEmoState - a handle returned by EE_EmoStateCreate
	
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful.

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EmoEngineEventGetEmoState(EmoEngineEventHandle hEvent, EmoStateHandle hEmoState);
	

	//! Retrieves the next EmoEngine event
	/*!
		Non-blocking call

		\param hEvent - a handle returned by EE_EmoEngineEventCreate()

		\return EDK_ERROR_CODE
                <ul>
		        <li> EDK_ERROR_CODE = EDK_OK if a new event has been retrieved
				<li> EDK_ERROR_CODE = EDK_NO_EVENT if no new events have been generated by EmoEngine
				</ul>
		
		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EngineGetNextEvent(EmoEngineEventHandle hEvent);

	
	//! Clear a specific EmoEngine event type or all events currently inside the event queue. Event flags can be combined together as one argument except EE_UnknownEvent and EE_EmulatorError.
	/*!
		\param eventTypes - EmoEngine event type (EE_Event_t), multiple events can be combined such as (EE_UserAdded | EE_UserRemoved)

		\return EDK_ERROR_CODE
                <ul>
		        <li> EDK_ERROR_CODE = EDK_OK if the events have been cleared from the queue
				<li> EDK_ERROR_CODE = EDK_INVALID_PARAMETER if input event types are invalid
				</ul>
		
		\sa EE_Event_t, edkErrorCode.h
	*/
	EDK_API int
		EE_EngineClearEventQueue(int eventTypes);

	
	//! Retrieves number of active users connected to the EmoEngine.
	/*!
		\param pNumUserOut - receives number of users

		\return EDK_ERROR_CODE
		        - EDK_ERROR_CODE = EDK_OK if successful.

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_EngineGetNumUser(unsigned int* pNumUserOut);


	//! Sets the player number displayed on the physical input device (currently the USB Dongle) that corresponds to the specified user
	/*!
		\param userId - EmoEngine user ID
		\param playerNum - application assigned player number displayed on input device hardware (must be in the range 1-4)
		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_SetHardwarePlayerDisplay(unsigned int userId, unsigned int playerNum);


	//! Loads an EmoEngine profile for the specified user.  
	/*!
		\param userId - user ID
		\param profileBuffer - pointer to buffer containing a serialized user profile previously returned from EmoEngine.
		\param length - buffer size (number of bytes)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE if the function succeeds in loading this profile

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_SetUserProfile(unsigned int userId, const unsigned char profileBuffer[], unsigned int length);


	//! Returns user profile data in a synchronous manner.
	/*!
	     Fills in the event referred to by hEvent with an EE_ProfileEvent event
		 that contains the profile data for the specified user.

		 \param userId - user ID
		 \param hEvent - a handle returned by EE_EmoEngineEventCreate()

		 \return EDK_ERROR_CODE
		         - EDK_ERROR_CODE = EDK_OK if successful

		 \sa edkErrorCode.h
	*/
	EDK_API int
		EE_GetUserProfile(unsigned int userId, EmoEngineEventHandle hEvent);


	//! Returns a serialized user profile for a default user in a synchronous manner.
	/*!
	    Fills in the event referred to by hEvent with an EE_ProfileEvent event
		that contains the profile data for the default user

		 \param hEvent - a handle returned by EE_EmoEngineEventCreate()

		 \return EDK_ERROR_CODE
		         - EDK_ERROR_CODE = EDK_OK if successful

		 \sa edkErrorCode.h
	*/
	EDK_API int
		EE_GetBaseProfile(EmoEngineEventHandle hEvent);
	

	//! Returns the number of bytes required to store a serialized version of the requested user profile.
	/*!	
		\param hEvt - an EmoEngineEventHandle of type EE_ProfileEvent
		\param pProfileSizeOut - receives number of bytes required by the profile

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_GetUserProfileSize(EmoEngineEventHandle hEvt, unsigned int* pProfileSizeOut);

	
	//! Copies a serialized version of the requested user profile into the caller's buffer.
	/*!		
		\param hEvt - an EmoEngineEventHandle returned in a EE_ProfileEvent event
		\param destBuffer - pointer to a destination buffer
		\param length - the size of the destination buffer in bytes

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_GetUserProfileBytes(EmoEngineEventHandle hEvt, unsigned char destBuffer[], unsigned int length);

	//!  Loads a user profile from disk and assigns it to the specified user
	/*!		
		\param userID - a valid user ID
		\param szInputFilename - platform-dependent filesystem path of saved user profile

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_LoadUserProfile(unsigned int userID, const char* szInputFilename);

	//!  Saves a user profile for specified user to disk
	/*!		
		\param userID - a valid user ID
		\param szOutputFilename - platform-dependent filesystem path for output file

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_SaveUserProfile(unsigned int userID, const char* szOutputFilename);

	//! Set threshold for Expressiv algorithms
	/*!
		\param userId - user ID
		\param algoName - Expressiv algorithm type
		\param thresholdName - Expressiv threshold type
		\param value - threshold value (min: 0 max: 1000)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_ExpressivAlgo_t, EE_ExpressivThreshold_t
	*/
	EDK_API int
		EE_ExpressivSetThreshold(unsigned int userId, EE_ExpressivAlgo_t algoName, EE_ExpressivThreshold_t thresholdName, int value);


	//! Get threshold from Expressiv algorithms
	/*!
		\param userId - user ID
		\param algoName - Expressiv algorithm type
		\param thresholdName - Expressiv threshold type
		\param pValueOut - receives threshold value

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_ExpressivAlgo_t, EE_ExpressivThreshold_t
	*/
	EDK_API int
		EE_ExpressivGetThreshold(unsigned int userId, EE_ExpressivAlgo_t algoName, EE_ExpressivThreshold_t thresholdName, int *pValueOut);


	//! Set the current facial expression for Expressiv training
	/*!
		Blocking call

		\param userId - user ID
		\param action - which facial expression would like to be trained

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa edkErrorCode.h, EE_ExpressivAlgo_t
	*/
	EDK_API int 
		EE_ExpressivSetTrainingAction(unsigned int userId, EE_ExpressivAlgo_t action);


	//! Set the control flag for Expressiv training
	/*!
		Blocking call

		\param userId - user ID
		\param control - pre-defined control command

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa edkErrorCode.h, EE_ExpressivTrainingControl_t
	*/
	EDK_API int 
		EE_ExpressivSetTrainingControl(unsigned int userId, EE_ExpressivTrainingControl_t control);


	//! Gets the facial expression currently selected for Expressiv training
	/*!
		Blocking call

		\param userId - user ID
		\param pActionOut - receives facial expression currently selected for training

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa EDK_ERROR_CODE, EE_ExpressivAlgo_t
	*/
	EDK_API int 
		EE_ExpressivGetTrainingAction(unsigned int userId, EE_ExpressivAlgo_t* pActionOut);

	
	//! Return the duration of a Expressiv training session
	/*!
		\param userId - user ID
		\param pTrainingTimeOut - receive the training time in ms

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_ExpressivGetTrainingTime(unsigned int userId, unsigned int* pTrainingTimeOut);


	//! Gets a list of the actions that have been trained by the user
	/*!
		Blocking call

		\param userId - user ID
		\param pTrainedActionsOut - receives a bit vector composed of EE_ExpressivAlgo_t contants

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa EDK_ERROR_CODE, EE_ExpressivAlgo_t
	*/
    EDK_API int
        EE_ExpressivGetTrainedSignatureActions(unsigned int userId, unsigned long* pTrainedActionsOut);


	//! Gets a flag indicating if the user has trained sufficient actions to activate a trained signature
	/*!
        *pfAvailableOut will be set to 1 if the user has trained EXP_NEUTRAL and at least
        one other Expressiv action.  Otherwise, *pfAvailableOut == 0.

		Blocking call

		\param userId - user ID
		\param pfAvailableOut - receives an int that is non-zero if a trained signature can be activated

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa EDK_ERROR_CODE
	*/
    EDK_API int
		EE_ExpressivGetTrainedSignatureAvailable(unsigned int userId, int* pfAvailableOut);

	//! Configures the Expressiv suite to use either the built-in, universal signature or a personal, trained signature
	/*!
        Note: Expressiv defaults to use its universal signature.  This function will fail if EE_ExpressivGetTrainedSignatureAvailable returns false.

		Blocking call

		\param userId - user ID
		\param sigType - signature type to use

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa EDK_ERROR_CODE, EE_ExpressivSignature_t
	*/
	EDK_API int 
		EE_ExpressivSetSignatureType(unsigned int userId, EE_ExpressivSignature_t sigType);

	//! Indicates whether the Expressiv suite is currently using either the built-in, universal signature or a trained signature
	/*!
		Blocking call

		\param userId - user ID
		\param pSigTypeOut - receives the signature type currently in use

		\return EDK_ERROR_CODE - current status of EmoEngine. If the query is successful, EDK_ERROR_CODE = OK.

		\sa EDK_ERROR_CODE, EE_ExpressivSignature_t
	*/
	EDK_API int 
		EE_ExpressivGetSignatureType(unsigned int userId, EE_ExpressivSignature_t* pSigTypeOut);


//DEPLOYMENT::STABLE_RELEASE::REMOVE_START

	//@@ These APIs have been obsoleted
	//@@ Use EE_CognitivSetActiveActions and EE_CognitivGetActiveActions instead

	//! Set the current Cognitiv level and corresponding action types
	/*!
		\param userId - user ID
		\param level - current level (min: 1, max: 4)
		\param level1Action - action type in level 1
		\param level2Action - action type in level 2
		\param level3Action - action type in level 3
		\param level4Action - action type in level 4

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivLevel_t, EE_CognitivAction_t
	*/
	EDK_API int 
		EE_CognitivSetCurrentLevel(unsigned int userId, EE_CognitivLevel_t level,
									EE_CognitivAction_t level1Action, EE_CognitivAction_t level2Action,
									EE_CognitivAction_t level3Action, EE_CognitivAction_t level4Action);

	
	//! Query the current Cognitiv level and corresponding action types
	/*!
		\param userId - user ID
		\param pLevelOut - current level (min: 1, max: 4)
		\param pLevel1ActionOut - action type in level 1
		\param pLevel2ActionOut - action type in level 2
		\param pLevel3ActionOut - action type in level 3
		\param pLevel4ActionOut - action type in level 4

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivLevel_t, EE_CognitivAction_t
	*/
	EDK_API int 
		EE_CognitivGetCurrentLevel(unsigned int userId, EE_CognitivLevel_t *pLevelOut,
									EE_CognitivAction_t* pLevel1ActionOut, EE_CognitivAction_t* pLevel2ActionOut,
									EE_CognitivAction_t* pLevel3ActionOut, EE_CognitivAction_t* pLevel4ActionOut);

//DEPLOYMENT::STABLE_RELEASE::REMOVE_END


	//! Set the current Cognitiv active action types
	/*!
		\param userId - user ID
		\param activeActions - a bit vector composed of EE_CognitivAction_t contants

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
	EDK_API int
		EE_CognitivSetActiveActions(unsigned int userId, unsigned long activeActions);

	
	//! Get the current Cognitiv active action types
	/*!
		\param userId - user ID
		\param pActiveActionsOut - receive a bit vector composed of EE_CognitivAction_t contants

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
	EDK_API int
		EE_CognitivGetActiveActions(unsigned int userId, unsigned long* pActiveActionsOut);

	
	//! Return the duration of a Cognitiv training session
	/*!
		\param userId - user ID
		\param pTrainingTimeOut - receive the training time in ms

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_CognitivGetTrainingTime(unsigned int userId, unsigned int* pTrainingTimeOut);

	
	//! Set the training control flag for Cognitiv training
	/*!
		\param userId - user ID
		\param control - pre-defined Cognitiv training control

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivTrainingControl_t
	*/
	EDK_API int 
		EE_CognitivSetTrainingControl(unsigned int userId, EE_CognitivTrainingControl_t control);

	
	//! Set the type of Cognitiv action to be trained
	/*!
		\param userId - user ID
		\param action - which action would like to be trained

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
	EDK_API int 
		EE_CognitivSetTrainingAction(unsigned int userId, EE_CognitivAction_t action);


	//! Get the type of Cognitiv action currently selected for training
	/*!
		\param userId - user ID
		\param pActionOut - action that is currently selected for training

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
	EDK_API int 
		EE_CognitivGetTrainingAction(unsigned int userId, EE_CognitivAction_t* pActionOut);


	//! Gets a list of the Cognitiv actions that have been trained by the user
	/*!
		Blocking call

		\param userId - user ID
		\param pTrainedActionsOut - receives a bit vector composed of EE_CognitivAction_t contants

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
    EDK_API int
        EE_CognitivGetTrainedSignatureActions(unsigned int userId, unsigned long* pTrainedActionsOut);
	
	
	//! Gets the current overall skill rating of the user in Cognitiv
	/*!
		Blocking call

		\param userId - user ID
		\param pOverallSkillRatingOut - receives the overall skill rating [from 0.0 to 1.0]

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
    EDK_API int
        EE_CognitivGetOverallSkillRating(unsigned int userId, float* pOverallSkillRatingOut);


	//! Gets the current skill rating for particular Cognitiv actions of the user
	/*!
		Blocking call

		\param userId - user ID
		\param action - a particular action of EE_CognitivAction_t contant
		\param pActionSkillRatingOut - receives the action skill rating [from 0.0 to 1.0]

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h, EE_CognitivAction_t
	*/
    EDK_API int
        EE_CognitivGetActionSkillRating(unsigned int userId, EE_CognitivAction_t action, float* pActionSkillRatingOut);

	
	//! Set the overall sensitivity for all Cognitiv actions
	/*!
		\param userId - user ID
		\param level - sensitivity level of all actions (lowest: 1, highest: 7)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_CognitivSetActivationLevel(unsigned int userId, int level);

	
	//! Set the sensitivity of Cognitiv actions
	/*!
		\param userId - user ID
		\param action1Sensitivity - sensitivity of action 1 (min: 1, max: 10)
		\param action2Sensitivity - sensitivity of action 2 (min: 1, max: 10)
		\param action3Sensitivity - sensitivity of action 3 (min: 1, max: 10)
		\param action4Sensitivity - sensitivity of action 4 (min: 1, max: 10)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_CognitivSetActionSensitivity(unsigned int userId,
										int action1Sensitivity, int action2Sensitivity,
										int action3Sensitivity, int action4Sensitivity);

	
	//! Get the overall sensitivity for all Cognitiv actions
	/*!
		\param userId - user ID
		\param pLevelOut - sensitivity level of all actions (min: 1, max: 10)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_CognitivGetActivationLevel(unsigned int userId, int *pLevelOut);

	
	//! Query the sensitivity of Cognitiv actions
	/*!
		\param userId - user ID
		\param pAction1SensitivityOut - sensitivity of action 1
		\param pAction2SensitivityOut - sensitivity of action 2
		\param pAction3SensitivityOut - sensitivity of action 3
		\param pAction4SensitivityOut - sensitivity of action 4

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int 
		EE_CognitivGetActionSensitivity(unsigned int userId,
										int* pAction1SensitivityOut, int* pAction2SensitivityOut,
										int* pAction3SensitivityOut, int* pAction4SensitivityOut);

	
	//! Start the sampling of Neutral state in Cognitiv
	/*!
		\param userId - user ID

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivStartSamplingNeutral(unsigned int userId);

	
	//! Stop the sampling of Neutral state in Cognitiv
	/*!
		\param userId - user ID

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivStopSamplingNeutral(unsigned int userId);

	
	//! Enable or disable signature caching in Cognitiv
	/*!
		\param userId  - user ID
		\param enabled - flag to set status of caching (1: enable, 0: disable)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivSetSignatureCaching(unsigned int userId, unsigned int enabled);


	//! Query the status of signature caching in Cognitiv
	/*!
		\param userId  - user ID
		\param pEnabledOut - flag to get status of caching (1: enable, 0: disable)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivGetSignatureCaching(unsigned int userId, unsigned int* pEnabledOut);


	//! Set the cache size for the signature caching in Cognitiv
	/*!
		\param userId - user ID
		\param size   - number of signatures to be kept in the cache (0: unlimited)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivSetSignatureCacheSize(unsigned int userId, unsigned int size);


	//! Get the current cache size for the signature caching in Cognitiv
	/*!
		\param userId - user ID
		\param pSizeOut - number of signatures to be kept in the cache (0: unlimited)

		\return EDK_ERROR_CODE
				- EDK_ERROR_CODE = EDK_OK if successful

		\sa edkErrorCode.h
	*/
	EDK_API int
		EE_CognitivGetSignatureCacheSize(unsigned int userId, unsigned int* pSizeOut);


	//! Returns a struct containing details about the specified EEG channel's headset 
    /*!
        \param channelId - channel identifier (see EmoStateDll.h)
        \param pDescriptorOut - provides detailed sensor location and other info

        \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful

        \sa EmoStateDll.h, edkErrorCode.h
	*/
	EDK_API int
		EE_HeadsetGetSensorDetails(EE_InputChannels_t channelId, InputSensorDescriptor_t* pDescriptorOut);


	//! Returns the current hardware version of the headset and dongle for a particular user
    /*!
        \param userId - user ID for query
		\param pHwVersionOut - hardware version for the user headset/dongle pair. hiword is headset version, loword is dongle version.

        \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful

        \sa EmoStateDll.h, edkErrorCode.h
	*/
	EDK_API int
		EE_HardwareGetVersion(unsigned int userId, unsigned long* pHwVersionOut);

	//! Returns the current version of the Emotiv SDK software
    /*!
		\param pszVersionOut - SDK software version in X.X.X.X format. Note: current beta releases have a major version of 0.
		\param nVersionChars - Length of char buffer pointed to by pszVersion argument.
		\param pBuildNumOut  - Build number.  Unique for each release.

        \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful

        \sa edkErrorCode.h
	*/
	EDK_API int
		EE_SoftwareGetVersion(char* pszVersionOut, unsigned int nVersionChars, unsigned long* pBuildNumOut);

	//! Returns the delta of the movement of the gyro since the previous call for a particular user
	/*!
		\param userId - user ID for query
		\param pXOut  - horizontal displacement
		\param pYOut  - vertical displacment

		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful

        \sa EmoStateDll.h
        \sa edkErrorCode.h
	*/
	EDK_API int
		EE_HeadsetGetGyroDelta(unsigned int userId, int* pXOut, int* pYOut);

	//! Re-zero the gyro for a particular user
	/*!
		\param userId - user ID for query

		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful

        \sa EmoStateDll.h
        \sa edkErrorCode.h
	*/
	EDK_API int
		EE_HeadsetGyroRezero(unsigned int userId);

	//! Returns a handle to memory that can hold an optimization paramaeter which is used to configure the behaviour of optimization
	/*!
		\return OptimizationParamHandle
	*/
	EDK_API OptimizationParamHandle
		EE_OptimizationParamCreate();

	//! Frees memory referenced by an optimization parameter handle
	/*!
		\param hParam - a handle returned by EE_OptimizationParamCreate()
	*/
	EDK_API void
		EE_OptimizationParamFree(OptimizationParamHandle hParam);

	//! Enable optimization. EmoEngine will try to optimize its performance according to the information passed in with optimization parameter. EmoEngine guarantees the correctness of the results of vital algorithms. For algorithms that are not vital, results are undefined.
	/*!
		\param hParam - a handle returned by EE_OptimizationParamCreate()
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationEnable(OptimizationParamHandle hParam);

	//! Determine whether optimization is on
	/*!
		\param pEnabledOut - receives information about whether optimization is on
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationIsEnabled(bool* pEnabledOut);

	//! Disable optimization
	/*!
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationDisable();

	//! Get optimization parameter.  If optimization is not enabled (this can be checked with EE_OptimmizationIsEnabled) then the results attached to the hParam parameter are undefined.
	/*!
		\param hParam - a handle returned by EE_OptimizationParamCreate()
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationGetParam(OptimizationParamHandle hParam);

	//! Get a list of vital algorithms of specific suite from optimization parameter
	/*!
		\param hParam - a handle returned by EE_OptimizationParamCreate()
		\param suite - suite that you are interested in
		\param pVitalAlgorithmBitVectorOut - receives a list of vital algorithm composed of EE_ExpressivAlgo_t, EE_AffectivAlgo_t or EE_CognitivAction_t depending on the suite parameter
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationGetVitalAlgorithm(OptimizationParamHandle hParam, EE_EmotivSuite_t suite, unsigned int* pVitalAlgorithmBitVectorOut);

	//! Set a list of vital algorithms of specific suite to optimization parameter
	/*!
		\param hParam - a handle returned by EE_OptimizationParamCreate()
		\param suite - suite that you are interested in
		\param vitalAlgorithmBitVector - a list of vital algorithm composed of EE_ExpressivAlgo_t, EE_AffectivAlgo_t or EE_CognitivAction_t depended on the suite parameter passed in
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
	EDK_API int
		EE_OptimizationSetVitalAlgorithm(OptimizationParamHandle hParam, EE_EmotivSuite_t suite, unsigned int vitalAlgorithmBitVector);

	//! Resets all settings and user-specific profile data for the specified detection suite
	/*!
		\param userId - user ID
		\param suite - detection suite (Expressiv, Affectiv, or Cognitiv)
		\param detectionBitVector - identifies specific detections.  Set to zero for all detections.
		\return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
	*/
    EDK_API int EE_ResetDetection(unsigned int userId, EE_EmotivSuite_t suite, unsigned int detectionBitVector);

//DEPLOYMENT::NON_PREMIUM_RELEASE::REMOVE_START
        //! Returns a handle to memory that can hold data. This handle can be reused by the caller to retrieve subsequent data.
        /*!
                \return DataHandle
        */
        EDK_API DataHandle EE_DataCreate();

        //! Frees memory referenced by a data handle.
        /*!
                \param hData - a handle returned by EE_DataCreate()
        */
        EDK_API void EE_DataFree(DataHandle hData);

        //! Updates the content of the data handle to point to new data since the last call
        /*!
                \param userId - user ID
                \param hData - a handle returned by EE_DataCreate()
                \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
        */
        EDK_API int EE_DataUpdateHandle(unsigned int userId, DataHandle hData);

        //! Extracts data from the data handle
        /*!
                \param hData - a handle returned by EE_DataCreate()
                \param channel - channel that you are interested in
                \param buffer - pre-allocated buffer
                \param bufferSizeInSample - size of the pre-allocated buffer
                \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
        */
        EDK_API int EE_DataGet(DataHandle hData, EE_DataChannel_t channel, double buffer[], unsigned int bufferSizeInSample);

        //! Returns number of sample of data stored in the data handle
        /*!
                \param hData - a handle returned by EE_DataCreate()
                \param nSampleOut - receives the number of sample of data stored in teh data handle
                \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
        */
        EDK_API int EE_DataGetNumberOfSample(DataHandle hData, unsigned int* nSampleOut);

        //! Sets the size of the data buffer. The size of the buffer affects how frequent EE_DataUpdateHandle() needs to be called to prevent data loss.
        /*!
                \param bufferSizeInSec - buffer size in second
                \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
        */
        EDK_API int EE_DataSetBufferSizeInSec(float bufferSizeInSec);

        //! Returns the size of the data buffer
        /*!
                \param pBufferSizeInSecOut - receives the size of the data buffer
                \return EDK_ERROR_CODE
                - EDK_ERROR_CODE = EDK_OK if successful
        */
        EDK_API int EE_DataGetBufferSizeInSec(float* pBufferSizeInSecOut);

        //! Controls acquisition of data from EmoEngine (which is off by default).
        /*!
            \param userId - user ID
                \param enable - If true, then enables data acquisition
                              - If false, then disables data acquisition
                \return EDK_ERROR_CODE
                        - EDK_ERROR_CODE = EDK_OK if the command succeeded
        */
        EDK_API int EE_DataAcquisitionEnable(unsigned int userId, bool enable);

        //! Returns whether data acquisition is enabled
        /*!
                \param userId - user ID
                \param pEnableOut - receives whether data acquisition is enabled
                \return EDK_ERROR_CODE
                        - EDK_ERROR_CODE = EDK_OK if the command succeeded
        */
        EDK_API int EE_DataAcquisitionIsEnabled(unsigned int userId, bool* pEnableOut);

        //! Sets sychronization signal
        /*!
                \param userId - user ID
                \param signal - value of the sychronization signal
                \return EDK_ERROR_CODE
                        - EDK_ERROR_CODE = EDK_OK if the command succeeded
        */
        EDK_API int EE_DataSetSychronizationSignal(unsigned int userId, int signal);

        //! Sets marker
        /*!
                \param userId - user ID
                \param marker - value of the marker
                \return EDK_ERROR_CODE
                        - EDK_ERROR_CODE = EDK_OK if the command succeeded
        */
        EDK_API int EE_DataSetMarker(unsigned int userId, int marker);
	
        //! Gets sampling rate
        /*!
                \param userId - user ID
                \param samplingRateOut - receives the sampling rate
                \return EDK_ERROR_CODE
                        - EDK_ERROR_CODE = EDK_OK if the command succeeded
        */
        EDK_API int EE_DataGetSamplingRate(unsigned int userId, unsigned int* samplingRateOut);
//DEPLOYMENT::NON_PREMIUM_RELEASE::REMOVE_END
};

#endif
