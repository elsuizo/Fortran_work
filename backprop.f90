
PROGRAM Neural_Simulator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Neural Network in Fortran90
!! Multilayer Perceptron trained with
!! the backpropagation learning algorithm
!! coded in Fortran90 by Phil Brierley
!! www.philbrierley.com
!! this code may be used and modified at will
!! compiled using Compaq Visual Fortran
!!-------------------------------------------------------
!! This code reads in data from a csv text file
!! For the neural network training process follow
!! the code in the subroutine 'sTrain_Net'
!! most of the other code is for the data handling
!!-------------------------------------------------------
!! modifications recommended:
!!
!! 1)
!! The data is split into a train,test & validation set.
!! The final model weights should be the best on the test
!! set.
!!
!! 2)
!! The reported errors are based on the normalised data
!! values. These could be scaled up to give actual errors
!!
!!-------------------------------------------------------
!! prefix logic
!! a = array
!! i = integer 
!! r = real
!! l = logical
!! c = character
!! f = function
!! s = subroutine
!! g = global variable
!!-------------------------------------------------------
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!we must declare all our variables
IMPLICIT NONE

!-------------------------------------
!    declarations
!-------------------------------------

!used for checking user input from keyboard
!the value may be compiler dependent so is set here
INTEGER :: giIOERR_OK   =   0

!a handle on the opened source data file
INTEGER, PARAMETER :: giUNIT = 10

  
!-------------------------------------
!    declare arrays
!-------------------------------------

!data arrays
REAL,ALLOCATABLE :: garDataArray(:,:)
REAL,ALLOCATABLE :: garTrainingInputs(:,:),garTrainingOutputs(:)
REAL,ALLOCATABLE :: garTestingInputs(:,:),garTestingOutputs(:)
REAL,ALLOCATABLE :: garValidationInputs(:,:),garValidationOutputs(:)
REAL,ALLOCATABLE :: garInputs_this_pattern(:)

!weight arrays
REAL,ALLOCATABLE :: garWIH(:,:),garWIH_Best(:,:)
REAL,ALLOCATABLE :: garWHO(:),garWHO_Best(:)
!you might also want to save the best test set weights!


!neuron outputs
REAL,ALLOCATABLE :: garHVAL(:)

!dummy arays used in matrix multiplication
REAL,ALLOCATABLE :: garDUMMY1(:,:),garDUMMY2(:,:)

!max and min values of each field
!use these to scale up the reported errors
REAL,ALLOCATABLE :: garMaxInp(:),garMinInp(:)
REAL             :: grMaxOut, grMinOut



!-------------------------------------
!    declare other system variables
!-------------------------------------

!network topolgy numbering (dependent on number of hidden neurons)
INTEGER :: giNHS			!Number Hidden Start
INTEGER :: giNHF			!Number Hidden Finish
INTEGER :: giNOS			!Number Output Start

!general network numbering (independent of number of hidden neurons)
INTEGER :: giINPPB		!INPputs Plus Bias
INTEGER :: giINPUTS
INTEGER :: iOUTPUTS
INTEGER :: giNDU		!Number Data Units

!information about the source data file
INTEGER :: giDATAFIELDS, giFILEROWS, giSKIPLINES

!number of patterns
INTEGER :: giPATS,giTRAINPATS,giTESTPATS,giVALIDPATS

!errors
REAL    :: grRMSE,grRMSEBEST,grRMSETEST


!-------------------------------------
!    user set variables
!-------------------------------------

!number of hidden neurons
INTEGER :: giHIDDDEN
	
!number of epochs to train for
INTEGER :: giEPOCHS

!learning rates
REAL    :: grALR
REAL    :: grBLR

!how often progress is output to screen
INTEGER :: giREDISPLAY



!-------------------------------------
!   main routine
!-------------------------------------
	
	! setup	
	CALL sSetUp

	! neural network training
	CALL sTrain_Net(giIOERR_OK)
	PRINT *,'BYE!'

!-------------------------------------
!    end of main routine
!-------------------------------------




!-------------------------------------
!   all the functions and subroutines
!-------------------------------------

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   THE MAIN NEURAL NETWORK LEARNING SUBS       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE sTrain_Net(iIOERR_OK)

INTEGER, INTENT(IN) ::	iIOERR_OK
INTEGER :: I,J,iPAT_NUM 
REAL    :: rOUTPUT_THIS_PAT, rOUTPRED, rER_THIS_PAT
REAL	:: rRAND


doLooping: DO  

	!get the number of epochs to train for or EXIT
	IF (.NOT. flGet_Number_Of_Epochs(giEPOCHS, iIOERR_OK)) EXIT

	!get the learning rates
	CALL fGet_Learning_Rates(grALR,grBLR,iIOERR_OK)

	!set how often the errors will be output to screen
	giREDISPLAY = fiGet_Screen_Output_Rate(iIOERR_OK)

	CALL sDisplay_Headers


	!do the required number of epochs
	doEpochs: DO J=1,giEPOCHS 

		!an epoch is when every pattern has been seen once
		doPats: DO I=1,giTRAINPATS  
			
			!! select a pattern at random
			CALL RANDOM_NUMBER(rRAND)
			iPAT_NUM=NINT(rRAND*(giTRAINPATS-1))+1

			!! set the data to this pattern
			garInputs_this_pattern(:)=garTrainingInputs(iPAT_NUM,:)
			rOUTPUT_THIS_PAT=garTrainingOutputs(iPAT_NUM)

			!! calculate the current error		
			rER_THIS_PAT= frCalc_Err_This_Pat &
			(garInputs_this_pattern,rOUTPUT_THIS_PAT,garWIH,garWHO)

			!! change weight hidden - output
			garWHO=garWHO-(grBLR*garHVAL*rER_THIS_PAT)

			!! change weight input - hidden 
			garDUMMY2(:,1)=garTrainingInputs(iPAT_NUM,:) 
			garDUMMY1(1,:)=rER_THIS_PAT*garWHO*(1-(garHVAL**2.00))
			garWIH=garWIH-(MATMUL(garDUMMY2,garDUMMY1)*grALR)

		ENDDO doPats ! one more epoch done


		!!  evaluate  'fitness' of of the network  after each epoch
		grRMSE = frCalculate_Error &
		(garTrainingInputs,garTrainingOutputs, garWIH, garWHO)

		!keep the new weights if an improvement has been made
		Call sKeep_Best_Weights(J)

		!! print errors to screen !!
		CALL sDisplay_Progress(J)


	ENDDO doEpochs
	
	! print final errors to screen
	CALL sDISPLAY_ERRORS
	
 ENDDO doLooping


END SUBROUTINE sTrain_Net



REAL FUNCTION frCalc_Err_This_Pat(arINPS_TP,rOUTPUT_TP,arWIHL,arWHOL)
! calculate the error on a specific pattern

	REAL,	DIMENSION(:),	INTENT(IN)	:: arINPS_TP
	REAL,	DIMENSION(:,:),	INTENT(IN)	:: arWIHL
	REAL,	DIMENSION(:),	INTENT(IN)	:: arWHOL
	REAL,	INTENT(IN)	:: rOUTPUT_TP

	REAL :: rOUTPREDL

	garHVAL=TANH(MATMUL(TRANSPOSE(arWIHL),arINPS_TP))
	garHVAL(UBOUND(garHVAL,1)) = 1
	rOUTPREDL=SUM(arWHOL*garHVAL)
	frCalc_Err_This_Pat =(rOUTPREDL-rOUTPUT_TP)


END FUNCTION frCalc_Err_This_Pat



REAL FUNCTION frCalculate_Error(arINPS,arOUT,arWIHL, arWHOL)
! calculate the overall error

	REAL, DIMENSION(:,:),	INTENT(IN)	:: arINPS
	REAL, DIMENSION(:),	INTENT(IN)	:: arOUT
	REAL, DIMENSION(:,:),	INTENT(INOUT)	:: arWIHL
	REAL, DIMENSION(:),	INTENT(INOUT)	:: arWHOL	

	REAL, DIMENSION(LBOUND(arINPS,2):UBOUND(arINPS,2))	:: arINPUTS_THIS_PAT
	REAL :: rSQERROR, rOUTPUT_THIS_PAT, rER_THIS_PAT

	INTEGER :: I,iLOWER,iUPPER

	iLOWER = LBOUND(arINPS,1)
	iUPPER = UBOUND(arINPS,1)

	! in this case the fitness function is the squared errors
	rSQERROR=0.0		

	DO I=iLOWER,iUPPER
		rOUTPUT_THIS_PAT = arOUT(I)
		arINPUTS_THIS_PAT(:)= arINPS(I,:)
		rER_THIS_PAT= frCalc_Err_This_Pat &
		(arINPUTS_THIS_PAT,rOUTPUT_THIS_PAT,arWIHL,arWHOL)
		rSQERROR=rSQERROR+(rER_THIS_PAT**2)
	ENDDO

	! root of the mean squared error
	frCalculate_Error=SQRT(rSQERROR/(iUPPER-iLOWER+1))	

END FUNCTION frCalculate_Error



SUBROUTINE sKEEP_BEST_WEIGHTS(iEpch)
! if the overall error has improved then keep the new weights

	INTEGER, INTENT(IN) :: iEpch

	!this will be on the first epoch
	IF (iEpch .EQ. 1) THEN
		grRMSEBEST = grRMSE
	ENDIF

	IF (grRMSE < grRMSEBEST) THEN
		garWIH_Best = garWIH
		garWHO_Best = garWHO
		grRMSEBEST = grRMSE
	ELSE
		garWIH = garWIH_Best
		garWHO = garWHO_Best
	ENDIF

END SUBROUTINE sKEEP_BEST_WEIGHTS





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  NON NEURAL SUBROUTINES AND FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE sSetUp

	CALL Random_Seed
	CALL sOpen_DataFile(giUNIT)
	CALL sScan_File(giUNIT,giDATAFIELDS, giFILEROWS, giSKIPLINES, giIOERR_OK)
	CALL sSet_Data_Constants(giPATS,giINPUTS,iOUTPUTS,giNDU,giINPPB, &
		giFILEROWS,giSKIPLINES,giDATAFIELDS)
	CALL sGet_Set_Sizes(giTRAINPATS,giTESTPATS,giVALIDPATS,giPATS, giIOERR_OK)
	CALL sAllocate_Data_Arrays
	CALL sRead_Data(giUNIT,giPATS,giNDU,garDataArray,giSKIPLINES,giIOERR_OK)
	CALL sCreate_Training_Data
	CALL sScale_Data
	CALL sSet_Weight_Constants(giIOERR_OK)
	CALL sAllocate_Weight_Arrays
	CALL sInitiate_Weights(garWIH,garWHO,garWIH_Best,garWHO_Best)

END SUBROUTINE sSetUp

!------------------
! Display
!------------------
SUBROUTINE sDisplay_Headers
!print to the screen

	IF (giTESTPATS > 0) THEN
		PRINT *,'epochs   TRAIN_error   TEST_error'
	ELSE
		PRINT *,'epochs   TRAIN_error'
	ENDIF

END SUBROUTINE sDisplay_Headers



SUBROUTINE sDisplay_Errors

	PRINT 100, frCalculate_Error &
	(garTrainingInputs,garTrainingOutputs, garWIH, garWHO)

	IF (giTESTPATS > 0) THEN 
		PRINT 110, frCalculate_Error &
		(garTestingInputs,garTestingOutputs,garWIH,garWHO)
	ENDIF

	IF (giVALIDPATS > 0) THEN
		PRINT 120, frCalculate_Error &
		(garValidationInputs,garValidationOutputs,garWIH,garWHO)
	ENDIF

100 FORMAT('TRAIN ERROR =',1X,F10.7)
110 FORMAT('TEST  ERROR =',1X,F10.7)
120 FORMAT('VAL   ERROR =',1X,F10.7)


END SUBROUTINE sDisplay_Errors



SUBROUTINE sDisplay_Progress(iEpch)

	INTEGER, INTENT(IN) :: iEpch

	IF ( (MODULO(iEpch,giREDISPLAY)==0) .OR. (iEpch==giEPOCHS) .OR. (iEpch==1) ) THEN
		IF (giTESTPATS > 0) THEN
			grRMSETEST = frCalculate_Error &
			(garTestingInputs,garTestingOutputs,garWIH, garWHO)
			PRINT 100,iEpch,grRMSEBEST,grRMSETEST
		ELSE
			PRINT 110,iEpch,grRMSEBEST
		ENDIF
	ENDIF

100 FORMAT(I5,4X,F10.7,4X,F10.7)
110 FORMAT(I5,4X,F10.7)

END SUBROUTINE sDisplay_Progress



!--------------------------
! input file handling
!--------------------------

INTEGER FUNCTION fiCount_Fields(iUNITNUMBER)
!count the number of fields by counting the delimiters	

	INTEGER, INTENT(IN)  :: iUNITNUMBER	
	CHARACTER (LEN = 20000):: miFIRSTLINE 

	READ(UNIT=iUNITNUMBER, FMT='(A)') miFIRSTLINE
	REWIND(UNIT=iUNITNUMBER)
	fiCount_Fields= fiCountF(',',TRIM(miFIRSTLINE)) + 1

END FUNCTION fiCount_Fields



INTEGER FUNCTION fiCountF(cLETTER, cSTRING)
! Count the number of occurrences of LETTER in STRING

	CHARACTER (1), INTENT(IN) :: cLETTER
	CHARACTER (*), INTENT(IN) :: cSTRING
	INTEGER :: I

	fiCountF = 0
	DO I = 1, LEN(cSTRING)
		IF (cSTRING(I:I) == cLETTER) fiCountF = fiCountF + 1
	END DO

END FUNCTION fiCountF



INTEGER FUNCTION fiCount_Rows(iUNITNUMBER)
!count the number of rows in the file

	INTEGER, INTENT(IN)  :: iUNITNUMBER	
 	CHARACTER (LEN = 3)  :: cALINE 

	REWIND(UNIT=iUNITNUMBER)

	fiCount_Rows= 0

	DO
	  READ(UNIT=iUNITNUMBER, FMT='(A)',END=100) cALINE

	  IF (TRIM(cALINE) .NE. '') THEN
		fiCount_Rows= fiCount_Rows + 1	
	  ELSE
		EXIT
	  ENDIF
		
	END DO
	
	100 CONTINUE
	  REWIND(UNIT=iUNITNUMBER)	

END FUNCTION fiCount_Rows



SUBROUTINE sOpen_DataFile(iUNITNUMBER)
! prompt user to enter data file name
	
	INTEGER, INTENT(IN) :: iUNITNUMBER
	CHARACTER*20000 cFILENAME

	CLOSE(iUNITNUMBER)

	PRINT *, 'Enter the datafile name'
	PRINT *, 'it must be comma delimited'
	PRINT *, 'with the the prediction variable in the last column.'
	PRINT *, 'If it is in the same folder as this executable just type the name'
	PRINT *, '?'

	DO 
		READ *, cFILENAME
		IF (cFILENAME == 'q' .or. cFILENAME == 'Q') STOP

		! open data file
		OPEN(UNIT=iUNITNUMBER, ERR=100, FILE=cFILENAME, STATUS='OLD')

		! File Opened OK so Exit
		EXIT

		! File not found
		100 PRINT *, 'Cannot find file ' // TRIM(cFILENAME) // &
		', re-enter or "Q" to quit'

	END DO

END SUBROUTINE sOpen_DataFile



SUBROUTINE sRead_Data &
(iUNITNUMBER, iPATTERNS, iFIELDS, arDATAARRAY, iHEADERROWS, iIOERR_OK)
! read in data	

	INTEGER, INTENT(IN)  :: iUNITNUMBER	
	INTEGER, INTENT(IN)  :: iPATTERNS
	INTEGER, INTENT(IN)  :: iFIELDS
	INTEGER, INTENT(IN)  :: iHEADERROWS
	INTEGER, INTENT(IN)  :: iIOERR_OK
	REAL, INTENT(OUT)    :: arDATAARRAY(iPATTERNS,iFIELDS)

	CHARACTER (LEN = 1)  :: cHEAD	
	
	INTEGER :: I,K
	INTEGER :: mis

	REWIND(UNIT=iUNITNUMBER)
	
	IF (iHEADERROWS > 0) THEN	
		DO I = 1,iHEADERROWS
			READ(UNIT=iUNITNUMBER, FMT='(A)',END=100) cHEAD
		ENDDO
		100 IF (I < iHEADERROWS) PRINT *,'Too many header rows...' 
	END IF	
	
	PRINT *, ' Reading data...'
	
	DO I=1,iPATTERNS
		READ(UNIT=iUNITNUMBER, FMT=*,iostat=mis, END=200)(arDATAARRAY(I,K),K=1,iFIELDS)
		IF(mis /= iIOERR_OK) then
		  PRINT *, 'Invalid data on line ', i + iHEADERROWS
		  STOP
		ENDIF
	ENDDO
	
	PRINT *, ' Data read OK!'

	200 IF (I < iPATTERNS) PRINT *,'Too many header rows...'

	CLOSE(iUNITNUMBER)

END SUBROUTINE sRead_Data



SUBROUTINE sScan_File(iUNITNUMBER, iDATAFIELDS, iFILEROWS, iSKIPLINES, iIOERR_OK)

	INTEGER, INTENT(IN)   :: iUNITNUMBER	
	INTEGER, INTENT(OUT)  :: iDATAFIELDS
	INTEGER, INTENT(OUT)  :: iFILEROWS
	INTEGER, INTENT(OUT)  :: iSKIPLINES
	INTEGER, INTENT(INOUT)   :: iIOERR_OK	

	PRINT *, 'Scanning file...'

	REWIND(UNIT=iUNITNUMBER)

	iDATAFIELDS = fiCount_Fields(iUNITNUMBER)
	iFILEROWS = fiCount_Rows(iUNITNUMBER)

	PRINT * , 'Fields = ',iDATAFIELDS
	PRINT * , 'Rows   = ',iFileRows

	DO
	  iSKIPLINES = fiGet_Header_Rows(iIOERR_OK)
		IF (iSKIPLINES < iFILEROWS) THEN
			EXIT
		ELSE
			PRINT *, 'You have more header rows than &
			there are rows in the data!'
		ENDIF
	END DO

END SUBROUTINE sScan_File



INTEGER FUNCTION fiGet_Header_Rows(iIOERR_OK)

	INTEGER, INTENT(IN)		:: iIOERR_OK
	INTEGER :: I

	DO
	 WRITE(*,'(A)',advance='no',iostat=I) 'How many header rows in the file?'
	 IF(I /= iIOERR_OK) EXIT
		READ(*,*,iostat=I) fiGet_Header_Rows
		IF(I /= iIOERR_OK) CYCLE
		IF (fiGet_Header_Rows< 0) THEN
		 PRINT *, 'OK - no header row...'
		 fiGet_Header_Rows= 0
		ENDIF
	 EXIT
	ENDDO

END FUNCTION fiGet_Header_Rows


!---------------------------
! getting user input
!---------------------------
INTEGER FUNCTION fiGet_Hidden_Neurons(iIOERR_OK)

	INTEGER, INTENT(IN)		:: iIOERR_OK
	INTEGER :: I

	DO
	 WRITE(*,'(A)',advance='no',iostat=I) 'How many hidden neurons?'
	 IF(I /= iIOERR_OK) EXIT
		READ(*,*,iostat=I) fiGet_Hidden_Neurons
		IF(I /= iIOERR_OK) CYCLE
		IF (fiGet_Hidden_Neurons< 1) THEN
		 PRINT *, 'Hidden neurons set to 1'
		 fiGet_Hidden_Neurons = 1
		ENDIF
	 EXIT
	ENDDO

END FUNCTION fiGet_Hidden_Neurons



SUBROUTINE sGet_Set_Sizes(iTRAIN,iTEST,iVALID,iTOTAL,iIOERR_OK)

	INTEGER, INTENT(IN)	:: iIOERR_OK, iTOTAL
	INTEGER, INTENT(OUT)	:: iTRAIN, iTEST, iVALID

	INTEGER	:: I, iPERCTR, iPERCTE

	iTRAIN = iTOTAL
	iTEST = 0
	iVALID = 0
	iPERCTR = 33
	iPERCTE = 0

	DO
		iPERCTR = 33
		WRITE(*,'(A)',advance='no',iostat=I) &
		'What percent do you want to use for the training sample? '

		IF(I /= iIOERR_OK) EXIT
		READ(*,*,iostat=I) iPERCTR
		IF(I /= iIOERR_OK) CYCLE

		IF (iPERCTR >= 100) THEN
			PRINT *, 'OK - all the data will be used for training...'
			iPERCTR = 100
		ENDIF

		IF (iPERCTR<33) THEN
			PRINT *, 'We will use 33%...'
			iPERCTR = 33
		ENDIF

		iTRAIN = iTOTAL * iPERCTR / 100

		EXIT

	ENDDO


	IF (iPERCTR < 100) THEN

	 DO

		WRITE(*,'(A)',advance='no',iostat=I) &
		'What percent do you want to use for the testing sample? '

		IF(I /= iIOERR_OK) EXIT
		READ(*,*,iostat=I) iPERCTE
		IF(I /= iIOERR_OK) CYCLE

		IF (iPERCTE>100-iPERCTR) CYCLE
		IF (iPERCTE<=0) CYCLE

		iTEST = iTOTAL * iPERCTE / 100
		IF (iTEST<1) iTEST = 1

		EXIT

	 ENDDO

	ENDIF


	iVALID = iTOTAL - iTRAIN - iTEST
	PRINT * , ' '
	PRINT * , ' ' 
	PRINT 100 ,iTOTAL
	PRINT 110 ,iPERCTR,iTRAIN
	PRINT 120 ,iPERCTE,iTEST
	PRINT 130 ,100-iPERCTR-iPERCTE,iVALID
	PRINT * , '' 

100 FORMAT ('Total patterns = ',I10)
110 FORMAT ('Train =',1X,I3,'%',I7,1X,'Patterns')
120 FORMAT ('Test  =',1X,I3,'%',I7,1X,'Patterns')
130 FORMAT ('Valid =',1X,I3,'%',I7,1X,'Patterns')

END SUBROUTINE sGet_Set_Sizes



LOGICAL FUNCTION flGet_Number_Of_Epochs(iEPOCHS, iIOERR_OK)

	INTEGER, INTENT(OUT)	:: iEPOCHS
	INTEGER, INTENT(IN)	:: iIOERR_OK
	INTEGER			:: i

	DO

	 WRITE(*,'(A)',advance='no',iostat=i) &
	 'How many epochs to train for (0 to exit) ?'

	 IF(i /= iIOERR_OK) EXIT
		READ(*,*,iostat=i) iEPOCHS
		IF(i /= iIOERR_OK) CYCLE

	 IF (iEPOCHS <= 0) THEN
		flGet_Number_Of_Epochs = .FALSE.
	 ELSE
		flGet_Number_Of_Epochs = .TRUE.
	 ENDIF	
	
	 EXIT

	ENDDO


END FUNCTION flGet_Number_Of_Epochs



SUBROUTINE fGet_Learning_Rates(ra,rb, iIOERR_OK)

	REAL, INTENT(OUT)	:: ra
	REAL, INTENT(OUT)	:: rb
	INTEGER, INTENT(IN)	:: iIOERR_OK
	INTEGER			:: I

	DO
	 WRITE(*,'(A)',advance='no',iostat=I) 'Learning Rate (>0 - 2) ?'
	 IF(I /= iIOERR_OK) EXIT
	 READ(*,*,iostat=I) ra
	 IF(i /= iIOERR_OK) CYCLE
	 rb = ra/10
	 EXIT
	ENDDO

END SUBROUTINE fGet_Learning_Rates



INTEGER FUNCTION fiGet_Screen_Output_Rate(iIOERR_OK) 

	INTEGER, INTENT(IN)	:: iIOERR_OK
	INTEGER			:: I

	DO
	 WRITE(*,'(A)',advance='no',iostat=I) 'Screen Output Rate (whole number)?'

	 IF(I /= iIOERR_OK) EXIT
		READ(*,*,iostat=I) fiGet_Screen_Output_Rate 
		IF(I /= iIOERR_OK) CYCLE

	 IF (fiGet_Screen_Output_Rate < 1) THEN
		fiGet_Screen_Output_Rate = 1
		PRINT *, 'OK - set it to 1!'
	 ENDIF	

	 EXIT
	ENDDO

END FUNCTION fiGet_Screen_Output_Rate



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DATA ALLOCATION SUBROUTINES AND FUNCTIONS   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE sAllocate_Data_Arrays()

	! set the array dimensions
	ALLOCATE(garDataArray(giPATS,giNDU)) !raw data read from file
	ALLOCATE(garTrainingInputs(giTRAINPATS,giINPPB)) !input patterns
	ALLOCATE(garTrainingOutputs(giTRAINPATS)) !output
	ALLOCATE(garTestingInputs(giTESTPATS,giINPPB)) !input patterns
	ALLOCATE(garTestingOutputs(giTESTPATS)) !output
	ALLOCATE(garValidationInputs(giVALIDPATS,giINPPB)) !input patterns
	ALLOCATE(garValidationOutputs(giVALIDPATS)) !output
	ALLOCATE(garInputs_this_pattern(giINPPB)) !pattern being presented
	ALLOCATE(garMaxInp(giINPPB))
	ALLOCATE(garMinInp(giINPPB))

END SUBROUTINE sAllocate_Data_Arrays



SUBROUTINE sAllocate_Weight_Arrays()

	ALLOCATE(garWIH(giINPPB,giNHS:giNHF))		!input-hidden weights
	ALLOCATE(garWIH_Best(giINPPB,giNHS:giNHF))	!best weights
	ALLOCATE(garWHO(giNHS:giNHF))			!hidden-output weights
	ALLOCATE(garWHO_Best(giNHS:giNHF))		!best weights
	ALLOCATE(garHVAL(giNHS:giNHF))			!hidden neuron outputs
	ALLOCATE(garDUMMY1(1,giNHS:giNHF))		!dummy matrix
	ALLOCATE(garDUMMY2(giINPPB,1))			!dummy matrix

END SUBROUTINE sAllocate_Weight_Arrays



SUBROUTINE sSet_Data_Constants &
	(iNPATS,iINPUTS,iNOUTPUTS,iNDU,iINPPB,iFILEROWS,iSKIPLINES,iDATAFIELDS)

	INTEGER, INTENT (OUT)	:: iNPATS
	INTEGER, INTENT (OUT)	:: iINPUTS
	INTEGER, INTENT (OUT)	:: iNOUTPUTS
	INTEGER, INTENT (OUT)	:: iNDU
	INTEGER, INTENT (OUT)	:: iINPPB
	
	INTEGER, INTENT (IN)	:: iFILEROWS
	INTEGER, INTENT (IN)	:: iSKIPLINES	
	INTEGER, INTENT (IN)	:: iDATAFIELDS 


	iNPATS = iFILEROWS - iSKIPLINES
	iINPUTS = iDATAFIELDS - 1
	iNOUTPUTS=1			! number of outputs (fixed)
	iNDU=iINPUTS+iNOUTPUTS		!Number Data Units
	iINPPB=iINPUTS+1		!INPut Plus Bias


END SUBROUTINE sSet_Data_Constants



SUBROUTINE sSet_Weight_Constants(iIOERR_OK)

	INTEGER, INTENT(IN) :: iIOERR_OK

	giHIDDDEN = fiGet_Hidden_Neurons(iIOERR_OK)  

	! number the neurons
	giHIDDDEN=giHIDDDEN+1			!accounts for bias to output
	giNHS=giINPPB+1				!Number Hidden Start
	giNHF=giINPPB+giHIDDDEN			!Number Hidden Finish
	giNOS=giNHF+1				!Number Output Start

END SUBROUTINE sSet_Weight_Constants



SUBROUTINE sCreate_Training_Data()
! create train, test and validation sets

	INTEGER :: iTRAINCOUNT,iTESTCOUNT,iVALIDCOUNT
	INTEGER :: I
	REAL :: rTRPC,rTEPC
	REAL :: rRAND	
	
	rTRPC = FLOAT(giTRAINPATS) / FLOAT(giPATS)
	rTEPC = rTRPC + (FLOAT(giTESTPATS) / FLOAT(giPATS))
	
	iTRAINCOUNT = 0
	iTESTCOUNT = 0
	iVALIDCOUNT = 0

	PRINT *, ' Allocating data...'

 DO I=1,giPATS
	CALL RANDOM_NUMBER(rRAND)

	IF ((rRAND <= rTRPC) .AND. (iTRAINCOUNT < giTRAINPATS)) THEN
		iTRAINCOUNT = iTRAINCOUNT + 1 
		garTrainingInputs(iTRAINCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garTrainingInputs(iTRAINCOUNT,giINPPB)=1
		garTrainingOutputs(iTRAINCOUNT)=garDataArray(I,giNDU)

	ELSEIF ((rRAND <= rTEPC) .AND. (iTESTCOUNT < giTESTPATS)) THEN
		iTESTCOUNT = iTESTCOUNT + 1
		garTestingInputs(iTESTCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garTestingInputs(iTESTCOUNT,giINPPB)=1
		garTestingOutputs(iTESTCOUNT)=garDataArray(I,giNDU)

	ELSEIF ((rRAND > rTEPC) .AND. (iVALIDCOUNT < giVALIDPATS)) THEN
		iVALIDCOUNT = iVALIDCOUNT + 1
		garValidationInputs(iVALIDCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garValidationInputs(iVALIDCOUNT,giINPPB)=1
		garValidationOutputs(iVALIDCOUNT)=garDataArray(I,giNDU)

	ELSEIF (iTRAINCOUNT < giTRAINPATS) THEN
		iTRAINCOUNT = iTRAINCOUNT + 1 
		garTrainingInputs(iTRAINCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garTrainingInputs(iTRAINCOUNT,giINPPB)=1
		garTrainingOutputs(iTRAINCOUNT)=garDataArray(I,giNDU)

	ELSEIF (iTESTCOUNT < giTESTPATS) THEN
		iTESTCOUNT = iTESTCOUNT + 1
		garTestingInputs(iTESTCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garTestingInputs(iTESTCOUNT,giINPPB)=1
		garTestingOutputs(iTESTCOUNT)=garDataArray(I,giNDU)

	ELSEIF (iVALIDCOUNT < giVALIDPATS) THEN
		iVALIDCOUNT = iVALIDCOUNT + 1
		garValidationInputs(iVALIDCOUNT,1:giINPUTS)=garDataArray(I,1:giINPUTS)
		garValidationInputs(iVALIDCOUNT,giINPPB)=1
		garValidationOutputs(iVALIDCOUNT)=garDataArray(I,giNDU)

	ENDIF

 ENDDO

	!the array 'data' is no longer required
	DEALLOCATE(garDataArray)				
	PRINT *, ' Data allocated OK!'

END SUBROUTINE sCreate_Training_Data




SUBROUTINE sScale_Data()
! scale the data

	INTEGER :: I

	PRINT *, ' Normalising data...'

	! find the max and min values
	garMaxInp(:) = MAXVAL(garTrainingInputs,1)
	garMinInp(:) = MINVAL(garTrainingInputs,1)
	grMaxOut = MAXVAL(garTrainingOutputs)
	grMinOut = MINVAL(garTrainingOutputs)


	! need to check if max = min

	! scale between -1 and 1
	DO i = 1,giTRAINPATS
	garTrainingInputs(i,1:giINPUTS) = &
	((garTrainingInputs(i,1:giINPUTS) - garMinInp(1:giINPUTS)) / &
	(garMaxInp(1:giINPUTS) - garMinInp(1:giINPUTS)) - 0.5) * 2
	ENDDO	

	garTrainingOutputs(:) = &
	((garTrainingOutputs(:) - grMinOut) / (grMaxOut - grMinOut) - 0.5) * 2
	
	IF (giTESTPATS > 0) THEN
	DO i = 1,giTESTPATS
	garTestingInputs(i,1:giINPUTS) = &
	((garTestingInputs(i,1:giINPUTS) - garMinInp(1:giINPUTS)) / &
	(garMaxInp(1:giINPUTS) - garMinInp(1:giINPUTS)) - 0.5) * 2
	ENDDO		
	garTestingOutputs(:) = &
	((garTestingOutputs(:) - grMinOut) / (grMaxOut - grMinOut) - 0.5) * 2
	ENDIF

	IF (giVALIDPATS > 0) THEN
	DO i = 1,giVALIDPATS
	garValidationInputs(i,1:giINPUTS) = &
	((garValidationInputs(i,1:giINPUTS) - garMinInp(1:giINPUTS)) / &
	(garMaxInp(1:giINPUTS) - garMinInp(1:giINPUTS)) - 0.5) * 2
	ENDDO	
	garValidationOutputs(:) = &
	((garValidationOutputs(:) - grMinOut) / (grMaxOut - grMinOut) - 0.5) * 2
	ENDIF

	PRINT *, ' Data normalised OK!'
	PRINT *, ''

END SUBROUTINE sScale_Data



SUBROUTINE sInitiate_Weights(arWIHL,arWHOL,arWIHBESTL,arWHOBESTL)
! generate initial random weights
	
	Real, Dimension (:,:),	Intent (INOUT) :: arWIHL
	Real, Dimension (:,:),	Intent (INOUT) :: arWIHBESTL	
	Real, Dimension (:), 	Intent (INOUT) :: arWHOL
	Real, Dimension (:),	Intent (INOUT) :: arWHOBESTL	

	INTEGER :: J, K
	REAL :: rRAND

	DO K=LBOUND(arWIHL,2),UBOUND(arWIHL,2)
		CALL RANDOM_NUMBER(rRAND)
		arWHOL(K)=((rRAND-0.5)*2)/10
		DO J=LBOUND(arWIHL,1),UBOUND(arWIHL,1)
			CALL RANDOM_NUMBER(rRAND)
			arWIHL(J,K)=((rRAND-0.5)*2)/10
		ENDDO
	ENDDO

	arWIHBESTL=arWIHL !record of best weights so far
	arWHOBESTL=arWHOL !record of best weights so far

END SUBROUTINE sInitiate_Weights


END PROGRAM Neural_Simulator
