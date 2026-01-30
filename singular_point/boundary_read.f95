 SUBROUTINE SOIL_READ
  USE BASIC 
  USE MATRIX
  IMPLICIT NONE 
  INTEGER :: I,J,K,N,SOILT
  INTEGER :: N1,N2,N3,N4
  REAL (KIND = 8) :: ALPHA,M2,A,B,THETAS,THETAR,HDRY,HWET,SS
  REAL (KIND = 8) :: K_S(2,2)
  CHARACTER (LEN = 256) :: PATH
  
  call getcwd(PATH)
  PATH = trim(PATH)//'/'
  DO N = 1,SOIL_TOTAL
  	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'soil'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'soil'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'soil'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'soil'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	IF (UIT.EQ.2) THEN 
 		READ (3,*) SS,ALPHA,M2,A,B,THETAS,THETAR
 		READ (3,*) SOILT
 		DO I = 1,2
 			READ (3,*) (K_S(I,J),J=1,2)
 		END DO
 		DO K = 1,CELL_TOT
 			IF (CELL_COE(K)%SOIL_UNIT.EQ.N) THEN 
 				CELL_COE(K)%SS = SS
 				CELL_COE(K)%ALPHA = ALPHA
 				CELL_COE(K)%M2 = M2
 				CELL_COE(K)%A = A; CELL_COE(K)%B = B
 				CELL_COE(K)%THETAS = THETAS
 				CELL_COE(K)%THETAR = THETAR
 				DO I = 1,2
 					DO J = 1,2
 						CELL_COE(K)%K_S(I,J) = K_S(I,J)
 					END DO
 				END DO
	 			CELL_COE(K)%SOIL_T = SOILT
 			END IF
 		END DO
 	ELSE   ! VG and Gardner models' coefficients
 		READ (3,*) SS,ALPHA,M2,THETAS,THETAR,HDRY,HWET
 		READ (3,*) SOILT
 		DO I = 1,2
 			READ (3,*) (K_S(I,J),J=1,2)
 		END DO
 		DO K = 1,CELL_TOT
 			IF (CELL_COE(K)%SOIL_UNIT.EQ.N) THEN 
 				CELL_COE(K)%SS = SS; CELL_COE(K)%ALPHA = ALPHA
 				CELL_COE(K)%M2 = M2
 				CELL_COE(K)%THETAS = THETAS
 				CELL_COE(K)%THETAR = THETAR
 				CELL_COE(K)%HDRY = HDRY; CELL_COE(K)%HWET = HWET
 				DO I = 1,2
 					DO J = 1,2
 						CELL_COE(K)%K_S(I,J) = K_S(I,J)
 					END DO
 				END DO
 				CELL_COE(K)%SOIL_T = SOILT
 			END IF
 		END DO
 	END IF
 	CLOSE (3)
  END DO

 END SUBROUTINE SOIL_READ
  
 SUBROUTINE HEAD_READ 
 USE BASIC
 IMPLICIT NONE
 INTEGER :: I,J,K,N,NK
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 INTEGER :: H_DEG
 REAL (KIND = 8) :: H1,H2,H_TIME1,H_TIME2,H_SET
 REAL (KIND = 8), ALLOCATABLE :: H_SET1(:)
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,HEAD_NUM      ! At this time, we only offer you less than 10000 boundary conditions. I think you might need to 
                        ! see a doctor (for mental problems) if you set up more than 10000 boundary conditions
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'head'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.10) THEN  ! PRESSURE HEAD 
 		READ (3,*) 
 	END IF
 	IF (ITYPE.EQ.0) THEN ! Constant
 		READ (3,*) H_DEG
 		IF (H_DEG.EQ.(-1)) THEN  ! Natural log distribution
 			ALLOCATE (H_SET1(3))
 			H_INPUT(N)%DEG = H_DEG
 			IF (NCYC.EQ.0) THEN 
 				ALLOCATE (H_INPUT(N)%H_COE(3))
 			END IF
 			READ (3,*) H_SET1(1),H_SET1(2),H_SET1(3)  ! HR,A,ALPHA
 			DO I = 1,3
 				H_INPUT(N)%H_COE(I) = H_SET1(I)
 			END DO
 			DEALLOCATE (H_SET1)
 		ELSE 
 			NK = 0
 			DO I = 0,H_DEG
 				DO J = 0,H_DEG-I
 					NK = NK+1
 				END DO
 			END DO
 			ALLOCATE (H_SET1(NK))
 			H_INPUT(N)%DEG = H_DEG
 			IF (NCYC.EQ.0) THEN 
 				ALLOCATE (H_INPUT(N)%H_COE(NK))
 			END IF
 			READ (3,*) (H_SET1(I),I=1,NK)
 			DO I = 1,NK
 				H_INPUT(N)%H_COE(I) = H_SET1(I)
 			END DO
 			DEALLOCATE (H_SET1)
 		END IF
 	ELSE ! Time-varying
 		READ (3,*) H_TIME1,H1
 		IF (NCYC.EQ.0) THEN 
 			H_SET = H1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.H_TIME1)
 			READ (3,*) H_TIME2,H2
 			IF (H_TIME2.GE.T) THEN 
 				H_SET = (H2*(T-H_TIME1)+H1*(H_TIME2-T))/(H_TIME2 - H_TIME1)
 			END IF
 			H_TIME1 = H_TIME2
 			H1 = H2
 		END DO
 		H_INPUT(N)%DEG = 0
 		ALLOCATE (H_INPUT(N)%H_COE(1))
 		H_INPUT(N)%H_COE(1) = H_SET
 	END IF
 	CLOSE (3)
 END DO
 
 
 END SUBROUTINE HEAD_READ
 
 SUBROUTINE FLUX_READ 
 USE BASIC
 IMPLICIT NONE
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: F_TIME1,F_TIME2
 REAL (KIND = 8) :: FX1,FX2,FY1,FY2,FZ1,FZ2
 REAL (KIND = 8) :: FX_SET,FY_SET,FZ_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,FLUX_NUM      ! At this time, we only offer you less than 10000 boundary conditions. I think you might need to 
                        ! see a doctor (for mental problems) if you set up more than 10000 boundary conditions
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'flux'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.0) THEN ! Constant
 		READ (3,*) FX_SET,FY_SET
 		QBX(N) = FX_SET; QBY(N) = FY_SET
 	ELSE  ! Time-varying
 		READ (3,*) F_TIME1,FX1,FY1
 		IF (NCYC.EQ.0) THEN 
 			QBX(N) = FX1; QBY(N) = FY1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.F_TIME1)
 			READ (3,*) F_TIME2,FX2,FY2
 			IF (F_TIME2.GE.T) THEN 
 				FX_SET = (FX2*(T-F_TIME1)+FX1*(F_TIME2-T))/(F_TIME2 - F_TIME1)
 				FY_SET = (FY2*(T-F_TIME1)+FY1*(F_TIME2-T))/(F_TIME2 - F_TIME1)
 			END IF
 			F_TIME1 = F_TIME2
 			FX1 = FX2; FY1 = FY2
 		END DO
 		QBX(N) = FX_SET; QBY(N) = FY_SET
 	END IF
 	CLOSE (3)
 END DO
 
 END SUBROUTINE FLUX_READ
 
 SUBROUTINE RIVER_READ 
 USE BASIC
 IMPLICIT NONE
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: H1,H2,H_TIME1,H_TIME2
 REAL (KIND = 8) :: B1,B2,THICK1,THICK2,AREA1,AREA2
 REAL (KIND = 8) :: H_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,RIVER_NUM      ! At this time, we only offer you less than 10000 boundary conditions. I think you might need to 
                        ! see a doctor (for mental problems) if you set up more than 10000 boundary conditions
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'river'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'river'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'river'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'river'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.0) THEN ! Constant
 		READ (3,*) RIVER_HEAD(N),RIVER_BED(N),RIVER_THICKNESS(N) !,RIVER_AREA(N)
 	ELSE  ! Time-varying
 		READ (3,*) H_TIME1,H1,B1,THICK1 !,AREA1
 		IF (NCYC.EQ.0) THEN 
 			RIVER_HEAD(N) = H1; RIVER_BED(N) = B1; RIVER_THICKNESS(N) = THICK1
 		!	RIVER_AREA(N) = AREA1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.H_TIME1)
 			READ (3,*) H_TIME2,H2,B2,THICK2
 			IF (H_TIME2.GE.T) THEN 
 				RIVER_HEAD(N) = (H2*(T-H_TIME1)+H1*(H_TIME2-T))/(H_TIME2 - H_TIME1)
 				RIVER_THICKNESS(N) = (THICK2*(T-H_TIME1)+THICK1*(H_TIME2-T))/(H_TIME2 - H_TIME1)
 				RIVER_BED(N) = (LOG(B2)*(T-H_TIME1)+LOG(B1)*(H_TIME2-T))/(H_TIME2-H_TIME1)
 				RIVER_BED(N) = EXP(RIVER_BED(N))
 			!	RIVER_AREA(N) = (AREA2*(T-H_TIME1)+AREA1*(H_TIME2-T))/(H_TIME2 - H_TIME1)
 			END IF
 			H_TIME1 = H_TIME2
 			H1 = H2; B1 = B2; THICK1 = THICK2 ! AREA1 = AREA2
 		END DO
 	END IF
 	CLOSE (3)
 END DO
  
 END SUBROUTINE RIVER_READ
 
 SUBROUTINE WELL_READ 
 USE BASIC
 IMPLICIT NONE
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: WELL1,WELL2,WELL_TIME1,WELL_TIME2
 REAL (KIND = 8) :: WELL_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,WELL_NUM      ! At this time, we only offer you less than 10000 boundary conditions. I think you might need to 
                        ! see a doctor (for mental problems) if you set up more than 10000 boundary conditions
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'well'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'well'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'well'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'well'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	read (3,*) WELL_X(N),WELL_Y(N),WELL_COE(N)%LOCAL_NUM
 !	read (3,*) WELL_LAYER(N)

 	IF (ITYPE.EQ.0) THEN ! Constant
 		READ (3,*) WELL_SET
 		QWELL(N) = WELL_SET
 	ELSE  ! Time-varying
 		READ (3,*) WELL_TIME1,WELL1
 		IF (NCYC.EQ.0) THEN 
 			WELL_SET = WELL1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.WELL_TIME1)
 			READ (3,*) WELL_TIME2,WELL2
 			IF (WELL_TIME2.GE.T) THEN 
 				WELL_SET = (WELL2*(T-WELL_TIME1)+WELL1*(WELL_TIME2-T))/(WELL_TIME2 - WELL_TIME1)
 			END IF
 			WELL_TIME1 = WELL_TIME2
 			WELL1 = WELL2
 		END DO
 		QWELL(N) = WELL_SET
 	END IF
 	CLOSE (3)
 END DO
 
 
 END SUBROUTINE WELL_READ
 
SUBROUTINE RAIN_READ 
 USE BASIC
 IMPLICIT NONE
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: F_TIME1,F_TIME2
 REAL (KIND = 8) :: FX1,FX2
 REAL (KIND = 8) :: FX_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,RAIN_NUM      ! At this time, we only offer you less than 10000 boundary conditions. I think you might need to 
                        ! see a doctor (for mental problems) if you set up more than 10000 boundary conditions
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'rain'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N - N1*10
 			NAME3 = trim(PATH)//'rain'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N - N1*100.0)/10.0)
 				N3 = N - N1*100 - N2*10
 				NAME3 = trim(PATH)//'rain'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N - N1*1000.0)/100.0)
 				N3 = FLOOR((N - N1*1000.0 - N2*100.0)/10.0)
 				N4 = N - N1*1000 - N2*100 - N3*10
 				NAME3 = trim(PATH)//'rain'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3, FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.0) THEN ! Constant
 		READ (3,*) FX_SET
 		Q_RAIN(N) = FX_SET
 	ELSE  ! Time-varying
 		READ (3,*) F_TIME1,FX1
 		IF (NCYC.EQ.0) THEN 
 			Q_RAIN(N) = FX1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.F_TIME1)
 			READ (3,*) F_TIME2,FX2
 			IF (F_TIME2.GE.T) THEN 
 				FX_SET = (FX2*(T-F_TIME1)+FX1*(F_TIME2-T))/(F_TIME2 - F_TIME1)
 			END IF
 			F_TIME1 = F_TIME2
 			FX1 = FX2
 		END DO
 		Q_RAIN(N) = FX_SET
 	END IF
 	CLOSE (3)
 END DO
 
 END SUBROUTINE RAIN_READ
 
 
 
