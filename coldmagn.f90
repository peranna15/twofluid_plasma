	PROGRAM COLDMAGN

	implicit none

        integer :: JR, NRAD, TNRAD, IREC, FLAGREC, ITREC, ITMREC, JK, JT, KFIN,NLIM,INDB

        logical :: isnan

        double precision, allocatable, dimension(:) :: DENE1, DENE2, VRE1, VRE2, DENI1, DENI2, VRI1, VRI2, Te1,Ti1,DLOGE1,DLOGE2
        double precision, allocatable, dimension(:) :: Te2, Ti2, ERAD, DEBYE, CHGSEP,GRADPE,GRADPI &
        ,DLOGI1,DLOGI2,GRADVE,VEMAX,VIMAX,VTHE1,VTHE2,VTHI1,VTHI2,BINT

	double precision :: R, R0, DEN0, TAU, VOL, EIN0, TFIN, DNICH, TITMREC, DR, W, PREC,SWTION,SWTELE,QTTIO0,QTTEL0,QTTIO,QTTEL

	double precision :: COLLEI, COUL, COLLHH, COLLVEL, COULI, DSEP, ECINT, RINT, EIN, EION, ENEJ, GRADVI,ELIM,DNAMB
	
	double precision :: TIME,TLOOP,SWTENE,BEXT,GRADVTHE,GRADVTHI
	
	OPEN(9,FILE='init_coldpla',STATUS='unknown')
	OPEN(30,FILE='input',STATUS='unknown')
!	OPEN(10,FILE='Profili.DAT',STATUS='unknown')
	OPEN(11,FILE='Timeprofs.DAT',STATUS='unknown')
	write(*,*) 'Densita gas unit 10^19 cm-3'
        read(9,*)
	read(9,*) DEN0
        write(*,*) DEN0
	write(*,*) 'tempo scarica'
        read(9,*)
	read(9,*) TAU	
	write(*,*) TAU	
	write(*,*) 'Volume scarica'
        read(9,*)
	read(9,*) VOL	
	write(*,*) VOL	
	write(*,*) 'Energia condensatori max J'
        read(9,*)
	read(9,*) EIN0	
	write(*,*) EIN0	
	write(*,*) 'passo temporale'
        read(9,*)
	read(9,*) W
        write(*,*) W
	write(*,*) 'numero steps radiali'
        read(9,*)
	read(9,*) TNRAD
	write(*,*) TNRAD
	NRAD=TNRAD
	write(*,*) 'tempo finale'
        read(9,*)
	read(9,*) TFIN
	write(*,*) TFIN
	
	write(*,*) 'raggio scarica iniziale'
        read(9,*)
	read(9,*) R0	
	write(*,*) R0	
	write(*,*) 'Posizione Nichel'
        read(9,*)
	read(9,*) DNICH
	write(*,*) DNICH
	write(*,*) 'Posizione record dati vs time'
        read(9,*)
	read(9,*) PREC
	write(*,*) PREC
	write(*,*) 'time steps record data'
        read(9,*)
	read(9,*) TITMREC
	write(*,*) TITMREC
	ITMREC=TITMREC
	write(*,*) 'dammi B esterno in Gauss'
        read(9,*)
	read(9,*) BEXT
	write(*,*) BEXT
        close(9)
!	--------------------------------
	write(30,*) DEN0
	write(30,*) TAU
	write(30,*) VOL
	write(30,*) EIN0
	write(30,*) W
	write(30,*) TNRAD
	write(30,*) TFIN
	write(30,*) R0
	write(30,*) DNICH
	write(30,*) PREC
	write(30,*) TITMREC
	write(30,*) BEXT
	CLOSE(30)
!	---------------------------

    allocate(DENE1(nrad))
	allocate(DENI1(nrad))
	allocate(CHGSEP(nrad))
	allocate(DEBYE(nrad))
	allocate(VRE1(nrad))
	allocate(VRE2(nrad))
	allocate(VRI1(nrad))
	allocate(VRI2(nrad))
	allocate(Te1(nrad))
	allocate(Ti1(nrad))
	allocate(Te2(nrad))
	allocate(Ti2(nrad))
	allocate(DENE2(nrad))
	allocate(DENI2(nrad))
	allocate(ERAD(nrad))
	allocate(GRADPE(nrad))
	allocate(GRADPI(nrad))
	allocate(DLOGE1(nrad))
	allocate(DLOGE2(nrad))
	allocate(DLOGI1(nrad))
	allocate(DLOGI2(nrad))
	allocate(GRADVE(nrad))
	allocate(VEMAX(nrad))
	allocate(VIMAX(nrad))
	allocate(VTHE1(nrad))
	allocate(VTHE2(nrad))
	allocate(VTHI1(nrad))
	allocate(VTHI2(nrad))
	allocate(BINT(nrad))
!       ------------------------	
	DR=DNICH/NRAD
	IREC=floor(PREC/DR)
!	------------------------
	EION=1.6*DEN0*VOL*13.6
	IF (EIN0.GT.EION) THEN
	EIN=EIN0-EION
		ELSE
	EIN=0
	ENDIF
	write(*,*) EIN,EIN0,EION
	QTTIO0=0
	QTTEL0=0
!	----------------------------------------
!	---initial conditions--------------
	R=0.
	DO JR=1,NRAD
		DENE1(JR)=DEN0*EXP(-((R)/R0)**2.)
	DENE2(JR)=DEN0*EXP(-((R)/R0)**2.)
	DENI1(JR)=DENE1(JR)
	DENI2(JR)=DENE2(JR)
	DLOGE1(JR)=-(R/R0)**2
	DLOGI1(JR)=-(R/R0)**2
	DLOGE2(JR)=-(R/R0)**2
	DLOGI2(JR)=-(R/R0)**2
	CHGSEP(JR)=0
	DEBYE(JR)=DR
	Te1(JR)=.1
	Ti1(JR)=.1
	Te2(JR)=.1
	Ti2(JR)=.1
	VRE1(JR)=5e7*SQRT(Te1(JR))
	VRI1(JR)=0
	VRE2(JR)=5e7*SQRT(Te1(JR))
	VRI2(JR)=0
	VTHE1(JR)=0
	VTHE2(JR)=0
	VTHI1(JR)=0
	VTHI2(JR)=0
	BINT(JR)=0
	ERAD(JR)=0
	VEMAX(JR)=0
	QTTIO0=QTTIO0+DENI1(JR)*(R+DR)*DR
	QTTEL0=QTTEL0+DENE1(JR)*(R+DR)*DR
	R=R+DR
	enddo
! 1	CONTINUE
	
	SWTION=1
	SWTELE=1
	ITREC=1
	TIME=0
	TLOOP=TFIN
	KFIN=floor(TLOOP/W)
!	----B along z and in Gauss--------
!	BEXT=1e3
	
!	----time evolution----------------
	DO JT=1,KFIN
		IF (JT*W.GT.5*TAU) THEN
	SWTENE=0
		ELSE
	SWTENE=1
		ENDIF
!	----------------------------------
!	------------------------------------------
!	---------------------------------------------
!	----calcolo E radiale-------------------
	R=0
	DO JR=1,NRAD
	BINT(JR)=0
	enddo
	
	DO JR=1,NRAD
	INDB=NRAD-JR+1
!	------------------------
	
!	---calcolo campo radiale-------------------
	RINT=0
	ECINT=0
	DO JK=1,JR
!	-----------------------------------
!	
			
	ECINT=ECINT+6e10*CHGSEP(JK)*RINT*MIN(DR,DEBYE(JK))
	RINT=RINT+DR
	enddo
! 4	CONTINUE
	IF (JR.EQ.1) THEN
	ERAD(JR)=0
		ELSE	
	ERAD(JR)=ECINT/R
	
		ENDIF
	BINT(INDB-1)=BINT(INDB) &
	+2*DR* &
	(DENI1(INDB)*VTHI1(INDB)-DENE1(INDB)*VTHE1(INDB))
	
	R=R+DR
	enddo
! 3	CONTINUE
!	----------end calcolo campo ele------------
!	-------loop radiale--------
 672	R=0
	DO JR=1,NRAD

!	------------------------
	ENEJ=EIN*(1-EXP(-(JT-1)*W/TAU))*EXP(-(R/R0)**2)
!	---termalizzazione--el ion--
	IF (DENE1(JR).GT.1e-20.AND.Te1(JR).GT.1) THEN
	COUL=26.-.5*LOG(DENE1(JR)/Te1(JR))
	COLLEI=3.2e10*DENE1(JR)*COUL/(Te1(JR)**1.5)
!	---collisione e-ion scambio momento-----
	COLLVEL=3e13*COUL*DENI1(JR)/Te1(JR)**1.5
		ELSE
	COUL=0
	COLLEI=0
	COLLVEL=0
		ENDIF
!	--collisione ion-ion scambio momento--
	IF (DENI1(JR).GT.1e-20.AND.Ti1(JR).GT.1) THEN
	COULI=26.-.5*LOG(DENI1(JR)/Ti1(JR))
	COLLHH=1.8e12*DENI1(JR)*COULI/Ti1(JR)**1.5
		ELSE
	COULI=0
	COLLHH=0
		ENDIF	
!	-------------------------
!	-----electrons------
	 IF (SWTELE.EQ.1) THEN
!	---------------------------------------

!	--------------------------------------------- 
		IF (JR.GT.1.AND.JR.LT.NRAD) THEN
	GRADPE(JR)=.5*(Te1(JR+1)-Te1(JR-1)+Te1(JR)* &
	(DLOGE1(JR+1)-DLOGE1(JR-1)))/DR
		
	GRADVE(JR)=.5*VRE1(JR)*(VRE1(JR+1)-VRE1(JR-1))
	
	GRADVTHE=.5*VRE1(JR)*(VTHE1(JR+1)-VTHE1(JR-1))
	
	ENDIF

	IF (JR.EQ.NRAD) THEN

	GRADPE(JR)=(Te1(JR)-Te1(JR-1)+Te1(JR)* &
	(DLOGE1(JR)-DLOGE1(JR-1)))/DR
	
	GRADVE(JR)=VRE1(JR)*(VRE1(JR)-VRE1(JR-1))
	
	GRADVTHE=VRE1(JR)*(VTHE1(JR)-VTHE1(JR-1))
	ENDIF
!	----------------------------------------------	
!	-------------------------------------------------
!	------------------------------------------   
    IF (JR.GT.1.AND.JR.LT.NRAD) THEN
    
	DLOGE2(JR)=DLOGE1(JR)-W*VRE1(JR)/R- &
        .5*W/DR*(VRE1(JR+1)-VRE1(JR-1)+ &
        VRE1(JR)*(DLOGE1(JR+1)-DLOGE1(JR-1)))	
    
!    TE2(JR)=TE1(JR)-(.666*W*VRI1(JR)*TE1(JR)/R+ &
!       .333*W/DR*TE1(JR)*(VRI1(JR+1)-VRI1(JR-1))+ &
!       .5*W/DR*VRI1(JR)*(TE1(JR+1)-TE1(JR-1))) &
!       +W*(TI1(JR)-TE1(JR))*COLLEI&
!       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR))&
!       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
!       *SWTENE

	TE2(JR)=TE1(JR) &
       +W*(TI1(JR)-TE1(JR))*COLLEI&
       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR))&
       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
       *SWTENE
        
		ENDIF
!	------------------------------------------		
		IF (JR.EQ.NRAD) THEN

    DLOGE2(JR)=DLOGE1(JR)-W*VRE1(JR)/R- &
        W/DR*(VRE1(JR)-VRE1(JR-1)+ &
        VRE1(JR)*(DLOGE1(JR)-DLOGE1(JR-1)))
        
!     TE2(JR)=TE1(JR)-(.666*W*VRI1(JR)*TE1(JR)/R+ &
!       .666*W/DR*TE1(JR)*(VRI1(JR)-VRI1(JR-1))+ &
!       W/DR*VRI1(JR)*(TE1(JR)-TE1(JR-1))) &
!       +W*(TI1(JR)-TE1(JR))*COLLEI&
!       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR))&
!       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
!       *SWTENE

	TE2(JR)=TE1(JR) &
       +W*(TI1(JR)-TE1(JR))*COLLEI&
       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR))&
       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
       *SWTENE
        		
		ENDIF
!	----------------------------------
		IF (JR.GT.1) THEN       
    VRE2(JR)=(1-W*COLLVEL)*VRE1(JR)-5e17*W*ERAD(JR) &
        -1.6e15*W*GRADPE(JR) &
        -W/DR*GRADVE(JR)-1.6e7*W*VTHE1(JR)*(BEXT+BINT(JR))
        
    VTHE2(JR)=(1-W*COLLVEL)*VTHE1(JR)+1.6e7*W*VRE1(JR)* &
    (BEXT+BINT(JR))-W/DR*GRADVTHE

	VEMAX(JR)=SQRT(ABS( &
	-1e18*ERAD(JR)*MIN(DR,DEBYE(JR)) &
	-3.2e15*GRADPE(JR)*DR))
!	-------------------------------------
    ENDIF
		
			IF (JR.EQ.1) THEN
	     
	 DLOGE2(JR)=DLOGE1(JR)-2*W*VRE1(JR+1)/DR
	 
!	
!	 TE2(JR)=TE1(JR)-1.33*W*VRI1(JR+1)*TE1(JR)/DR &
!       +W*(TI1(JR)-TE1(JR))*COLLEI&
!       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR)) &
!       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
!       *SWTENE

	TE2(JR)=TE1(JR) &
       +W*(TI1(JR)-TE1(JR))*COLLEI&
       -7e5*W*SQRT(Te1(JR))*DENE1(JR)*(1+13.6/Te1(JR)) &
       +W*ENEJ/(2.4*TAU*DEN0*VOL)*EXP(-(JT-1)*W/TAU) &
       *SWTENE
        			 
		ENDIF

!	----end electrons--------------------
		ENDIF
!	--------------ions--------------------

!	-----------------------------------------------
	IF (JR.GT.1.AND.JR.LT.NRAD) THEN
	GRADVI=.5*VRI1(JR)*(VRI1(JR+1)-VRI1(JR-1))
	
	GRADPI(JR)=.5*(TI1(JR+1)-TI1(JR-1)+TI1(JR)* &
	(DLOGI1(JR+1)-DLOGI1(JR-1)))/DR
	
	GRADVTHI=.5*VRI1(JR)*(VTHI1(JR+1)-VTHI1(JR-1))
	ENDIF

	IF (JR.EQ.NRAD) THEN

	GRADPI(JR)=(TI1(JR)-TI1(JR-1)+TI1(JR)* &
	(DLOGI1(JR)-DLOGI1(JR-1)))/DR
	GRADVI=VRI1(JR)*(VRI1(JR)-VRI1(JR-1))
	GRADVTHI=VRI1(JR)*(VTHI1(JR)-VTHI1(JR-1))
	ENDIF


	IF (JR.GT.1.AND.JR.LT.NRAD) THEN

	DLOGI2(JR)=DLOGI1(JR)-W*VRI1(JR)/R- &
        .5*W/DR*(VRI1(JR+1)-VRI1(JR-1)+ &
        VRI1(JR)*(DLOGI1(JR+1)-DLOGI1(JR-1)))
        
!    TI2(JR)=TI1(JR)-.666*W*VRI1(JR)*TI1(JR)/R- &
!       .333*W/DR*TI1(JR)*(VRI1(JR+1)-VRI1(JR-1))- &
!       .5*W/DR*VRI1(JR)*(TI1(JR+1)-TI1(JR-1)) &
!       -W*(TI1(JR)-TE1(JR))*COLLEI

	TI2(JR)=TI1(JR) &
       -W*(TI1(JR)-TE1(JR))*COLLEI
              
		ENDIF
		

	IF (JR.EQ.NRAD) THEN

	DLOGI2(JR)=DLOGI1(JR)-W*VRI1(JR)/R- &
        W/DR*(VRI1(JR)-VRI1(JR-1)+ &
        VRI1(JR)*(DLOGI1(JR)-DLOGI1(JR-1)))
        
!         TI2(JR)=TI1(JR)-.666*W*VRI1(JR)*TI1(JR)/R- &
!       .666*W/DR*TI1(JR)*(VRI1(JR)-VRI1(JR-1))- &
!       W/DR*VRI1(JR)*(TI1(JR)-TI1(JR-1)) & 
!       -W*(TI1(JR)-TE1(JR))*COLLEI

	TI2(JR)=TI1(JR) &
       -W*(TI1(JR)-TE1(JR))*COLLEI

		ENDIF

!	-------------------------------
		
	IF (JR.EQ.1) THEN

	DLOGI2(JR)=DLOGI1(JR)-2*W*VRI1(JR+1)/DR 
	
!	TI2(JR)=TI1(JR)-1.332*W*VRI1(JR+1)*TI1(JR)/DR &
!       -W*(TI1(JR)-TE1(JR))*COLLEI

	TI2(JR)=TI1(JR) &
       -W*(TI1(JR)-TE1(JR))*COLLEI
	 
		ENDIF
			IF (JR.GT.1) THEN
	
	VRI2(JR)=(1-W*COLLHH)*VRI1(JR)+3e14*W*ERAD(JR) &
        -1e12*W*GRADPI(JR) &
        -W/DR*GRADVI+1e4*W*VTHI1(JR)*(BEXT+BINT(JR))
    
    VTHI2(JR)=(1-W*COLLHH)*VTHI1(JR)-1e4*W*VRI1(JR)* &
    (BEXT+BINT(JR))-W/DR*GRADVTHI
    
    VIMAX(JR)=(SQRT(ABS( &
	+6e14*ERAD(JR)*DR &
	-2e12*GRADPI(JR)*DR)))
    
    	ENDIF
!	-------end ions--------------------------
!   	  ENDIF
!	------------------------------------
	R=R+DR
	enddo

! 5	CONTINUE


!	--------nuovi valori----------------
	
	QTTIO=0
	QTTEL=0
	R=0
	DO JR=1,NRAD
!	---------------------------------
!	-------------------------------
	IF (JR.EQ.1) THEN
	VRI1(1)=0
	VTHI1(1)=0
	ENDIF
	
	DLOGI1(JR)=DLOGI2(JR)
	
	DENI1(JR)=DEN0*EXP(DLOGI1(JR))
		
		
		IF (VRI2(JR).LT.0.AND.ABS(VRI2(JR)).GE.VIMAX(JR)) THEN
	VRI1(JR)=-VIMAX(JR)
		ENDIF
	IF (VRI2(JR).GT.0.AND.ABS(VRI2(JR)).GE.VIMAX(JR)) THEN
	VRI1(JR)=+VIMAX(JR)
		ENDIF
	IF (VRI2(JR).GT.-VIMAX(JR).AND.VRI2(JR).LT.VIMAX(JR)) THEN
	VRI1(JR)=VRI2(JR)
	ENDIF	
	
	VTHI1(JR)=VTHI2(JR)
	
	IF (DENI2(JR).LT.1e-4*DENI1(1)) THEN
		TI1(JR)=.1
			ELSE
		TI1(JR)=TI2(JR)
			ENDIF
!	--------------------------------
	IF (SWTELE.EQ.1) THEN
!	-----------------------------------------
	IF (JR.EQ.1) THEN
	VRE1(1)=0
	VTHE1(1)=0
	ENDIF
		
		DLOGE1(JR)=DLOGE2(JR)
		
	IF (VRE2(JR).LT.0.AND.ABS(VRE2(JR)).GE.VEMAX(JR)) THEN
	VRE1(JR)=-VEMAX(JR)
		ENDIF
	IF (VRE2(JR).GT.0.AND.ABS(VRE2(JR)).GE.VEMAX(JR)) THEN
	VRE1(JR)=+VEMAX(JR)
		ENDIF
	IF (VRE2(JR).GT.-VEMAX(JR).AND.VRE2(JR).LT.VEMAX(JR)) THEN
	VRE1(JR)=VRE2(JR)
	ENDIF
			
!	-----------------------------

	DENE1(JR)=DEN0*EXP(DLOGE1(JR))
	
	IF (DENE2(JR).LT.1e-4*DENE1(1)) THEN
		TE1(JR)=.1
		ELSE
	TE1(JR)=TE2(JR)
		ENDIF
!	VTHE1(JR)=VTHE2(JR)
	
	IF (VTHE2(JR).LT.0.AND.ABS(VTHE2(JR)).GE.VEMAX(JR)) THEN
	VTHE1(JR)=-VEMAX(JR)
		ENDIF
	IF (VTHE2(JR).GT.0.AND.ABS(VTHE2(JR)).GE.VEMAX(JR)) THEN
	VTHE1(JR)=+VEMAX(JR)
		ENDIF
	IF (VTHE2(JR).GT.-VEMAX(JR).AND.VTHE2(JR).LT.VEMAX(JR)) THEN
	VTHE1(JR)=VTHE2(JR)
	ENDIF
			
!	----------------------------------------
!	---------------------------------
	IF (DENE1(JR)+DENI1(JR).GT.0) THEN
            DEBYE(JR)=2.5e-7*SQRT((Te1(JR)+Ti1(JR)) &
            /(DENE1(JR)+DENI1(JR)))
        ELSE
            DEBYE(JR)=DR
        ENDIF
        
     IF (DEBYE(JR).GT.DR) THEN
     DEBYE(JR)=DR
     ENDIF

	IF (DEBYE(JR).LT.DR) THEN
	
	DSEP=DENI1(JR) &
        *DEBYE(JR)/(DR+DEBYE(JR))
     
     IF (DENE1(JR).GT.DENI1(JR)+DSEP.AND. &
     DEBYE(JR).LT.DR) THEN
     DENE1(JR)=DENI1(JR)+DSEP
     DLOGE1(JR)=LOG(DENE1(JR)/DEN0)
     	ENDIF
     IF (DENE1(JR).LT.DENI1(JR)-DSEP.AND. &
     DEBYE(JR).LT.DR) THEN
     DENE1(JR)=DENI1(JR)-DSEP
     DLOGE1(JR)=LOG(DENE1(JR)/DEN0)
     	ENDIF
     IF (DENE1(JR).LE.DENI1(JR)+DSEP.AND. &
        DENE1(JR).GE.DENI1(JR)-DSEP) THEN
	DENE1(JR)=DEN0*EXP(DLOGE1(JR))
	ENDIF

			ENDIF
        CHGSEP(JR)=(DENI1(JR)-DENE1(JR)*SWTELE)
!	-----------------------------------
		ENDIF
!	------------------------------	
	
	QTTIO=QTTIO+DENI1(JR)*(R+DR)*DR
	
	QTTEL=QTTEL+DENE1(JR)*(R+DR)*DR


	R=R+DR
	enddo
!	---------------------------------------
	
!	-------------------------------------------
	IF (DENE1(1).LE..01*DEN0.AND.SWTELE.EQ.1) THEN
	SWTELE=0
	W=40*W
	ENDIF
! 6	CONTINUE
!	---------------------------------
!	IF (JT*W.GE.2e-9) THEN
!		ITMREC=1
!		ENDIF
!	-----------------------------------------------
!	IF (ITREC.GT.6e-8/W) ITMREC=1
	IF (JT.EQ.ITREC) THEN
	FLAGREC=1
	ITREC=ITREC+ITMREC
		ELSE
	FLAGREC=0
	ITREC=ITREC
	ENDIF
	
	IF (ABS(VRE1(2)).GT.5e9) THEN
	FLAGREC=1
	ENDIF
	IF (FLAGREC.EQ.1) THEN
!	---write profili-----------------
	OPEN(10,FILE='Profili.DAT',STATUS='unknown')
	R=0
	DO JR=1,NRAD
!	write(*,*) 'Sto scrivendo i profili allo step', JT
	write(10,7) R,Te1(JR),TI1(JR),DENE1(JR),DENI1(JR), &
	VRE1(JR),VRI1(JR),VTHE1(JR),VTHI1(JR),ERAD(JR),&
	CHGSEP(JR),BEXT+BINT(JR),&
	DENI1(JR)*VTHI1(JR)-DENE1(JR)*VTHE1(JR),TIME
	R=R+DR
	enddo
        flush(10)
        close(10)
	write(11,8) TIME,TE1(IREC),TI1(IREC),DENE1(IREC),VRE1(IREC), &
        DENI1(IREC),VRI1(IREC),BINT(IREC),QTTIO,QTTEL,SWTELE, &
        DENI1(NRAD-5)*VRI1(NRAD-5)-DENE1(NRAD-5)*VRE1(NRAD-5)
        flush(11)
       
		ENDIF
!	--------------------------------------------
	TIME=TIME+W
	IF (TIME.GE.TFIN) THEN
	GOTO 342
	ENDIF
	     enddo
!	------------------------------------
!	---end loop temporale------2----------
 342	write(*,*) 'ho finito, culo'
 		write(*,*) TIME,W,TFIN,KFIN,DENI1(1)
!	-------------------------------
!8	FORMAT(e10.4,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4, &
!        1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4)
7      FORMAT(15(e10.4,1X))
8      FORMAT(15(e10.4,1X))
 	STOP

END PROGRAM COLDMAGN
