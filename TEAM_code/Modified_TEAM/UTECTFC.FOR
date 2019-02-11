C-----* EJECTION IN A 2D-DUCT (UNSTEADY COMPRESSIBLE)*---------

C*****************************************************************
C*******************                          ********************
C*******************         UTECTFC          ********************
C*******************                          ********************
C*****************************************************************
C********   This code is modified off of the orginal   ***********
C********   TEAM code for a 2D jet.                    ***********
C********   Questions: Please contact hparish.         ***********
C*****************************************************************
C*****************************************************************
      INCLUDE 'UCOMMON.INC'
      CALL SET
      CALL GRID
      CALL INLET
C********** START ITERATION
      TIME=0.0D0
    9 TIME=TIME+DT
   10 NITER=NITER+1
      CALL CALCU
      CALL CALCV
      CALL CALCP
      IF(NE.EQ.1) CALL CALCE
      IF(NTF.EQ.1) CALL CALCTE
      IF(NTF.EQ.1) CALL CALCED
      CALL PROPS
      CALL OMB
      RESORM=RESORM/FLOWIN
      RESORU=RESORU/XMONIN
      RESORV=RESORV/XMONIN
      RESORT=RESORT/XTMIN
      RESORK=RESORK/XKMIN
      RESORE=RESORE/XEMIN
      SORCE=DMAX1(RESORU,RESORV,RESORM,RESORT)
      WRITE(6,31)NITER,TIME,FLOWE,RESORU,RESORV,RESORM,RESORT,RESORK
      IF(SORCE.GE.1.E+10.AND.NITER.GT.30)STOP
      IF(SORCE.LE.SORMAX)GO TO 20
      GO TO 10
   20 CONTINUE
C************NEXT TIME STEP      
      DO 1 I=1,NI
      DO 1 J=1,NJ
      UO(I,J)=U(I,J)
      VO(I,J)=V(I,J)
      TMO(I,J)=TM(I,J)
      PO(I,J)=P(I,J)
      DENO(I,J)=DEN(I,J)
      TEO(I,J)=TE(I,J)
      EDO(I,J)=ED(I,J)
   1  CONTINUE
C************ENTRA COND. CORR. 
      DO 2 J=1,NJ
      TMO(1,J)=TMAE
      PO(2,J)=PAE
      PO(1,J)=PO(2,J)
      P(2,J)=PO(2,J)
      P(1,J)=P(2,J)
      TM(1,J)=TMO(1,J)
   2  CONTINUE
C********CHECK TIME
      IF(TIME.LT.TMAX) GOTO 9
C********** END ITERATION
   31 FORMAT(1X,I4,1P10E10.3)
C***************************DATA FILES********************************
      PI=ATAN(1.0)*4.0
      DEN1=PATM/(GC*273.15D0)
      XK=2.926D-4/DEN1
      XF=XL/(NIM1-1)
      DEN0=PSTGE/(GC*TSTGE)
      OPEN(9,FILE='UPAR2.DAT')
      OPEN(11,FILE='UE2.DAT')
      WRITE(11,*) 'TITLE="X" '
      WRITE(11,*) 'VARIABLES=X,Y,U,V,VEC,MACH,P,TM,DEN,TE,ED,TI,LS,VIS,
     +DEF,DELA,DELAX,DEL2A,DEL2X'
      WRITE(11,*) 'GEOMETRY F=POINT,M=GRID,C=WHITE,S=LOCAL'
      WRITE(11,*) '1'
      WRITE(11,*) NI
      DO 5 I=1,NI
      WRITE(11,899) X(I),Y(1)
    5 CONTINUE
      WRITE(11,*) 'GEOMETRY F=POINT,M=GRID,C=WHITE,S=LOCAL'
      WRITE(11,*) '1'
      WRITE(11,*) IJ1P1
      DO 55 I=1,IJET1
      WRITE(11,899) XU(I),Y(NJ)
   55 CONTINUE
      XU1=XU(IJET1)-YL*DCOS(AL*PI/180)/4.D0
      YN1=Y(NJ)+YL*DSIN(AL*PI/180)/4.D0
      WRITE(11,899) XU1,YN1
      WRITE(11,*) 'GEOMETRY F=POINT,M=GRID,C=WHITE,S=LOCAL'
      WRITE(11,*) '1'
      WRITE(11,*) NI-IJ2M1
      XU2=DJ*DSIN(AL*PI/180)+XU1
      YN2=DJ*DCOS(AL*PI/180)+YN1
      WRITE(11,899) XU2,YN2
      DO 555 I=IJET2,NIM1        
      WRITE(11,899) XU(I),Y(NJ)
  555 CONTINUE
  899 FORMAT(1X,2(D17.8,2X))
      WRITE(11,*) 'ZONE t="1" '
      WRITE(11,1008) NI,NJ
 1008 FORMAT(1X,'I=',I4,'J=',I4,'F=POINT')
      DO 113 J=2,NJM1
      DO 113 I=1,NI
      DENN=DEN(I,J)+FY(J)*(DEN(I,J+1)-DEN(I,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      DELA(I,J)=-1.D1*XK*(DENN-DENS)/(SYCV(J)*(1+XK*DEN(I,J)))
      D2=-1.D0*XK*(DEN(I,J+1)-DEN(I,J))/(VYCV(J)*(1+XK*DENN))
      D1=-1.D0*XK*(DEN(I,J)-DEN(I,J-1))/(VYCV(J-1)*(1+XK*DENS))
      DEL2A(I,J)=(D2-D1)/SYCV(J)
      IF(J.EQ.2) THEN
      DELA(I,1)=DELA(I,2)
      DEL2A(I,1)=DEL2A(I,2)
      END IF
      IF(J.EQ.NJM1) THEN
      DELA(I,NJ)=DELA(I,NJM1)
      DEL2A(I,NJ)=DEL2A(I,NJM1)
      END IF
  113 CONTINUE
      DO 111 J=1,NJ
      DO 112 I=1,NI
      IF(I.EQ.NI) THEN
      DELAX(NI,J)=DELAX(NIM1,J)
      DEL2X(NI,J)=DEL2X(NIM1,J)
      GO TO 898
      END IF
      IF(I.EQ.1) GO TO 898
      DENE=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DELAX(I,J)=-1.D1*XK*(DENE-DENW)/(SXCV(I)*(1+XK*DEN(I,J)))
      D2=-1.D0*XK*(DEN(I+1,J)-DEN(I,J))/(UXCV(I)*(1+XK*DENE))
      D1=-1.D0*XK*(DEN(I,J)-DEN(I-1,J))/(UXCV(I-1)*(1+XK*DENW))
  898 IF(I.EQ.2) THEN
      DELAX(1,J)=DELAX(2,J)
      DEL2X(1,J)=DEL2X(2,J)
      END IF
C     XC=X(I)/XF
C     IC=X(I)/XF
C     IB=IC+1
      IF(I.EQ.2*(I/2)) THEN
C     XG(I)=IB-XC
      XG(I)=1.0D0
      ELSE
C     XG(I)=XC-IB+1
      XG(I)=0.0D0
      END IF
      IF(DABS(V(I,J)).LT.1.E-33)V(I,J)=0.0D0
      V1=V(I,J)
      IF(J.EQ.1)GOTO 118
      V1=(V(I,J)+V(I,J-1))/2
  118 IF(DABS(U(I,J)).LT.1.E-33)U(I,J)=0.0D0
      U1=U(I,J)
      IF(I.EQ.1)GOTO 119
      IF(J.EQ.NJ) THEN
      IF(I.GT.IJET1.AND.I.LE.IJET2) THEN
      U1=DABS(VIN)*DCOS(AL*PI/180)
      GOTO 119
      END IF
      END IF
      U1=(U(I,J)+U(I-1,J))/2
  119 VEC=SQRT(U1**2+V1**2)
      DEX=1.D0*XK*(DEN(I,J)-DEN0)/(1.+XK*DEN(I,J))
      JD=DEX/XF
      DEX=DEX/XF
      DEX=DEX-JD
      XG(I)=DABS(DABS(XG(I)-2*(DEX-0.5D0))-1.0D0)
      TI(I,J)=100*SQRT(2*TE(I,J))/VIN
      P(I,J)=P(I,J)+PATM
      VMCH=VEC/SQRT(GC*TM(I,J)*SG)
      VMJ=VIN/SQRT(SG*GC*TMAJ)
      VME=UAV/SQRT(SG*GC*TMAE)
      SL=0.09*TE(I,J)**(1.5)/(ED(I,J)+TINY)
      WRITE(11,1000) X(I),Y(J),U1,V1,VEC,VMCH,P(I,J),TM(I,J),DEN(I,J)
     +,TE(I,J),ED(I,J),TI(I,J),SL,VIS(I,J),XG(I),DELA(I,J),DELAX(I,J)
     +,DEL2A(I,J),DEL2X(I,J)
 1000 FORMAT(1X,23(D17.8,2X)) 
  112 CONTINUE
  111 CONTINUE
      END FILE 11
      PMJ=P(IJ1P1,NJ)
      DO 1035 I=IJ1P1,IJET2
      IF(PMJ.LT.P(I,NJ)) PMJ=P(I,NJ)
 1035 CONTINUE
      FTE=0.0D0
      DO 1036 J=2,NJM1
      FTE=FTE+DEN(1,J)*SYCV(J)*U(1,J)*TE(1,J)
 1036 CONTINUE
      TKEEN=FTE/FLOWE
      PAE=PAE+PATM
      PSTGJ=PMJ*(1.+(SG-1.)*(VMJ**2)/2.)**(SG/(SG-1.))
      TIJ=100*SQRT(2*TKEIN)/VIN
      TIAE=100*SQRT(2*TKEEN)/VIN
      WRITE(9,1010) NI
      WRITE(9,1011) NJ
      WRITE(9,1012) TIME
      WRITE(9,1013) SORCE
      WRITE(9,1014) XL
      WRITE(9,1015) YL
      WRITE(9,1016) XJET
      WRITE(9,1017) DJ
      WRITE(9,1018) AL
      WRITE(9,1019) VIN
      WRITE(9,1020) UAV
      WRITE(9,1021) VMJ
      WRITE(9,1022) VME
      WRITE(9,1023) TSTGJ
      WRITE(9,1024) TSTGE
      WRITE(9,1025) PSTGJ
      WRITE(9,1026) PSTGE
      WRITE(9,1027) TMAJ
      WRITE(9,1028) TMAE
      WRITE(9,1029) PMJ
      WRITE(9,1030) PAE
      WRITE(9,1031) FLOWIN
      WRITE(9,1032) FLOWE
      WRITE(9,1033) TIJ
      WRITE(9,1034) TIAE
 1010 FORMAT(1X,'NO.GRID X-DIR.=',I5) 
 1011 FORMAT(1X,'NO.GRID Y-DIR.=',I5) 
 1012 FORMAT(1X,'TIME=',D17.7) 
 1013 FORMAT(1X,'MAX. ERROR=',D17.7) 
 1014 FORMAT(1X,'DUCT LENGTH=',D17.7) 
 1015 FORMAT(1X,'DUCT WIDTH=',D17.7) 
 1016 FORMAT(1X,'JET POS.=',D17.7) 
 1017 FORMAT(1X,'JET DIA.=',D17.7) 
 1018 FORMAT(1X,'JET ANGLE=',D17.7) 
 1019 FORMAT(1X,'JET VEL.=',D17.7) 
 1020 FORMAT(1X,'ENTRAIN AVE.VEL.=',D17.7) 
 1021 FORMAT(1X,'JET MACH=',D17.7) 
 1022 FORMAT(1X,'ENTRAIN AVE.MACH=',D17.7) 
 1023 FORMAT(1X,'JET STAG.TEMP.=',D17.7) 
 1024 FORMAT(1X,'ENTRAIN STAG.TEMP.=',D17.7) 
 1025 FORMAT(1X,'JET STAG.PRES.=',D17.7) 
 1026 FORMAT(1X,'ENTRAIN STAG.PRES.=',D17.7) 
 1027 FORMAT(1X,'JET AVE.TEMP.=',D17.7) 
 1028 FORMAT(1X,'ENTRAIN AVE.TEMP.=',D17.7) 
 1029 FORMAT(1X,'JET MAX.PRES.=',D17.7) 
 1030 FORMAT(1X,'ENTRAIN AVE.PRES.=',D17.7) 
 1031 FORMAT(1X,'JET M.FLOW RATE=',D17.7) 
 1032 FORMAT(1X,'ENTRAIN M.FLOW RATE=',D17.7) 
 1033 FORMAT(1X,'JET TURB.INTENSITY(%)=',D17.7) 
 1034 FORMAT(1X,'ENTRAIN AVE.TURB.INTENSITY(%)=',D17.7) 
      END FILE 9
      STOP
      END
C*****************************PARAMETERS******************************
      SUBROUTINE USER
      INCLUDE 'UCOMMON.INC'
C********** SET DEFAULT VALUES
      ENTRY SET
      PI=ATAN(1.0)*4.0
      URFU=0.4D0
      URFV=0.4D0
      URFP=0.4D0
      URFK=0.2D0
      URFGEN=1.0D0
      URFGTM=1.0D0
      URFE=0.2D0
      URFTM=1.0D0
      URFVIS=1.0D0
      URFDEN=1.0D0
      NSWPU=1
      NSWPV=1
      NSWPP=2
      NSWPT=1
      NSWPK=1
      NSWPE=1
      NCSA=1
      NCF=1
      NTF=1
      NE=1
      NQWT=1
      NQET=0
      NQNT=0
      NQST=0
      NQUV=0
      NQKE=0
      NITER=0
      MAXIT=5000
      IPREF=1
      JPREF=1
      IMON=7
      JMON=7
      SORMAX=1.D-3
      SG=1.4D0
      GC=287.0D0
      SPHP=1007.0D0
      SPHV=SPHP-GC
      PRANL=0.7D0
      PRANT=0.9D0
      VISCOS=1.8D-5
      DENSIT=1.18473D0
      AQC=26.0D0
      DT=1.0D-1
      TMAX=0.2D0
      XJET=0.2D0
      DJ=.01D0
      AL=10.0D0
      VIN=526.1136336D0
      UIN=0.0D0
      PATM=101.325D3
      TSTGE=298.0D0
      TSTGJ=310.0D0
      PSTGE=PATM
      IF(NTF.EQ.1) THEN
      TKEIN=2.0D-3*VIN*VIN
      TEDIN=TKEIN**(1.5)/0.3       
      ELSE
      TKEIN=0.0D0
      TEDIN=0.0D0
      END IF
      CMU=0.09D0
      C1=1.44D0
      C2=1.92D0
      CAPPA=.4187D0
      ELOG=9.793D0
      PRED=1.3D0
C     PRED=CAPPA**2/(C2-C1)/SQRT(CMU)
      PRTE=1.0D0
      CMU25=CMU**0.25
      CMU75=CMU**0.75
      RETURN
C********** GRID SPECIFICATION & CALCULATE GEOMETRICAL QUANTITIES
      ENTRY GRID
      R(1)=1.D0
      IND=1
      NIM1=M-1
      NJM1=N-1
      XL=4.0D0
      YL=0.1D0
      X0=XJET-DJ/(2.*DSIN(AL*PI/180))
      X1=XJET+DJ/(2.*DSIN(AL*PI/180))
      HSTEP=0.0D0
C     RAJET=1.0D0
      RACL=1.0D0
      RAY=1.0/RACL
      RAXL=1.0D0
      RAXR=1.0D0
      RAXJ=1.0D0
      IJET1=15
      IJET2=20
      IJ1P1=IJET1+1
      IJ1M1=IJET1-1
      IJ2P1=IJET2+1
      IJ2M1=IJET2-1
      JSTEP=1
      JSTP1=JSTEP+1
      JSTM1=JSTEP-1
      IL=IJET1
      ILP1=IL+1
      XU(1)=0.0D0
      YV(1)=0.0D0
      XU(IJET1)=X0
      XU(IJET2)=X1
      XU(NIM1)=XL
      YV(NJM1)=YL
      TOY=0.0D0
C     DO 103 J=2,JSTEP
C 103 TOY=TOY+RAJET**(J-2)
C     DY=HSTEP/TOY
C     DO 104 J=2,JSTEP
C     JJ=JSTP1-J
C 104 YV(JJ)=YV(JJ+1)-DY*RAJET**(J-2)
      TOY=0.0D0
      JCL=NJM1
C     JCL=(NJM1-JSTEP)/2+JSTEP
      JCLP1=JCL+1
      TOY=0.0D0
      DO 82 J=JSTP1,JCL
   82 TOY=TOY+RACL**(J-JSTP1)
      DY=YL/TOY
      DO 124 J=JSTP1,JCL
  124 YV(J)=YV(J-1)+DY*RACL**(J-JSTP1)
      IF(JCL.EQ.NJM1) GOTO 923
      TOY=0.0D0
      DO 105 J=JCLP1,NJM1
  105 TOY=TOY+RAY**(J-JCLP1)
      DY=(YL-0.5*YL)/TOY
      DO 106 J=JCLP1,NJM1
  106 YV(J)=YV(J-1)+DY*RAY**(J-JCLP1)
  923 TOX=0.0D0
      DO 107 I=2,IL
  107 TOX=TOX+RAXL**(I-2)
      DX=X0/TOX
      DO 108 I=2,IL
  108 XU(I)=XU(I-1)+DX*RAXL**(I-2)
      TOX=0.0D0
      DO 111 I=ILP1,IJET2
  111 TOX=TOX+RAXJ**(I-ILP1)
      DX=(X1-X0)/TOX
      DO 132 I=ILP1,IJET2
  132 XU(I)=XU(I-1)+DX*RAXJ**(I-ILP1)
      TOX=0.0D0
      DO 109 I=IJ2P1,NIM1
  109 TOX=TOX+RAXR**(I-IJ2P1)
      DX=(XL-X1)/TOX
      DO 110 I=IJ2P1,NIM1
  110 XU(I)=XU(I-1)+DX*RAXR**(I-IJ2P1)
      XU(1)=0.0D0
      YV(1)=0.0D0
      CALL INIT
      RETURN
C********** INLET CONDITION & INITIAL GUESS
      ENTRY INLET
      PI=ATAN(1.0)*4.0
      PREF=0.0D0
      DO 2213 J=1,NJ
      IF(NQWT.EQ.1) THEN
      TM(1,J)=298.0D0
      ELSE 
      QUW(J)=0.0D0
      END IF
      MQW(J)=NQWT
 2213 CONTINUE
      DO 2214 J=1,NJ
      IF(NQET.EQ.1) THEN
      TM(NI,J)=298.D0
      ELSE      
      QUE(J)=0.0D0
      END IF
      MQE(J)=NQET
 2214 CONTINUE
      DO 213 I=1,NI
      IF(NQST.EQ.1) THEN
      TM(I,1)=298.0D0
      ELSE 
      QUS(I)=0.0D0
      END IF
      MQS(I)=NQST
      U(I,1)=0.0D0
      V(I,1)=0.0D0
  213 U(I,NJ)=0.0D0
      DO 235 J=2,NJM1
      P(2,J)=PREF
      P(1,J)=P(2,J)
      TM(1,J)=TSTGE
      P(NI,J)=PREF
  235 CONTINUE
      DO 234 I=1,IJET1
      IF(NQNT.EQ.1) THEN
      TM(I,NJ)=298.0D0
      ELSE 
      QUN(I)=0.0D0
      END IF
      MQN(I)=NQNT
      V(I,NJM1)=0.0D0
  234 V(I,NJ)=0.0D0
      DO 2434 I=IJ2P1,NI
      IF(NQNT.EQ.1) THEN
      TM(I,NJ)=298.0D0
      ELSE 
      QUN(I)=0.0D0
      END IF
      MQN(I)=NQNT
      V(I,NJM1)=0.0D0
 2434 V(I,NJ)=0.0D0
      DO 216 I=IJ1P1,IJET2
      IF(XU(IJET2)-XU(I).LT.TINY) GOTO 246
      U(I,NJ)=VIN*DCOS(AL*PI/180)
  246 CONTINUE
      V(I,NJM1)=-VIN*DSIN(AL*PI/180)
      V(I,NJ)=-VIN*DSIN(AL*PI/180)
C     TMAJ=298.0D0
      TMAJ=TSTGJ-(SG-1.)*(VIN**2)/(2.*SG*GC)
      MQN(I)=1
      TM(I,NJ)=TMAJ
  216 CONTINUE   
      IF(NTF.EQ.1) THEN
      TLS=0.6*DJ
      DO 145 I=IJ1P1,IJET2
      U1=VIN*DCOS(AL*PI/180)
      V1=-VIN*DSIN(AL*PI/180)
      TE(I,NJ)=2.0D-03*(V1**2+U1**2)
      ED(I,NJ)=TE(I,NJ)**(1.5)/TLS
  145 CONTINUE
      END IF
      JSTP1=JSTEP+1
C     DELTA=0.50*HSTEP
C     TKEIN1=0.8E-02*UIN*UIN
C     A=(TKEIN-TKEIN1)/DELTA
C     DO 144 J=JSTP1,NJM1
C     YY=Y(J)-HSTEP
C     IF(YY.GT.HSTEP) YY=2.0*HSTEP-YY
C     U(1,J)=UIN*(YY/DELTA)**(1.0/5.6)
C     U(1,J)=UIN
C     IF(YY.LE.HSTEP.AND.YY.GE.DELTA) U(1,J)=UIN
C     TE(1,J)=A*YY+TKEIN1
C     TE(1,J)=TKEIN1
C     IF(YY.LE.HSTEP.AND.YY.GE.DELTA) TE(1,J)=TKEIN
C     TLS1=CAPPA*YL
C     ED(1,J)=TE(1,J)**(1.5)/TLS1
C     ED(1,J)=TE(1,J)**(1.5)/(CMU25*.012*.051)
C     IF(YY.LT.DELTA) ED(1,J)=TE(1,J)**(1.5)/TLS1
C 144 CONTINUE
      DO 202 J=2,NJM1
      DO 202 I=2,NI
      U(I,J)=UIN
      TE(I,J)=TKEIN
      ED(I,J)=TEDIN
  202 CONTINUE
      DO 503 I=2,NIM1
      YPLUSS(I)=11.0D0
  503 YPLUSN(I)=11.0D0
C     DO 204 J=1,JSTEP
C 204 XPLUSW(J)=11.0D0
C----CALCULATE INLET MASS FLOWRATE--------------------------------------
      FLOWIN=0.0D0
      XMONIN=0.0D0
C     DO 657 J=2,NJM1
C     ARDEN=DEN(1,J)*SYCV(J)
C     FLOWIN=FLOWIN+ARDEN*U(1,J)
C     XMONIN=XMONIN+ARDEN*U(1,J)**2
C 657 CONTINUE
      DO 350 I=1,NI
      DO 350 J=1,NJ
      IF(NCF.EQ.1) THEN
      DEN(I,J)=(P(I,J)+PATM)/(GC*TM(I,J))
      ELSE
      DEN(I,J)=DENSIT
      END IF
      IF(NE.EQ.1) THEN
      VISCO(I,J)=1.458D-6*TM(I,J)**(1.5)/(TM(I,J)+110.4D0)
      ELSE
      VISCO(I,J)=VISCOS
      END IF
      TVIS(I,J)=CMU*DEN(I,J)*TE(I,J)**2/(ED(I,J)+TINY)
      IF(TVIS(I,J).LT.TINY) TVIS(I,J)=TINY
      VIS(I,J)=TVIS(I,J)+VISCO(I,J)
      GAT(I,J)=TVIS(I,J)/PRANT+VISCO(I,J)/PRANL
      UO(I,J)=U(I,J)
      VO(I,J)=V(I,J)
      DENO(I,J)=DEN(I,J)
      PO(I,J)=P(I,J)
      TMO(I,J)=TM(I,J)
      TEO(I,J)=TE(I,J)
      EDO(I,J)=ED(I,J)
  350 CONTINUE
      DO 658 I=IJ1P1,IJET2
      ARDEN=DEN(I,NJ)*SXCV(I)
      FLOWIN=FLOWIN+ARDEN*DABS(V(I,NJ))
      XMONIN=XMONIN+ARDEN*DABS(V(I,NJ))*VIN
  658 CONTINUE
      RETURN
C********** UPDATE PROPERTIES
      ENTRY PROPS
      DO 900 I=1,NI
      DO 900 J=1,NJ
      IF(NCF.EQ.1) THEN
      DENNEW=(P(I,J)+PATM)/(GC*TM(I,J))
      DEN(I,J)=DEN(I,J)+URFDEN*(DENNEW-DEN(I,J))
      END IF
      IF(NE.EQ.1) THEN
      VISCONW=1.458D-6*TM(I,J)**(1.5)/(TM(I,J)+110.4D0)
      VISCO(I,J)=VISCO(I,J)+URFVIS*(VISCONW-VISCO(I,J))
      END IF
      TVISNW=CMU*DEN(I,J)*TE(I,J)**2/(ED(I,J)+TINY)
      IF(TVISNW.LT.TINY) TVISNW=TINY
      TVIS(I,J)=TVIS(I,J)+URFVIS*(TVISNW-TVIS(I,J))
      VIS(I,J)=TVIS(I,J)+VISCO(I,J)
      IF(NE.EQ.1) THEN
      GAT(I,J)=TVIS(I,J)/PRANT+VISCO(I,J)/PRANL
      END IF
  900 CONTINUE
      RETURN
C********** U-VELOCITY B C
      ENTRY MODU
C
C-----TOP WALL
      J=NJM1
      DO 210 I=2,IJ1M1
      YP=0.5*SYCV(J)
      YPLUSA=YPLUSN(I)+FX(I)*(YPLUSN(I+1)-YPLUSN(I))
      SQRTK=SQRT(TE(I,J)+FX(I)*(TE(I+1,J)-TE(I,J)))
      DENU=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      VISC=VISCO(I,J)+FX(I)*(VISCO(I+1,J)-VISCO(I,J))
      IF (YPLUSA.LE.11.63) GO TO 211
      TMULT=DENU*CMU25*SQRTK*CAPPA/DLOG(ELOG*YPLUSA)
      GO TO 212
  211 TMULT=VISC/YP
  212 TAUN(I)=-TMULT*U(I,J)
      SP(I,J)=SP(I,J)-TMULT*UXCV(I)
  210 AN(I,J)=0.0D0
      TAUN(IJET1)=TAUN(IJ1M1)
      J=NJM1
      DO 223 I=IJ2P1,NIM1
      YP=0.5*SYCV(J)
      SQRTK=SQRT(TE(I,J)+FX(I)*(TE(I+1,J)-TE(I,J)))
      YPLUSA=YPLUSN(I)+FX(I)*(YPLUSN(I+1)-YPLUSN(I))
      VISC=VISCO(I,J)+FX(I)*(VISCO(I+1,J)-VISCO(I,J))
      DENU=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      IF (YPLUSA.LE.11.63) GO TO 221
      TMULT=DENU*CMU25*SQRTK*CAPPA/DLOG(ELOG*YPLUSA)
      GO TO 222
  221 TMULT=VISC/YP
  222 TAUN(I)=-TMULT*U(I,J)
      SP(I,J)=SP(I,J)-TMULT*UXCV(I)
  223 AN(I,J)=0.0D0
      TAUN(IJET2)=TAUN(IJ2P1)
      TAUN(1)=TAUN(2)
      TAUN(NI)=TAUN(NIM1)
C------BOTTOM WALL
      J=2
      DO 219 I=2,NIM1
      YP=0.5*SYCV(J)
      SQRTK=SQRT(TE(I,J)+FX(I)*(TE(I+1,J)-TE(I,J)))
      YPLUSA=YPLUSS(I)+FX(I)*(YPLUSS(I+1)-YPLUSS(I))
      DENU=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      VISC=VISCO(I,J)+FX(I)*(VISCO(I+1,J)-VISCO(I,J))
      IF (YPLUSA.LE.11.63) GO TO 214
      TMULT=DENU*CMU25*SQRTK*CAPPA/DLOG(ELOG*YPLUSA)
      GO TO 215
  214 TMULT=VISC/YP
  215 TAUS(I)=-TMULT*U(I,J)
      SP(I,J)=SP(I,J)-TMULT*UXCV(I)
  219 AS(I,J)=0.0D0
      TAUS(1)=TAUS(2)
      TAUS(NI)=TAUS(NIM1)
C-----WEST BOUNDARY(ENTRAINMENT)
      DO 256 J=2,NJM1
      IF(NCF.EQ.1) THEN
      DENWE=DEN(1,J)
      DENWW=DEN(2,J)+FX(2)*(DEN(3,J)-DEN(2,J))
      SP(2,J)=SP(2,J)-AW(2,J)*(1.D0-DENWW/DENWE)
      END IF
      AW(2,J)=0.0D0
  256 CONTINUE
C-----OUTLET
      DO 206 J=2,NJM1
      U(NI,J)=U(NIM1,J)
      AE(NIM1,J)=0.0D0
  206 CONTINUE
      RETURN
C
C
      ENTRY MODV
      RETURN
C
C
      ENTRY CORTE
      DO 307 J=2,NJM1
      TE(NI,J)=TE(NIM1,J)
  307 CONTINUE
      DO 362 J=2,NJM1
      TE(1,J)=TE(2,J)
  362 CONTINUE   
C     DO 363 I=2,NIM1
C     TE(I,1)=TE(I,2)
C 363 CONTINUE   
      RETURN

      ENTRY CORED
      DO 304 J=2,NJM1
      ED(NI,J)=ED(NIM1,J)
  304 CONTINUE
      DO 364 J=2,NJM1
      ED(1,J)=ED(2,J)
  364 CONTINUE   
C     DO 365 I=2,NIM1
C     ED(I,1)=ED(I,2)
C 365 CONTINUE   
      RETURN
      
      ENTRY CORE 
      DO 845 J=2,NJM1
      IF(MQE(J).EQ.0) THEN
      TM(NI,J)=TM(NIM1,J)
      END IF
  845 CONTINUE   
      DO 836 J=2,NJM1
      IF(MQW(J).EQ.0) THEN 
      TM(1,J)=TM(2,J)
      END IF
  836 CONTINUE
      DO 837 I=1,NI
      IF(MQS(I).EQ.0) THEN 
      TM(I,1)=TM(I,2)
      END IF
  837 CONTINUE
      DO 838 I=1,IJET1
      IF(MQN(I).EQ.0) THEN 
      TM(I,NJ)=TM(I,NJM1)
      END IF
  838 CONTINUE
      DO 839 I=IJ2P1,NI
      IF(MQN(I).EQ.0) THEN 
      TM(I,NJ)=TM(I,NJM1)
      END IF
  839 CONTINUE
      RETURN

      ENTRY MODP
      DO 911 I=2,NIM1
      DO 911 J=2,NJM1
      AREAW=RSYCV(J)
      AREAS=SXCV(I)*RV(J-1)
      UKW=1.0/(SG*GC*TM(I-1,J))
      UKP=1.0/(SG*GC*TM(I,J))
      UKE=1.0/(SG*GC*TM(I+1,J))
      UKS=1.0/(SG*GC*TM(I,J-1))
      UKN=1.0/(SG*GC*TM(I,J+1))
      SP(I,J)=(UKE*DMAX1(-U(I,J),0.0D0)-UKP*DMAX1(U(I,J),0.0D0))*AREAW
     ++(UKW*DMAX1(U(I-1,J),0.0D0)-UKP*DMAX1(-U(I-1,J),0.0D0))*AREAW
     ++(UKN*DMAX1(-V(I,J),0.0D0)-UKP*DMAX1(V(I,J),0.0D0))*AREAS
     ++(UKS*DMAX1(V(I,J-1),0.0D0)-UKP*DMAX1(-V(I,J-1),0.0D0))*AREAS
  911 CONTINUE  
      DO 912 I=1,NI
      UKN=1.0/(SG*GC*TM(I,NJ))
      AREAS=SXCV(I)*RV(NJM1-1)
      SP(I,NJM1)=SP(I,NJM1)-UKN*DMAX1(-V(I,NJM1),0.0D0)*AREAS
      AN(I,NJM1)=0.0D0
      AS(I,2)=0.0D0
  912 CONTINUE
C     DO 913 J=1,NJ
C     AE(NIM1,J)=0.0D0
C 913 CONTINUE
      RETURN
      
      ENTRY MODE
C-----TOP WALL
      PI=ATAN(1.0D0)*4.0D0
      J=NJM1
      DO 810 I=2,IJET1
      YP=0.5*SYCV(J)
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GENT(I,J)*VOL/SPHP
      GENCOU=(DABS(TAUN(I-1)*U(I-1,J))+DABS(TAUN(I)*U(I,J)))/(2*YP)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GENT(I,J)-VIS(I,J)*(DVDX+DUDY)**2    
      GENT(I,J)=GENCOU+GENRES
      SU(I,J)=SU(I,J)+GENT(I,J)*VOL/SPHP
      SQRTK=SQRT(TE(I,J))
      YPLUSA=YPLUSN(I)
      VISC=VISCO(I,J)
      IF (YPLUSA.LE.11.63) GO TO 811
      TPLUS=PRANT*DLOG(ELOG*YPLUSA)/CAPPA
     ++PRANT*PI*SQRT(AQC/CAPPA)*(PRANL/PRANT-1.)
     +*(PRANT/PRANL)**(0.25)/(4.*DSIN(PI/4.0))
      GO TO 812
  811 TPLUS=PRANL*YPLUSA
  812 IF(MQN(I).NE.1) THEN 
      TM(I,NJ)=-(TPLUS/YPLUSA)*QUN(I)*YP/(SPHP*VISC)+TM(I,NJM1)
      ELSE
      QUN(I)=(YPLUSA/TPLUS)*SPHP*VISC*(TM(I,NJM1)-TM(I,NJ))/YP
      END IF
      SU(I,J)=SU(I,J)+(YPLUSA/TPLUS)*TM(I,NJ)*VISC*SXCV(I)/YP
      SP(I,J)=SP(I,J)-(YPLUSA/TPLUS)*VISC*SXCV(I)/YP
      AN(I,J)=0.0D0
  810 CONTINUE  
      J=NJM1
      DO 823 I=IJ2P1,NIM1
      YP=0.5*SYCV(J)
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GENT(I,J)*VOL/SPHP
      GENCOU=(DABS(TAUN(I-1)*U(I-1,J))+DABS(TAUN(I)*U(I,J)))/(2*YP)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GENT(I,J)-VIS(I,J)*(DVDX+DUDY)**2
      GENT(I,J)=GENCOU+GENRES
      SU(I,J)=SU(I,J)+GENT(I,J)*VOL/SPHP
      SQRTK=SQRT(TE(I,J))
      YPLUSA=YPLUSN(I)
      VISC=VISCO(I,J)
      IF (YPLUSA.LE.11.63) GO TO 821
      TPLUS=PRANT*DLOG(ELOG*YPLUSA)/CAPPA
     ++PRANT*PI*SQRT(AQC/CAPPA)*(PRANL/PRANT-1.)
     +*(PRANT/PRANL)**(0.25)/(4.*DSIN(PI/4.0))
      GO TO 822
  821 TPLUS=PRANL*YPLUSA
  822 IF(MQN(I).NE.1) THEN 
      TM(I,NJ)=-(TPLUS/YPLUSA)*QUN(I)*YP/(SPHP*VISC)+TM(I,NJM1)
      ELSE
      QUN(I)=(YPLUSA/TPLUS)*SPHP*VISC*(TM(I,NJM1)-TM(I,NJ))/YP
      END IF
      SU(I,J)=SU(I,J)+(YPLUSA/TPLUS)*TM(I,NJ)*VISC*SXCV(I)/YP
      SP(I,J)=SP(I,J)-(YPLUSA/TPLUS)*VISC*SXCV(I)/YP
      AN(I,J)=0.0D0
  823 CONTINUE  
C------BOTTOM WALL
      J=2
      DO 819 I=2,NIM1
      YP=0.5*SYCV(J)
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GENT(I,J)*VOL/SPHP
      GENCOU=(DABS(TAUS(I-1)*U(I-1,J))+DABS(TAUS(I)*U(I,J)))/(2*YP)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GENT(I,J)-VIS(I,J)*(DVDX+DUDY)**2
      GENT(I,J)=GENCOU+GENRES
      SU(I,J)=SU(I,J)+GENT(I,J)*VOL/SPHP
      SQRTK=SQRT(TE(I,J))
      YPLUSA=YPLUSS(I)
      VISC=VISCO(I,J)
      IF (YPLUSA.LE.11.63) GO TO 814
      TPLUS=PRANT*DLOG(ELOG*YPLUSA)/CAPPA
     ++PRANT*PI*SQRT(AQC/CAPPA)*(PRANL/PRANT-1.)
     +*(PRANT/PRANL)**(0.25)/(4.*DSIN(PI/4.0))
      GO TO 815
  814 TPLUS=PRANL*YPLUSA
  815 IF(MQS(I).NE.1) THEN 
      TM(I,1)=-(TPLUS/YPLUSA)*QUS(I)*YP/(SPHP*VISC)+TM(I,2)
      ELSE
      QUS(I)=(YPLUSA/TPLUS)*SPHP*VISC*(TM(I,2)-TM(I,1))/YP
      END IF
      SU(I,J)=SU(I,J)+(YPLUSA/TPLUS)*TM(I,1)*VISC*SXCV(I)/YP
      SP(I,J)=SP(I,J)-(YPLUSA/TPLUS)*VISC*SXCV(I)/YP
      AS(I,J)=0.0D0
  819 CONTINUE  
C-----WEST BOUNDARY(ENTRAINMENT)
      IF(NQWT.EQ.0) THEN 
      DO 816 J=1,NJ
      AW(2,J)=0.0D0
  816 CONTINUE
      END IF
C-----OUTLET
      IF(NQET.EQ.0) THEN
      DO 844 J=1,NJ
  844 AE(NIM1,J)=0.0D0
      END IF
      RETURN
      
      ENTRY OMB
      FLOWE=0.
      DO 263 J=2,NJM1
      DENU=DEN(2,J)+FX(2)*(DEN(3,J)-DEN(2,J))
      ARDEN=DENU*SYCV(J)
      FLOWE=FLOWE+ARDEN*U(2,J)
      U(1,J)=DENU*U(2,J)/DEN(1,J)
  263 CONTINUE   
      FLOWIN=0.0D0
      XMONIN=0.0D0
      XKMIN=0.0D0
      XEMIN=0.0D0
      XTMIN=0.0D0
      DO 904 I=IJ1P1,IJET2
      ARDEN=DEN(I,NJ)*SXCV(I)
      FLOWIN=FLOWIN+ARDEN*DABS(V(I,NJ))
      XMONIN=XMONIN+ARDEN*DABS(V(I,NJ))*VIN
      XKMIN=XKMIN+SXCV(I)*DABS(V(I,NJ))*TE(I,NJ)
      XEMIN=XEMIN+SXCV(I)*DABS(V(I,NJ))*ED(I,NJ)
      XTMIN=XTMIN+SXCV(I)*DABS(V(I,NJ))*TM(I,NJ)
  904 CONTINUE
C     DO 906 J=2,NJM1
C     ARDEN=DEN(1,J)*SYCV(J)
C     FLOWIN=FLOWIN+ARDEN*DABS(U(1,J))
C     XMONIN=XMONIN+ARDEN*U(1,J)**2
C     XKMIN=XKMIN+SYCV(J)*DABS(U(1,J)*TE(1,J))
C     XEMIN=XEMIN+SYCV(J)*DABS(U(1,J)*ED(1,J))
C     XTMIN=XTMIN+SYCV(J)*DABS(U(1,J)*TM(1,J))
C 906 CONTINUE
C     ARDENT=0.
c     FLOW=0.
c     DO 264 J=2,NJM1
c     DENU=DEN(NI,J)
c     ARDEN=DENU*SYCV(J)
c     ARDENT=ARDENT+ARDEN
c 264 FLOW=FLOW+ARDEN*U(NIM1,J)
c     UINC=(FLOWIN+FLOWE-FLOW)/ARDENT
      DO 205 J=2,NJM1
      U(NI,J)=U(NIM1,J)
  205 CONTINUE
      RETURN

      ENTRY MODTE
      DO 44 J=1,NJ
   44 AE(NIM1,J)=0.0D0
C-----TOP WALL
      J=NJM1
      DO 610 I=2,IJET1
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GEN(I,J)*VOL
      SP(I,J)=SP(I,J)+(CMU*DEN(I,J)**2*DABS(TE(I,J))/TVIS(I,J))*VOL
      YP=0.5*SYCV(J)
      SQRTK=SQRT(TE(I,J))
      GENCOU=0.5*(DABS(TAUN(I-1)*U(I-1,J))+DABS(TAUN(I)*U(I,J)))/YP
      YPLUSN(I)=DEN(I,J)*SQRTK*CMU25*YP/VISCO(I,J)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GEN(I,J)-TVIS(I,J)*(DVDX+DUDY)**2
      GEN(I,J)=GENCOU+GENRES
      IF (YPLUSN(I).LE.11.63) GO TO 611
      DITERM=DEN(I,J)*CMU75*SQRTK*DLOG(ELOG*YPLUSN(I))/(CAPPA*YP)
      GO TO 612
  611 DITERM=DEN(I,J)*CMU75*SQRTK*YPLUSN(I)/YP
  612 SU(I,J)=GEN(I,J)*VOL+SU(I,J)
      SP(I,J)=-DITERM*VOL+SP(I,J)
  610 AN(I,J)=0.0D0
      J=NJM1
      DO 623 I=IJ2P1,NI
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GEN(I,J)*VOL
      SP(I,J)=SP(I,J)+(CMU*DEN(I,J)**2*DABS(TE(I,J))/TVIS(I,J))*VOL
      YP=0.5*SYCV(J)
      SQRTK=SQRT(TE(I,J))
      GENCOU=0.5*(DABS(TAUN(I-1)*U(I-1,J))+DABS(TAUN(I)*U(I,J)))/YP
      YPLUSN(I)=DEN(I,J)*SQRTK*CMU25*YP/VISCO(I,J)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GEN(I,J)-TVIS(I,J)*(DVDX+DUDY)**2
      GEN(I,J)=GENCOU+GENRES
      IF (YPLUSN(I).LE.11.63) GO TO 624
      DITERM=DEN(I,J)*CMU75*SQRTK*DLOG(ELOG*YPLUSN(I))/(CAPPA*YP)
      GO TO 625
  624 DITERM=DEN(I,J)*CMU75*SQRTK*YPLUSN(I)/YP
  625 SU(I,J)=GEN(I,J)*VOL+SU(I,J)
      SP(I,J)=-DITERM*VOL+SP(I,J)
  623 AN(I,J)=0.0D0
C-----BOTTOM WALL
      J=2
      DO 613 I=2,NI
      VOL=RSYCV(J)*SXCV(I)
      SU(I,J)=SU(I,J)-GEN(I,J)*VOL
      SP(I,J)=SP(I,J)+(CMU*DEN(I,J)**2*DABS(TE(I,J))/TVIS(I,J))*VOL
      YP=0.5*SYCV(J)
      SQRTK=SQRT(TE(I,J))
      GENCOU=0.5*(DABS(TAUS(I-1)*U(I-1,J))+DABS(TAUS(I)*U(I,J)))/YP
      YPLUSS(I)=DEN(I,J)*SQRTK*CMU25*YP/VISCO(I,J)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENRES=GEN(I,J)-TVIS(I,J)*(DVDX+DUDY)**2
      GEN(I,J)=GENCOU+GENRES
      IF (YPLUSS(I).LE.11.63) GO TO 615
      DITERM=DEN(I,J)*CMU75*SQRTK*DLOG(ELOG*YPLUSS(I))/(CAPPA*YP)
      GO TO 614
  615 DITERM=DEN(I,J)*CMU75*SQRTK*YPLUSS(I)/YP
  614 SU(I,J)=GEN(I,J)*VOL+SU(I,J)
      SP(I,J)=-DITERM*VOL+SP(I,J)
  613 AS(I,J)=0.0D0
C-----WEST BOUNDARY(ENTRAINMENT)      
      DO 616 J=2,NJM1
  616 AW(2,J)=0.0D0
      RETURN
C
C
      ENTRY MODED
      DO 23 J=1,NJ
   23 AE(NIM1,J)=0.0D0
C------TOP WALL
      J=NJM1
      DO 710 I=2,IJET1
      YP=0.5*SYCV(J)
      TERM=CMU75/(CAPPA*YP)
      SU(I,J)=GREAT*TERM*(TE(I,J))**1.5
  710 SP(I,J)=-GREAT
      DO 711 I=IJ2P1,NIM1
      YP=0.5*SYCV(J)
      TERM=CMU75/(CAPPA*YP)
      SU(I,J)=GREAT*TERM*(TE(I,J))**1.5
  711 SP(I,J)=-GREAT
C-----BOTTOM WALL
      J=2
      DO 730 I=2,NIM1
      YP=0.5*SYCV(J)
      TERM=CMU75/(CAPPA*YP)
      SU(I,J)=GREAT*TERM*(TE(I,J))**1.5
      SP(I,J)=-GREAT
  730 CONTINUE
C-----WEST BOUNDARY(ENTRAINMENT)      
      DO 731 J=2,NJM1 
  731 AW(2,J)=0.0D0
      RETURN
      END
C***************************INITIALIZE********************************      
      SUBROUTINE INIT
      INCLUDE 'UCOMMON.INC'
C********** DAFAULT VALUES
      GREAT=1.D+30
      TINY=1.D-30
C********** CALCULATE GEOMETRICAL QUANTITIES
      NI=NIM1+1
      NJ=NJM1+1
      NIM2=NIM1-1
      NJM2=NJM1-1
      X(1)=XU(1)
      X(NI)=XU(NIM1)
      Y(1)=YV(1)
      Y(NJ)=YV(NJM1)
      DO 100 I=2,NIM1
  100 X(I)=0.5*(XU(I-1)+XU(I))
      DO 150 J=2,NJM1
  150 Y(J)=0.5*(YV(J-1)+YV(J))
      SXCV(1)=TINY
      SXCV(NI)=TINY
      SYCV(1)=TINY
      SYCV(NJ)=TINY
      DO 200 I=2,NIM1
  200 SXCV(I)=XU(I)-XU(I-1)
      DO 250 J=2,NJM1
  250 SYCV(J)=YV(J)-YV(J-1)
      DO 300 I=1,NIM1
  300 UXCV(I)=X(I+1)-X(I)
      DO 350 J=1,NJM1
  350 VYCV(J)=Y(J+1)-Y(J)
      UCVI(1)=0.0D0
      UCVIP(NIM1)=0.0D0
      VCVJ(1)=0.0D0
      VCVJP(NJM1)=0.0D0
      DO 400 I=2,NIM1
      UCVI(I)=0.5*SXCV(I)
  400 UCVIP(I-1)=UCVI(I)
      DO 450 J=2,NJM1
      VCVJ(J)=0.5*SYCV(J)
  450 VCVJP(J-1)=VCVJ(J)
      DO 500 I=1,NIM1
  500 FX(I)=UCVI(I)/UXCV(I)
      DO 550 J=1,NJM1
  550 FY(J)=VCVJ(J)/VYCV(J)
      IF(IND.EQ.2)GO TO 600
      DO 650 J=1,NJ
      R(J)=1.D0
  650 RV(J)=1.D0
      GO TO 700
  600 DO 750 J=2,NJ
  750 R(J)=R(J-1)+VYCV(J-1)
      RV(1)=R(1)
      DO 800 J=2,NJM1
  800 RV(J)=RV(J-1)+SYCV(J)
  700 DO 850 J=1,NJ
  850 RSYCV(J)=R(J)*SYCV(J)
      DO 900 J=1,NJM1
  900 RVYCV(J)=0.5*(R(J+1)+R(J))*VYCV(J)
      DO 925 I=1,NIM1
  925 B3XS(I)=UCVI(I)*UCVIP(I)*UXCV(I)
      DO 926 I=1,NI
  926 B3XU(I)=SXCV(I)**3*0.25
      DO 927 J=1,NJM1
  927 B3YS(J)=VCVJ(J)*VCVJP(J)*VYCV(J)
      DO 928 J=1,NJ
  928 B3YV(J)=SYCV(J)**3*0.25
C********** SET VARIABLES TO ZERO
      DO 950 I=1,NI
      DO 950 J=1,NJ
      U(I,J)=0.0D0
      V(I,J)=0.0D0
      P(I,J)=0.0D0
      TM(I,J)=(TSTGE+TSTGJ)/2
      PP(I,J)=0.0D0
      DU(I,J)=0.0D0
      DV(I,J)=0.0D0
      SU(I,J)=0.0D0
      SP(I,J)=0.0D0
      TE(I,J)=0.D0
      ED(I,J)=0.D0
      GEN(I,J)=0.0D0
      GENT(I,J)=0.0D0
      AE(I,J)=0.0D0
      AW(I,J)=0.0D0
      AN(I,J)=0.0D0
      AS(I,J)=0.0D0
      SMPW(J)=0.0D0
  950 CONTINUE
      RETURN
      END
C***************************U-MOMENTUM*********************************
      SUBROUTINE CALCU
      INCLUDE 'UCOMMON.INC'
      LOGICAL LX,LY 
      NQUICK=NQUV
C********** ASSEMBLY OF COEFFICIENTS
      DO 100 I=2,NIM1
      DO 100 J=2,NJ
C********** AREAS,VOLUME,DELX & DELY
      AREAW=RSYCV(J)
      AREAS=UXCV(I)*RV(J-1)
      DELX1=UCVIP(I-1)
      DELX2=UCVI(I)
      DELY1=VCVJ(J-1)
      DELY2=VCVJP(J-1)
      B3X=B3XU(I)
      B3Y=B3YS(J-1)
      LX=I.GT.2.AND.I.LT.NIM1
      LY=J.GT.2.AND.J.LT.NJ
      VOL=RSYCV(J)*UXCV(I)
      SP(I,J)=0.0D0
C********** CALCULATE CONVECTION
      DENWE=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      DENWW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      CW=0.5*(DENWE*U(I,J)+DENWW*U(I-1,J))*AREAW
      DENSE=DEN(I+1,J-1)+FY(J-1)*(DEN(I+1,J)-DEN(I+1,J-1))
      DENSW=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      CS=0.5*(DENSE*V(I+1,J-1)+DENSW*V(I,J-1))*AREAS
C********** CALCULATE DIFFUSION
      DW=VIS(I,J)*AREAW/SXCV(I)
      GAMSE=VIS(I+1,J-1)+FY(J-1)*(VIS(I+1,J)-VIS(I+1,J-1))
      GAMSW=VIS(I,J-1)+FY(J-1)*(VIS(I,J)-VIS(I,J-1))
      GAMS=0.5*(GAMSE+GAMSW)
      DS=GAMS*AREAS/VYCV(J-1)
C********** SCHEMES
      NSUP=-1
      IF(CW.GT.0.0D0)NSUP=1
      UP=0.5*FLOAT(1+NSUP)
      UM=0.5*FLOAT(1-NSUP)
      NSVP=-1
      IF(CS.GT.0.0D0)NSVP=1
      VP=0.5*FLOAT(1+NSVP)
      VM=0.5*FLOAT(1-NSVP)
      IF(LX.AND.NQUICK.EQ.1)GO TO 111
      TEMPW=DW-DABS(CW)*0.1
      COEW=0.0D0
      TEMPW=TEMPW/(DW+TINY)
      IF(TEMPW.GT.TINY)COEW=DW*TEMPW**5
C     IF(.NOT.LX.AND.NQUICK.EQ.1)COEW=DW-0.5*DABS(CW)
C     COEW=DW
      SPPX=0.0D0
      SUPX=0.0D0
      AE(I,J)=0.0D0
      GO TO 300
  111 KXP=I-NSUP
      LXP=KXP-(1+NSUP)/2
      DELX3=UM*DELX2-UP*DELX1-SXCV(KXP)*NSUP
      B1X=DELX2*DELX3*(DELX2-DELX3)
      B2X=DELX1*DELX3*(DELX3+DELX1)
      BX=B1X-B2X+B3X
      TEMP=B3X/BX*CW
      COEW=DW+(UM*B1X+UP*B2X)/BX*CW
      CW=CW-TEMP
      TEMPP=TEMP*UP
      TEMPM=TEMP*UM
      SPPX=-TEMPP
      SUPX=TEMPP*U(LXP,J)
      AW(I-1,J)=AW(I-1,J)-TEMPP
      AE(I,J)=TEMPM
      SP(I-1,J)=SP(I-1,J)+TEMPM
      SU(I-1,J)=SU(I-1,J)-TEMPM*U(LXP,J)
  300 CONTINUE
      IF(LY.AND.NQUICK.EQ.1)GO TO 222
      TEMPS=DS-DABS(CS)*0.1
      COES=0.0D0
      TEMPS=TEMPS/(DS+TINY)
      IF(TEMPS.GT.TINY)COES=DS*TEMPS**5
C     IF(.NOT.LY)COES=DS-DMAX1(FY(J-1)*CS,-(1.-FY(J-1))*CS)
C     COES=DS
      SPPY=0.0D0
      SUPY=0.0D0
      AN(I,J)=0.0D0
      GO TO 400
  222 KYP=J-1-NSVP
      LYP=KYP+(1-NSVP)/2
      DELY3=VM*DELY2-VP*DELY1-VYCV(KYP)*NSVP
      B1Y=DELY2*DELY3*(DELY2-DELY3)
      B2Y=DELY1*DELY3*(DELY1+DELY3)
      BY=B1Y-B2Y+B3Y
      TEMP=B3Y/BY*CS
      COES=DS+(VM*B1Y+VP*B2Y)/BY*CS
      CS=CS-TEMP
      TEMPP=TEMP*VP
      TEMPM=TEMP*VM
      SPPY=-TEMPP
      SUPY=TEMPP*U(I,LYP)
      AS(I,J-1)=AS(I,J-1)-TEMPP
      AN(I,J)=TEMPM
      SP(I,J-1)=SP(I,J-1)+TEMPM
      SU(I,J-1)=SU(I,J-1)-TEMPM*U(I,LYP)
  400 CONTINUE
C********** ASSEMBLE COEFFICIENTS FOR E,W,N & S
      AW(I,J)=COEW+CW*UP
      AE(I-1,J)=AW(I,J)-CW+AE(I-1,J)
      AS(I,J)=COES+CS*VP
      AN(I,J-1)=AS(I,J)-CS+AN(I,J-1)
C********** CALCULATE SOURCES
      DUDXE=(U(I+1,J)-U(I,J))/SXCV(I+1)
      DUDXW=(U(I,J)-U(I-1,J))/SXCV(I)
      DVDYE=R(J)*(V(I+1,J)-V(I+1,J-1))/SYCV(J)
      DVDYW=R(J)*(V(I,J)-V(I,J-1))/SYCV(J)
      IF(J.NE.NJ)DVDXN=RV(J)*(V(I+1,J)-V(I,J))/UXCV(I)
      DVDXS=RV(J-1)*(V(I+1,J-1)-V(I,J-1))/UXCV(I)
      IF(J.NE.NJ)
     +GAMN=0.5*(VIS(I,J)+VIS(I+1,J)+FY(J)*(VIS(I,J+1)+VIS(I+1,J+1)
     +-VIS(I,J)-VIS(I+1,J)))
      SUT=(2.*DEN(I,J)*TE(I,J)/3.-2.*DEN(I+1,J)*TE(I+1,J)/3.)/UXCV(I)
      SUC=-2.*(VIS(I+1,J)*(DUDXE+DVDYE)-VIS(I,J)*(DUDXW+DVDYW))
     +/(3.*UXCV(I)) 
      SU(I,J)=(DUDXE*VIS(I+1,J)-DUDXW*VIS(I,J))/UXCV(I)
     ++(DVDXN*GAMN-DVDXS*GAMS)/RSYCV(J)+(P(I,J)-P(I+1,J))/UXCV(I)
     ++SUT*NTF+SUC*NCF
      SU(I,J)=SU(I,J)*VOL+SUPX+SUPY
      SP(I,J)=SP(I,J)*VOL+SPPX+SPPY
      DENOE=DENO(I,J)+FX(I)*(DENO(I+1,J)-DENO(I,J))
      APO(I,J)=DENOE*VOL/DT
      SU(I,J)=SU(I,J)+APO(I,J)*UO(I,J)
      SP(I,J)=SP(I,J)-APO(I,J)
  100 CONTINUE
C********** PROBLEM MODIFICATION
      CALL MODU
C********** FINAL COEFF. ASSEMBLY AND RESIDUAL CALCULATION
      RESORU=0.0D0
      DO 200 I=2,NIM1
      DO 200 J=2,NJM1
      AP(I,J)=AW(I,J)+AE(I,J)+AS(I,J)+AN(I,J)-SP(I,J)
      DU(I,J)=RSYCV(J)/AP(I,J)
      RESOR=AN(I,J)*U(I,J+1)+AS(I,J)*U(I,J-1)+AE(I,J)*U(I+1,J)
     ++AW(I,J)*U(I-1,J)-AP(I,J)*U(I,J)+SU(I,J)
      IF(-SP(I,J).GE.GREAT)RESOR=0.0D0
      RESORU=RESORU+DABS(RESOR)
  200 CONTINUE
C********** SOLUTION OF DIFFERENCE EQUATIONS
      CALL LISOLV(2,2,NIM1,NJM1,U,NSWPU,URFU)
      RETURN
      END
C************************V-MOMENTUM***********************************      
      SUBROUTINE CALCV
      INCLUDE 'UCOMMON.INC'
      LOGICAL LX,LY 
C********** ASSEMBLY OF COEFFICIENTS
      DO 100 I=2,NI
      DO 100 J=2,NJM1
C********** AREAS,VOLUME,DELX & DELY
      AREAS=SXCV(I)*R(J)
      AREAW=RVYCV(J)
      DELX1=UCVI(I-1)
      DELX2=UCVIP(I-1)
      DELY1=VCVJP(J-1)
      DELY2=VCVJ(J)
      B3X=B3XS(I-1)
      B3Y=B3YV(J)
      LX=I.GT.2.AND.I.LT.NI
      LY=J.GT.2.AND.J.LT.NJM1
      VOL=RVYCV(J)*SXCV(I)
      SP(I,J)=0.0D0
C********** CALCULATE CONVECTION
      DENWN=DEN(I-1,J+1)+FX(I-1)*(DEN(I,J+1)-DEN(I-1,J+1))
      DENWS=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      CW=0.5*(DENWN*U(I-1,J+1)*R(J+1)+DENWS*U(I-1,J)*R(J))*VYCV(J)
      DENSN=DEN(I,J)+FY(J)*(DEN(I,J+1)-DEN(I,J))
      DENSS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      CS=0.5*(DENSN*V(I,J)*RV(J)+DENSS*V(I,J-1)*RV(J-1))*SXCV(I)
C********** CALCULATE DIFFUSION
      GAMWN=VIS(I-1,J+1)+FX(I-1)*(VIS(I,J+1)-VIS(I-1,J+1))
      GAMWS=VIS(I-1,J)+FX(I-1)*(VIS(I,J)-VIS(I-1,J))
      GAMW=0.5*(GAMWN+GAMWS)
      DW=GAMW*AREAW/UXCV(I-1)
      DS=VIS(I,J)*AREAS/SYCV(J)
C********** SCHEMES
      NSUP=-1
      IF(CW.GT.0.0D0)NSUP=1
      UP=0.5*FLOAT(1+NSUP)
      UM=0.5*FLOAT(1-NSUP)
      NSVP=-1
      IF(CS.GT.0.0D0)NSVP=1
      VP=0.5*FLOAT(1+NSVP)
      VM=0.5*FLOAT(1-NSVP)
      IF(LX.AND.NQUICK.EQ.1)GO TO 111
      TEMPW=DW-DABS(CW)*0.1
      COEW=0.0D0
      TEMPW=TEMPW/(DW+TINY)
      IF(TEMPW.GT.TINY)COEW=DW*TEMPW**5
C     IF(.NOT.LX)COEW=DW-DMAX1(FX(I-1)*CW,-(1.D0-FX(I-1))*CW)
C     COEW=DW
      SPPX=0.0D0
      SUPX=0.0D0
      AE(I,J)=0.0D0
      GO TO 300
  111 KXP=I-1-NSUP
      LXP=KXP+(1-NSUP)/2
      DELX3=UM*DELX2-UP*DELX1-UXCV(KXP)*NSUP
      B1X=DELX2*DELX3*(DELX2-DELX3)
      B2X=DELX1*DELX3*(DELX1+DELX3)
      BX=B1X-B2X+B3X
      TEMP=B3X/BX*CW
      COEW=DW+(UM*B1X+UP*B2X)/BX*CW
      CW=CW-TEMP
      TEMPP=TEMP*UP
      TEMPM=TEMP*UM
      SPPX=-TEMPP
      SUPX=TEMPP*V(LXP,J)
      AW(I-1,J)=AW(I-1,J)-TEMPP
      AE(I,J)=TEMPM
      SP(I-1,J)=SP(I-1,J)+TEMPM
      SU(I-1,J)=SU(I-1,J)-TEMPM*V(LXP,J)
  300 CONTINUE
      IF(LY.AND.NQUICK.EQ.1)GO TO 222
      TEMPS=DS-DABS(CS)*0.1
      COES=0.0D0
      TEMPS=TEMPS/(DS+TINY)
      IF(TEMPS.GT.TINY)COES=DS*TEMPS**5
C     IF(.NOT.LY.AND.NQUICK.EQ.1)COES=DS-0.5*DABS(CS)
C     COES=DS
      SPPY=0.0D0
      SUPY=0.0D0
      AN(I,J)=0.0D0
      GO TO 400
  222 KYP=J-NSVP
      LYP=KYP-(1+NSVP)/2
      DELY3=VM*DELY2-VP*DELY1-SYCV(KYP)*NSVP
      B1Y=DELY2*DELY3*(DELY2-DELY3)
      B2Y=DELY1*DELY3*(DELY1+DELY3)
      BY=B1Y-B2Y+B3Y
      TEMP=B3Y/BY*CS
      COES=DS+(VM*B1Y+VP*B2Y)/BY*CS
      CS=CS-TEMP
      TEMPP=TEMP*VP
      TEMPM=TEMP*VM
      SPPY=-TEMPP
      SUPY=TEMPP*V(I,LYP)
      AS(I,J-1)=AS(I,J-1)-TEMPP
      AN(I,J)=TEMPM
      SP(I,J-1)=SP(I,J-1)+TEMPM
      SU(I,J-1)=SU(I,J-1)-TEMPM*V(I,LYP)
  400 CONTINUE
C********** ASSEMBLE COEFFICIENTS FOR E,W,N & S
      AW(I,J)=COEW+CW*UP
      AE(I-1,J)=AW(I,J)-CW+AE(I-1,J)
      AS(I,J)=COES+CS*VP
      AN(I,J-1)=AS(I,J)-CS+AN(I,J-1)
C********** CALCULATE SOURCES
      DUDYE=(U(I,J+1)-U(I,J))/VYCV(J)
      DUDYW=(U(I-1,J+1)-U(I-1,J))/VYCV(J)
      DUDXN=(U(I,J+1)-U(I-1,J+1))/SXCV(I)
      DUDXS=(U(I,J)-U(I-1,J))/SXCV(I)
      DVDYN=R(J+1)*(V(I,J+1)-V(I,J))/SYCV(J+1)
      DVDYS=R(J)*(V(I,J)-V(I,J-1))/SYCV(J)
      IF(I.NE.NI)
     +GAME=0.5*(VIS(I,J+1)+VIS(I,J)+FX(I)*(VIS(I+1,J+1)+VIS(I+1,J)
     +-VIS(I,J+1)-VIS(I,J)))
      SUT=(2.*DEN(I,J)*TE(I,J)/3.-2.*DEN(I,J+1)*TE(I,J+1)/3.)/VYCV(J)
      SUC=-2.*(VIS(I,J+1)*(DUDXN+DVDYN)-VIS(I,J)*(DUDXS+DVDYS))
     +/(3.*VYCV(J)) 
      SU(I,J)=(GAME*DUDYE-GAMW*DUDYW)/SXCV(I)
     ++(DVDYN*VIS(I,J+1)-DVDYS*VIS(I,J))/RVYCV(J)
     ++(P(I,J)-P(I,J+1))/VYCV(J)
     ++SUT*NTF+SUC*NCF 
      IF(IND.EQ.2)SP(I,J)=-(VIS(I,J)+VIS(I,J+1))/((R(J)+R(J+1))*0.5)**2
      SU(I,J)=SU(I,J)*VOL+SUPX+SUPY
      SP(I,J)=SP(I,J)*VOL+SPPX+SPPY
      DENON=DENO(I,J)+FY(J)*(DENO(I,J+1)-DENO(I,J))
      APO(I,J)=DENON*VOL/DT
      SU(I,J)=SU(I,J)+APO(I,J)*VO(I,J)
      SP(I,J)=SP(I,J)-APO(I,J)
  100 CONTINUE
C********** PROBLEM MODIFICATION
      CALL MODV
C********** FINAL COEFF. ASSEMBLY AND RESIDUAL CALCULATION
      RESORV=0.0D0
      DO 200 I=3,NIM1
      DO 200 J=2,NJM2
      AP(I,J)=AW(I,J)+AE(I,J)+AS(I,J)+AN(I,J)-SP(I,J)
      DV(I,J)=0.5*(R(J)+R(J+1))*SXCV(I)/AP(I,J)
      RESOR=AN(I,J)*V(I,J+1)+AS(I,J)*V(I,J-1)+AE(I,J)*V(I+1,J)
     ++AW(I,J)*V(I-1,J)-AP(I,J)*V(I,J)+SU(I,J)
      IF(-SP(I,J).GE.GREAT)RESOR=0.0D0
      RESORV=RESORV+DABS(RESOR)
  200 CONTINUE
C********* SOLUTION OF DIFFERENCE EQUATIONS
      CALL LISOLV(3,2,NIM1,NJM2,V,NSWPV,URFV)
      RETURN
      END
C************************PRESSURE FIELD*******************************
      SUBROUTINE CALCP
      INCLUDE 'UCOMMON.INC'
      DO 100 I=2,NI
      DO 100 J=2,NJ
C********** AREAS
      AREAW=RSYCV(J)
      AREAS=SXCV(I)*RV(J-1)
C********** CALCULATE COEFFICIENTS
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      UKW=1.0/(SG*GC*TM(I-1,J))
      UKS=1.0/(SG*GC*TM(I,J-1))
      UKP=1.0/(SG*GC*TM(I,J))
      IF(U(I-1,J).GE.0.0D0) THEN
      KW=1
      CWK=U(I-1,J)*UKW*AREAW
      ELSE
      KW=0
      CWK=U(I-1,J)*UKP*AREAW
      END IF
      IF(V(I,J-1).GE.0.0D0) THEN
      KS=1
      CSK=V(I,J-1)*UKS*AREAS
      ELSE
      KS=0
      CSK=V(I,J-1)*UKP*AREAS
      END IF
      AW(I,J)=DENW*AREAW*DU(I-1,J)+CWK*KW*NCSA
      AE(I-1,J)=AW(I,J)-CWK*NCSA
      AS(I,J)=DENS*AREAS*DV(I,J-1)+CSK*KS*NCSA
      AN(I,J-1)=AS(I,J)-CSK*NCSA
C********** CALCULATE SOURCES
      CW=DENW*U(I-1,J)*AREAW
      CS=DENS*V(I,J-1)*AREAS
      SU(I,J)=CW+CS
      IF(J.NE.2)SU(I,J-1)=SU(I,J-1)-CS
      IF(I.NE.2)SU(I-1,J)=SU(I-1,J)-CW
      SP(I,J)=0.0D0
  100 CONTINUE
C********** PROBLEM MODIFICATION
      IF(NCSA.EQ.1) CALL MODP
C********** FINAL COEFFICIENT ASSEMBLY AND RESIDUAL SOURCE CALCULATION
      RESORM=0.0D0
      DO 200 I=3,NIM1
      DO 200 J=2,NJM1
      SU(I,J)=SU(I,J)+(DENO(I,J)-DEN(I,J))*SXCV(I)*SYCV(J)/DT
      AP(I,J)=AW(I,J)+AE(I,J)+AS(I,J)+AN(I,J)-SP(I,J)
  200 RESORM=RESORM+DABS(SU(I,J))
C********** SOLUTION OF DIFFERENCE EQUATION
      CALL LISOLV(3,2,NIM1,NJM1,PP,NSWPP,1.0D0)
C********** CORRECT VELOCITY & PRESSURE
C     DO 205 J=2,NJM1
C     PP(NI,J)=PP(NIM1,J)
C     PP(1,J)=PP(2,J)
C 205 CONTINUE
      DO 201 I=1,NI
      PP(I,1)=PP(I,2)
      PP(I,NJ)=PP(I,NJM1)
  201 CONTINUE
C**********VELOCITIES
      DO 400 I=2,NIM1
      DO 400 J=2,NJM1
      U(I,J)=U(I,J)+DU(I,J)*(PP(I,J)-PP(I+1,J))
      V(I,J)=V(I,J)+DV(I,J)*(PP(I,J)-PP(I,J+1))
  400 CONTINUE
C********** PRESSURE
      IF(NCSA.EQ.1) THEN
      DO 902 I=1,NI
      DO 902 J=1,NJ
      UKP=1.0/(SG*GC*TM(I,J))
      DENNEW=DEN(I,J)+UKP*PP(I,J)
      DEN(I,J)=DEN(I,J)+URFDEN*(DENNEW-DEN(I,J))
  902 CONTINUE
      END IF
      IF(NCF.EQ.0) THEN
      PPREF=PP(IPREF,JPREF)
      ELSE
      PPREF=0.0D0
      END IF
      DO 500 I=1,NI
      DO 500 J=1,NJ
      P(I,J)=P(I,J)+URFP*(PP(I,J)-PPREF)
      PP(I,J)=0.0D0
c	print*,'i:',i
  500 CONTINUE
	i=i-1
c	print*,'i:',i
c     pause
C********** ENTRAINMENT ASSOCIATED WITH CONSTANT PRESSURE B C
      FLOWE=0.0D0
      ARDENT=0.0D0
      DO 263 J=2,NJM1
      DENU=DEN(2,J)+FX(2)*(DEN(3,J)-DEN(2,J))
      ARDEN=DENU*SYCV(J)
      ARDENT=ARDENT+ARDEN
      FLOWE=FLOWE+ARDEN*U(2,J)
      U(1,J)=DENU*U(2,J)/DEN(1,J)
  263 CONTINUE   
      UAV=FLOWE/ARDENT
      TMAE=TSTGE
      IF(UAV.LT.TINY) THEN  
      PAE=0.0D0
      ELSE
      TMAE=TSTGE-(SG-1.)*(UAV**2)/(2.*SG*GC)
      PAE=PSTGE*(TMAE/TSTGE)**(SG/(SG-1.0D0))-PATM
      END IF
      DO 204 J=2,NJM1
      TM(1,J)=TMO(I,J)
C     TM(1,J)=TM(1,J)+0.1*(TMAE-TM(1,J))
      P(2,J)=PO(2,J)
C     P(2,J)=P(2,J)+0.1*(PAE-P(2,J))
      P(1,J)=P(2,J)
  204 CONTINUE
      DO 203 I=1,NI
      P(I,1)=P(I,2)
      P(I,NJ)=P(I,NJM1)
  203 CONTINUE
      RETURN
      END
C*************************TEMPERATURE FIELD****************************
      SUBROUTINE CALCE
      INCLUDE 'UCOMMON.INC'
      LOGICAL LX,LY 
      NQUICK=NQKE
C********** ASSEMBLY OF COEFFICIENTS
      DO 100 I=2,NI
      DO 100 J=2,NJ
C********** AREAS,VOLUME,DELX,DELY & TVIS
      AREAW=RSYCV(J)
      AREAS=SXCV(I)*RV(J-1)
      VOL=RSYCV(J)*SXCV(I)
      DELX1=UCVI(I-1)
      DELX2=UCVIP(I-1)
      DELY1=VCVJ(J-1)
      DELY2=VCVJP(J-1)
      B3X=B3XS(I-1)
      B3Y=B3YS(J-1)
      LX=I.GT.2.AND.I.LT.NI
      LY=J.GT.2.AND.J.LT.NJ
C********** CALCULATE CONVECTION
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      CW=DENW*U(I-1,J)*AREAW
      CS=DENS*V(I,J-1)*AREAS
C********** CALCULATE DIFFUSION
      GAMW=GAT(I-1,J)+FX(I-1)*(GAT(I,J)-GAT(I-1,J))
      GAMS=GAT(I,J-1)+FY(J-1)*(GAT(I,J)-GAT(I,J-1))
      DW=GAMW*AREAW/UXCV(I-1)
      DS=GAMS*AREAS/VYCV(J-1)
C********** SCHEMES
      NSUP=-1
      IF(CW.GT.0.0D0)NSUP=1
      UP=0.5*FLOAT(1+NSUP)
      UM=0.5*FLOAT(1-NSUP)
      NSVP=-1
      IF(CS.GT.0.0D0)NSVP=1
      VP=0.5*FLOAT(1+NSVP)
      VM=0.5*FLOAT(1-NSVP)
      IF(LX.AND.NQUICK.EQ.1)GO TO 111
      TEMPW=DW-DABS(CW)*0.1
      COEW=0.0D0
      TEMPW=TEMPW/(DW+TINY)
      IF(TEMPW.GT.TINY)COEW=DW*TEMPW**5
C     IF(.NOT.LX)COEW=DW-DMAX1(FX(I-1)*CW,-(1.-FX(I-1))*CW)
C     COEW=DW
      SPPX=0.0D0
      SUPX=0.0D0
      AE(I,J)=0.0D0
      GO TO 300
  111 KXP=I-1-NSUP
      LXP=KXP+(1-NSUP)/2
      DELX3=UM*DELX2-UP*DELX1-UXCV(KXP)*NSUP
      B1X=DELX2*DELX3*(DELX2-DELX3)
      B2X=DELX1*DELX3*(DELX1+DELX3)
      BX=B1X-B2X+B3X
      TEMP=B3X/BX*CW
      COEW=DW+(UM*B1X+UP*B2X)/BX*CW
      CW=CW-TEMP
      TEMPP=TEMP*UP
      TEMPM=TEMP*UM
      SPPX=-TEMPP
      SUPX=TEMPP*TM(LXP,J)
      AW(I-1,J)=AW(I-1,J)-TEMPP
      AE(I,J)=TEMPM
      SP(I-1,J)=SP(I-1,J)+TEMPM
      SU(I-1,J)=SU(I-1,J)-TEMPM*TM(LXP,J)
  300 CONTINUE
      IF(LY.AND.NQUICK.EQ.1)GO TO 222
      TEMPS=DS-DABS(CS)*0.1
      COES=0.0D0
      TEMPS=TEMPS/(DS+TINY)
      IF(TEMPS.GT.TINY)COES=DS*TEMPS**5
C     IF(.NOT.LY)COES=DS-DMAX1(FY(J-1)*CS,-(1.-FY(J-1))*CS)
C     COES=DS
      SPPY=0.0D0
      SUPY=0.0D0
      AN(I,J)=0.0D0
      GO TO 400
  222 KYP=J-1-NSVP
      LYP=KYP+(1-NSVP)/2
      DELY3=VM*DELY2-VP*DELY1-VYCV(KYP)*NSVP
      B1Y=DELY2*DELY3*(DELY2-DELY3)
      B2Y=DELY1*DELY3*(DELY1+DELY3)
      BY=B1Y-B2Y+B3Y
      TEMP=B3Y/BY*CS
      COES=DS+(VM*B1Y+VP*B2Y)/BY*CS
      CS=CS-TEMP
      TEMPP=TEMP*VP
      TEMPM=TEMP*VM
      SPPY=-TEMPP
      SUPY=TEMPP*TM(I,LYP)
      AS(I,J-1)=AS(I,J-1)-TEMPP
      AN(I,J)=TEMPM
      SP(I,J-1)=SP(I,J-1)+TEMPM
      SU(I,J-1)=SU(I,J-1)-TEMPM*TM(I,LYP)
  400 CONTINUE
C********** ASSEMBLE COEFFICIENTS FOR E,W,N & S
      AW(I,J)=COEW+CW*UP
      AE(I-1,J)=AW(I,J)-CW+AE(I-1,J)
      AS(I,J)=COES+CS*VP
      AN(I,J-1)=AS(I,J)-CS+AN(I,J-1)
C********** CALCULATE SOURCES
      IF(I.EQ.NI.OR.J.EQ.NJ) GOTO 556
      PE=P(I,J)+FX(I)*(P(I+1,J)-P(I,J))
      PW=P(I-1,J)+FX(I-1)*(P(I,J)-P(I-1,J))
      PN=P(I,J)+FY(J)*(P(I,J+1)-P(I,J))
      PS=P(I,J-1)+FY(J-1)*(P(I,J)-P(I,J-1))
      DPDX=(PE-PW)/SXCV(I)
      DPDY=(PN-PS)/SYCV(J)
      DPDT=(P(I,J)-PO(I,J))/DT
      DUDX=(U(I,J)-U(I-1,J))/SXCV(I)
      DVDY=(V(I,J)-V(I,J-1))/SYCV(J)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GENT(I,J)=VIS(I,J)*(2.*(DUDX**2+DVDY**2)+(DUDY+DVDX)**2)
      GT1=DPDX*(U(I,J)+U(I-1,J))/2.+DPDY*(V(I,J)+V(I,J-1))/2.+DPDT 
      IF(GT1.GT.0.0D0) THEN
      K1=1 
      K2=0
      ELSE 
      K1=0
      K2=1
      END IF
      IF(IND.EQ.2)GENT(I,J)=GENT(I,J)+2.*(0.5*(V(I,J)+V(I,J-1))/R(J))**2
      GE1=-2.*VIS(I,J)*(DUDX+DVDY)**2/(3.*TM(I,J))
      SU(I,J)=GENT(I,J)/SPHP+GT1*K1/SPHP
      SP(I,J)=GE1*NCF/SPHP+GT1*K2/(SPHP*TM(I,J))
  556 SU(I,J)=SU(I,J)*VOL+SUPX+SUPY
      SP(I,J)=SP(I,J)*VOL+SPPX+SPPY
      APO(I,J)=DEN(I,J)*VOL/DT
      SU(I,J)=SU(I,J)+APO(I,J)*TMO(I,J)
      SP(I,J)=SP(I,J)-APO(I,J)
  100 CONTINUE
C********** PROBLEM MODIFICATION
      CALL MODE
C***** FINAL COEFFICIENT ASSEMBLY AND RESIDUAL SOURCE CALCULATION
      RESORT=0.0D0
      DO 200 I=2,NIM1
      DO 200 J=2,NJM1
      SU(I,J)=SUTM(I,J)+URFGTM*(SU(I,J)-SUTM(I,J))
      AP(I,J)=AN(I,J)+AS(I,J)+AE(I,J)+AW(I,J)-SP(I,J)
      RESOR=AN(I,J)*TM(I,J+1)+AS(I,J)*TM(I,J-1)+AE(I,J)*TM(I+1,J)
     ++AW(I,J)*TM(I-1,J)-AP(I,J)*TM(I,J)+SU(I,J)
      IF(-SP(I,J).GE.GREAT)RESOR=0.0D0
      RESORT=RESORT+DABS(RESOR)
C********** UNDER-RELAXATION
      SUTM(I,J)=SU(I,J)
  200 CONTINUE
C********* SOLUTION OF DIFFERENCE EQUATIONS
      CALL LISOLV(2,2,NIM1,NJM1,TM,NSWPT,URFTM)
      CALL CORE
      RETURN
      END
C***********************TURB. ENERGY FIELD****************************
      SUBROUTINE CALCTE
      INCLUDE 'UCOMMON.INC'
      LOGICAL LX,LY 
      NQUICK=NQKE
C********** ASSEMBLY OF COEFFICIENTS
      DO 100 I=2,NI
      DO 100 J=2,NJ
C********** AREAS,VOLUME,DELX,DELY & TVIS
      AREAW=RSYCV(J)
      AREAS=SXCV(I)*RV(J-1)
      VOL=RSYCV(J)*SXCV(I)
      DELX1=UCVI(I-1)
      DELX2=UCVIP(I-1)
      DELY1=VCVJ(J-1)
      DELY2=VCVJP(J-1)
      B3X=B3XS(I-1)
      B3Y=B3YS(J-1)
      LX=I.GT.2.AND.I.LT.NI
      LY=J.GT.2.AND.J.LT.NJ
C********** CALCULATE CONVECTION
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      CW=DENW*U(I-1,J)*AREAW
      CS=DENS*V(I,J-1)*AREAS
C********** CALCULATE DIFFUSION
      GMW=VISCO(I-1,J)+FX(I-1)*(VISCO(I,J)-VISCO(I-1,J))
      GMS=VISCO(I,J-1)+FY(J-1)*(VISCO(I,J)-VISCO(I,J-1))
      GAMW=TVIS(I-1,J)+FX(I-1)*(TVIS(I,J)-TVIS(I-1,J))
      GAMS=TVIS(I,J-1)+FY(J-1)*(TVIS(I,J)-TVIS(I,J-1))
      DW=(GAMW/PRTE+GMW)*AREAW/UXCV(I-1)
      DS=(GAMS/PRTE+GMS)*AREAS/VYCV(J-1)
C********** SCHEMES
      NSUP=-1
      IF(CW.GT.0.0D0)NSUP=1
      UP=0.5*FLOAT(1+NSUP)
      UM=0.5*FLOAT(1-NSUP)
      NSVP=-1
      IF(CS.GT.0.0D0)NSVP=1
      VP=0.5*FLOAT(1+NSVP)
      VM=0.5*FLOAT(1-NSVP)
      IF(LX.AND.NQUICK.EQ.1)GO TO 111
      TEMPW=DW-DABS(CW)*0.1
      COEW=0.0D0
      TEMPW=TEMPW/(DW+TINY)
      IF(TEMPW.GT.TINY)COEW=DW*TEMPW**5
C     IF(.NOT.LX)COEW=DW-DMAX1(FX(I-1)*CW,-(1.-FX(I-1))*CW)
C     COEW=DW
      SPPX=0.0D0
      SUPX=0.0D0
      AE(I,J)=0.0D0
      GO TO 300
  111 KXP=I-1-NSUP
      LXP=KXP+(1-NSUP)/2
      DELX3=UM*DELX2-UP*DELX1-UXCV(KXP)*NSUP
      B1X=DELX2*DELX3*(DELX2-DELX3)
      B2X=DELX1*DELX3*(DELX1+DELX3)
      BX=B1X-B2X+B3X
      TEMP=B3X/BX*CW
      COEW=DW+(UM*B1X+UP*B2X)/BX*CW
      CW=CW-TEMP
      TEMPP=TEMP*UP
      TEMPM=TEMP*UM
      SPPX=-TEMPP
      SUPX=TEMPP*TE(LXP,J)
      AW(I-1,J)=AW(I-1,J)-TEMPP
      AE(I,J)=TEMPM
      SP(I-1,J)=SP(I-1,J)+TEMPM
      SU(I-1,J)=SU(I-1,J)-TEMPM*TE(LXP,J)
  300 CONTINUE
      IF(LY.AND.NQUICK.EQ.1)GO TO 222
      TEMPS=DS-DABS(CS)*0.1
      COES=0.0D0
      TEMPS=TEMPS/(DS+TINY)
      IF(TEMPS.GT.TINY)COES=DS*TEMPS**5
C     IF(.NOT.LY)COES=DS-DMAX1(FY(J-1)*CS,-(1.-FY(J-1))*CS)
C     COES=DS
      SPPY=0.0D0
      SUPY=0.0D0
      AN(I,J)=0.0D0
      GO TO 400
  222 KYP=J-1-NSVP
      LYP=KYP+(1-NSVP)/2
      DELY3=VM*DELY2-VP*DELY1-VYCV(KYP)*NSVP
      B1Y=DELY2*DELY3*(DELY2-DELY3)
      B2Y=DELY1*DELY3*(DELY1+DELY3)
      BY=B1Y-B2Y+B3Y
      TEMP=B3Y/BY*CS
      COES=DS+(VM*B1Y+VP*B2Y)/BY*CS
      CS=CS-TEMP
      TEMPP=TEMP*VP
      TEMPM=TEMP*VM
      SPPY=-TEMPP
      SUPY=TEMPP*TE(I,LYP)
      AS(I,J-1)=AS(I,J-1)-TEMPP
      AN(I,J)=TEMPM
      SP(I,J-1)=SP(I,J-1)+TEMPM
      SU(I,J-1)=SU(I,J-1)-TEMPM*TE(I,LYP)
  400 CONTINUE
C********** ASSEMBLE COEFFICIENTS FOR E,W,N & S
      AW(I,J)=COEW+CW*UP
      AE(I-1,J)=AW(I,J)-CW+AE(I-1,J)
      AS(I,J)=COES+CS*VP
      AN(I,J-1)=AS(I,J)-CS+AN(I,J-1)
C********** CALCULATE SOURCES
      IF(I.EQ.NI.OR.J.EQ.NJ) GO TO 555
      DENE=DEN(I,J)+FX(I)*(DEN(I+1,J)-DEN(I,J))
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DENN=DEN(I,J)+FY(J)*(DEN(I,J+1)-DEN(I,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      PE=P(I,J)+FX(I)*(P(I+1,J)-P(I,J))
      PW=P(I-1,J)+FX(I-1)*(P(I,J)-P(I-1,J))
      PN=P(I,J)+FY(J)*(P(I,J+1)-P(I,J))
      PS=P(I,J-1)+FY(J-1)*(P(I,J)-P(I,J-1))
      DUDX=(U(I,J)-U(I-1,J))/SXCV(I)
      DVDY=(V(I,J)-V(I,J-1))/SYCV(J)
      IF(J.NE.NJ)
     +UN=0.5*(U(I,J)+U(I-1,J)+FY(J)*(U(I,J+1)+U(I-1,J+1)-U(I,J)-U(I-1,J)
     +))
      US=0.5*(U(I,J-1)+U(I-1,J-1)+FY(J-1)*(U(I,J)+U(I-1,J)-U(I,J-1)
     +-U(I-1,J-1)))
      DUDY=(UN-US)/SYCV(J)
      IF(I.NE.NI)
     +VE=0.5*(V(I,J)+V(I,J-1)+FX(I)*(V(I+1,J)+V(I+1,J-1)-V(I,J)-V(I,J-1)
     +))
      VW=0.5*(V(I-1,J)+V(I-1,J-1)+FX(I-1)*(V(I,J)+V(I,J-1)-V(I-1,J)
     +-V(I-1,J-1)))
      DVDX=(VE-VW)/SXCV(I)
      GEN(I,J)=TVIS(I,J)*(2.*(DUDX**2+DVDY**2)+(DUDY+DVDX)**2)
      IF(IND.EQ.2)GEN(I,J)=GEN(I,J)+2.*TVIS(I,J)*(0.5*(V(I,J)+V(I,J-1))
     +/R(J))**2
      GENC(I,J)=-2.*TVIS(I,J)*((DUDX+DVDY)**2)/3.
     +-2.*DEN(I,J)*TE(I,J)*(DUDX+DVDY)/3. 
      GE2(I,J)=((DENE-DENW)*(PE-PW)/(SXCV(I)**2) 
     ++(DENN-DENS)*(PN-PS)/(SYCV(J)**2))*TVIS(I,J)/(DEN(I,J)**2)
      GE3=GENC(I,J)+GE2(I,J)
      IF(GE3.GT.0.0D0) THEN
      K1=1 
      K2=0
      ELSE 
      K1=0
      K2=1
      END IF
      SU(I,J)=GEN(I,J)+K1*GE3*NCF
      SP(I,J)=-CMU*DEN(I,J)**2*TE(I,J)/TVIS(I,J)
     ++NCF*K2*GE3/TE(I,J)
  555 SU(I,J)=SU(I,J)*VOL+SUPX+SUPY
      SP(I,J)=SP(I,J)*VOL+SPPX+SPPY
      APO(I,J)=DEN(I,J)*VOL/DT
      SU(I,J)=SU(I,J)+APO(I,J)*TEO(I,J)
      SP(I,J)=SP(I,J)-APO(I,J)
  100 CONTINUE
C********** PROBLEM MODIFICATION
      CALL MODTE
C***** FINAL COEFFICIENT ASSEMBLY AND RESIDUAL SOURCE CALCULATION
      RESORK=0.0D0
      DO 200 I=2,NIM1
      DO 200 J=2,NJM1
      SU(I,J)=SUTE(I,J)+URFGEN*(SU(I,J)-SUTE(I,J))
      AP(I,J)=AN(I,J)+AS(I,J)+AE(I,J)+AW(I,J)-SP(I,J)
      RESOR=AN(I,J)*TE(I,J+1)+AS(I,J)*TE(I,J-1)+AE(I,J)*TE(I+1,J)
     ++AW(I,J)*TE(I-1,J)-AP(I,J)*TE(I,J)+SU(I,J)
      IF(-SP(I,J).GE.GREAT)RESOR=0.0D0
      RESORK=RESORK+DABS(RESOR)
C********** UNDER-RELAXATION
      SUTE(I,J)=SU(I,J)
  200 CONTINUE
C********* SOLUTION OF DIFFERENCE EQUATIONS
      CALL LISOLV(2,2,NIM1,NJM1,TE,NSWPK,URFK)
      CALL CORTE
      RETURN
      END
C**************************ENERGY DISSIPATION FIELD*******************      
      SUBROUTINE CALCED
      INCLUDE 'UCOMMON.INC'
      LOGICAL LX,LY 
C********** ASSEMBLY OF COEFFICIENTS
      DO 100 I=2,NI
      DO 100 J=2,NJ
C********** AREAS,VOLUME,DELX & DELY
      AREAW=RSYCV(J)
      AREAS=SXCV(I)*RV(J-1)
      VOL=RSYCV(J)*SXCV(I)
      DELX1=UCVI(I-1)
      DELX2=UCVIP(I-1)
      DELY1=VCVJ(J-1)
      DELY2=VCVJP(J-1)
      B3X=B3XS(I-1)
      B3Y=B3YS(J-1)
      LX=I.GT.2.AND.I.LT.NI
      LY=J.GT.2.AND.J.LT.NJ
C********** CALCULATE CONVECTION
      DENW=DEN(I-1,J)+FX(I-1)*(DEN(I,J)-DEN(I-1,J))
      DENS=DEN(I,J-1)+FY(J-1)*(DEN(I,J)-DEN(I,J-1))
      CW=DENW*U(I-1,J)*AREAW
      CS=DENS*V(I,J-1)*AREAS
C********** CALCULATE DIFFUSION
      GMW=VISCO(I-1,J)+FX(I-1)*(VISCO(I,J)-VISCO(I-1,J))
      GMS=VISCO(I,J-1)+FY(J-1)*(VISCO(I,J)-VISCO(I,J-1))
      GAMW=TVIS(I-1,J)+FX(I-1)*(TVIS(I,J)-TVIS(I-1,J))
      GAMS=TVIS(I,J-1)+FY(J-1)*(TVIS(I,J)-TVIS(I,J-1))
      DW=(GAMW/PRED+GMW)*AREAW/UXCV(I-1)
      DS=(GAMS/PRED+GMS)*AREAS/VYCV(J-1)
C********** SCHEMES
      NSUP=-1
      IF(CW.GT.0.0D0)NSUP=1
      UP=0.5*FLOAT(1+NSUP)
      UM=0.5*FLOAT(1-NSUP)
      NSVP=-1
      IF(CS.GT.0.0D0)NSVP=1
      VP=0.5*FLOAT(1+NSVP)
      VM=0.5*FLOAT(1-NSVP)
      IF(LX.AND.NQUICK.EQ.1)GO TO 111
      TEMPW=DW-DABS(CW)*0.1
      COEW=0.0D0
      TEMPW=TEMPW/(DW+TINY)
      IF(TEMPW.GT.TINY)COEW=DW*TEMPW**5
C     IF(.NOT.LX)COEW=DW-DMAX1(FX(I-1)*CW,-(1.-FX(I-1))*CW)
C     COEW=DW
      SPPX=0.0D0
      SUPX=0.0D0
      AE(I,J)=0.0D0
      GO TO 300
  111 KXP=I-1-NSUP
      LXP=KXP+(1-NSUP)/2
      DELX3=UM*DELX2-UP*DELX1-UXCV(KXP)*NSUP
      B1X=DELX2*DELX3*(DELX2-DELX3)
      B2X=DELX1*DELX3*(DELX1+DELX3)
      BX=B1X-B2X+B3X
      TEMP=B3X/BX*CW
      COEW=DW+(UM*B1X+UP*B2X)/BX*CW
      CW=CW-TEMP
      TEMPP=TEMP*UP
      TEMPM=TEMP*UM
      SPPX=-TEMPP
      SUPX=TEMPP*ED(LXP,J)
      AW(I-1,J)=AW(I-1,J)-TEMPP
      AE(I,J)=TEMPM
      SP(I-1,J)=SP(I-1,J)+TEMPM
      SU(I-1,J)=SU(I-1,J)-TEMPM*ED(LXP,J)
  300 CONTINUE
      IF(LY.AND.NQUICK.EQ.1)GO TO 222
      TEMPS=DS-DABS(CS)*0.1
      COES=0.0D0
      TEMPS=TEMPS/(DS+TINY)
      IF(TEMPS.GT.TINY)COES=DS*TEMPS**5
C     IF(.NOT.LY)COES=DS-DMAX1(FY(J-1)*CS,-(1.-FY(J-1))*CS)
C     COES=DS
      SPPY=0.0D0
      SUPY=0.0D0
      AN(I,J)=0.0D0
      GO TO 400
  222 KYP=J-1-NSVP
      LYP=KYP+(1-NSVP)/2
      DELY3=VM*DELY2-VP*DELY1-VYCV(KYP)*NSVP
      B1Y=DELY2*DELY3*(DELY2-DELY3)
      B2Y=DELY1*DELY3*(DELY1+DELY3)
      BY=B1Y-B2Y+B3Y
      TEMP=B3Y/BY*CS
      COES=DS+(VM*B1Y+VP*B2Y)/BY*CS
      CS=CS-TEMP
      TEMPP=TEMP*VP
      TEMPM=TEMP*VM
      SPPY=-TEMPP
      SUPY=TEMPP*ED(I,LYP)
      AS(I,J-1)=AS(I,J-1)-TEMPP
      AN(I,J)=TEMPM
      SP(I,J-1)=SP(I,J-1)+TEMPM
      SU(I,J-1)=SU(I,J-1)-TEMPM*ED(I,LYP)
  400 CONTINUE
C********** ASSEMBLE COEFFICIENTS FOR E,W,N & S
      AW(I,J)=COEW+CW*UP
      AE(I-1,J)=AW(I,J)-CW+AE(I-1,J)
      AS(I,J)=COES+CS*VP
      AN(I,J-1)=AS(I,J)-CS+AN(I,J-1)
C********** CALCULATE SOURCES
      IF(I.EQ.NI.OR.J.EQ.NJ) GO TO 554
      GE3=GENC(I,J)-GE2(I,J)
      IF(GE3.GT.0.0D0) THEN
      K1=1 
      K2=0
      ELSE 
      K1=0
      K2=1
      END IF
      SU(I,J)=C1*CMU*DEN(I,J)*TE(I,J)*(GEN(I,J)
     ++GE3*K1*NCF)/TVIS(I,J)
      SP(I,J)=-C2*CMU*(DEN(I,J)**2)*TE(I,J)/TVIS(I,J)
     ++K2*GE3*NCF*C1/TE(I,J)
  554 SU(I,J)=SU(I,J)*VOL+SUPX+SUPY
      SP(I,J)=SP(I,J)*VOL+SPPX+SPPY
      APO(I,J)=DEN(I,J)*VOL/DT
      SU(I,J)=SU(I,J)+APO(I,J)*EDO(I,J)
      SP(I,J)=SP(I,J)-APO(I,J)
  100 CONTINUE
C********** PROBLEM MODIFICATION
      CALL MODED
C***** FINAL COEFFICIENT ASSEMBLY AND RESIDUAL SOURCE CALCULATION
      RESORE=0.0D0
      DO 200 I=2,NIM1
      DO 200 J=2,NJM1
      SU(I,J)=SUED(I,J)+URFGEN*(SU(I,J)-SUED(I,J))
      AP(I,J)=AN(I,J)+AS(I,J)+AE(I,J)+AW(I,J)-SP(I,J)
      RESOR=AN(I,J)*ED(I,J+1)+AS(I,J)*ED(I,J-1)+AE(I,J)*ED(I+1,J)
     ++AW(I,J)*ED(I-1,J)-AP(I,J)*ED(I,J)+SU(I,J)
      IF(-SP(I,J).GE.GREAT)RESOR=0.0D0
      RESORE=RESORE+DABS(RESOR)
C********** UNDER-RELAXATION
      SUED(I,J)=SU(I,J)
  200 CONTINUE
C********* SOLUTION OF DIFFERENCE EQUATIONS
      CALL LISOLV(2,2,NIM1,NJM1,ED,NSWPE,URFE)
      CALL CORED
      RETURN
      END
C************************EQUATION SOLVER*******************************      
      SUBROUTINE LISOLV(ISTART,JSTART,NI,NJ,PHI,NSW,URF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (M=82,N=10)
      DIMENSION PHI(M,N),A(M),B(M),C(M),D(M)
      COMMON/COEF/AP(M,N),AN(M,N),AS(M,N),AE(M,N),AW(M,N)
     +,SU(M,N),SP(M,N)
      NIM1=NI
      NJM1=NJ
      FAC=(1.0D0-URF)/URF
      JSTM1=JSTART-1
      ISTM1=ISTART-1
      DO 2000 IT=1,NSW
      A(JSTM1)=0.0D0
C-----COMMENCE W-E SWEEP
      DO 100 I=ISTART,NIM1
      C(JSTM1)=PHI(I,JSTM1)
C-----COMMENCE S-N TRAVERSE
      DO 101 J=JSTART,NJM1
C-----ASSEMBLE TDMA COEFFICIENTS
      A(J)=AN(I,J)
      B(J)=AS(I,J)
      C(J)=AE(I,J)*PHI(I+1,J)+AW(I,J)*PHI(I-1,J)+SU(I,J)
     ++AP(I,J)*PHI(I,J)*FAC 
      D(J)=AP(I,J)/URF
C-----CALCULATE COEFFICIENTS OF RECURRENCE FORMULA
      TERM=1./(D(J)-B(J)*A(J-1))
      A(J)=A(J)*TERM
  101 C(J)=(C(J)+B(J)*C(J-1))*TERM
C-----OBTAIN NEW PHI@S
      DO 102 JJ=JSTART,NJM1
      J=NJ+JSTART-JJ
  102 PHI(I,J)=A(J)*PHI(I,J+1)+C(J)
  100 CONTINUE
      A(ISTM1)=0.0D0
C-----COMMENCE S-N SWEEP
      DO 1000 J=JSTART,NJM1
      C(ISTM1)=PHI(ISTM1,J)
C-----COMMENCE W-E TRAVERSE
      DO 1010 I=ISTART,NIM1
C-----ASSEMBLE TDMA COEFFICIENTS
      A(I)=AE(I,J)
      B(I)=AW(I,J)
      C(I)=AN(I,J)*PHI(I,J+1)+AS(I,J)*PHI(I,J-1)+SU(I,J)
     ++AP(I,J)*PHI(I,J)*FAC 
      D(I)=AP(I,J)/URF
C-----CALCULATE COEFFICIENTS OF RECURRENCE FORMULA
      TERM=1./(D(I)-B(I)*A(I-1))
      A(I)=A(I)*TERM
 1010 C(I)=(C(I)+B(I)*C(I-1))*TERM
C-----OBTAIN NEW PHI@S
      DO 1020 II=ISTART,NIM1
      I=NI+ISTART-II
 1020 PHI(I,J)=A(I)*PHI(I+1,J)+C(I)
 1000 CONTINUE
 2000 CONTINUE
      RETURN
      END
C***************************END OF CODE*****************************
