        PROGRAM GMA
C
C   DOUBLE PREC. VERSION
C
C   G M A  GAUSS-MARKOV-AITKEN LEAST-SQUARES NUCLEAR DATA EVALUATION PROGRAM
C
C          W.P. POENITZ,ANL
C
C*******************************************************************************
CVP 20 Aug 2003: size of parameter vector was increased from 995 to 1200
CVP and B(X) from X=500000 to X=720600 
CVP 2001:Minor correction was introduced by V.G. Pronyaev in the version
CVP of GMA submitted to the Nuclear Data Section of the IAEA in
CVP 1998. These corrections are explained by comments with CVP
CVP and their begin and end is marked by comment lines with CVPBEG**********
CVP and CVPEND**************************************************************  
CVP The references, major formulaes and variables used are added
CVP to keep tracks and for convinience
CVP The file with this version of GMA was renamed in GMAR.FOR 
CVP
CVP The general description of the system of the codes for data
CVP preparation and processing, database for standard neutron reactions
CVP and Gauss-Markov-Aitken least-square procedure used for
CVP fitting is given in [1] Report ANL/NDM-139, W.P. Poenitz and S.E. Aumeier
CVP "The Simultaneous Evaluation of the Standards and Other Cross sections 
CVP of Imporatance for Technology", (1997).
CVP
CVP The basic formulae of the GMA formalism can be taken from e.g. book
CVP [2] S.L. Meyer "Data Analysis for Scientists and Engineeres", 
CVP John Wiley&Sons (1975), formulae (33.15), page 389.
CVP 
CVP Detailed deduction of the formulaes used in GMA code including for 
CVP different type of cross section measurements and their combinations
CVP is given in "Proceedings of the Conference on Nuclear Data Evaluation 
CVP Methods and Procedures", held at BNL, USA, September 22-25, 1980, 
CVP [3] Report BNL-NSC-51363, Vol. 1 (1981), p. 249, W.P. Poenitz," Data 
CVP Interpretation, Objective Evaluation Procedures and Mathematical 
CVP Techniques for the Evaluation of Energy Dependent Ratio, Shape and 
CVP Cross Section Data". The combining of the correlated experimental 
CVP data in sub-blocks substentially reduced the time of matrices inversion. 
CVP 
CVP Notations for most important variables are given by authors in comments.
CVP Few were added under CVP comments 
C
      IMPLICIT REAL*8 (A-H,O-Z)
      character*4 AKON,ACON,DUM
      character*8 CLABL,BREF,CLAB
      character*16 TYPE
      INTEGER*4 ih,im,is,i100,iy,imo,iday
C
C      KAS        indexes of experimental cross sections
C      NT         id of cross sections involved in measured quantity
C      IDEN       data set info (see below)
C      NSETN      shape data set numbers
C      IPP        i/o choices
C      NENF       tags of norm. uncertainty components
C      NETG       tags of energy dependent uncertainty components
C      NCSST      data set Nr.s for cross correlations
C      NEC        error component pairs for cross correlations
C      NRED       data set Nr.s for downweighting
C      KA         indexes
C      NELIM      data set Nr.s to exclude from evaluation
C
      INTEGER*2 KAS(250,5),NT(5),IDEN(30,8),NSETN(200),IPP(8),
     1 NENF(40,10),NETG(11,40),NCSST(10),NEC(2,10,10),NRED(160)
     2 ,KA(1200,250),NELIM(40)
C
C
C      MINIMUM COMPUTER SPACE REQUIREMENTS
C
C      LDA   DIMENSION OF MAX SIZE OF DATA BLOCK,ADJUST E,CSS,AA,DCS,
C            ECOR,KAS,KA,CO,AM
C      LDB   DIMENSION OF PARAMETER SPACE,ADJUST B,BM,DE,AA,EN,CS
C            BM,DE,AA,EN,
C
C   Parameters/apriori
C
C      EN    ENERGY GRID
C      CS    APRIORI CROSS SECTIONS
C      MCS   INDEXES
C
      COMMON /APR/ EN(1200),CS(1200),MCS(35,3)
C
C   Data block / data set   
C
C      E       ENERGIES OF EXPERIMENTAL DATA SET
C      CSS     MEASUREMENT VALUES OF EXPERIMENTAL DATA SET
C      DCS     TOTAL UNCERTAINTIES OF EXPERIMENTAL VALUES
C      CO      ENERGY DEPENDENT UNCERTAINTY COMPONENTS
C      ECOR    CORRELATION MATRIX OF EXPERIMENTS IN DATA BLOCK
C      ENFIS   ENERGIES OF FISSION SPECTRUM
C      FIS     FISSION SPECTRUM*BINWIDTH
C      ENFF    NORMALIZATION UNCERTAINTIES COMPONENTS
C      EPAF    UNCERTAINTY COMPONENT PARAMETERS
C      FCFC    CROSS CORRELATION FACTORS
C
      COMMON /DATA/ E(250),CSS(250),DCS(250),FIS(250),ECOR(250,250)
     1 ,ENFIS(250),CO(12,250),ENFF(30,10),EPAF(3,11,30),FCFC(10,10)
C
C      AA      COEFFICIENT MATRIX
C      AM      MEASUREMENT VECTOR
CVP    ECOR    inverse of correlation matrix of measurement vector AM
CVP            or relative weight matrix
CVP    BM      vector accumulating (in data block cycle) sum of 
CVP            AA(transpose)*ECOR*AM
CVP    B       matrix accumulating (in data block cycle) sum of 
CVP            AA(transpose)*ECOR*AA
C      DE      ADJUSTMENT VECTOR
CVP            equals B(inverse)*BM
C
      COMMON /GAUSS/ AA(1200,250),AM(250),DE(1200),BM(1200),B(720600)
C
C  IDEN(K,I)  K DATA SET SEQUENCE
C
C             I=1 NO OF DATA POINTS IN SET
C             I=2 INDEX OF FIRST VALUE IN ECOR
C             I=3 YEAR
C             I=4 DATA SET TAG
C             I=5 NO OF OTHER SETS WITH WHICH CORRELATIONS ARE GIVEN
C             I=6 DATA SET NR
C             I=7 MT-TYPE OF DATA
C             I=8 1 ABSOLUTE, 2 SHAPE
C
C  CONTROLS/LABELS
C
      COMMON /LABL/ AKON(12),CLAB(35,2),CLABL(4),TYPE(10),BREF(4)
C
      open(3,file='data.gma')
      open(4,file='gma.res')
      OPEN(5,FILE='plot.dta')
C
C      CONTROL
C      MODE    MODE OF EXP. CORRELATION MATRIX CONSTRUCTION
C      I/OC    I/O CONTROL
C      APRI    INPUT OF APRIORI CROSS SECTIONS
C      FIS*    FISSION FLUX FACTORS
C      BLCK    START OF DATA BLOCK
C      DATA    INPUT OF DATA
C      END*    END OF EXP  DATA, START OF EVALUATION
C      EDBL    END OF EXP. DATA SET BLOCK
C      ELIM    DATA SETS TO BE ELIMINATED
C
      AKON( 1)= 'APRI'
      AKON( 2)= 'DATA'
      AKON( 3)= 'END*'
      AKON( 4)= 'BLCK'
      AKON( 5)= 'I/OC'
      AKON( 7)= 'EDBL'
      AKON( 8)= 'FIS*'
      AKON( 9)= 'MODE'
      AKON(10)= 'STOP'
      AKON(11)= 'ELIM'
C
C      MEASUREMENT TYPE IDENTIFICATION
C
      TYPE(1)= 'CROSS SECTION   '
      TYPE(2)= 'CS-SHAPE        '
      TYPE(3)= 'RATIO           '
      TYPE(4)= 'RATIO SHAPE     '
      TYPE(5)= 'TOTAL  CS       ' 
      TYPE(6)= 'FISSION AVERAGE '
      TYPE(7)= 'ABS. S1/(S2+S3) '
      TYPE(8)= 'SHAPE OF SUM    '
      TYPE(9)= 'SHP. S1/(S2+S3) '
C
C      INITIALIZE PARAMETERS
C
C      NTOT TOTAL NO OF DATA POINTS
C      LDA   MAX NO IN DATA BLOCK
C      LDB NO OF UNKNOWNS
C
C      DATA DE/1200*0.D0/,BM/1200*0.D0/,B/720600*0.D0/,NRED/160*0/
C
      MODREP=0
      NTOT=0
      SIGMA2=0.D0
      LDF=200
      LDA=250
      LDB=1200
      LDBB2=LDB*(LDB+1)/2
      MODC=3
      NSHP=0
      NFIS=0
      IELIM=0
      LLL=0
CVP nullning of the vector BM and vector (matrix) B which keep
CVP accumilate sum in cycle on experimental data
CVPBEG********************************************************
      DO 777 I=1,1200
  777 BM(I)=0.D0
      DO 778 I=1,720600
  778 B(I)=0.D0
CVPEND******************************************************** 
C
C      CONTROL
C
C      AKON     CONTROL WORD  A4
C      MC1-8    PARAMETERS   8I5
C
      WRITE(4,110)
  110 FORMAT(1H1,' PROGRAM  GMA'//,'    GAUSS-MARKOV-AITKEN LEAST SQUARE
     1S NUCLEAR DATA EVALUATION'//,'    W.P.POENITZ,ANL'///)
C
C      Control parameter input
C
   50 READ(3,100) ACON,MC1,MC2,MC3,MC4,MC5,MC6,MC7,MC8
  100 FORMAT(A4,1X,8I5)
C
      DO 10 K=1,11
      IF(ACON .EQ. AKON(K)) GO TO (1,2,3,4,5,6,7,8,9,51,170),K
   10 CONTINUE
    6 WRITE(4,104) ACON
  104 FORMAT(A4,2X,'  CONTROL CODE UNKNOWN')
      stop
C
C      test option:  forced stop for testing purpose
C
   51 WRITE(4,107)
  107 FORMAT( '  REQUESTED STOP ' )
      STOP
C
C      input:  data set numbers which are to be excluded from the evaluation
C
  170 IELIM=MC1
      READ(3,171)(NELIM(K),K=1,IELIM)
  171 FORMAT(16I5)
      GO TO 50
C
C      I/O CONTROL
C
    5 IPP(1)=MC1
      IPP(2)=MC2
      IPP(3)=MC3
      IPP(4)=MC4
      IPP(5)=MC5
      IPP(6)=MC6
      IPP(7)=MC7
      IPP(8)=MC8
      GO TO 50
C
C      INPUT OF CROSS SECTIONS TO BE EVALUATED,ENERGY GRID AND APRIORI CS
C
C      MC1=NE      NO OF ENERGIES/PARAMETERS (total)
C      MC2=NC      NO OF CROSS SECTION TYPES
C      NR          NO OF PARAMETERS (cross sections)
C
    1 NE=MC1
      NC=MC2
      NR=0
      NE1=NE+1
      DO 33 K=1,NC
      READ(3,120) (CLAB(K,L),L=1,2)
  120 FORMAT(2A8)
      DO 34 L=1,NE1
      READ(3,103) EX1,CSX1
  103 FORMAT(2E10.4)
      IF(EX1 .EQ. 0.D0) GO TO 35
      NR=NR+1
      EN(NR)=EX1
      CS(NR)=CSX1
   34 CONTINUE
   35 MCS(K,1)=L-1
   33 CONTINUE
C
C  MCS(k,l) CONTAINES IN l= 1 TOTAL NUMBER EACH CS
C                           2 START BOUNDARY
C                           3 END OF BOUNDARY
C
      MCS(1,2)=1
      MCS(1,3 )=MCS(1,1)
      DO 30 K=2,NC
      MCS(K,2)=MCS(K-1,2)+MCS(K-1,1)
   30 MCS(K,3)=MCS(K-1,3)+MCS(K,1)
      WRITE(4,134)
  134 FORMAT(//2X,36HCROSS SECTIONS OF PRESENT EVALUATION//)
      DO 22 K=1,NC
   22 WRITE(4,135) K,(CLAB(K,L),L=1,2)
  135 FORMAT(10X,I3,5X,2A8)
      IF(IPP(1) .EQ. 0) GO TO  665
      WRITE(4,136)
  136 FORMAT(1H1//,2X,35HENERGIES AND APRIORI CROSS SECTIONS//)
      DO 24 K=1,NC
      WRITE(4,137) (CLAB(K,L),L=1,2)
  137 FORMAT(/ '     INDEX     E/MEV   ',7X,2A8 /)
      JC1=MCS(K,2)
      JC2=MCS(K,3)
      LQ=0
      DO 23 L=JC1,JC2
      LQ=LQ+1
  138 FORMAT(2X,2I4,3X,E10.4,3X,F15.8)
   23 WRITE(4,138) LQ,L,EN(L),CS(L)
   24 CONTINUE
  665 CONTINUE
      WRITE(4,113) NR
  113 FORMAT(/,' TOTAL NO OF PARAMETERS ',I4/)
C
C      for checking
C
      if(ipp(7) .eq. 0) go to 50
      write(4,4390)
 4390 format(' No of Parameters per Cross Section '/)
      WRITE(4,154) (MCS(K,1),K=1,NC)
      write(4,4391)
 4391 format(/' Start Address '/)
      WRITE(4,154) (MCS(K,2),K=1,NC)
      write(4,4392)
 4392 format(/' End Address '/)
  154 FORMAT(3X,3HMCS,10I5)
      WRITE(4,154) (MCS(K,3),K=1,NC)
      GO TO 50
C
C      BLOCK INPUT
C
C      N     TOTAL NO OF DATA POINTS IN BLOCK
C      ID    TOTAL NO OF DATA SETS IN BLOCK
C
    4 N=0
      ID=0
      DO 32 K=1,LDA
      AM(K)=0.D0
      DO 81 I=1,LDB
CVP line with BM(I)=0.D0  was commented because BM(I) cleaning should 
CVP done outside of the cycle on measured data 
CVPBEG*****************************************************************
CVP      BM(I)=0.D0
CVPEND*****************************************************************
      KA(I,K)=0
   81 AA(I,K)=0.D0
      DO 165 J=1,5
  165 KAS(K,J)=0
      DO 32 L=1,LDA
   32 ECOR(K,L)=0.D0
      NADD=1
      if(modrep .ne. 0) go to 50
      WRITE(4,108)
  108 FORMAT(/' DATABLOCK************************DATABLOCK**************
     1******************************************DATABLOCK '/)
      GO TO 50
C
C      DATA SET INPUT
C
C      MC1 NS      DATA SET NO
C      MC2 MT      TYPE OF MEASUREMENT
C      MC3 NCOX    CORRELATION MATRIX GIVEN IF .NE. 0
C      MC4 NCT     NO OF CROSS SECTIONS INVOLVED
C      MC5 NT(1)   CROSS SECTION IDENTIFICATION
C      MC6 NT(2)   SAME
C      MC7 NT(3)   SAME
C      MC8 NNCOX   DIVIDE UNCERTAINTIES BY 10.
C
    2 NS=MC1
      MT=MC2
      NCOX=MC3
      NCT=MC4
      NT(1)=MC5
      NT(2)=MC6
      NT(3)=MC7
      NNCOX=MC8
      IF(NCT .GT. 3) READ(3,123) (NT(L),L=4,NCT)
  123 FORMAT(16I5)
      ID=ID+1
      IDEN(ID,2)=N+1
      IDEN(ID,6)=NS
      IDEN(ID,7)=MT
C
C       identify absolute or shape data
C
      MTTP=1
      IDEN(ID,8)=1
      IF(MT .EQ. 2 .OR. MT .EQ. 4) GO TO 510
      IF(MT .EQ. 8 .OR. MT .EQ. 9) GO TO 510
      GO TO 511
  510 MTTP=2
      IDEN(ID,8)=2
  511 CONTINUE
CVP      if(modrep .ne. 0) go to 140
      WRITE(4,142)
  142 FORMAT(//, ' ***********DATASET**************************** '/)
      NU=NCT
      IF(NCT .GT. 4) NU=4
      NCT2=NCT-NU
      NU1=NU+1
      WRITE(4,139)MC1,TYPE(MT),((CLAB(NT(K),L),L=1,2),K=1,NU)
  139 FORMAT(2X,8HDATA SET,I5,2X,A16,4(2X,2A8))
      IF(NCT2 .LE. 0) GO TO 140
      WRITE(4,149)((CLAB(NT(K),L),L=1,2),K=NU1,NCT2)
  149 FORMAT(2X,6(2X,2A8))
  140 CONTINUE
C
C      NAME ID AND REFERENCE I/O
C
      READ(3,131) (IDEN(ID,L),L=3,5),(CLABL(L),L=1,4),(BREF(L),L=1,4)
      WRITE(*,132) (IDEN(ID,L),L=3,4),(CLABL(L),L=1,4),(BREF(L),L=1,4)
CVP      if(modrep .ne. 0) go to 183
      WRITE(4,132) (IDEN(ID,L),L=3,4),(CLABL(L),L=1,4),(BREF(L),L=1,4)
  131 FORMAT(3I5,4A8,4A8)
  132 FORMAT(/5H YEAR,I5,4H TAG,I3,10H AUTHOR:  ,4A8,4A8/)
  183 NCCS=IDEN(ID,5)
C
C      READ(3,    ) NORMALIZATION UNCERTAINTIES
C
      XNORU=0.D0
      IF(MTTP .EQ. 2 ) GO TO 200
      READ(3,201) (ENFF(ID,L),L=1,10),(NENF(ID,L),L=1,10)
  201 FORMAT(10F5.1,10I3)
C
C      CALCULATE TOTAL NORMALIZATION UNCERTAINTY
C
      DO 208 L=1,10
  208 XNORU=XNORU+(ENFF(ID,L))**2
  200 CONTINUE
C
C      READ(3,    ) ENERGY DEPENDENT UNCERTAINTY PARAMETERS
C
      READ(3,202) ((EPAF(L,K,ID),L=1,3),NETG(K,ID),K=1,11)
  202 FORMAT(3F5.2,I3)
C
C      READ(3,    ) CORRELATIONS INFORMATION
C
      IF(NCCS .EQ. 0) GO TO 203
      DO 204 K=1,NCCS
      READ(3,205) NCSST(K),((NEC(M,L,K),M=1,2),L=1,10)
      READ(3,841) (FCFC(L,K),L=1,10)
  841 FORMAT(10F5.1)
  204 CONTINUE
  205 FORMAT(I5,20I3)
  203 CONTINUE
C
C      ACCOUNTING
C
C      N,NADD      NO OF TOTAL DATA POINTS SO FAR IN BLOCK
C      ID          NO OF EXPERIMENTAL DATA SETS
C      NP          NO OF DATA POINTS IN THIS SET
C
      NALT=NADD
      DO 21 KS=1,LDA
      READ(3,109) E(NADD),CSS(NADD),(CO(L,NADD),L=1,12)
  109 FORMAT(2E10.4,12F5.1)
      IF(E(NADD) .EQ. 0.D0) GO TO 95
C
C      SORT EXP ENERGIES  TO FIND CORRESPONDING INDEX OF EVALUATION EN
C
C      KAS(I,L)   GIVES INDEX OF EVALUATION ENERGY FOR I.TH EXP POINT
C                 AND L.TH CROSS SECTION
C
      IF(MT .EQ. 6) GO TO 70
C
C      NCT is the number of cross sections involved
C
      DO 48 L=1,NCT
      JE=MCS(NT(L),2)
      JI=MCS(NT(L),3)
      DO 12 K=JE,JI
      E1=.999*EN(K)
      E2=1.001D0*EN(K)
      IF(E(NADD) .GT. E1 .AND. E(NADD) .LT. E2) GO TO 75
   12 CONTINUE
      GO TO 15
   75 CONTINUE
      KAS(NADD,L)=K
C
C      Exception for dummy data sets
C
      IF(NS .GE. 900 .AND. NS .LE. 909) CSS(NADD)=CS(K)
   48 CONTINUE
C
C      this is the Axton special (uncertainties have been multiplied by 10
C         in order to preserve precision beyond 0.1%)
C
   70 IF(NNCOX .EQ. 0) GO TO 59
      DO 57 LZ=1,11
   57 CO(LZ,NADD)=CO(LZ,NADD)/10.
   59 CONTINUE
C
C      test option:  as set with mode control
C
C      changing weights of data based on year or data set tag
C
      IF(MOD2 .EQ. 0) GO TO 320
      IF(MOD2 .GT. 1000) GO TO 321
      IF(MOD2 .EQ. 10) GO TO 322
      if(mod2 .gt. 10) go to 320
      GO TO (331,336,336,336,336,336,336,336,336),MOD2
  331 CONTINUE
C
C      downweighting data sets with tags .NE. 1
C
      IF(IDEN(ID,4) .EQ. 1) GO TO 320
  342 DO 340 I=3,11
  340 CO(I,NADD)=AMO3*CO(I,NADD)
      GO TO 320
C
C      downweighting based on year of measurement
C
  321 CONTINUE
      IF(IDEN(ID,3) .LT. MOD2) GO TO 342
      GO TO 320
  322 CONTINUE
C
C      downweighting of specified data sets
C
      DO 391 IST=1,NELI
      IF(IDEN(ID,6) .EQ. NRED(IST)) GO TO 342
  391 CONTINUE
      GO TO 320
  336 WRITE(4,339) NS
  339 FORMAT('  WEIGHTING OPTION NOT IMPLEMENTED, DATA SET  ',I5/)
  320 CONTINUE
C
C      CALCULATE TOTAL UNCERTAINTY  DCS
C
      RELU=0.
      DO 207 L=3,11
  207 RELU=RELU+(CO(L,NADD))**2
      DCS(NADD)=DSQRT(XNORU+RELU)
      NADD=NADD+1
   15 CONTINUE
   21 CONTINUE
C
C      all data of set have been red
C
   95 NP=NADD-NALT
      IF(IELIM .EQ. 0) GO TO 173
C
C      data set excluded ?
C
      DO 172 K=1,IELIM
      IF(NS .EQ. NELIM(K)) GO TO 517
  172 CONTINUE
C
C      NO VALID DATA POINTS OR SET REQUESTED TO BE EXCLUDED
C
  173 CONTINUE
      if(np .eq. 1 .and. mttp .eq. 2) go to 517
      IF(NP .NE. 0) GO TO 167
  517 WRITE(4,168) NS
  168 FORMAT(' SET ',I5,' W/O VALID POINTS OR ELIMINATED'/)
      ID=ID-1
      NADD=NALT
      GO TO 50
C
C      continue for valid data
C
  167 IDEN(ID,1)=NP
      NADD1=NADD-1
  460 CONTINUE
C
C      CONSTRUCT ECOR
C
C         MODE  1   INPUT OF ECOR
C               2   UNCORRELATED
C               3   ALL CORRELATED ERRORS GIVEN
C
      IF(NCOX .EQ. 0) GO TO 5001
      MODAL=MODC
      MODC=1
 5001 CONTINUE
      GO TO (54,1779,56,56,56,56),MODC
C
C      INPUT OF ECOR
C
   54 CONTINUE
      DO 61 KS=1,NCOX
  161 FORMAT(10F8.5)
   61 READ(3,161) (ECOR(KS,L),L=1,KS)
      MODC=MODAL
      GO TO 79
C
C       CONSTRUCT ECOR FROM UNCERTAINTY COMPONENTS
C
   56 CONTINUE
      ECOR(NALT,NALT)=1.0D0
      IF(NP .EQ. 1) GO TO 1789
      NALT1=NALT+1
      DO 62 KS=NALT1,NADD1
      C1=DCS(KS)
      KS1=KS-1
      DO 162 KT=NALT,KS1
      Q1=0.
      C2=DCS(KT)
      DO 215 L=3,11
      IF(NETG(L,ID) .EQ. 9) GO TO 214
      IF(NETG(L,ID) .EQ. 0) GO TO 214
      FKS=EPAF(1,L,ID)+EPAF(2,L,ID)
      XYY=EPAF(2,L,ID)-(E(KS)-E(KT))/(EPAF(3,L,ID)*E(KS))
      IF(XYY .LT. 0.D0) XYY=0.
      FKT=EPAF(1,L,ID)+XYY
      Q1=Q1+CO(L,KS)*CO(L,KT)*FKS*FKT
  214 CONTINUE
  215 CONTINUE
      CERR=(Q1+XNORU)/(C1*C2)
C     FIXED THE SINGLE PRECISION ASSIGNMENT
      IF(CERR .GT. .99) CERR= .99D0
      ECOR(KS,KT)=CERR
  162 CONTINUE
      ECOR(KS,KS)=1.
   62 CONTINUE
C
C   ADD CROSS CORRELATIONS OF EXPERIMENTAL DATA BLOCK
C
 1789 CONTINUE
      IF(ID .EQ. 1) GO TO 79
      IF(NCCS .EQ. 0) GO TO 79
      ID1=ID-1
      DO 271 I=1,NCCS
      NSET=NCSST(I)
      DO 272 II=1,ID1
      IF(IDEN(II,6) .EQ. NSET) GO TO 273
  272 CONTINUE
C
C   CORRELATED DATA SET NOT FOUND AHEAD OF PRESENT DATA
C   SET WITHIN DATA BLOCK
C
      WRITE(4,274) NSET,NS
  274 FORMAT('CORRELATED DATA SET  ',I5,' NOT FOUND FOR SET ',I5)
      GO TO 275
  273 NCPP=IDEN(II,1)
C
C      cross correlation
C
      MTT=IDEN(II,8)
      IF(MTT .EQ. 2 .AND. NP .EQ. 1) GO TO 275
      IF(MTTP .EQ. 2 .AND. NCPP .EQ. 1) GO TO 275
  469 NCST=IDEN(II,2)
      NCED=NCPP+NCST-1
      DO 278 K=NALT,NADD1
      C1=DCS(K)
      DO 279 KK=NCST,NCED
      C2=DCS(KK)
      Q1=0.
      DO 281 KKK=1,10
      NC1=NEC(1,KKK,I)
      NC2=NEC(2,KKK,I)
      IF(NC1 .GT. 21 .OR. NC2 .GT. 21) GO TO 2811
      IF(NC1 .EQ. 0) GO TO 2753
      IF(NC2 .EQ. 0) GO TO 2753
      AMUFA=FCFC(KKK,I)
      IF(NC1 .GT. 10) GO TO 310
      C11=ENFF(ID,NC1)
      GO TO 311
  310 NC1=NC1-10
      IF(NETG(NC1,ID) .EQ. 9) GO TO 2800
      FKT=EPAF(1,NC1,ID)+EPAF(2,NC1,ID)
      GO TO 2801
 2800 FKT=1.
 2801 CONTINUE
      C11=FKT*CO(NC1,K)
  311 CONTINUE
      IF(NC2 .GT. 10) GO TO 312
      C22=ENFF(II,NC2)
      GO TO 313
  312 NC2=NC2-10
      IF(NETG(NC2,II) .EQ. 9) GO TO 2802
      XYY=EPAF(2,NC2,II)-ABS(E(K)-E(KK))/(EPAF(3,NC2,II)*E(KK))
      IF(XYY .LT. 0.D0) XYY=0.D0
      FKS=EPAF(1,NC2,II)+XYY
      GO TO 2803
 2802 FKS=1.
 2803 CONTINUE
      C22=FKS*CO(NC2,KK)
  313 CONTINUE
      Q1=Q1+AMUFA*C11*C22
 2811 CONTINUE
  281 CONTINUE
 2753 ECOR(K,KK)=Q1/(C1*C2)
  279 CONTINUE
  278 CONTINUE
  275 CONTINUE
  271 CONTINUE
      GO TO 79
C
C       UNCORRELATED OR SINGLE VALUE
C
 1779 DO 74 KLK=NALT,NADD1
   74 ECOR(KLK,KLK)=1.
   79 CONTINUE
      IF(MT .EQ. 6) GO TO 28
C
C   output of KAS for checking
C
      IF(IPP(7) .EQ. 0) GO TO 2309
      DO 47 K=1,NCT
   47 WRITE(4,702) K,(KAS(L,K),L=NALT,NADD1)
  702 FORMAT(20I5)
C
C      DETERMINE APRIORI NORMALIZATION FOR SHAPE MEASUREMENTS
C
 2309 CONTINUE
      IF(MTTP .EQ. 1 ) GO TO 2828
      NSHP=NSHP+1
      NSETN(NSHP)=NS
      L=NR+NSHP
      IF(L .GT. LDB) GO TO 76
 2828 CONTINUE
      write(4,5173)
CVP   PRIOR/EXP column is added
 5173 format(/'  ENERGY/MEV   VALUE    ABS. UNCERT. ',
     1 ' PRIOR/EXP UNCERT./%    DIFF./%',
     2 '  VAL.*SQRT(E)'/)
      AP=0.
      WWT=0.
      DO 29 K=NALT,NADD1
      CSSK=CSS(K)
      DCSK=DCS(K)
      WXX=1.D0/(DCSK*DCSK)
      WWT=WWT+WXX
   97 KX=KAS(K,1)
      KY=KAS(K,2)
      KZ=KAS(K,3)
      AX=CS(KX)
      IF(MT .EQ. 4 .OR. MT .EQ. 3) AX=AX/CS(KY)
      IF(MT .EQ. 8 .OR. MT .EQ. 5) AX=AX+CS(KY)
      IF(MT .EQ. 5 .AND. NCT .EQ. 3) AX=AX+CS(KZ)
      IF(MT .EQ. 8 .AND. NCT .EQ. 3) AX=AX+CS(KZ)
      IF(MT .EQ. 9 .OR. MT .EQ. 7) AX=AX/(CS(KY)+CS(KZ))
      AZ=AX/CSSK
CVPBEG Assigning uncertainties as % error relative the prior
      IF(MPPP .EQ. 1) DCS(K)=AZ*DCS(K)
CVPEND 
C
C      DATA OUTPUT
C
      IF(IPP(2) .EQ. 0) GO TO 99
      SECS=DSQRT(E(K))*CSSK
      FDQ=DCSK*CSSK/100.D0
      DIFF=(CSSK-AX)*100.D0/AX
CVP   AZ print out was added
      WRITE(4,133) E(K),CSSK,FDQ,AZ,DCSK,DIFF,SECS
CVP   Print out for Ratio of pior/exp value is added
  133 FORMAT(2X,E10.4,2X,E10.4,2X,E10.4,3X,F6.4,3X,F6.2,
     1 3X,F10.1,3X,F10.4)
   99 CONTINUE
C
   29 AP=AP+AZ*WXX
      AP=AP/WWT
CVP      if(modrep .ne. 0) go to 2627
      WRITE(4,111) L,AP,NS,(CLABL(LL),LL=1,4)
  111 FORMAT(/' APRIORI NORM ',I4,F10.4,I5,2X,4A8)
 2627 continue
      IF(MTTP .EQ. 2 ) GO TO 2826
      GO TO 28
 2826 CONTINUE
      IF(MODREP .NE. 0) GO TO 63
C     FIXED THE SINGLE PRECISION ASSIGNMENT
      AP=1.0D0/AP
      CS(L)=AP
      GO TO 28
   63 AP=CS(L)
   28 CONTINUE
C
C      FILL AA,AM,AND COV
C
      DO 18 KS=NALT,NADD1
C     FIXED THE SINGLE PRECISION ASSIGNMENT
      DQQQ=DCS(KS)*CSS(KS)*0.01D0
      J =KAS(KS,1)
      I =KAS(KS,2)
      I8=KAS(KS,3)
      IF(J .EQ. 0 .AND. MT .NE. 6) GO TO 89
      IF(I .NE. 0) GO TO 148
      IF(MT .EQ. 3) GO TO 89
      IF(MT .EQ. 5) GO TO 89
      IF(MT .EQ. 4) GO TO 89
      IF(MT .EQ. 7) GO TO 89
      IF(MT .EQ. 8) GO TO 89
      IF(MT .EQ. 9) GO TO 89
  148 CONTINUE
      IF(I8 .NE. 0) GO TO 147
      IF(MT .EQ. 7) GO TO 89
      IF(MT .EQ. 9) GO TO 89
  147 N=N+1
      IF(MT .EQ. 6) GO TO 46
      IF(MT .EQ. 5) GO TO 45
      IF(MT .EQ. 8) GO TO 248
      KA(J,1)=KA(J,1)+1
      KR=KA(J,1)
      KA(J,KR+1)=N
      GO TO (41,42,43,44,446,446,247,446,249),MT
  446 WRITE(4,447)
  447 FORMAT(10H ERROR 446)
      STOP
C
C      CROSS SECTION
C
   41 CX=CS(J)
      AA(J,KR)=CX/DQQQ
      GO TO 36
C
C      CROSS SECTION SHAPE    L is shape data norm. const. index
C
   42 continue
      CX=CS(J)*AP
      CXX=CX/DQQQ
      AA(J,KR)=CXX
      KA(L,1)=KA(L,1)+1
      KR=KA(L,1)
      KA(L,KR+1)=N
      AA(L,KR)= CXX
      GO TO 36
C
C      RATIO
C
   43 CX=CS(J)/CS(I)
      CCX=CX/DQQQ
      AA(J,KR)=CCX
      KA(I,1)=KA(I,1)+1
      KR=KA(I,1)
      KA(I,KR+1)=N
      AA(I,KR)=-CCX
      GO TO 36
C
C      RATIO SHAPE
C
   44 continue
      CX=CS(J)*AP/CS(I)
      CXX=CX/DQQQ
      AA(J,KR)=CXX
      KA(I,1)=KA(I,1)+1
      KR=KA(I,1)
      KA(I,KR+1)=N
      AA(I,KR)=-CXX
      KA(L,1)=KA(L,1)+1
      KR=KA(L,1)
      KA(L,KR+1)=N
      AA(L,KR)= CXX
      GO TO 36
C
C      TOTAL CROSS SECTION
C
   45 CX=0.
      DO 49 I=1,NCT
      II=KAS(KS,I)
   49 CX=CX+CS(II)
      DO 60 I=1,NCT
      J =KAS(KS,I)
      KA(J,1)=KA(J,1)+1
      KR=KA(J,1)
      KA(J,KR+1)=N
   60 AA(J,KR)=CS(J)/DQQQ
      GO TO 36
C
C      FISSION AVERAGE
C
   46 K=0
      IF(NT(1) .EQ. 9) K=1
      JA=MCS(NT(1),2)
      JE=MCS(NT(1),3)
      NW1=1
      IF(NT(1) .EQ. 9) NW1=2
      NW=NW1
      FL=0.
      SFL=0.
      J1=JA+1
      J2=JE-1
      DO 53 LI=J1,J2
      NW=NW+1
      FL=FL+FIS(NW)
      EL1=(EN(LI)+EN(LI-1))*0.5D0
      EL2=(EN(LI)+EN(LI+1))*0.5D0
      DE1=(EN(LI)-EL1)*0.5D0
      DE2=(EL2-EN(LI))*0.5D0
      SS1=.5*(CS(LI)+0.5D0*(CS(LI)+CS(LI-1)))
      SS2=.5*(CS(LI)+0.5D0*(CS(LI)+CS(LI+1)))
      CSSLI=(SS1*DE1+SS2*DE2)/(DE1+DE2)
      SFL=SFL+CSSLI*FIS(NW)
   53 CONTINUE
      FL=FL+FIS(1)+FIS(NW+1)
      SFL=SFL+FIS(1)*CS(JA)+FIS(NW+1)*CS(JE)
      SFIS=SFL/FL
      WRITE(4,156) EAVR,SFIS,FL,CSS(KS),DCS(KS)
  156 FORMAT( 'AP FISSION AVERAGE ',3F10.4,'  EXP. VAL. ',2F10.4)
      CX=SFIS
      DO 39 J=JA,JE
      K=K+1
      KA(J,1)=KA(J,1)+1
      KR=KA(J,1)
      KA(J,KR+1)=N
      IF(J .EQ. JA .OR. J .EQ. JE) GO TO 195
      EL1=(EN(J)+EN(J-1))*0.5D0
      EL2=(EN(J)+EN(J+1))*0.5D0
      DE1=(EN(J)-EL1)*0.5D0
      DE2=(EL2-EN(J))*0.5D0
      SS1=.5*(CS(J)+0.5D0*(CS(J)+CS(J-1)))
      SS2=.5*(CS(J)+0.5D0*(CS(J)+CS(J+1)))
      CSSJ=(SS1*DE1+SS2*DE2)/(DE1+DE2)
      GO TO 196
  195 CSSJ=CS(J)
  196 CONTINUE
   39 AA(J,KR)=CSSJ*FIS(K)/DQQQ
      GO TO 36
C
C   ABSOLUTE RATIO S1/(S2+S3)
C
  247 CX=CS(J)/(CS(I)+CS(I8))
      IF(I .EQ. J) GO TO 251
      CBX=CX/DQQQ
      AA(J,KR)=CBX
      KA(I,1)=KA(I,1)+1
      KR=KA(I,1)
      KA(I,KR+1)=N
      CBX2=CBX*CBX*DQQQ/CS(J)
      CCX=CBX2*CS(I)
      AA(I,KR)=-CCX
      KA(I8,1)=KA(I8,1)+1
      KR=KA(I8,1)
      KA(I8,KR+1)=N
      CCX=CBX2*CS(I8)
      AA(I8,KR)=-CCX
      GO TO 36
  251 CBX=CX*CX*CS(I8)/(CS(J)*DQQQ)
      AA(J,KR)=CBX
      KA(I8,1)=KA(I8,1)+1
      KR=KA(I8,1)
      KA(I8,KR+1)=N
      AA(I8,KR)=-CBX
      GO TO 36
C
C   SHAPE OF SUM
C
  248 CX=0.
      DO 253 I=1,NCT
      II=KAS(KS,I)
  253 CX=CX+CS(II)*AP
      APDQ=AP/DQQQ
      DO 254 I=1,NCT
      J=KAS(KS,I)
      KA(J,1)=KA(J,1)+1
      KR=KA(J,1)
      KA(J,KR+1)=N
  254 AA(J,KR)=CS(J)*APDQ
      KA(L,1)=KA(L,1)+1
      KR=KA(L,1)
      KA(L,KR+1)=N
      AA(L,KR)=CX/DQQQ
      GO TO 36
C
C   SHAPE OF RATIO S1/(S2+S3)
C
  249 CONTINUE
      CII8=CS(I)+CS(I8)
      CX=AP*CS(J)/CII8
      CBX=CX/DQQQ
      IF(I .EQ. J) GO TO 390
      AA(J,KR)=CBX
      KA(I,1)=KA(I,1)+1
      KR=KA(I,1)
      KA(I,KR+1)=N
      CDX=CBX*CS(I)/CII8
      AA(I,KR)=-CDX
      KA(I8,1)=KA(I8,1)+1
      KR=KA(I8,1)
      KA(I8,KR+1)=N
      CDX=CBX*CS(I8)/CII8
      AA(I8,KR)=-CDX
      KA(L,1)=KA(L,1)+1
      KR=KA(L,1)
      KA(L,KR+1)=N
      AA(L,KR)=CBX
      GO TO 36
  390 CCX=CBX*CS(I8)/CII8
      AA(J,KR)=CCX
      KA(I8,1)=KA(I8,1)+1
      KR=KA(I8,1)
      KA(I8,KR+1)=N
      AA(I8,KR)=-CCX
      KA(L,1)=KA(L,1)+1
      KR=KA(L,1)
      KA(L,KR+1)=N
      AA(L,KR)=CBX
      GO TO 36
   36 AM(N)=(CSS(KS)-CX)/DQQQ
      GO TO 667
   89 WRITE(4,704) MC1,KS
  704 FORMAT( '  DATA POINT BUT NOT AN AP FOR SET ',I5,' NO ',I4)
  667 CONTINUE
   18 CONTINUE
      GO TO 50
C
C    Data BLOCK complete
C
C      FILL IN SYMMETRIC TERM
C

    7 N1=N-1
      IF(ID .EQ. 0) GO TO 50
      WRITE(4,2830) N
 2830 FORMAT(80X,4HN = ,I5)
      IREP=0
      DO 25 K=1,N1
      K1=K+1
      DO 25 L=K1,N
      IF(MODC .EQ. 2) ECOR(L,K)=0.
   25 ECOR(K,L)=ECOR(L,K)
C
C      output of correlation matrix of data block
C
      IF(IPP(3) .EQ. 0) GO TO 67
      IF(N .EQ. 1) GO TO 67
      IF(MODC .EQ. 2) GO TO 67
      WRITE(4,101)
  101 FORMAT(1H*//,'   CORRELATION MATRIX OF DATA BLOCK'/)
      DO 65 K=1,N
   65 WRITE(4,151) (ECOR(K,L),L=1,K)
  151 FORMAT(1X,24F5.2)
   67 CONTINUE
C
C      INVERT ECOR
C
      IF(MODC .EQ. 2) GO TO 19
      IF(N .EQ. 1) GO TO 19
      CALL DPOFA(ECOR,LDA,N,INFO)
      IF(INFO .NE. 0) GO TO 31
      JOB=1
      CALL DPODI(ECOR,LDA,N,DET,JOB)
      DO 17 K=2,N
      L1=K-1
      DO 17 L=1,L1
   17 ECOR(K,L)=ECOR(L,K)
C
C      output of inverted correlation matrix of data block
C
      IF (IPP(5) .EQ. 0) GO TO 19
      DO 20 K=1,N
   20 WRITE(4,151) (ECOR(K,L),L=1,K)
   19 CONTINUE
C
C      GET MATRIX PRODUCTS
C
      NRS=NR+NSHP
      DO 90 I=1,NRS
      NI=KA(I,1)
      IF(NI .EQ. 0) GO TO 92
      DO 83 J=I,NRS
      NJ=KA(J,1)
      IF(NJ .EQ. 0) GO TO 84
      IJ=J*(J-1)/2+I
      DO 85 MI=1,NI
      MIX=KA(I,MI+1)
      DO 85 MJ=1,NJ
      MJX=KA(J,MJ+1)
      B(IJ)=B(IJ)+AA(I,MI)*AA(J,MJ)*ECOR(MIX,MJX)
   85 CONTINUE
   84 CONTINUE
   83 CONTINUE
   92 CONTINUE
   90 CONTINUE
      DO 91 I=1,NRS
      NI=KA(I,1)
      IF(NI .EQ. 0) GO TO 93
      DO 86 MI=1,NI
      MIX=KA(I,MI+1)
      DO 86 MJ=1,N
      BM(I)=BM(I)+AA(I,MI)*ECOR(MIX,MJ)*AM(MJ)
   86 CONTINUE
   93 CONTINUE
   91 CONTINUE
      DO 26 I=1,N
      SUX=0.
      DO 52 J=1,N
   52 SUX=SUX+ECOR(I,J)*AM(J)
   26 SIGMA2=SIGMA2+AM(I)*SUX
      NTOT=NTOT+N
      SIGL=SIGMA2/DBLE(NTOT)
      WRITE(*,476) N,NTOT,NSHP,NRS,SIGL
      if(n .gt. lda) stop
      if(nrs .gt. ldb) stop
      if(modrep .ne. 0)go to 50
      WRITE(4,476) N,NTOT,NSHP,NRS,SIGL
  476 FORMAT(/' ADDED ',I5,' TO GIVE ',I5,' TOTAL',2I5,F10.2/)
      GO TO 50
C
C      GETTING THE RESULT
C
    3 CONTINUE
      REWIND 3
      write(*,6919)
 6919 format(' start getting the result ')
      SIGMAA=SIGMA2/REAL(NTOT-NRS)
      WRITE(4,9679)SIGMAA
 9679 FORMAT(/' UNCERTENTY SCALING   ',E12.4/)
      NRST=NRS*(NRS+1)/2
      IF(IPP(8) .EQ. 0) GO TO 51
      IF(IPP(4) .EQ. 0) GO TO 68
      WRITE(4,116)
  116 FORMAT(1H*//,'  MATRIX PRODUCT'//)
      WRITE(4,152) (B(I),I=1,NRST)
  152 FORMAT(2X,10E10.4)
   68 CONTINUE
      WRITE(4,2840) LDB,NRS,NTOT
 2840 FORMAT(80X,9HLDB,NRS= ,2I6,6H  NTOT,I8)
 7103 FORMAT(2E16.8)
      write(*,6918)
 6918 format(' start on matrix inversion ')
      CALL DPPFA(B,NRS,INFO)
      IF(INFO .NE. 0) GO TO 40
      WRITE(4,9171)
      WRITE(*,9171)
 9171 FORMAT(' INVERT SOLUTION MATRIX')
      JOB=1
      CALL DPPDI(B,NRS,DET,JOB)
      write(*,6917)
 6917 format(' completed inversion of matrix')
CC
      DO 13 I=1,NRS
      DE(I)=0.
      DO 13 K=1,NRS
      IK=K*(K-1)/2+I
      IF(K .LT. I)IK=I*(I-1)/2+K
   13 DE(I)=DE(I)+B(IK)*BM(K)
C
C      output of the result
C
      DO 14 L=1,NC
      WRITE(4,117) (CLAB(L,I),I=1,2)
      WRITE(5,117) (CLAB(L,I),I=1,2)
  117 FORMAT(1H1,'   RESULT',5X,2A8//)
      WRITE(4,112)
  112 FORMAT('   E/MEV         CS/B            DCS/B       DCS/%',
     1'     DIF/%    CS*SQRT(E)'/)
      JA=MCS(L,2)
      JI=MCS(L,3)
      FLX=0.
      DO 77 K=JA,JI
      KBK=K*(K-1)/2+K
      DDX=CS(K)*DSQRT(B(KBK))
      CXX=CS(K)*(1.D0+DE(K))
      CXXD=100.D0*(CXX-CS(K))/CXX
      DO 705 KK=1,NFIS
      IF(ENFIS(KK).GT. .999D0*EN(K).AND.ENFIS(KK).LT. 1.001D0*EN(K))
     1 GO TO 703
  705 CONTINUE
      GO TO 706
  703 IF(K .EQ. JA .OR. K .EQ. JI) GO TO 295
      EL1=(EN(K)+EN(K-1))*0.5D0
      EL2=(EN(K)+EN(K+1))*0.5D0
      DE1=(EN(K)-EL1)*0.5D0
      DE2=(EL2-EN(K))*0.5D0
      SS1=.5*(CXX+0.5D0*(CXX+(1.D0+DE(K-1))*CS(K-1)))
      SS2=.5*(CXX+0.5D0*(CXX+(1.D0+DE(K+1))*CS(K+1)))
      CSSK=(SS1*DE1+SS2*DE2)/(DE1+DE2)
      GO TO 296
  295 CSSK=CXX
  296 CONTINUE
      FLX=FLX+FIS(KK)*CSSK
  706 FQW=DDX*100.D0/CXX
      SECS=DSQRT(EN(K))*CXX
      WRITE(4,153) EN(K),CXX,DDX,FQW,CXXD,SECS
      WRITE(5,153) EN(K),CXX,DDX,FQW,CXXD,SECS
      IF(MODAP .EQ. 0) GO TO 58
      IF(MODAP .EQ. 2 .AND. K .LE. MCS(5,3)) GO TO 58
      CS(K)=CXX
   58 CONTINUE
  153 FORMAT(1X,E10.4,2F15.8,2X,F6.1,3X,F7.2,3X,F10.5)
   77 CONTINUE
      WRITE(4,158) FLX
  158 FORMAT(1H*//,'  FISSION AVERAGE ' ,F8.4//)
   14 CONTINUE
C
C   OUTPUT OF NORM. FACTORS FOR SHAPE DATA
C
      WRITE(4,114)
  114 FORMAT(1H*///, '  NORMALIZATION  OF SHAPE DATA '///)
      NR1=NR+1
      LLX=0
      IF(NSHP .EQ. 0) GO TO 292
      DO 82 K=NR1,NRS
      LLX=LLX+1
      KK=K*(K-1)/2+K
      ZCS=CS(K)
      DDX=CS(K)*DSQRT(B(KK))
      CXX=CS(K)*(1.D0+DE(K))
      DDXD=DDX*100.D0/CXX
      WRITE(4,115) K,NSETN(LLX),CXX,DDX,DDXD,ZCS
      CS(K)=CXX
   82 CONTINUE
  115 FORMAT(2I6,4F10.4)
 5115 FORMAT(2I6,2D20.12,F10.4)
  292 CONTINUE
C
C     reset for repeat of fit with replaced apriori from first fit
C
      IF(MODAP .EQ. 0 .or. modrep .eq. MODAP) GO TO 64
      MODREP=MODREP+1
      NTOT=0
      SIGMA2=0.D0
      NSHP=0
      DO 72 L=1,LDB
      DE(L)=0.D0
   72 BM(L)=0.D0
      DO 73 L=1,LDBB2
   73 B(L)=0.D0
CC
      DO 69 L=1,2000
      READ(3,130)DUM   
  130 FORMAT(A4)
      IF(DUM .EQ. AKON(4)) GO TO 4
   69 CONTINUE
   64 CONTINUE
C
C   OUTPUT OF CORRELATION MATRIX OF THE RESULT
C
      IF(IPP(6) .eq. 0) go to 184
      DO 78  K=1,NC
      J1=MCS(K,2)
      J2=MCS(K,3)
      DO 80 L=1,K
      WRITE(4,122) (CLAB(K,J),J=1,2),(CLAB(L,J),J=1,2)
  122 FORMAT(1H1, '  CORRELATION MATRIX OF THE RESULT   ',2A8,2A8///)
      J3=MCS(L,2)
      J4=MCS(L,3)
      IF(K .EQ. L) GO TO 87
      DO 88 I=J1,J2
      II=I*(I-1)/2+I
      DO 16 J=J3,J4
      IJ=I*(I-1)/2+J
      JJ=J*(J-1)/2+J
   16 BM(J)=B(IJ)/DSQRT(B(II)*B(JJ))
   88 WRITE(4,151) (BM( M),M=J3,J4)
      GO TO 300
   87 DO 55 I=J1,J2
      II=I*(I-1)/2+I
      DO 27 J=J1,I
      IJ=I*(I-1)/2+J
      JJ=J*(J-1)/2+J
   27 BM(J)=B(IJ)/DSQRT(B(II)*B(JJ))
   55  WRITE(4,151) (BM( M),M=J1,I)
  300 CONTINUE
   80 CONTINUE
   78 CONTINUE
  184 continue
      STOP
C
C      INPUT OF FISSION SPECTRUM
C
    8 CONTINUE
      IF(MC1 .EQ. 0) GO TO 692
      DO 690 K=1,LDF
      READ(3,119) ENFIS(K),FIS(K)
  119 FORMAT(2E13.5)
      IF(ENFIS(K) .EQ. 0.) GO TO 691
  690 CONTINUE
  691 NFIS=K-1
      GO TO 38
CC
CC      MAXWELLIAN SPECTRUM
CC
  692 EAVR=DBLE(MC2)/1000.D0
      NFIS=MCS(MC3,1)
      JA=MCS(MC3,2)
      JE=MCS(MC3,3)
      LL=0
      DO 693 K=JA,JE
      LL=LL+1
  693 ENFIS(LL)=EN(K)
      NFIS1=NFIS-1
      FISUM=0.D0
      DO 695 K=2,NFIS1
      E1=(ENFIS(K-1)+ENFIS(K))/2.D0
      E2=(ENFIS(K+1)+ENFIS(K))/2.D0
      DE12=E2-E1
      F1=DSQRT(E1)*DEXP(-1.5D0*E1/EAVR)
      F2=DSQRT(E2)*DEXP(-1.5D0*E2/EAVR)
      E12=(E1+E2)/2.D0
      F3=DSQRT(E12)*DEXP(-1.5D0*E12/EAVR)
      FIS(K)=((F1+F2)*.5D0+F3)*.5D0
  695 FISUM=FISUM+DE12*FIS(K)
      FIS(1)=DSQRT(EN(JA))*DEXP(-1.5D0*EN(JA)/EAVR)
      FIS(NFIS)=DSQRT(EN(JE))*DEXP(-1.5D0*EN(JE)/EAVR)
      DE12=(ENFIS(2)-ENFIS(1))/2.D0
      DE13=ENFIS(1)+DE12
      FISUM=FISUM+FIS(1)*DE13
      DE14=(ENFIS(NFIS)-ENFIS(NFIS1))/2.D0
      FISUM=FISUM+FIS(NFIS)*2.D0*DE14
   38 WRITE(4,800)
  800 FORMAT(/' FISSION SPECTRUM * BIN WIDTH'/)
      IF(MC1 .NE. 0) GO TO 189
      DO 696 K=2,NFIS1
      E1=(ENFIS(K-1)+ENFIS(K))/2.D0
      E2=(ENFIS(K+1)+ENFIS(K))/2.D0
      DE12=E2-E1
  696 FIS(K)=FIS(K)*DE12/FISUM
      FIS(NFIS)=FIS(NFIS)*DE14/FISUM
      FIS(1)=FIS(1)*DE13/FISUM
  189 continue
      do 694 kq=1,nfis
  694 WRITE(4,157) ENFIS(KQ),FIS(KQ)
  157 FORMAT(2F10.6)
      GO TO 50
C
C   MODE DEFINITION
C
    9 MODC=MC1
      MOD2=MC2
CVPBEG MPPP=1 allows to use anti-PPP option, when errors of exp data 
CVP   are taken as % uncertainties from true (posterior) evaluation 
      MPPP=MC5
CVPEND
      AMO3=DBLE(MC3)/10.
      MODAP=MC4
      IF(MC2 .NE. 10) GO TO 50
C
C      test option:  input of data set numbers which are to be downweighted
C
      K1=1
      K2=16
      DO 678 K=1,10
      READ(3,677) (NRED(I),I=K1,K2)
      IF(NRED(K2) .EQ. 0) GO TO 679
      K1=K1+16
      K2=K2+16
  678 CONTINUE
  677 FORMAT(16I5)
  679 DO 680 K=K1,K2
      IF(NRED(K) .EQ. 0) GO TO 681
  680 CONTINUE
  681 NELI=K-1
      GO TO 50
C
C      MESSAGE
C
   76 WRITE(4,701) MC1
  701 FORMAT( '   OVERFLOW OF UNKNOWN-VECTOR SPACE WITH SET  ',I3)
      STOP
C
C      ATTEMPT TO MAKE CORR. MATRIX POSITIVE DEFINITE
C
   31 WRITE(4,105)
      IREP=IREP+1
      N1=N-1
      DO 2211 K=1,N1
      K1=K+1
      DO 2211 L=K1,N
      IF(MODC .EQ. 2)ECOR(L,K)=0.
 2211 ECOR(K,L)=ECOR(L,K)
      DO 2212 K=1,N
 2212 ECOR(K,K)=1.
      CXZ=0.10D0
      DO 37 K=1,N
      DO 37 L=1,N
      ECOR(K,L)=ECOR(K,L)/(1.D0+CXZ)
      IF(K .EQ. L) ECOR(K,L)=1.
   37 CONTINUE
      IF(IREP .LT. 15) GO TO 67
      GO TO 50
C
C      MESSAGES
C
   40 WRITE(4,106)
      STOP
  105 FORMAT(/' EXP BLOCK CORREL. MATRIX NOT PD',20X,'***** WARNING *')
  106 FORMAT( '  SOLUTION  CORREL. MATRIX NOT PD ' )
      END
      SUBROUTINE DPODI(A,LDA,N,DET,JOB)                                 
C
C    this file contains all linpack and blas routines required by GMA
C
      INTEGER LDA,N,JOB                                                 
      DOUBLE PRECISION A(LDA,1)                                         
      DOUBLE PRECISION DET(2)                                           
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION S                                                
      INTEGER I,J,JM1,K,KP1                                             
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         S = 10.0D0                                                     
         DO 50 I = 1, N                                                 
            DET(1) = A(I,I)**2*DET(1)                                   
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DET(1) .GE. 1.0D0) GO TO 20                             
               DET(1) = S*DET(1)                                        
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (DET(1) .LT. S) GO TO 40                                 
               DET(1) = DET(1)/S                                        
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
      IF (MOD(JOB,10) .EQ. 0) GO TO 140                                 
         DO 100 K = 1, N                                                
            A(K,K) = 1.0D0/A(K,K)                                       
            T = -A(K,K)                                                 
            CALL DSCAL(K-1,T,A(1,K),1)                                  
            KP1 = K + 1                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = A(K,J)                                               
               A(K,J) = 0.0D0                                           
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
         DO 130 J = 1, N                                                
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 120                                   
            DO 110 K = 1, JM1                                           
               T = A(K,J)                                               
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)                        
  110       CONTINUE                                                    
  120       CONTINUE                                                    
            T = A(J,J)                                                  
            CALL DSCAL(J,T,A(1,J),1)                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPOFA(A,LDA,N,INFO)                                    
      INTEGER LDA,N,INFO                                                
      DOUBLE PRECISION A(LDA,1)                                         
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER J,JM1,K                                                   
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0                                                   
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 20                                    
            DO 10 K = 1, JM1                                            
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)                 
               T = T/A(K,K)                                             
               A(K,J) = T                                               
               S = S + T*T                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            S = A(J,J) - S                                              
            IF (S .LE. 0.0D0) GO TO 40                                  
            A(J,J) = DSQRT(S)                                           
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPDI(AP,N,DET,JOB)                                    
      INTEGER N,JOB                                                     
      DOUBLE PRECISION AP(1)                                            
      DOUBLE PRECISION DET(2)                                           
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION S                                                
      INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1                           
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         S = 10.0D0                                                     
         II = 0                                                         
         DO 50 I = 1, N                                                 
            II = II + I                                                 
            DET(1) = AP(II)**2*DET(1)                                   
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DET(1) .GE. 1.0D0) GO TO 20                             
               DET(1) = S*DET(1)                                        
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (DET(1) .LT. S) GO TO 40                                 
               DET(1) = DET(1)/S                                        
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
      IF (MOD(JOB,10) .EQ. 0) GO TO 140                                 
         KK = 0                                                         
         DO 100 K = 1, N                                                
            K1 = KK + 1                                                 
            KK = KK + K                                                 
            AP(KK) = 1.0D0/AP(KK)                                       
            T = -AP(KK)                                                 
            CALL DSCAL(K-1,T,AP(K1),1)                                  
            KP1 = K + 1                                                 
            J1 = KK + 1                                                 
            KJ = KK + K                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = AP(KJ)                                               
               AP(KJ) = 0.0D0                                           
               CALL DAXPY(K,T,AP(K1),1,AP(J1),1)                        
               J1 = J1 + J                                              
               KJ = KJ + J                                              
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
         JJ = 0                                                         
         DO 130 J = 1, N                                                
            J1 = JJ + 1                                                 
            JJ = JJ + J                                                 
            JM1 = J - 1                                                 
            K1 = 1                                                      
            KJ = J1                                                     
            IF (JM1 .LT. 1) GO TO 120                                   
            DO 110 K = 1, JM1                                           
               T = AP(KJ)                                               
               CALL DAXPY(K,T,AP(J1),1,AP(K1),1)                        
               K1 = K1 + K                                              
               KJ = KJ + 1                                              
  110       CONTINUE                                                    
  120       CONTINUE                                                    
            T = AP(JJ)                                                  
            CALL DSCAL(J,T,AP(J1),1)                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPFA(AP,N,INFO)                                       
      INTEGER N,INFO                                                    
      DOUBLE PRECISION AP(1)                                            
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER J,JJ,JM1,K,KJ,KK                                          
         JJ = 0                                                         
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0                                                   
            JM1 = J - 1                                                 
            KJ = JJ                                                     
            KK = 0                                                      
            IF (JM1 .LT. 1) GO TO 20                                    
            DO 10 K = 1, JM1                                            
               KJ = KJ + 1                                              
               T = AP(KJ) - DDOT(K-1,AP(KK+1),1,AP(JJ+1),1)             
               KK = KK + K                                              
               T = T/AP(KK)                                             
               AP(KJ) = T                                               
               S = S + T*T                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            JJ = JJ + J                                                 
            S = AP(JJ) - S                                              
            IF (S .LE. 0.0D0) GO TO 40                                  
            AP(JJ) = DSQRT(S)                                           
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPSL(AP,N,B)                                          
      INTEGER N                                                         
      DOUBLE PRECISION AP(1),B(1)                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,KK                                                   
      KK = 0                                                            
      DO 10 K = 1, N                                                    
         T = DDOT(K-1,AP(KK+1),1,B(1),1)                                
         KK = KK + K                                                    
         B(K) = (B(K) - T)/AP(KK)                                       
   10 CONTINUE                                                          
      DO 20 KB = 1, N                                                   
         K = N + 1 - KB                                                 
         B(K) = B(K)/AP(KK)                                             
         KK = KK - K                                                    
         T = -B(K)                                                      
         CALL DAXPY(K-1,T,AP(KK+1),1,B(1),1)                            
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      subroutine daxpy(n,da,dx,incx,dy,incy)                            
c                                                                       
c     constant times a vector plus a vector.                            
c     uses unrolled loops for increments equal to one.                  
c     jack dongarra, linpack, 3/11/78.                                  
c                                                                       
      double precision dx(1),dy(1),da                                   
      integer i,incx,incy,ix,iy,m,mp1,n                                 
c                                                                       
      if(n.le.0D0)return                                                  
      if (da .eq. 0.0d0) return                                         
      if(incx.eq.1.D0 .and. incy.eq.1D0)go to 20                               
c                                                                       
c        code for unequal increments or equal increments                
c          not equal to 1                                               
c                                                                       
      ix = 1                                                            
      iy = 1                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                 
      do 10 i = 1,n                                                     
        dy(iy) = dy(iy) + da*dx(ix)                                     
        ix = ix + incx                                                  
        iy = iy + incy                                                  
   10 continue                                                          
      return                                                            
c                                                                       
c        code for both increments equal to 1                            
c                                                                       
c                                                                       
c        clean-up loop                                                  
c                                                                       
   20 m = mod(n,4)                                                      
      if( m .eq. 0 ) go to 40                                           
      do 30 i = 1,m                                                     
        dy(i) = dy(i) + da*dx(i)                                        
   30 continue                                                          
      if( n .lt. 4 ) return                                             
   40 mp1 = m + 1                                                       
      do 50 i = mp1,n,4                                                 
        dy(i) = dy(i) + da*dx(i)                                        
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)                            
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)                            
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)                            
   50 continue                                                          
      return                                                            
      end                                                               
      double precision function ddot(n,dx,incx,dy,incy)                 
c                                                                       
c     forms the dot product of two vectors.                             
c     uses unrolled loops for increments equal to one.                  
c     jack dongarra, linpack, 3/11/78.                                  
c                                                                       
      double precision dx(1),dy(1),dtemp                                
      integer i,incx,incy,ix,iy,m,mp1,n                                 
c                                                                       
      ddot = 0.0d0                                                      
      dtemp = 0.0d0                                                     
      if(n.le.0)return                                                  
      if(incx.eq.1.D0 .and.incy.eq. 1.D0) go to 20
c                                                                       
c        code for unequal increments or equal increments                
c          not equal to 1                                               
c                                                                       
      ix = 1                                                            
      iy = 1                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                 
      do 10 i = 1,n                                                     
        dtemp = dtemp + dx(ix)*dy(iy)                                   
        ix = ix + incx                                                  
        iy = iy + incy                                                  
   10 continue                                                          
      ddot = dtemp                                                      
      return                                                            
c                                                                       
c        code for both increments equal to 1                            
c                                                                       
c                                                                       
c        clean-up loop                                                  
c                                                                       
   20 m = mod(n,5)                                                      
      if( m .eq. 0 ) go to 40                                           
      do 30 i = 1,m                                                     
        dtemp = dtemp + dx(i)*dy(i)                                     
   30 continue                                                          
      if( n .lt. 5 ) go to 60                                           
   40 mp1 = m + 1                                                       
      do 50 i = mp1,n,5                                                 
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +             
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue                                                          
   60 ddot = dtemp                                                      
      return                                                            
      end                                                               
      subroutine  dscal(n,da,dx,incx)                                   
c                                                                       
c     scales a vector by a constant.                                    
c     uses unrolled loops for increment equal to one.                   
c     jack dongarra, linpack, 3/11/78.                                  
c     modified to correct problem with negative increment, 8/21/90.     
c                                                                       
      double precision da,dx(1)                                         
      integer i,incx,ix,m,mp1,n                                         
c                                                                       
      if(n.le.0)return                                                  
      if(incx.eq.1)go to 20                                             
c                                                                       
c        code for increment not equal to 1                              
c                                                                       
      ix = 1                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                 
      do 10 i = 1,n                                                     
        dx(ix) = da*dx(ix)                                              
        ix = ix + incx                                                  
   10 continue                                                          
      return                                                            
c                                                                       
c        code for increment equal to 1                                  
c                                                                       
c                                                                       
c        clean-up loop                                                  
c                                                                       
   20 m = mod(n,5)                                                      
      if( m .eq. 0 ) go to 40                                           
      do 30 i = 1,m                                                     
        dx(i) = da*dx(i)                                                
   30 continue                                                          
      if( n .lt. 5 ) return                                             
   40 mp1 = m + 1                                                       
      do 50 i = mp1,n,5                                                 
        dx(i) = da*dx(i)                                                
        dx(i + 1) = da*dx(i + 1)                                        
        dx(i + 2) = da*dx(i + 2)                                        
        dx(i + 3) = da*dx(i + 3)                                        
        dx(i + 4) = da*dx(i + 4)                                        
   50 continue                                                          
      return                                                            
      end                                                               
                                                                        
