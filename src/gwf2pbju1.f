C *****************************************************************************
C -- Polyline Boundary Junction (PBJ) Package
C ---- by Leland Scantlebury and James R. Craig
C -----------------------------------------------------------------------------
      module pbjmodule
        integer :: nsegments
        integer :: IPBJCB                                   ! Write CBC flag
        integer :: IPRPBJ=1                                 ! LST Printout flag (hardcoded, for now)
        integer :: pbjmode                                  ! 0-spec head, 1-drain, 2-stage dep
        integer :: condtype                                 ! 0-conductivity, 1-unit conductivity (per len) 2-leakance coeff
                                                            ! (-1 == not read, head spec mode)
        integer, allocatable, dimension(:,:)   :: segnodes  ! Nodes to which flow is interpolated, for each segment
        integer, allocatable, dimension(:,:,:) :: segIA     ! Off-diagonal location in AMAT of nodes listed in segia
        integer, allocatable, dimension(:)     :: pbjact    ! Active nodes for SP (1-active, 0-inactive)
        real,    allocatable, dimension(:,:)   :: bw        ! Barycentric coordinates for segment start/end points
        real,    allocatable, dimension(:,:)   :: segelevs  ! Segment start/end point streambed elevations
        real,    allocatable, dimension(:)     :: seglens   ! Segment lengths (for unit conductivity, etc)
        real*8,  allocatable, dimension(:,:)   :: cond      ! Segment conductivity
        real*8,  allocatable, dimension(:,:)   :: seghead   ! Specified heads (pbjmode=0) or external stage (pbjmode=2) at each end
        double precision, parameter :: zero=0.0D0
        CHARACTER*100 PBJ_VERSION
        DATA PBJ_VERSION /'PBJ -- POLYLINE BOUNDARY JUNCTION PACKAGE,
     1 VERSION 1, 7/14/2020 INPUT READ FROM UNIT'/
        
      contains
C -----------------------------------------------------------------------------
        integer function find_n2m_inJA(n, m, seg) result(ija)
          use GLOBAL,     ONLY:IOUT,IA,JA
          integer, intent(in)       :: n,m,seg
          integer                   :: i
      
          ija = -1
          if (m==n) then  ! Small cheat
            ija = IA(n)
            return
          end if
          do i=1, (IA(n+1)-IA(n)-1)    ! num nodes attached to n
            if (JA(IA(n)+i) == m) then
              ija = IA(n)+i
              exit                     ! Leave loop if number found
            end if
          end do
          if (ija < 0) then            ! Error in number never found
            write(IOUT,*) ' UNCONNECTED NODES LISTED FOR PBJ SEGMENT'
            write(IOUT,*) 'SEGMENT: ', seg, 'NODES: ', m, n
            call USTOP(' ')
          end if
          return
        end function find_n2m_inJA
C -----------------------------------------------------------------------------
      
      end module pbjmodule
      
      
C *****************************************************************************
C -- Standard MODFLOW Package Routines
C -----------------------------------------------------------------------------
      
      subroutine GWF2PBJU1AR(IN)
C     ******************************************************************
C     Read PBJ input file and allocate arrays
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,IA,JA,AMAT
      use pbjmodule
      integer  :: i, j, k, m, n
      integer, allocatable, dimension(:) :: nodetemp
      real, allocatable, dimension(:) :: weighttemp, elevtemp
      character*200 line
      
      CHARACTER*24 ANAME(4)
      DATA ANAME(1) /'           SEGMENT NODES'/
      DATA ANAME(2) /'     BARYCENTRIC WEIGHTS'/
      DATA ANAME(3) /'      SEGMENT ELEVATIONS'/
      DATA ANAME(4) /'          SEGMENT LENGTH'/
C     ------------------------------------------------------------------      
      
C-----Identify package, initialize variables
      write(IOUT,1) trim(PBJ_VERSION), in
    1 format(1X,/1X,A,I5)
      nsegments = 0
      
C-----Read in number of segments and IPBJCB flag
      call URDCOM(IN,IOUT,line)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,nsegments,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPBJCB,R,IOUT,IN)
      
C-----Read in Mode and Conductance Type
      LLOC=1
      call URDCOM(IN,IOUT,line)
      call URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
C-----Identify Mode
      if (LINE(ISTART:ISTOP)=='HEADSPEC') then
        pbjmode = 0
        condtype = -1
      else if (LINE(ISTART:ISTOP)=='DRAIN') then
        pbjmode = 1
      else if (LINE(ISTART:ISTOP)=='EXTSTAGE') then
        pbjmode = 2
      else
        write(IOUT,*) ' INVALID PBJ PACKAGE MODE'
        call USTOP(' ')
      end if
C-----Identify Conductivity Type
      if (pbjmode /= 0) then
        call URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
        if (LINE(ISTART:ISTOP)=='CONDUCTANCE') then
          condtype = 0
        else if (LINE(ISTART:ISTOP)=='UNITCOND') then
          condtype = 1
        else if (LINE(ISTART:ISTOP)=='LEAKCOEF') then
          condtype = 2
        else
          write(IOUT,*) ' INVALID PBJ PACKAGE CONDUCTIVITY TYPE'
          call USTOP(' ')
        end if
      end if
      
C-----Initialize Arrays based upon nsegments
      call initialize_pbj_arrays()
      
C-----Read segment nodes
C     There's not an existing routine for reading a 3D array of integers, so
C     read as a 1D array and then convert
      K = 0
      allocate(nodetemp(nsegments*3))           ! 3 Nodes per segment (triangles!)
      call U1DINT(nodetemp,ANAME(1),nsegments*3,K,IN,IOUT)
      
C-----Read segment barycentric weights
      allocate(weighttemp(nsegments*6))         ! 3 weights for every segment start/end
      call U1DREL(weighttemp,ANAME(2),nsegments*6,K,IN,IOUT)
      
C-----Read segment start/end elevations
      if (pbjmode /= 0) then
        allocate(elevtemp(nsegments*2))           ! 1 elevation for every segment start/end
        call U1DREL(elevtemp,ANAME(3),nsegments*2,K,IN,IOUT)
      end if
      
C-----Read segment lengths (if required by condtype)
      if (condtype > 0) then
        call U1DREL(seglens,ANAME(4),nsegments,K,IN,IOUT)
      end if

C-----Move temporary arrays to their final resting places
      do i=1, nsegments
        do j=1, 3                               ! Loop over nodes
          segnodes(j,i) = nodetemp((i-1)*3 + j)
          n = segnodes(j,i)                     ! current diagonal node
          bw(j,i)   = weighttemp((i-1)*6 + j)
          bw(j+3,i) = weighttemp((i-1)*6 + j+3)
C-----Find connected nodes in JA (needed later for entering into AMAT)
          do k=1, 3
            m = nodetemp((i-1)*3 + k)           ! current off-diag node
            segIA(k,j,i) = find_n2m_inJA(n,m,i)
          end do
        end do
        if (pbjmode /= 0) then                    ! Elev not needed for spec head
          do j=1, 2                               ! Loop over segment ends
            segelevs(j,i) = elevtemp((i-1)*2 + j)
          end do
        end if
      end do
      
C-----Deallocate
      deallocate(nodetemp, weighttemp)
      if (pbjmode /= 0) deallocate(elevtemp)
      
C-----Return
      return
      end subroutine GWF2PBJU1AR

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1RP(IN)
C     ******************************************************************
C     Read in stress period (time-varying) values
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      use pbjmodule
      USE GLOBAL,       ONLY:IOUT,IFREFM,NEQS
      integer :: itmp
      CHARACTER*20 ANAME(0:2)
      DATA ANAME(0) /'        CONDUCTIVITY'/
      DATA ANAME(1) /'   UNIT CONDUCTIVITY'/
      DATA ANAME(2) /' LEAKANCE COEFFICENT'/
      CHARACTER*20 HNAME(0:2)
      DATA HNAME(0) /'      SPECIFIED HEAD'/
      DATA HNAME(1) /''/                     !PBJmode 1 (Drain) only requires time variant conductivities
      DATA HNAME(2) /'      EXTERNAL STAGE'/
C     ------------------------------------------------------------------
C-----Identify package
      write(IOUT,6) trim(PBJ_VERSION), in
    6 format(1X,/1X,A,I5)
      
C-----SPECIFIED HEADS or EXTERNAL STAGE
      if (pbjmode /= 1) then
        qspec = 0.0          ! Reset
C-----Read ITMP (flag to re-use data or not)
        if (IFREFM.EQ.0) then
          read(IN,'(I10)') itmp
        else
          read(IN,*) itmp
        end if
C-----Process itmp
        IF(itmp < 0) THEN
          WRITE(IOUT,7) 'REUSING PBJ ',trim(HNAME(pbjmode))
    7     FORMAT(1X,/1X,2A,' VALUES FROM LAST STRESS PERIOD')
C-----Read in values
        else
          pbjact = 0
          call GWF2PBJU1R(seghead,2,itmp,HNAME(pbjmode),IN,IOUT,IFREFM)
        end if
      end if
      
C-----CONDUCTANCES
      if (pbjmode > 0) then
C-----Read ITMP (flag to re-use data or not)
        if (IFREFM.EQ.0) then
          read(IN,'(I10)') itmp
        else
          read(IN,*) itmp
        end if
C-----Process itmp
        if(itmp < 0) THEN
          WRITE(IOUT,8) trim(ANAME(condtype))
    8     FORMAT(1X,/1X, 
     1         'REUSING PBJ ',A,' VALUES FROM LAST STRESS PERIOD')
        else
          pbjact = 0
          call GWF2PBJU1R(cond,2,itmp,ANAME(condtype),IN,IOUT,IFREFM)
C-----Where possible, pre-process conductances for stress period
          if (condtype > 0) then
C-----Unit conductivities & leakance coefficients can be multiplied by segment lengths
            do i=1, nsegments
              cond(:,i) = cond(:,i)*seglens(i)
            end do
          end if
        end if !end of itmp check
      end if !end of pbjmode check
        
      return
      end subroutine GWF2PBJU1RP

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1FM
C     ******************************************************************
C     ADD SEGMENT FLOW TO SOURCE TERM
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      use pbjmodule
      use GLOBAL,     ONLY:IBOUND,HNEW,RHS,AMAT,IA
      integer            :: h, i, j, k, N, ija
      double precision   :: heads(2), bwe(2), segfunc, debug, estflow
      double precision, parameter :: bignum = 10.0**10.0
C     ------------------------------------------------------------------
      
C-----If no segments, return
      if (nsegments <= 0) return

C-----------------------------------------------------------------------
C-----PBJ Mode 0 - Specified heads
      if (pbjmode == 0) then
        do i=1, nsegments
          if (pbjact(i)==0) cycle                       ! Skip inactive segments
C-----For each segment, add to AMAT and RHS for each triangle connection
          do j=1,3
            bwe(1) = bw(j  ,i)
            bwe(2) = bw(j+3,i)            
            n = segnodes(j,i)                        ! Current node
            do k=1,3
C-----Add giant conductance to AMAT
              segfunc = (bignum*bwe(1)*bw(k,i)+      ! Segment function evaluation
     1                   bignum*bwe(2)*bw(k+3,i))/2  ! numerically integrated
              ija = segIA(k,j,i)
              AMAT(ija) = AMAT(ija) - segfunc
            end do
C-----Add giant conductance to RHS, with specified head
            RHS(n) = RHS(n)-(bignum*bwe(1)*seghead(1,i)
     1                     + bignum*bwe(2)*seghead(2,i))/2
          end do
        end do ! end segment loop
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C-----PBJ Mode 1 - Drain or External Heads
      else if (pbjmode == 1) then
C-----Loop over segments
        do i=1, nsegments          
          if (pbjact(i)==0) cycle                       ! Skip inactive segments
C-----Interpolate to get current segment heads (start/end)
          heads = zero
          do j=1,3
            heads(1) = heads(1) + HNEW(segnodes(j,i)) * bw(j  ,i)
            heads(2) = heads(2) + HNEW(segnodes(j,i)) * bw(j+3,i)
          end do
C-----Move to next segment if heads are both below elevation
          if ((heads(1) <= segelevs(1,i)).and. 
     1        (heads(2) <= segelevs(2,i))) cycle
C-----Calculate "effective" weights (0 if head below elevation)
          do j=1,3
            debug = zero
            bwe(1) = bw(j  ,i)
            bwe(2) = bw(j+3,i)
            if (heads(1) <= segelevs(1,i)) bwe(1) = zero
            if (heads(2) <= segelevs(2,i)) bwe(2) = zero
C-----For each barycentric coordinate add conductances to nodes in AMAT
            n = segnodes(j,i)                           ! Current node
            do k=1,3
              segfunc = (cond(1,i)*bwe(1)*bw(k,i)+      ! Segment function evaluation
     1                   cond(2,i)*bwe(2)*bw(k+3,i))/2  ! numerically integrated
              ija = segIA(k,j,i)
              AMAT(ija) = AMAT(ija) - segfunc
              !debug = debug + segfunc * HNEW(segnodes(k,i))  ! FIX COND
            end do
C-----Add terms to RHS
              RHS(N) = RHS(N)-(cond(1,i)*bwe(1)*segelevs(1,i)
     1                       + cond(2,i)*bwe(2)*segelevs(2,i))/2
          end do
        end do
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C-----PBJ Mode 2 - External Heads (River condition)
      else
C-----Loop over segments
        do i=1, nsegments          
          if (pbjact(i)==0) cycle                       ! Skip inactive segments
C-----Interpolate to get current segment heads (start/end)
          heads = zero
          do j=1,3
            heads(1) = heads(1) + HNEW(segnodes(j,i)) * bw(j  ,i)
            heads(2) = heads(2) + HNEW(segnodes(j,i)) * bw(j+3,i)
          end do
          do j=1,3
            bwe(1) = bw(j  ,i)
            bwe(2) = bw(j+3,i)
C-----For each barycentric coordinate add conductances to nodes in AMAT
            n = segnodes(j,i)                           ! Current node
C-----Compare average gw head to average stream elevation, if lower then percolation conditions
            if (sum(heads(1:2))/2 >= sum(segelevs(1:2,i))/2) then
              do k=1,3
                segfunc = (cond(1,i)*bwe(1)*bw(k,i)+      ! Segment function evaluation
     1                     cond(2,i)*bwe(2)*bw(k+3,i))/2  ! numerically integrated
                ija = segIA(k,j,i)
                AMAT(ija) = AMAT(ija) - segfunc
              end do
C-----Add terms to RHS
              RHS(N) = RHS(N)-((cond(1,i)*bwe(1)*seghead(1,i)
     1                        + cond(2,i)*bwe(2)*seghead(2,i)))/2
            else
C-----Percolation - only add to RHS
              RHS(N) = RHS(N)-(cond(1,i)*bwe(1)*
     1                                    (seghead(1,i)-segelevs(1,i))
     2                        + cond(2,i)*bwe(2)*
     3                                    (seghead(2,i)-segelevs(2,i))
     4                        )/2
            end if
          end do
        end do
C-----------------------------------------------------------------------
      end if
      
      return
      end subroutine GWF2PBJU1FM

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1BD(KSTP,KPER)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR POLYLINE BOUNDARY JUNCTIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,IUNSTR,
     *                 NODES,NEQS,INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      use pbjmodule
C
      character*16          :: TEXT
      integer               :: i, j
      double precision      :: Q, RATOUT, ROUT, heads(2), avgover
C
      DATA TEXT /'    PBJ SEGMENTS'/
C     ------------------------------------------------------------------
      
C-----TODO: Handle cell-by-cell flow calculations/output
      
C-----No segments, no service
      if (nsegments <= 0) return
      
      RATOUT = zero
C-----------------------------------------------------------------------
C-----PBJ Mode 0 - Specified heads
      if (pbjmode == 0) then
        ! Calculated in GWF2PBJU1BDCHWR using FLOWJA
        ! No need to add to budget terms here
        return
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C-----PBJ Mode 1 - Drain or External Heads
      else if (pbjmode == 1) then
C-----Loop over segments calculate total segment flow
      do i=1, nsegments
C-----Interpolate heads to segment start/end
        heads = zero
        do j=1,3
          heads(1) = heads(1) + HNEW(segnodes(j,i)) * bw(j  ,i)
          heads(2) = heads(2) + HNEW(segnodes(j,i)) * bw(j+3,i)
        end do
C-----Move to next segment if heads are both below elevation
        if ((heads(1) <= segelevs(1,i)).and. 
     1      (heads(2) <= segelevs(2,i))) cycle
C-----Calculate flow out of segment, add to cumulative count
        heads(1) = max(heads(1)-segelevs(1,i), zero)
        heads(2) = max(heads(2)-segelevs(2,i), zero)
        Q = (cond(1,i) * heads(1) + cond(2,i) * heads(2))/2        ! Sign convention: positive (+) = leaving aquifer
        RATOUT = RATOUT + Q
      end do
C-----------------------------------------------------------------------      

C-----------------------------------------------------------------------      
C-----PBJ Mode 2 - External Heads (River condition)
      else
C-----Loop over segments calculate total segment flow
      do i=1, nsegments
C-----Interpolate heads to segment start/end
        heads = zero
        do j=1,3
          heads(1) = heads(1) + HNEW(segnodes(j,i)) * bw(j  ,i)
          heads(2) = heads(2) + HNEW(segnodes(j,i)) * bw(j+3,i)
        end do
C-----Compare average gw head to average stream elevation, if lower then percolation conditions
        if (sum(heads(1:2))/2 >= sum(segelevs(1:2,i))/2) then  ! GW Head dependent, RATE=CRIV*(HRIV-HNEW)
          heads(1) = seghead(1,i) - heads(1)
          heads(2) = seghead(2,i) - heads(2)
        else                                                   ! Percolating, RATE=CRIV*(HRIV-RBOT)
          heads(1) = seghead(1,i) - segelevs(1,i)
          heads(2) = seghead(2,i) - segelevs(2,i)
        end if
        Q = (cond(1,i) * heads(1) + cond(2,i) * heads(2))/2
        RATOUT = RATOUT + Q
      end do
      
C-----------------------------------------------------------------------
      end if

C-----MOVE RATES,VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
      ROUT=RATOUT
      VBVL(3,MSUM)=zero
      VBVL(4,MSUM)=ROUT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
      
C-----INCREMENT BUDGET TERM COUNTER.
      MSUM=MSUM+1
      return
      end subroutine GWF2PBJU1BD

C -----------------------------------------------------------------------------

      SUBROUTINE GWF2PBJU1DA
C-----Deallocate PBJ MEMORY
      USE pbjmodule
C
        deallocate(segnodes, bw, segelevs, cond)
C
      RETURN
      END

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1R(rlist,ncol,nlst,label,INPACK,IOUT,IFREFM)
C     ******************************************************************
C     Read in a list of values to an array, intended for
C     reading stress period arrays in the PBJ package.
C     Arguments:
C      - rlist array to be read into (ncol,nlst) (real*8)
C      - ncol  how many "columns" of data to be read in
C      - nlst  number of values to be read in (e.g. itmp)
C      - label name of the data type being read in, be to printed in LST
C      - IN, IOUT, and IFREFM globals
C     Note: write statements assume ncol < 10
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE pbjmodule
      integer,intent(in)               :: nlst,ncol,INPACK,IOUT,IFREFM
      integer                          :: i, j, iseg, LLOC, ICLOSE
      real*8,intent(inout)             :: rlist(ncol,nlst)
      real                             :: istream(ncol), SFAC
      character*(*),intent(in)         :: label
      character*200                    :: line, outline, FNAME
      INCLUDE 'openspec.inc'
C     ------------------------------------------------------------------
      
C-----Check for and decode EXTERNAL and OPEN/CLOSE records.
      IN=INPACK
      ICLOSE=0
      READ(IN,'(A)') LINE
      SFAC=1.
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
         IN=I
         IF(IPRFLG.EQ.1)WRITE(IOUT,111) IN
  111    FORMAT(1X,'Reading list on unit ',I4)
         READ(IN,'(A)') LINE
      ELSE IF(LINE(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=LINE(ISTART:ISTOP)
         IN=NUNOPN
         IF(IPRFLG.EQ.1)WRITE(IOUT,115) IN,trim(FNAME)
  115    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         OPEN(UNIT=IN,FILE=trim(FNAME),ACTION=ACTION(1))
         ICLOSE=1
         READ(IN,'(A)') LINE
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CONSTANT') THEN
        if (nlst /= nsegments) then            ! ERROR
          write(IOUT,*) ' Cannot set values to CONSTANT unless all
     1 segments are active'
          call USTOP(' ')
        end if
C-----If constant, read, save, LEAVE
        do j=1, ncol
          call URWORD(LINE,LLOC,ISTART,ISTOP,3,I,istream(j),IOUT,IN)
          rlist(j,:) = istream(j)
        end do
        if(IPRPBJ==1) then
          write(IOUT,117) trim(label),(istream(j)*SFAC, j=1, ncol)
 117      format(1x,a,' CONSTANT ',9F10.0)
        end if
        pbjact = 1
        return
      END IF
C
C3------Check for SFAC record.
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SFAC') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SFAC,IOUT,IN)
         IF(IPRFLG.EQ.1) THEN
           WRITE(IOUT,116) SFAC
  116      FORMAT(1X,'LIST SCALING FACTOR=',1PG12.5)
         ENDIF
         READ(IN,'(A)') LINE
      END IF
      
C-----Write out label to list file
      if (IPRPBJ==1) then
        outline = ''
        write(IOUT,1) trim(label), ' PER SEGMENT'
        write(IOUT,2) 'SEGMENT NO. ', (i, i=1, ncol)
 1      format(/,16X,2A)
 2      format(4X,A,I14,8I16)
      end if
      
C-----Read in values
      do i=1, nlst
        if (i>1) read(IN,'(A)') line            ! First line was read & scanned for flags
        if(IFREFM.EQ.0) then
          read(line,'(1I10,9F10.0)') iseg, (istream(j), j=1, ncol)
        else
          LLOC=1
          call URWORD(LINE,LLOC,ISTART,ISTOP,2,iseg,R,IOUT,IN)
          do j=1, ncol
            call URWORD(LINE,LLOC,ISTART,ISTOP,3,I,istream(j),IOUT,IN)
          end do
        end if
C-----Check for issues
        if ((iseg < 0).or.(iseg > nsegments)) then
          write(IOUT,*) ' Segment no. higher than maximum'
          call USTOP(' ')
        end if
C-----Put into proper array location
        do j=1, ncol
          pbjact(iseg) = 1
          rlist(j,iseg) = istream(j)*SFAC
        end do
C-----Write to LST, if writing
        if (IPRPBJ==1) then
          write(IOUT,3) iseg, (istream(j)*SFAC, j=1, ncol)
 3        format(4X,1I10,9E16.4)
        end if
      end do

C-----Done reading the list.  If file is open/close, close it.
      IF(ICLOSE.NE.0) CLOSE(UNIT=IN)
      
      return
      end subroutine GWF2PBJU1R

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1BDCHWR(KKSTP,KKPER)
C     ******************************************************************
C     Computes specified-head ("constant head") budget termf
C     LS - Not currently run in mfusg.f
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE pbjmodule
      USE GLOBAL, only:IBOUND,IA,FLOWJA
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      integer                          :: i, j, n, nconn
      double precision                 :: flow, rate, pbjin, pbjout,
     1                                    cin, cout, nsum
      character*16          :: TEXT
      DATA TEXT /'  PBJ SPEC HEADS'/
C     ------------------------------------------------------------------
      
C-----Only run if pbj is in specified head mode
      if (pbjmode /= 0) return
      
      pbjin  = zero
      pbjout = zero
      rate   = zero
      
C-----TODO: Handle cell-by-cell flow calculations/output      
      
C-----Loop over PBJ segments, calculate node water balances
C-----Assumption: all gains/losses due to PBJ
      do i=1, nsegments
        do j=1, 3
          n = segnodes(j,i)
          nsum = 0
          if (IBOUND(n) < 1) cycle
C-----Loop over all node connections
          do nconn = IA(N)+1,IA(N+1)-1
            flow = FLOWJA(nconn)
            nsum = nsum + flow
            rate = rate + flow
          end do
          if (nsum < 0) then
            pbjout = pbjout - nsum
          else
            pbjin = pbjin + nsum
          end if
        end do
      end do
      
C-----Add to GWF Budget
      CIN=pbjin
      COUT=pbjout
      VBVL(1,MSUM)=VBVL(1,MSUM)+CIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+COUT*DELT
      VBVL(3,MSUM)=CIN
      VBVL(4,MSUM)=COUT
      VBNM(MSUM)=TEXT
      MSUM=MSUM+1

      return
      end subroutine GWF2PBJU1BDCHWR

C -----------------------------------------------------------------------------

C *****************************************************************************
C -- PBJ Helper Routines
C -----------------------------------------------------------------------------
      
      subroutine initialize_pbj_arrays()
        use pbjmodule
        allocate(segnodes   (3, nsegments),
     1           segIA      (3, 3, nsegments),
     2           bw         (6, nsegments),
     3           segelevs   (2, nsegments),
     4           seglens    (nsegments),
     5           cond       (2,nsegments),
     6           pbjact     (nsegments))
        if (pbjmode /= 1) then
          allocate(seghead(2, nsegments))
        end if
      end subroutine initialize_pbj_arrays

C -----------------------------------------------------------------------------
