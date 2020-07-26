C *****************************************************************************
C -- Polyline Boundary Junction (PBJ) Package
C ---- by Leland Scantlebury and James R. Craig
C -----------------------------------------------------------------------------
      module pbjmodule
        integer :: nsegments
        integer :: IPBJCB                                   ! Write CBC flag
        integer :: IPRPBJ=1                                 ! LST Printout flag (hardcoded, for now)
        integer :: pbjmode                                  ! 0-spec head, 1-drain, 2-head dep
        integer :: condtype                                 ! 0-conductivity, 1-unit conductivity (per len) 2-leakance coeff
                                                            ! (-1 == not read, head spec mode)
        integer, allocatable, dimension(:,:)   :: segnodes  ! Nodes to which flow is interpolated, for each segment
        integer, allocatable, dimension(:,:,:) :: segIA     ! Off-diagonal location in AMAT of nodes listed in segia
        real*8,  allocatable, dimension(:,:)   :: bw        ! Barycentric coordinates for segment start/end points
        real*8,  allocatable, dimension(:,:)   :: segelevs  ! Segment start/end point streambed elevations
        real  ,  allocatable, dimension(:)     :: seglens   ! Segment lengths (for unit conductivity, etc)
        real,    allocatable, dimension(:,:)   :: cond      ! Segment conductivity
        character(len=16),    dimension(0)    :: PBJAUX    ! Junk
        double precision, parameter :: zero=0.0D0
        
      contains
C -----------------------------------------------------------------------------
        integer function find_n2m_inJA(n, m) result(ija)
          use GLOBAL,     ONLY:IA,JA
          integer, intent(in)       :: n,m
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
            write(IOUT,*) 'SEGMENT: ', i, 'NODES: ', m, n
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
      DATA ANAME(4) /'     SEGMENT CONDUCTANCE'/
C     ------------------------------------------------------------------      
      
C-----Identify package, initialize variables
      write(IOUT,1)in
    1 format(1X,/1X,'PBJ -- POLYLINE BOUNDARY JUNCTION PACKAGE,',
     1' VERSION 1, 7/14/2020 INPUT READ FROM UNIT ',I4)
      nsegments = 0
      
C-----Read in number of segments and IPBJCB flag
      call URDCOM(IN,IOUT,line)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,nsegments,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IPBJCB,R,IOUT,IN)
      
C-----Read in Mode and Conductance Type
      call URDCOM(IN,IOUT,line)
      call URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
C-----Identify Mode
      if (LINE(ISTART:ISTOP)=='HEADSPEC') then
        pbjmode = 0
        condtype = -1
      else if (LINE(ISTART:ISTOP)=='DRAIN') then
        pbjmode = 1
      else if (LINE(ISTART:ISTOP)=='EXTHEAD') then
        pbjmode = 2
      else
        write(IOUT,*) ' INVALID PBJ PACKAGE MODE'
        call USTOP(' ')
      end if
C-----Identify Conductivity Type
      if (pbjmode /= 0) then
        call URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
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
      allocate(elevtemp(nsegments*2))           ! 1 elevation for every segment start/end
      call U1DREL(elevtemp,ANAME(3),nsegments*2,K,IN,IOUT)
      
C-----Read conductances (currently do not vary by stress period)
CLS      call U1DREL(cond,ANAME(4),nsegments,K,IN,IOUT)

C-----Move temporary arrays to their final resting places
      do i=1, nsegments
        do j=1, 3                               ! Loop over nodes
          segnodes(i,j) = nodetemp((i-1)*3 + j)
          n = segnodes(i,j)                     ! current diagonal node
          bw(i,j)   = weighttemp((i-1)*6 + j)
          bw(i,j+3) = weighttemp((i-1)*6 + j+3)
C-----Find connected nodes in JA (needed later for entering into AMAT)
          do k=1, 3
            m = nodetemp((i-1)*3 + k)           ! current off-diag node
            segIA(i,j,k) = find_n2m_inJA(n,m)
          end do
        end do
        do j=1, 2                               ! Loop over segment ends
          segelevs(i,j) = elevtemp((i-1)*2 + j)
        end do
      end do
      
C-----Deallocate
      deallocate(nodetemp, weighttemp, elevtemp)
      
C-----Return
      return
      end subroutine GWF2PBJU1AR

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1RP(IN)  ! CURRENTLY NOT CALLED IN MFUSG.F (Because it doesn't do anything)
C     ******************************************************************
C     Read segment conductances
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      use pbjmodule
      USE GLOBAL,       ONLY:IOUT,IFREFM,NEQS
      integer :: itmp
C     ------------------------------------------------------------------
C-----Identify package
      write(IOUT,1)in
    1 format(1X,/1X,'PBJ -- POLYLINE BOUNDARY JUNCTION PACKAGE,',
     1' VERSION 1, 7/14/2020 INPUT READ FROM UNIT ',I4)
      
C-----Read ITMP (flag to re-use data or not)
      if (IFREFM.EQ.0) then
        read(IN,'(I10)') itmp
      else
        read(IN,*) itmp
      end if

C-----Process itmp
      IF(itmp < 0) THEN
        WRITE(IOUT,7)
    7   FORMAT(1X,/1X, 'REUSING PBJ VALUES FROM LAST STRESS PERIOD')
C-----Nothing else to do!
        return
      END IF
      
C-----If using conductances, read in conductances
      if ((itmp >= 0).and.(condtype == 0)) then
!      SUBROUTINE ULSTRDU(NLIST,RLIST,LSTBEG,LDIM,MXLIST,IAL,INPACK,IOUT,
!     1     LABEL,CAUX,NCAUX,NAUX,IFREFM,NODES,ISCLOC1,ISCLOC2,
!     2     IPRFLG)
        CALL ULSTRDU(itmp,cond,1,nsegments,2,1,IN,IOUT,
     1     'SEGMENT        NODE         CONDUCTANCE1  CONDUCTANCE2',
     2     PBJAUX,0,0,IFREFM,NEQS,3,4,IPRPBJ)
      else
        write(*,*) "OH GOD IMPLEMENT MORE"
      end if
      
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
      double precision   :: heads(2), bwe(2), hw, debug, estflow
C     ------------------------------------------------------------------
      
C-----If no segments, return
      if (nsegments <= 0) return
      
C-----Loop over segments
      do i=1, nsegments
C-----Interpolate to get current segment heads (start/end)
        heads = zero
        do j=1,3
          heads(1) = heads(1) + HNEW(segnodes(i,j)) * bw(i,j  )
          heads(2) = heads(2) + HNEW(segnodes(i,j)) * bw(i,j+3)
        end do
C-----Move to next segment if heads are both below elevation
        if ((heads(1) <= segelevs(i,1)).and. 
     1      (heads(2) <= segelevs(i,2))) cycle
C-----Calculate "effective" weights (0 if head below elevation)
        do j=1,3
          debug = zero
          bwe(1) = bw(i,j)
          bwe(2) = bw(i,j+3)
          if (heads(1) <= segelevs(i,1)) bwe(1) = zero
          if (heads(2) <= segelevs(i,2)) bwe(2) = zero
C-----For each barycentric coordinate add conductances to nodes in AMAT
          n = segnodes(i,j)                           ! Current node
          do k=1,3
            hw = (bwe(1)*bw(i,k)+bwe(2)*bw(i,k+3))/2  ! head "weight"
            ija = segIA(i,j,k)
            AMAT(ija) = AMAT(ija) - cond(i,1)*hw  ! FIX COND
            debug = debug + cond(i,1) * hw * HNEW(segnodes(i,k))  ! FIX COND
          end do
C-----Add non-head related terms to RHS
          estflow = -cond(i,1)*((bwe(1)*segelevs(i,1)  ! FIX COND
     1                            + bwe(2)*segelevs(i,2)))/2
          RHS(N) = RHS(N)-cond(i,1)*((bwe(1)*segelevs(i,1)  ! FIX COND
     1                            + bwe(2)*segelevs(i,2)))/2
!          write(*,*) '     Segment ', i
!          write(*,*) '     Heads = ', heads(1), heads(2)
!          write(*,*) 'Total AMAT = ', debug
!          write(*,*) '       RHS = ', estflow
!          write(*,*) '      Flow = ', debug + estflow
        end do
      end do      
      
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
      
C-----Loop over segments calculate total segment flow
      RATOUT = zero
      do i=1, nsegments
C-----Interpolate heads to segment start/end
        heads = zero
        do j=1,3
          heads(1) = heads(1) + HNEW(segnodes(i,j)) * bw(i,j  )
          heads(2) = heads(2) + HNEW(segnodes(i,j)) * bw(i,j+3)
        end do
C-----Move to next segment if heads are both below elevation
        if ((heads(1) <= segelevs(i,1)).and. 
     1      (heads(2) <= segelevs(i,2))) cycle
C-----Calculate flow out of segment, add to cumulative count
        heads(1) = max(heads(1)-segelevs(i,1), zero)
        heads(2) = max(heads(2)-segelevs(i,2), zero)
        avgover = (heads(1)+heads(2))/2
        Q = cond(i,1) * avgover               ! Sign convention: positive (+) = leaving aquifer  ! FIX COND
        RATOUT = RATOUT + Q
      end do
      
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
      
C *****************************************************************************
C -- PBJ Helper Routines
C -----------------------------------------------------------------------------
      
      subroutine initialize_pbj_arrays()
        use pbjmodule
        allocate(segnodes   (nsegments, 3),
     1           segIA      (nsegments, 3, 3),
     2           bw         (nsegments, 6),
     3           segelevs   (nsegments, 2),
     4           cond       (nsegments, 2))
      end subroutine initialize_pbj_arrays

C -----------------------------------------------------------------------------
