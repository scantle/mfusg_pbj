C *****************************************************************************
C -- Polyline Boundary Junction (PBJ) Package
C ---- by Leland Scantlebury and James R. Craig
C -----------------------------------------------------------------------------
      module pbjmodule
        integer :: nsegments
        integer, allocatable, dimension(:,:)   :: segnodes  ! Nodes to which flow is interpolated, for each segment
        integer, allocatable, dimension(:,:,:) :: segIA     ! Off-diagonal location in AMAT of nodes listed in segia
        real*8,  allocatable, dimension(:,:)   :: bw        ! Barycentric coordinates for segment start/end points
        real*8,  allocatable, dimension(:,:)   :: segelevs  ! Segment start/end point streambed elevations
        real*8,  allocatable, dimension(:)     :: cond      ! Segment conductivity
        
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
      
C-----Read in number of segments and options (TODO: Have options)
      call URDCOM(IN,IOUT,line)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,nsegments,R,IOUT,IN)
      
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
      call U1DREL(cond,ANAME(4),nsegments,K,IN,IOUT)

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
      USE GLOBAL,       ONLY:IOUT,IFREFM
      integer :: itmp
C     ------------------------------------------------------------------
cLS No reason to have any of the parameters vary by SP currently      
cC-----Identify package
c      write(IOUT,1)in
c    1 format(1X,/1X,'PBJ -- POLYLINE BOUNDARY JUNCTION PACKAGE,',
c     1' VERSION 1, 7/14/2020 INPUT READ FROM UNIT ',I4)
c      
cC-----Read ITMP (flag to re-use data or not)
c      if (IFREFM.EQ.0) then
c          read(IN,'(I10)') itmp
c      else
c          read(IN,*) itmp
c      end if
c      
cC-----
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
      double precision   :: heads(2), bwe(2), hw
      double precision, parameter :: zero=0.0D0
C     ------------------------------------------------------------------
      
C-----If no segments, return
      if (nsegments <= 0) return
      
C-----Loop over segments
      do i=1, nsegments
C-----Interpolate to get segment heads (start/end)
        heads = 0
        do j=1,3
          heads(1) = heads(1) + HNEW(segnodes(i,j)) * bw(i,j  )
          heads(2) = heads(2) + HNEW(segnodes(i,j)) * bw(i,j+3)
C-----Calculate "effective" weights (0 if head below elevation)
          bwe(1)   = bw(i,j)
          bwe(2) = bw(i,j+3)
          if (heads(1) <= segelevs(i,1)) bwe(1) = 0
          if (heads(2) <= segelevs(i,2)) bwe(2) = 0
C-----For each barycentric coordinate add conductances to nodes in AMAT
          n = segnodes(i,j)                     ! Current node
          do k=1,3
            hw = (bwe(1)*bw(i,k)+bwe(2)*bw(i,k+3))/2
          ! Figure out where to put it
            ija = segIA(i,j,k)
            AMAT(ija) = AMAT(ija) - cond(i)*hw
          end do
          
C-----Add to RHS
          RHS(N) = RHS(N)-cond(i)*(bwe(1)*segelevs(i,1)
     1                           + bwe(2)*segelevs(i,2))/2
        end do
      end do      
      
      end subroutine GWF2PBJU1FM
      
C *****************************************************************************
C -- PBJ Helper Routines
C -----------------------------------------------------------------------------
      
      subroutine initialize_pbj_arrays()
        use pbjmodule
C     Package currently assumes voronoi grid with triangular interpolation
C     Thus, three nodes per segment
        allocate(segnodes   (nsegments, 3),
     1           segIA      (nsegments, 3, 3),
     2           bw         (nsegments, 6),
     3           segelevs   (nsegments, 2),
     4           cond       (nsegments))
      end subroutine initialize_pbj_arrays

C -----------------------------------------------------------------------------

      subroutine destory_pbj_arrays()
        use pbjmodule
      
        deallocate(segnodes, bw, segelevs, cond)
      
      end subroutine destory_pbj_arrays

C -----------------------------------------------------------------------------
