C *****************************************************************************
C -- Polyline Boundary Junction (PBJ) Package
C ---- by Leland Scantlebury and James R. Craig
C -----------------------------------------------------------------------------
      module pbjmodule
        integer :: nsegments
        integer, allocatable, dimension(:,:) :: segnodes    ! Nodes to which flow is interpolated, for each segment
        real,    allocatable, dimension(:,:) :: baryweights ! Barycentric coordinates for segment start/end points
        real,    allocatable, dimension(:,:) :: segelevs    ! Segment start/end point streambed elevations
        real,    allocatable, dimension(:)   :: cond        ! Segment conductivity
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
      USE GLOBAL,       ONLY:IOUT
      use pbjmodule
      integer  :: i, j, k
      integer, allocatable, dimension(:) :: nodetemp
      real, allocatable, dimension(:) :: weighttemp, elevtemp
      character*200 line
      
      CHARACTER*24 ANAME(3)
      DATA ANAME(1) /'           SEGMENT NODES'/
      DATA ANAME(2) /'     BARYCENTRIC WEIGHTS'/
      DATA ANAME(3) /'      SEGMENT ELEVATIONS'/
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
      call initialize_arrays()
      
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
      call U1DREL(elevtemp,ANAME(2),nsegments*2,K,IN,IOUT)

C-----Move temporary arrays to their final resting places
      do i=1, nsegments
          do j=1, 3
              segnodes(i,j) = nodetemp((i-1)*3 + j)
              do k=1, 2
                  baryweights(i,j*k) = weighttemp((i-1)*6 + j*k)
              end do
          end do
          do j=1, 2
              segelevs(i,j) = elevtemp((i-1)*2 + j)
          end do
      end do
      
C-----Deallocate
      deallocate(nodetemp, weighttemp, elevtemp)
      
C-----Return
      return
      end subroutine GWF2PBJU1AR

C -----------------------------------------------------------------------------

      subroutine GWF2PBJU1RP(IN)
C     ******************************************************************
C     Read segment conductances
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      use pbjmodule
      USE GLOBAL,       ONLY:IOUT
C     ------------------------------------------------------------------
      
      write(*,*) 'YES'
      
      end subroutine GWF2PBJU1RP

      
C *****************************************************************************
C -- PBJ Helper Routines
C -----------------------------------------------------------------------------
      
      subroutine initialize_arrays()
      use pbjmodule
C     Package currently assumes voronoi grid with triangular interpolation
C     Thus, three nodes per segment
        allocate(segnodes   (nsegments, 3),
     1           baryweights(nsegments, 6),
     2           segelevs   (nsegments, 2),
     3           cond       (nsegments))
      end subroutine initialize_arrays

C -----------------------------------------------------------------------------

