! PDAF in AWI-CM2 / Fesom 2.0

! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

module g_comm_pdaf

  use, intrinsic :: ISO_FORTRAN_ENV
  use mod_parallel_pdaf, ONLY: comm_filter, mype_filter

  implicit none

contains

!#ifdef DEBUG
!  ! Only needed in debug mode
!  subroutine check_mpi_comm(rn, sn, r_mpitype, s_mpitype, rPE, sPE)
!    USE g_PARSUP
!    IMPLICIT NONE
!
!    ! General version of the communication routine for 2D nodal fields
!
!    integer, intent(in) :: sn, rn, r_mpitype(:), s_mpitype(:), rPE(:), sPE(:)
!    integer  :: n, sdebug, rdebug, status(MPI_STATUS_SIZE), request
!
!    DO n=1,rn    
!       call MPI_TYPE_SIZE(r_mpitype(n), rdebug, MPIerr)
!       CALL MPI_ISEND(rdebug, 1, MPI_INTEGER, rPE(n), 10, COMM_FILTER, request, MPIerr)
!    END DO
!
!    DO n=1, sn
!       call MPI_RECV(sdebug, 1, MPI_INTEGER, sPE(n), 10, COMM_FILTER,    &
!            status, MPIerr)
!       call MPI_TYPE_SIZE(s_mpitype(n), rdebug, MPIerr)
!       if (sdebug /= rdebug) then
!          print *, "Mismatching MPI send/recieve message lengths."
!          print *,"Send/receive process numbers: ", mype, '/', sPE(n)
!          print *,"Number of send/receive bytes: ", sdebug, '/', rdebug
!          call MPI_ABORT( COMM_FILTER, 1 )
!       end if
!    END DO
!    CALL MPI_BARRIER(COMM_FILTER,MPIerr)
!
!  END SUBROUTINE check_mpi_comm
!#endif

!============================================================================
!
subroutine gather_nod3D_pdaf(arr3D, arr3D_global)

! Make nodal information available to master PE 
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

integer      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real64) ::  arr3D_global(:,:)
real(real64), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to receive directly into arr3D_global

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_DOUBLE_PRECISION, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

end if
end subroutine gather_nod3D_pdaf
!_pdaf
!============================================================================
!
subroutine gather_real4_nod3D_pdaf(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real32)  ::  arr3D(:,:)
real(real32)  ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

end if
end subroutine gather_real4_nod3D_pdaf
!=======================================================

subroutine gather_int2_nod3D_pdaf(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

integer(int16) ::  arr3D(:,:)
integer(int16) ::  arr3D_global(:,:)
integer(int16), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_SHORT, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_nod2D*nl1, MPI_SHORT, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

end if
end subroutine gather_int2_nod3D_pdaf


!==============================================
subroutine gather_nod2D_pdaf(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real64) ::  arr2D_global(:)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))

      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_DOUBLE_PRECISION, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_nod2D_pdaf
!==============================================
subroutine gather_real4_nod2D_pdaf(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_real4_nod2D_pdaf

!==============================================
subroutine gather_int2_nod2D_pdaf(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

integer(int16) ::  arr2D(:)
integer(int16) ::  arr2D_global(:)
integer(int16), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_SHORT, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_SHORT, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_int2_nod2D_pdaf

!============================================================================
subroutine gather_elem3D_pdaf(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real64) ::  arr3D_global(:,:)
real(real64), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_DOUBLE_PRECISION, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_DOUBLE_PRECISION, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_elem3D_pdaf

!===================================================================

subroutine gather_real4_elem3D_pdaf(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER        :: nl1
integer        :: n

real(real32)   ::  arr3D(:,:)
real(real32)   ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_real4_elem3D_pdaf


!===================================================================

subroutine gather_int2_elem3D_pdaf(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

integer(int16) ::  arr3D(:,:)
integer(int16) ::  arr3D_global(:,:)
integer(int16), allocatable :: recvbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D, ende, err_alloc
integer        :: max_loc_Dim, i, status(MPI_STATUS_SIZE)

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global
! (Carefull with duplicate interface elements, coming from two
!  PEs at once!)

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
! 
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_SHORT, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   

      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                        = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   
   call MPI_SEND( arr3D, myDim_elem2D*nl1, MPI_SHORT, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_int2_elem3D_pdaf


!==============================================
subroutine gather_elem2D_pdaf(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real64) ::  arr2D_global(:)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_DOUBLE_PRECISION, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_DOUBLE_PRECISION, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF
end if

end subroutine gather_elem2D_pdaf

!==============================================
subroutine gather_real4_elem2D_pdaf(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF
end if

end subroutine gather_real4_elem2D_pdaf

!==============================================
subroutine gather_int2_elem2D_pdaf(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

integer(int16) ::  arr2D(:)
integer(int16) ::  arr2D_global(:)
integer(int16), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_SHORT, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_SHORT, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF
end if

end subroutine gather_int2_elem2D_pdaf


!============================================================================
subroutine gather_real8to4_nod3D_pdaf(arr3D, arr3D_global)

! Make nodal information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64)  ::  arr3D(:,:)
real(real32)   ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
real(real32), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, n3D, ierr

 if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1, ubound(arr3D_global,2)))

      do  n = 1, npes-1
         n3D = (remPtr_nod2D(n+1) - remPtr_nod2D(n))*nl1
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(1,start), n3D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_nod2D(1:myDim_nod2D)) = arr3D(1:nl1,1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE

   allocate(sendbuf(nl1,myDim_nod2D))
   sendbuf(1:nl1,1:myDim_nod2D) = arr3D(1:nl1,1:myDim_nod2D)
   
   call MPI_SEND(sendbuf, myDim_nod2D*nl1, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   deallocate(sendbuf)
   
ENDIF

end if

end subroutine gather_real8to4_nod3D_pdaf
!==============================================
subroutine gather_real8to4_nod2D_pdaf(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64)  ::  arr2D(:)
real(real32)   :: arr2D_global(:)
real(real32)   :: sendbuf(myDim_nod2D)
real(real64), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

! Consider MPI-datatypes to recv directly into arr2D_global!

 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)
IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))

      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   sendbuf(1:myDim_nod2D) = real(arr2D(1:myDim_nod2D),real32)

   call MPI_SEND(sendbuf, myDim_nod2D, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

end if
end subroutine gather_real8to4_nod2D_pdaf
!==============================================
!============================================================================
subroutine gather_real8to4_elem3D_pdaf(arr3D, arr3D_global)

! Make element information available to master PE 
!
! Use only with 3D arrays stored in (vertical, horizontal) way

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH


IMPLICIT NONE

INTEGER      :: nl1
integer      :: n

real(real64) ::  arr3D(:,:)
real(real32)  ::  arr3D_global(:,:)
real(real32), allocatable :: recvbuf(:,:)
real(real32), allocatable :: sendbuf(:,:)
integer        :: req(npes-1)
integer        :: start, e3D


 if (npes> 1) then
CALL MPI_BARRIER(COMM_FILTER,MPIerr)
nl1=ubound(arr3D,1)

! Consider MPI-datatypes to recv directly into arr3D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(nl1,remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e3D = (remPtr_elem2D(n+1) - remPtr_elem2D(n))*nl1
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(1,start), e3D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr3D_global(1:nl1,myList_elem2D(1:myDim_elem2D)) = arr3D(1:nl1,1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr3D_global(1:nl1, remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                     = recvbuf(1:nl1, 1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else
      arr3D_global(:,:) = arr3D(:,:)
   endif

ELSE
   allocate(sendbuf(nl1,myDim_elem2D))
   sendbuf(1:nl1,1:myDim_elem2D) = arr3D(1:nl1,1:myDim_elem2D)
   
   call MPI_SEND(sendbuf, myDim_elem2D*nl1, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   deallocate(sendbuf)
ENDIF

end if
end subroutine gather_real8to4_elem3D_pdaf
!==============================================
subroutine gather_real8to4_elem2D_pdaf(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
real(real32)  :: sendbuf(myDim_elem2D)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)
! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   sendbuf(1:myDim_elem2D) = real(arr2D(1:myDim_elem2D),real32)
   call MPI_SEND(sendbuf, myDim_elem2D, MPI_REAL, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

end if
end subroutine gather_real8to4_elem2D_pdaf
!==============================================
subroutine gather_elem2D_i_pdaf(arr2D, arr2D_global)
! Make element information available to master PE 
  use g_PARSUP
  use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
  use o_MESH
  IMPLICIT NONE

  integer                       :: n
  integer                       :: arr2D(:)
  integer                       :: arr2D_global(:)
  integer, allocatable          :: recvbuf(:)
  integer                       :: req(npes-1)
  integer                       :: start, e2D
  CALL MPI_BARRIER(COMM_FILTER,MPIerr)
  ! Consider MPI-datatypes to recv directly into arr2D_global!
  IF ( mype_filter == 0 ) THEN
     if (npes > 1) then
        allocate(recvbuf(remPtr_elem2D(npes)))
        do  n = 1, npes-1
            e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
            start = remPtr_elem2D(n)
            call MPI_IRECV(recvbuf(start), e2D, MPI_INTEGER, n, 2, COMM_FILTER, req(n), MPIerr)
        enddo      
        arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
        call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
        arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                          = recvbuf(1 : remPtr_elem2D(npes)-1)
        deallocate(recvbuf)
     else
        arr2D_global(:) = arr2D(:)
     endif
  ELSE
     call MPI_SEND(arr2D, myDim_elem2D, MPI_INTEGER, 0, 2, COMM_FILTER, MPIerr )
  ENDIF
end subroutine gather_elem2D_i_pdaf
!============================================================================
subroutine gather_nod2D_i_pdaf(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH

IMPLICIT NONE

integer              :: n
integer              :: arr2D(:)
integer              :: arr2D_global(:)
integer, allocatable :: recvbuf(:)
integer              :: req(npes-1)
integer              :: start, n2D

if (npes> 1) then

CALL MPI_BARRIER(COMM_FILTER,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype_filter == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global, 1)))

      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_INTEGER, n, 2, COMM_FILTER, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)

      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_INTEGER, 0, 2, COMM_FILTER, MPIerr )
   
ENDIF

endif
end subroutine gather_nod2D_i_pdaf
!============================================================================
!
subroutine gather_edg2D_pdaf(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH
IMPLICIT NONE

real(real64) ::  arr2D(:)
real(real64) ::  arr2Dglobal(:)

integer                                  ::  i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  ibuf
REAL(real64), ALLOCATABLE, DIMENSION(:) ::  rbuf

IF ( mype_filter == 0 ) THEN
    arr2Dglobal(myList_edge2D(1:myDim_edge2D))=arr2D(1:myDim_edge2D)
    DO  n = 1, npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, COMM_FILTER, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(rbuf(buf_size), ibuf(buf_size))

       CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
                      1, COMM_FILTER, status, MPIerr )

       CALL MPI_RECV(rbuf(1), buf_size, MPI_DOUBLE_PRECISION, sender, &
                      2, COMM_FILTER, status, MPIerr )
       arr2Dglobal(ibuf)=rbuf
       DEALLOCATE(ibuf, rbuf)
    ENDDO
ELSE
    CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, COMM_FILTER, MPIerr )
    CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
                   COMM_FILTER, MPIerr )
    CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_DOUBLE_PRECISION, 0, 2,&
                   COMM_FILTER, MPIerr )
ENDIF
CALL MPI_BARRIER(COMM_FILTER,MPIerr)
end subroutine gather_edg2D_pdaf
!
!============================================================================
!
subroutine gather_edg2D_i_pdaf(arr2D, arr2Dglobal)
! A 2D version of the previous routine
use g_PARSUP
use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
USE o_MESH
IMPLICIT NONE

integer  ::  arr2D(:)
integer  ::  arr2Dglobal(:)

integer                                  ::  i, n, buf_size, sender, status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:)       ::  ibuf, vbuf

IF ( mype_filter == 0 ) THEN
    arr2Dglobal(myList_edge2D(1:myDim_edge2D))=arr2D(1:myDim_edge2D)
    DO  n = 1, npes-1
       CALL MPI_RECV( buf_size, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                      0, COMM_FILTER, status, MPIerr )
       sender = status(MPI_SOURCE)
       ALLOCATE(ibuf(buf_size), vbuf(buf_size))

       CALL MPI_RECV(ibuf(1), buf_size, MPI_INTEGER, sender, &
                      1, COMM_FILTER, status, MPIerr )

       CALL MPI_RECV(vbuf(1), buf_size, MPI_INTEGER, sender, &
                      2, COMM_FILTER, status, MPIerr )
       arr2Dglobal(ibuf)=vbuf
       DEALLOCATE(ibuf, vbuf)
    ENDDO
ELSE
    CALL MPI_SEND( myDim_edge2D, 1, MPI_INTEGER, 0, 0, COMM_FILTER, MPIerr )
    CALL MPI_SEND( myList_edge2D(1), myDim_edge2D, MPI_INTEGER, 0, 1, &
                   COMM_FILTER, MPIerr )
    CALL MPI_SEND( arr2D(1), myDim_edge2D, MPI_INTEGER, 0, 2,&
                   COMM_FILTER, MPIerr )
ENDIF
CALL MPI_BARRIER(COMM_FILTER,MPIerr)
end subroutine gather_edg2D_i_pdaf
!==============================================

end module g_comm_pdaf



module g_comm_auto_pdaf
use g_comm_pdaf
implicit none

interface gather_nod_pdaf
      module procedure gather_nod3D_pdaf
      module procedure gather_nod2D_pdaf
      module procedure gather_real4_nod3D_pdaf
      module procedure gather_real4_nod2D_pdaf
      module procedure gather_int2_nod3D_pdaf
      module procedure gather_int2_nod2D_pdaf
      module procedure gather_real8to4_nod3D_pdaf
      module procedure gather_real8to4_nod2D_pdaf
      module procedure gather_nod2D_i_pdaf
end interface gather_nod_pdaf

interface gather_elem_pdaf
      module procedure gather_elem3D_pdaf
      module procedure gather_elem2D_pdaf
      module procedure gather_real4_elem3D_pdaf
      module procedure gather_real4_elem2D_pdaf
      module procedure gather_int2_elem3D_pdaf
      module procedure gather_int2_elem2D_pdaf
      module procedure gather_real8to4_elem3D_pdaf
      module procedure gather_real8to4_elem2D_pdaf
      module procedure gather_elem2D_i_pdaf
end interface gather_elem_pdaf

interface gather_edge_pdaf
      module procedure gather_edg2D_pdaf
      module procedure gather_edg2D_i_pdaf
end interface gather_edge_pdaf


private  ! hides items not listed on public statement 
public :: gather_nod_pdaf, gather_elem_pdaf, gather_edge_pdaf
end module g_comm_auto_pdaf
