! PDAF in AWI-CM2 / Fesom 2.0

! 2020-04 - Longjiang Mu - Initial commit for AWI-CM3

SUBROUTINE ignore_nod_pdaf()
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE mod_assim_pdaf, & 
      ONLY: depth_excl, depth_excl_no, mesh_fesom
  USE g_PARSUP, &
      ONLY: myDim_nod2D, MPI_COMM_FESOM, myList_edge2D, myDim_edge2D, myDim_elem2D

  INTEGER :: i, j, el, ed, elnodes(3) ! Counter
  CHARACTER(4) :: mype_string
  INTEGER, ALLOCATABLE :: depth_excl_tmp(:), depth_excl_tmp2(:) 

! nodes excluding
  ! get dimension
  depth_excl_no=0
! boundary
  DO ed=1, myDim_edge2D
    if(myList_edge2D(ed) > mesh_fesom%edge2D_in) then
      do i = 1, 2
        if (mesh_fesom%edges(i,ed)<=myDim_nod2D) then
          depth_excl_no = depth_excl_no + 1  
        endif
      enddo
    endif
  END DO
!  WRITE(mype_string,'(i4.4)') mype_filter
!  open(mype_filter,file='nodes_excl.'//mype_string)
!  write(mype_filter,*)mesh_fesom%bc_index_nod2D 
! connected to boundary
  do el = 1, myDim_elem2D
    elnodes = mesh_fesom%elem2D_nodes(:,el)
  
!  write(mype_filter,*) el, elnodes
    ! find nodes on mesh_fesom%edges that have nodes on boundary
    if (sum(mesh_fesom%bc_index_nod2D(elnodes)) <= 2 ) then
        do i = 1, 3
            if (elnodes(i)<=myDim_nod2D .and. mesh_fesom%bc_index_nod2D(elnodes(i))==1) then
                depth_excl_no = depth_excl_no + 1
            endif
        end do
    end if
  end do

!  close(mype_filter)
  !DO i=1, myDim_nod2d
  !  IF (mesh_fesom%geo_coord_nod2D(2,i)<-60. .and. mesh_fesom%depth(i) > -100.) THEN
  !    depth_excl_no = depth_excl_no + 1
  !  ENDIF
  !END DO
! mesh problem in mediterranean sea
  !DO i=1, myDim_nod2d
  !  IF (mesh_fesom%geo_coord_nod2D(2,i)>36. .and. mesh_fesom%geo_coord_nod2D(2,i)<40. .and.  &
  !      mesh_fesom%geo_coord_nod2D(1,i)>13. .and. mesh_fesom%geo_coord_nod2D(1,i)<18. ) THEN
  !    depth_excl_no = depth_excl_no + 1
  !  ENDIF
  !END DO

  IF (depth_excl_no==0) RETURN

  ALLOCATE(depth_excl_tmp2(depth_excl_no))
  ALLOCATE(depth_excl_tmp(depth_excl_no))
  depth_excl_tmp = 0
  depth_excl_tmp2 = 0

  ! get data
  depth_excl_no=0
  DO ed=1,myDim_edge2D
    if(myList_edge2D(ed) > mesh_fesom%edge2D_in) then
       do i = 1, 2
          if (mesh_fesom%edges(i,ed)<=myDim_nod2D) then
            depth_excl_no = depth_excl_no + 1              
            depth_excl_tmp(depth_excl_no) = mesh_fesom%edges(i,ed)
          endif
       end do
    endif
  END DO

  do el = 1, myDim_elem2D
    elnodes = mesh_fesom%elem2D_nodes(:,el)
    if (sum(mesh_fesom%bc_index_nod2D(elnodes)) <= 2 ) then
        do i = 1, 3
            if (elnodes(i)<=myDim_nod2D .and. mesh_fesom%bc_index_nod2D(elnodes(i))==1) then
                depth_excl_no = depth_excl_no + 1
                depth_excl_tmp(depth_excl_no)=elnodes(i)
            endif
        end do
    end if
  end do

  !DO i=1, myDim_nod2d
  !  IF (mesh_fesom%geo_coord_nod2D(2,i)<-60. .and. mesh_fesom%depth(i) > -100.) THEN
  !    depth_excl_no =  depth_excl_no + 1
  !    depth_excl_tmp(depth_excl_no) = i
  !  ENDIF
  !END DO

! mesh problem in mediterranean sea
  !DO i=1, myDim_nod2d
  !  IF (mesh_fesom%geo_coord_nod2D(2,i)>36. .and. mesh_fesom%geo_coord_nod2D(2,i)<40. .and.  &
  !      mesh_fesom%geo_coord_nod2D(1,i)>13. .and. mesh_fesom%geo_coord_nod2D(1,i)<18.) THEN
  !    depth_excl_no = depth_excl_no + 1
  !    depth_excl_tmp(depth_excl_no) = i
  !  ENDIF
  !END DO  
  
! remove duplicates
  IF (depth_excl_no>1) THEN 
      depth_excl_no=1
      depth_excl_tmp2(1) = depth_excl_tmp(1)
      outer: do i=2,size(depth_excl_tmp)
         do j=1,depth_excl_no
            if (depth_excl_tmp2(j) == depth_excl_tmp(i)) then
               cycle outer
            end if
         end do
         depth_excl_no = depth_excl_no + 1
         depth_excl_tmp2(depth_excl_no) = depth_excl_tmp(i)
      end do outer

      ! save to depth_excl
      ALLOCATE(depth_excl(depth_excl_no))
      do j=1,depth_excl_no
        depth_excl(j) = depth_excl_tmp2(j)
      end do

      call mysort
  ELSE 
      ALLOCATE(depth_excl(1))
      depth_excl(1)=depth_excl_tmp(1)
  ENDIF

  IF (ALLOCATED(depth_excl_tmp)) DEALLOCATE(depth_excl_tmp)
  IF (ALLOCATED(depth_excl_tmp2)) DEALLOCATE(depth_excl_tmp2)

!  WRITE(mype_string,'(i4.4)') mype_filter
!  open(mype_filter,file='nodes_excl.'//mype_string)
!  write(mype_filter,*) depth_excl
!  close(mype_filter)  
END SUBROUTINE ignore_nod_pdaf

subroutine mysort()
  USE mod_assim_pdaf, ONLY: depth_excl, depth_excl_no
  integer::  min,i,j,tmp

  do i=1,depth_excl_no-1
    min=i
    do j=i+1,depth_excl_no
        if (depth_excl(j) < depth_excl(min)) then
            min=j
        end if
    end do
    tmp = depth_excl(min)
    depth_excl(min) = depth_excl(i)
    depth_excl(i) = tmp 
  end do

end subroutine

