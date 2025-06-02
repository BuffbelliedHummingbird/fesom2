subroutine compute_tracer_elems(tr_nod,tr_elem)

    ! USE o_ARRAYS, &
    !  ONLY: UV, Unode
    USE g_PARSUP, &
       ONLY: myDim_elem2D, eDim_elem2D, &
             myDim_nod2D, eDim_nod2D
    USE mod_assim_pdaf, &
       ONLY: mesh_fesom
    USE g_comm_auto                    ! contains: interface exchange_nod()

    IMPLICIT NONE

    ! *** local variables ***
    integer          :: elem, nz, &    ! counters
                        ule, &         ! upper vertical boundary index at element, default==1
                        nle            ! number of levels as element considering bottom topography
                        
    ! *** arguments ***
    real, intent(inout) :: tr_nod (mesh_fesom%nl-1, myDim_nod2D)   ! tracer at nodes
    real, intent(out)   :: tr_elem(mesh_fesom%nl-1,myDim_elem2D)   ! tracer at elements
    
    ! *** variables imported from modules ***
    ! ulevels(:)               ! upper boundary index of all vertical element loops, default==1
    ! nlevels(:)               ! number of levels at elements considering bottom topography
    ! elem2D_nodes(:,e)        ! 3 nodes of element e
    
    ! #include "associate_mesh.h"
    
    ! In FESOM2.1, the velocities are defined at the centroids of the triangular elements.
    ! Note on trigonometry:
    ! 1. The centroid of a triangle is at the intersection of its medians.
    ! 2. The medians of a triangle divide the triangle into 6 small triangles with equal area.
    ! Thus, each of the 3 vertices is associated with an equal area share of the element
    ! and weighted equally for interpolation.
    
    elements: DO elem = 1, myDim_elem2D
    
        ule = mesh_fesom% ulevels(elem)
        nle = mesh_fesom% nlevels(elem)
        
        watercolumn: DO nz = ule, nle-1

        tr_elem(nz,elem) = (  tr_nod(nz,mesh_fesom% elem2D_nodes(1,elem)) &
                            + tr_nod(nz,mesh_fesom% elem2D_nodes(2,elem)) &
                            + tr_nod(nz,mesh_fesom% elem2D_nodes(3,elem)))  /3
        
        END DO watercolumn
    END DO elements
    
end subroutine compute_tracer_elems
