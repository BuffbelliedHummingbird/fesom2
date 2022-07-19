subroutine compute_vel_elems()

    ! USE MOD_MESH
    ! USE o_PARAM
    USE o_ARRAYS, &
       ONLY: UV, Unode
    USE g_PARSUP, &
       ONLY: myDim_elem2D
    USE mod_assim_pdaf, &
       ONLY: mesh_fesom

    IMPLICIT NONE

    ! *** local variables ***
    integer          :: elem, nz, &    ! counters
                        ule, &         ! upper vertical boundary index at element, default==1
                        nle            ! number of levels as element considering bottom topography
    
    ! *** variables imported from modules ***
    ! ulevels(:)               ! upper boundary index of all vertical element loops, default==1
    ! nlevels(:)               ! number of levels at elements considering bottom topography
    ! elem2D_nodes(:,e)        ! 3 nodes of element e
    ! UV (2,nz,elem)           ! u and v velocities defined at element centroids
    ! Unode(2,nz,node)         ! u and v velocities defined at nodes
    
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

        UV(1,nz,elem) = (  Unode(1,nz,mesh_fesom% elem2D_nodes(1,elem)) &
                         + Unode(1,nz,mesh_fesom% elem2D_nodes(2,elem)) &
                         + Unode(1,nz,mesh_fesom% elem2D_nodes(3,elem)))  /3
        UV(2,nz,elem) = (  Unode(2,nz,mesh_fesom% elem2D_nodes(1,elem)) &
                         + Unode(2,nz,mesh_fesom% elem2D_nodes(2,elem)) &
                         + Unode(2,nz,mesh_fesom% elem2D_nodes(3,elem)))  /3
        
        END DO watercolumn
    END DO elements
    
end subroutine compute_vel_elems
