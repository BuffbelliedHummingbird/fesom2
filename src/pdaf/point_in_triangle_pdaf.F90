SUBROUTINE point_in_triangle_pdaf(el2D,   pt)
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

!//TODO: needs refinement
  use g_PARSUP
  use o_param
  use g_rotate_grid
  USE mod_assim_pdaf, ONLY: mesh_fesom
  implicit none

  INTEGER, INTENT(OUT)                   :: el2D
  REAL(kind=8), DIMENSION(2), INTENT(IN) :: pt

  integer                                :: elem, elnodes(3), q
  real(kind=8)                           :: alph, mean_lon, mean_lat, rlon, rlat
  real(kind=8)                           :: xe(4), ye(4), xt1, xt2, yt1, yt2, x1, x2, y1, y2
  real(kind=8)                           :: s1, s2, angle

  el2D=0

  DO elem=1, myDim_elem2D

     elnodes(1:3)=mesh_fesom%elem2D_nodes(:, elem)
     ! ========
     ! Very rough criteria to reduce the work
     ! ========
     if (rotated_grid) then
        do q=1,3
           rlon=mesh_fesom%coord_nod2D(1, elnodes(q))
           rlat=mesh_fesom%coord_nod2D(2, elnodes(q))
           call r2g(xe(q), ye(q), rlon, rlat)
        end do
     else
        write(*,*) 'point_in_triangle, something is wrong'
        call par_ex
        stop
        xe(1:3)=mesh_fesom%coord_nod2D(1, elnodes(1:3))
        ye(1:3)=mesh_fesom%coord_nod2D(2, elnodes(1:3))
     endif

     !//TODO: replace by below??
 !    do q=1,3
 !     xe(q)=geolon2d(mesh_fesom%elem2D_nodes(q, elem))
 !     ye(q)=geolat2d(mesh_fesom%elem2D_nodes(q, elem))
 !    end do

            !=====
            ! Cyclicity:
            ! Remove 2*pi jumps in x-coordinate
            ! of nodes in triangle  
            ! =====
            xe(2)=xe(2)-xe(1)
            xe(3)=xe(3)-xe(1)
            DO q=2,3
             if(xe(q)> pi) xe(q)=xe(q)-2*pi
             if(xe(q)<-pi) xe(q)=xe(q)+2*pi
            END DO
            xe(2)=xe(2)+xe(1)
            xe(3)=xe(3)+xe(1)
            ! =====
            ! Now x(2) and x(3) are measured consistent to x(1) 
            mean_lon=sum(xe(1:3))/3.0_8
            mean_lat=sum(ye(1:3))/3.0_8
            xt1=maxval(abs(xe(1:3)-mean_lon))
            yt1=maxval(abs(ye(1:3)-mean_lat))
            ! =====
            ! Cyclicity
            ! =====
            if(xt1> pi) xt1=xt1-2*pi
            if(xt1<-pi) xt1=xt1+2*pi

            x1=pt(1)-mean_lon
            if(x1>pi) x1=x1-2*pi
            if(x1<-pi) x1=x1+2*pi

            if((abs(x1)>xt1).or.(abs(pt(2)-mean_lat)>yt1)) cycle

            ! ========
            ! Find if the regular mesh node is within elem
            ! ========
            xe(4)=xe(1)
            ye(4)=ye(1)
            angle=0.
            DO q=1,3
               xt1=xe(q)
               xt2=xe(q+1)
               yt1=ye(q)
               yt2=ye(q+1)
               x1=xt1-pt(1)
               x2=xt2-pt(1)
               y1=yt1-pt(2)
               y2=yt2-pt(2)
               ! =====
               ! Cyclicity check
               ! =====
               if(x1>pi) x1=x1-2*pi
               if(x1<-pi) x1=x1+2*pi
               if(x2>pi) x2=x2-2*pi
               if(x2<-pi) x2=x2+2*pi

               alph=x1*y2-x2*y1
               s1=x1*x1+y1*y1
               s2=x2*x2+y2*y2
               alph=alph/sqrt(s1*s2)
               alph=asin(alph)
               IF ((xt1-xt2)**2+(yt1-yt2)**2 > max(s1,s2)) THEN
               if (alph>0) alph=pi-alph
               if (alph<=0) alph=-pi-alph
               END IF
               angle=angle+alph
             END DO
             IF (abs(angle)>pi) THEN
             el2D=elem
             exit
             END IF
         ENDDO
 ! ============
 ! A triangle does not belong to a single PE, thus some i,j will be 
 ! updated on several PEs. It is assumed that elem belongs to a PE 
 ! that own the first node.
 ! ============
  if(el2D>0) then
    if (mesh_fesom%elem2D_nodes(1, el2D) > myDim_nod2D)  then
            el2D=0
    end if
  end if
END SUBROUTINE point_in_triangle_pdaf
