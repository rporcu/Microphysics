!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_functions_des                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE stl_functions_des

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ClosestPtPointTriangle                                  !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine ClosestPtPointTriangle(pointp, points, closest_point)
      USE param1, only: zero, one

      IMPLICIT NONE

! points are the three nodes of the triangle
! point p is the sphere center
      real(c_real), intent(in), dimension(3,3) :: points
      real(c_real), intent(in), dimension(3) :: pointp
      real(c_real), intent(out), dimension(3) ::  closest_point
! Local variables
      real(c_real), dimension(3) :: pointa, pointb, pointc
      real(c_real), dimension(3) :: ab, ac, ap, bp,cp
      real(c_real) :: d1, d2, d3, d4, vc, v, d5, d6, vb, w, va, denom

      pointa = points(1,:)
      pointb = points(2,:)
      pointc = points(3,:)

      ab = pointb - pointa
      ac = pointc - pointa
      ap = pointp - pointa
      d1 = DOT_PRODUCT(ab, ap)
      d2 = DOT_PRODUCT(ac, ap)

      IF(d1 <= Zero .AND. d2 <= zero) then
         closest_point = pointa
         return
      end if

! Check if P in vertex region outside B
      bp = pointp - pointb
      d3 = DOT_PRODUCT(ab, bp);
      d4 = DOT_PRODUCT(ac, bp);
      if (d3 >= zero .and. d4 <= d3) then
         closest_point = pointb
         return
      endif

! Check if P in edge region of AB, if so return projection of P onto AB
      vc = d1*d4 - d3*d2;
      if (vc <= zero .and. d1 >= zero .and. d3 <= zero) then
         v = d1 / (d1 - d3);
         closest_point =  pointa + v * ab;
         return
      end if

! Check if P in vertex region outside C
      cp = pointp - pointc
      d5 = DOT_PRODUCT(ab, cp)
      d6 = DOT_PRODUCT(ac, cp)
      if (d6 >= zero .and. d5 <= d6) then
         closest_point  = pointc
         return
      endif

! Check if P in edge region of AC, if so return projection of P onto AC
      vb = d5*d2 - d1*d6

      if (vb <= zero .and. d2 >= zero .and. d6 <= zero) then
         w = d2 / (d2 - d6)
         closest_point = pointa + w * ac
         return
      end if

! Check if P in edge region of BC, if so return projection of P onto BC
      va = d3*d6 - d5*d4
      if (va <= zero .and.(d4 - d3) >= zero .and. (d5 - d6) >= zero) then
         w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
         closest_point = pointb + w * (pointc - pointb)
         return
      end if


! P inside face region. Compute Q through its barycentric coordinates (u,v,w)
      denom = one / (va + vb + vc)
      v = vb * denom
      w = vc * denom
      closest_point = pointa + ab * v + ac * w;
      return
      end Subroutine ClosestPtPointTriangle



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: intersectLnPlane                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine intersectLnPlane(ref_line, dir_line, ref_plane,       &
         norm_plane, line_param)

      USE param1, only: zero

      IMPLICIT NONE

! Reference point and direction of the line
      real(c_real), INTENT(IN) :: REF_LINE(3),  DIR_LINE(3)
! reference point and normal of the plane
      real(c_real), INTENT(IN) :: REF_PLANE(3), NORM_PLANE(3)

! line is parameterized as p = p_ref + t * dir_line, t is line_param
      real(c_real), intent(out) :: line_param

      !local vars
      real(c_real) :: denom

      denom = DOT_PRODUCT(dir_line, norm_plane)

      if(denom*denom.gt.zero) then
         line_param = DOT_PRODUCT(ref_plane-ref_line, norm_plane)
         line_param = line_param/denom
      endif

      return
      end subroutine intersectLnPlane

!......................................................................!
! Subroutine TRI_BOX_OVERLAP                                           !
! Author: J.Musser                                   Date: 10-22-2015  !
!                                                                      !
! Purpose: Determine if a box (DES grid cell) intersects the triangle  !
!    (SLT). Note that the DES grid size is slightly increased to       !
!    capture STLs near the boarder of the cell. Otherwise, some        !
!    collisions could ve over looked.                                  !
!                                                                      !
! Author: Tomas Akenine-Moller                   Accessed: 10-22-2015  !
! REF: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/    !
!         code/tribox2.txt                                             !
!......................................................................!
      SUBROUTINE TRI_BOX_OVERLAP(pCENTER, pHALFSIZE, pVERTS, pOVERLAP)

      IMPLICIT NONE

      real(c_real), INTENT(IN) :: pCENTER(3), pHALFSIZE(3)
      real(c_real), INTENT(IN) :: pVERTS(3,3)
      LOGICAL, INTENT(OUT) :: pOVERLAP

      real(c_real) :: v0(3), v1(3), v2(3)
      real(c_real) :: fex, fey, fez
      real(c_real) :: normal(3), e0(3), e1(3), e2(3)

      pOVERLAP = .FALSE.

      v0 = pVERTS(1,:) - pCENTER
      v1 = pVERTS(2,:) - pCENTER
      v2 = pVERTS(3,:) - pCENTER

      e0 = v1-v0
      e1 = v2-v1
      e2 = v0-v2

      fex = abs(e0(1))
      fey = abs(e0(2))
      fez = abs(e0(3))

      if(ATEST_X01(e0(3),e0(2),fez,fey)) return
      if(ATEST_Y02(e0(3),e0(1),fez,fex)) return
      if(ATEST_Z12(e0(2),e0(1),fey,fex)) return

      fex = abs(e1(1))
      fey = abs(e1(2))
      fez = abs(e1(3))

      if(ATEST_X01(e1(3),e1(2),fez,fey)) return
      if(ATEST_Y02(e1(3),e1(1),fez,fex)) return
      if(ATEST_Z0 (e1(2),e1(1),fey,fex)) return

      fex = abs(e2(1))
      fey = abs(e2(2))
      fez = abs(e2(3))

      if(ATEST_X2 (e2(3),e2(2),fez,fey)) return
      if(ATEST_Y1 (e2(3),e2(1),fez,fex)) return
      if(ATEST_Z12(e2(2),e2(1),fey,fex)) return

      if(findMin(v0(1),v1(1),v2(1)) > phalfsize(1) .OR. &
         findMax(v0(1),v1(1),v2(1)) <-phalfsize(1)) return

      if(findMin(v0(2),v1(2),v2(2)) > phalfsize(2) .OR. &
         findMax(v0(2),v1(2),v2(2)) <-phalfsize(2)) return

      if(findMin(v0(3),v1(3),v2(3)) > phalfsize(3) .OR. &
         findMax(v0(3),v1(3),v2(3)) <-phalfsize(3)) return


      normal(1) = e0(2)*e1(3)-e0(3)*e1(2)
      normal(2) = e0(3)*e1(1)-e0(1)*e1(3)
      normal(3) = e0(1)*e1(2)-e0(2)*e1(1)

      if(.NOT.planeBoxOverlap(normal,v0,phalfsize)) return

      pOVERLAP = .TRUE.

      RETURN

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Function: planeBoxOverlap                                            !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION planeBoxOverlap(norm, vert, maxbox)

      real(c_real) :: norm(3), vert(3), maxbox(3)

      integer :: lc
      real(c_real) :: vmin(3), vmax(3), v

      do lc=1,3
         v=vert(lc)
         if(norm(lc) > 0.0d0) then
            vmin(lc) = -maxbox(lc) - v
            vmax(lc) =  maxbox(lc) - v
         else
            vmin(lc) = maxbox(lc) - v
            vmax(lc) =-maxbox(lc) - v
         endif
      enddo

      if(dot_product(norm,vmin) > 0.0d0) then
         planeBoxOverlap = .false.
         return
      elseif(dot_product(norm,vmax) >= 0.0d0) then
         planeBoxOverlap = .true.
         return
      endif

      planeBoxOverlap = .false.
      return

      RETURN
      END FUNCTION planeBoxOverlap

!``````````````````````````````````````````````````````````````````````!
! Function: findMin                                                    !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMin(x0,x1,x2)

      real(c_real) :: x0,x1,x2

      findMin = x0

      if(x1<findMin) findMin=x1
      if(x2<findMin) findMin=x2

      RETURN
      END FUNCTION findMin

!``````````````````````````````````````````````````````````````````````!
! Function: findMax                                                    !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMax(x0,x1,x2)

      real(c_real) :: x0,x1,x2

      findMax = x0

      if(x1>findMax) findMax=x1
      if(x2>findMax) findMax=x2

      RETURN
      END FUNCTION findMax

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_X01                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X01(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p0, p2, rad

      p0 = a*v0(2) - b*v0(3)
      p2 = a*v2(2) - b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X01=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X01

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_X2                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X2(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p0, p1, rad

      p0 = a*v0(2) - b*v0(3)
      p1 = a*v1(2) - b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X2=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X2

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Y02                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y02(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p0, p2, rad

      p0 = -a*v0(1) + b*v0(3)
      p2 = -a*v2(1) + b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y02=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y02

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Y1                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y1(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p0, p1, rad

      p0 = -a*v0(1) + b*v0(3)
      p1 = -a*v1(1) + b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y1=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y1

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Z12                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z12(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p1, p2, rad

      p1 = a*v1(1) - b*v1(2)
      p2 = a*v2(1) - b*v2(2)

      if(p2<p1) then; lMIN=p2; lMAX=p1
      else; lMIN=p1; lMAX=p2; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z12=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z12

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Z0                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z0(a,b,fa,fb)

      real(c_real) :: a, b, fa, fb
      real(c_real) :: lMin, lMax, p0, p1, rad

      p0 = a*v0(1) - b*v0(2)
      p1 = a*v1(1) - b*v1(2)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z0=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z0

      END SUBROUTINE TRI_BOX_OVERLAP

      END MODULE STL_FUNCTIONS_DES
