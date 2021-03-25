      module overlapcodes
        ! Overlap codes  from 
        !Quintela Matos, I. & Escobedo, F. 
        !Congruent phase behavior of a binary compound crystal of colloidal spheres 
        !and dimpled cubes.
        ! The Journal of Chemical Physics 153, 214503 (2020).
        !General flag for overlap
        !0 overlap
        !1 no overlap
        subroutine check(rsq,i,j,itypei,itypej,flag,dscmin)
        include 'parameters.h'
        include 'coords.h'
        integer i,j,flag,p,s,itypei,itypej,flagc,flag1
        double precision sepaxis(3),dscmin,rsq
        double precision vectora(3,6),vectorb(3,6)
    
        !overlap check between cube-sphere and cube-cube   
        dscmin = 0.0
        flag=0
    
        if (itypei.eq.2.and.itypej.eq.2) return ! sphere-sphere
        ! polyhedron-polyhedron
        if (itypei.eq.1.and.itypej.eq.1) then  
    !  These 6 checks are not needed for octahedra
    !check for contact between cubes
        do p =1,3
        sepaxis(:)=axes(:,p,i)
        call test(sepaxis,i,j,flag)
        if (flag.eq.1) goto 999
        sepaxis(:)=axes(:,p,j)
        call test(sepaxis,i,j,flag)
        if (flag.eq.1) goto 999
        enddo
  
        do p=1,3
            do s=1,3
                call crossp(axes(:,p,i),axes(:,s,j),sepaxis)
                if (sepaxis(1).ne.0.0.or.sepaxis(2).ne.0.0.or.
     &   sepaxis(3).ne.0.0) then
                    call test(sepaxis,i,j,flag)
                    if (flag.eq.1) goto 999
                endif              
            enddo
        enddo
        !check if contact is within indentation area
        if (flag.eq.0) then
            call cubeoverlap(i,j,flag) 
            if (flag.eq.1) goto 999
        endif
            ! end polyhedron-polyhedron
            !sphere-polyhedron
        elseif (itypei.eq.1) then
            call sctest(rsq,i,j,flag,dscmin,itypei,itypej)
            if (flag.eq.1) goto 999
        elseif (itypej.eq.1) then
            call sctest(rsq,j,i,flag,dscmin,itypei,itypej)
            if (flag.eq.1) goto 999
        endif
  999     return
        end
        subroutine sctest(rsqsc,i,j,flag,dscmin,itypei,itypej)
        ! Based on the Arvo algorithm for overlap between cubes and
        ! spheres from:
        ! J. Arvo, Graphic Gems (Academic Press Professional, Inc., San Diego, CA,1990).
        !Adapted to acccount for spherical indentation on cube facet
        !
        include 'parameters.h'
        include 'coords.h'
        integer k,flag,flag2,i,j,itypei,itypej
        double precision dsq,dscmin,rsqsc
        double precision ald(3)
        !i cube index
        !itype(i) =1
        !j sphere index
        !itype(j) =2
        flag=0
        dsq =0.0d0
        !dsq=squared distance between cube and sphere
        do k =1,3
           ald(k) = abs(axes(1,k,i)*rij(k)+axes(2,k,i)*rij(k)
     &     +axes(3,k,i)*rij(k))
           if (ald(k).gt.1.0)then
            dsq = dsq + (ald(k)-1.0)*(ald(k)-1.0)
           endif
        enddo
        !projection of rij to the vector perpendicular to the cube facet 
        dscmin = (max(ald(1),ald(2),ald(3)))**2
        !check for contact between cube and sphere
        if (dsq>roversqh22) then
            flag = 1
            return
        !check for minimum distance between indented facet and sphere 
        else if ((dscmin>dminsq).and.(dscmin<dmax2)) then
            ! check for contance outside the indentation area
            !squared lateral distance
            dlat2 = rsqsc - dscmin
            s1 = sqrt(roversqh22 -(dscmin-roverh11)**2.0)
            !maximum squared lateral distance
            dlatm2 = (dh - s1)**2
            if(dlat2.lt.dlatm2) then
                flag = 1
                return
            endif
        endif
        return
        end


        subroutine sotest(rsqsc,ipart,jpart,flag,dscmin)
        include 'parameters.h'
        include 'coords.h'
        !Based on point to triangle distnace algorithm from:
        !P. J. Schneider and D. H. Eberly, Geometric Tools for Computer Graphics 2003
        !Adapted to acccount for spherical indentation on octahedra facet
        real*8 P(3)
        real*8 dist2,dscmin,rsqsc,PP0(3)
        real*8 V(6,3),ri(3),dist(6),rt(3)
        integer nlist(6),flag,ipart,jpart,i,k,ii   
        flag = 0

        ! ipart is always polyhedron
        ! jpart is always sphere
        P(:) = q(:,jpart)
        dx = DMIN3(P(1) - q(1,ipart))
        dy = DMIN3(P(2) - q(2,ipart))
        dz = DMIN3(P(3) - q(3,ipart))
        do k=1,3
         rt(k)=HN(k,1)*dx+HN(k,2)*dy+HN(k,3)*dz
        enddo
        P(:) = rt(:)
        ri(:) = 0.0d0
        call calcvertice(ri,ipart,V)
        

        do i=1,6
         dx = P(1)-V(i,1)
         dy = P(2)-V(i,2)
         dz = P(3)-V(i,3)
         dist(i) = dx**2 + dy**2 + dz**2
         nlist(i)=i
        end do

!    Finding 3 vertices closest to point
!  sorting neighbors
!Outuput dist is sorted
        call sortindex(6,dist,nlist)
        
        ii = nlist(1)

         !check if closest vertice overlap
        if(dist(1).lt.roversqh22) return

        call pointTriangleDistance(V(nlist(1),:),V(nlist(2),:)
     &   ,V(nlist(3),:) ,P,dist2,PP0,dscmin,dist(1))

        if(dist2.gt.roversqh22) then
                flag = 1
                return
        else 
       
                if ((dscmin>dminsq)) then
                dlat2 = rsqsc - dscmin
                s1 = roversqh22 -(sqrt(dscmin)-rindeh)**2.0
                dlatm2 = (dh - sqrt(s1))*(dh - sqrt(s1))
                write(*,*) dlat2,dlatm2,rindeh
                write(*,*) s1,rindeh
                write(*,*) dscmin,dminsq
                pause
                                
                               
                    
                      if(dlat2.lt.dlatm2) then
            
          flag = 1 

                                    return
                            else

                                    return
                            endif
           else
                   return
           ENDIF
        endif
cc          write(6,*) 'dist=',dist2,roversqh22
         
         END
       
		subroutine cubeoverlap(i,j,flag)
		include 'parameters.h'
		include 'coords.h'
		integer flag, i, j,c1,c2,floc,ii,jj,ll,k
		double precision vert2(3),vertji(3),vertij(3)
		double precision dij,dji,distf(6),f1(4,3),qcf(3)
		double precision verts1(8,3),faces1(6,4,3),cfaces1(6,3)
		double precision vertsi(8,3),facesi(6,4,3),cfacesi(6,3)
        double precision vertsj(8,3),facesj(6,4,3),cfacesj(6,3)
		double precision lineoi(3,3),lineoj(3,3),distest(3),divt(3)
		double precision lineorigin(3,3),dins(3),dins2(3)
		double precision midv2,madv2,midv3,madv3,midv1,madv1
		flag = 0
		!distance between vertice of i closest to the center of j and center oj j

        call vertice(j,i,dji,vertji,vertsj,facesj,cfacesj,lineoj)
		call vertice(i,j,dij,vertij,vertsi,facesi,cfacesi,lineoi)

        if (dij<dji) then
            cvmin = dij
			c1 = j
			c2 = i
			verts1 = vertsj
			faces1 = facesj
			cfaces1 = cfacesj
			vert2(:) = vertij
            lineorigin = lineoi
		else
			cvmin = dji
			c1 = i !C1 CENTER
			c2 = j !C2 VETICES
			verts1 = vertsi
			faces1 = facesi
			cfaces1 = cfacesi
			vert2(:) = vertji
			lineorigin = lineoj
		endif

		dins = vert2 - q(:,c1)

		do jj=1,3
			dins(jj) = DMIN3(dins(jj))
		enddo
		do k=1,3
          dins2(k)=HN(k,1)*dins(1)+HN(k,2)*dins(2)+HN(k,3)*dins(3)
        enddo

		d1l=abs(axes(1,1,c1)*dins2(1)+axes(2,1,c1)*dins2(2)
     &   + axes(3,1,c1)*dins2(3))
        d2l=abs(axes(1,2,c1)*dins2(1)+axes(2,2,c1)*dins2(2)
     &   + axes(3,2,c1)*dins2(3))
        d3l=abs(axes(1,3,c1)*dins2(1)+axes(2,3,c1)*dins2(2)
     &   + axes(3,3,c1)*dins2(3))
		ddl = max(d1l,d2l,d3l)

       if (ddl<(a6/2)) then
        do ii =1,6
            qcf = cfaces1(ii,:)- q(:,c2)
            do jj=1,3
                qcf(jj) = DMIN3(qcf(jj))
		 	enddo
		 	distf(ii) = qcf(1)**2+ qcf(2)**2+qcf(3)**2
		enddo
		floc = minloc(distf,integer)
		f1 = faces1(floc,:,:)
		distest = q(:,c1) - vert2(:)
        do jj=1,3
			distest(jj) = DMIN3(distest(jj))
		enddo
		do k=1,3
            divt(k)=HN(k,1)*distest(1)+HN(k,2)*distest(2)
     &        +HN(k,3)*distest(3)
        enddo
		c1l=axes(1,1,c1)*divt(1)+axes(2,1,c1)*divt(2)
     &   + axes(3,1,c1)*divt(3)
		c2l=axes(1,2,c1)*divt(1)+axes(2,2,c1)*divt(2)
     &   + axes(3,2,c1)*divt(3)
		c3l=axes(1,3,c1)*divt(1)+axes(2,3,c1)*divt(2)
     &   + axes(3,3,c1)*divt(3)
		c1l = abs(c1l)
		c2l = abs(c2l)
		c3l = abs(c3l)
		drsqd = max(c1l,c2l,c3l)


		if (drsqd>(a6/2-dp)) then
			call pointrectangle(c1,f1,vert2,dprsq)
            if(dprsq<dpsq) then
				call planeline(c1,f1,cfaces1(floc,:),lineorigin,vert2,flag)
				if (flag.eq.1) return
			endif
		endif
	   endif
	   return
	   end
       subroutine vertice(i,j,dij,minvert,verts,faces,cfaces,lineorigin)
       include 'parameters.h'
       include 'coords.h'
       integer  i, j,ploc,ll,k
       double precision verts(8,3),dist(8),qv(3)
       double precision vt(3),rsqv(8),minvert(3)
       double precision faces(6,4,3),cfaces(6,3)
       double precision uv(8,3),lineorigin(3,3)
        !This subroutine calculate the 8 vertices of a cube and then finds
        !which of these vertices is closer to the center of another cube j
       do k=1,8
        uv(k,1)=daxis(k,1)*axes(1,1,i)+daxis(k,2)*axes(1,2,i)+
     &   daxis(k,3)*axes(1,3,i) + q(1,i)*HN(1,1)
        uv(k,2) = daxis(k,1)*axes(2,1,i)+daxis(k,2)*axes(2,2,i)+
     &   daxis(k,3)*axes(2,3,i) + q(2,i)*HN(1,1)
        uv(k,3) = daxis(k,1)*axes(3,1,i)+daxis(k,2)*axes(3,2,i)+
     &   daxis(k,3)*axes(3,3,i) + q(3,i)*HN(1,1)
      enddo
      
      call issquare(uv,faces,cfaces)
      do ll=1,8
      do k=1,3
        verts(ll,k)=HNI(k,1)*uv(ll,1)+HNI(k,2)*uv(ll,2)
     &  +HNI(k,3)*uv(ll,3)
      enddo
      enddo


      do k =1,6
        do ll =1,3
      cfaces(k,ll)=HNI(ll,1)*cfaces(k,1)+HNI(ll,2)*cfaces(k,2)
     &  + HNI(ll,3)*cfaces(k,3)
        enddo
       do l =1,4
        do ll =1,3
        faces(k,l,ll)=HNI(ll,1)*faces(k,l,1)
     &   + HNI(ll,2)*faces(k,l,2)
     &   + HNI(ll,3)*faces(k,l,3)
        enddo
       enddo
      enddo
      
      do k=1,8
        qv(:) = verts(k,:) - q(:,j)

        do l=1,3
            qv(l) = DMIN3(qv(l))
        enddo
        rsqv(k) = qv(1)**2+qv(2)**2+qv(3)**2
      enddo
      dij = minval(rsqv)
      ploc = minloc(rsqv,integer)
      minvert(:) = verts(ploc,:)
      ll = 1
      do k=1,8
        distmin = distver(uv(ploc,:),uv(k,:))
        if (abs(distmin-4.0)<(1e-4)) then
            lineorigin(ll,:) = verts(k,:)
            ll = ll +1
        endif
      enddo
      return
      end

      subroutine issquare(verts,faces,cfaces)
       include 'parameters.h'
       include 'coords.h'
       integer  i, k,ploc
       double precision verts(8,3),dist(8),distx(3,3)
       double precision faces(6,4,3),cfaces(6,3)
       double precision xx(3,3),qx(3),yy(3,3),diag2
       double precision y12(3),y13(3),y23(3),sid2
       k = 1
       m = 1
       sid2= a6**2
       diag2 = 3*sid2
       do i =2,8
        dist(i) = distver(verts(1,:),verts(i,:))
        if (abs(dist(i)-sid2)<(1e-4)) then
            xx(k,:) = verts(i,:)
            k = k +1
        elseif (abs(dist(i)-diag2)<(1e-4)) then
            qx(:) = verts(i,:)
        elseif (abs(dist(i)-8.0)<(1e-4)) then
            yy(m,:) = verts(i,:)
            m = m +1
        endif
       enddo
       do i = 1,3
           do j =1,3
            distx(j,i) = distver(xx(j,:),yy(i,:))
           enddo
           if (abs(distx(1,i)-distx(2,i))<(1e-4)) then
                y12(:) = yy(i,:)
           endif
           if (abs(distx(1,i)-distx(3,i))<(1e-4)) then
                y13(:) = yy(i,:)
           endif
           if (abs(distx(2,i)-distx(3,i))<(1e-4)) then
                y23(:) = yy(i,:)
           endif
       enddo

       do i =1,3
           faces(i,1,:) = verts(1,:)
           faces(i+3,1,:) = qx(:)
           faces(i,2,:) = xx(i,:)
           faces(i+3,4,:) = xx(i,:)
        enddo
        faces(1,4,:) = xx(3,:)
        faces(1,3,:) = y13(:)
        faces(2,4,:) = xx(1,:)
        faces(2,3,:) = y12(:)
        faces(3,4,:) = xx(2,:)
        faces(3,3,:) = y23(:)
        faces(4,2,:) = y12(:)
        faces(4,3,:) = y13(:)
        faces(5,2,:) = y12(:)
        faces(5,3,:) = y23(:)
        faces(6,2,:) = y13(:)
        faces(6,3,:) = y23(:)
        do i =1,6
          cfaces(i,:) = sum(faces(i,:,:),1)/4
        enddo
        return
       end

      subroutine crossprod(cross,a, b)
            double precision cross(3)
            double precision a(3), b(2)
          
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
      return
      END 

      subroutine pointTriangleDistance(B0,B1,B2,P,sqrdistance
        & ,PP0,dscmin,minv)
          include 'parameters.h'
          include 'coords.h'
          integer iI
          double precision TRI(3,3),P(3),ri(3)
          double precision sqrdistance, PP0(3),dscmin
          double precision B0(3),E0(3),E1(3),D0(3)
          double precision B1(3),B2(3),rt(3)
        
          DO II=1,3
          E0(II) = B1(II) - B0(II)
          E1(II)= B2(II) - B0(II)
          D0(II) = B0(II) - P(II)
          ENDDO
          a = E0(1)*E0(1)+E0(2)*E0(2)+E0(3)*E0(3)        
          b = E0(1)*E1(1)+E0(2)*E1(2)+E0(3)*E1(3)        
          c = E1(1)*E1(1)+E1(2)*E1(2)+E1(3)*E1(3)               
          d = E0(1)*D0(1)+E0(2)*D0(2)+E0(3)*D0(3)
          e = E1(1)*D0(1)+E1(2)*D0(2)+E1(3)*D0(3)
          f = D0(1)*D0(1)+D0(2)*D0(2)+D0(3)*D0(3)
   !        if (f==0) then
   !         WRITE(*,*) DMIN3(B1 - B0)*HN(1,1)
   !         WRITE(*,*) B1 - B0
   !         WRITE(*,*) B0-P
   !         WRITE(*,*) HN(1,1)
   !         PAUSE
   !         dscmin=minv
   !         sqrdistance=minv
   !         PP0 = BO
   !         RETURN
   !        endif
          det = a * c - b * b
          s = b * e - c * d
          t = b * d - a * e
          !write(*,*) a,s,b,s,t,c,e,f,d,det
          if ((s + t) <= det) then
           if (s < 0.0) then
               if (t < 0.0) then
                   ! region4
                   if (d < 0) then
                       t = 0.0
                       if (-d >= a) then
                           s = 1.0
                       else
                           s = -d / a
                       endif
                   else
                       s = 0.0
                       if (e >= 0.0) then
                           t = 0.0
                       else
                           if (-e >= c) then
                               t = 1.0
                           else
                               t = -e / c
                               ! of region 4
                           endif
                        endif
                   endif      
               else
                   ! region 3
                   s = 0
                   if (e >= 0) then
                       t = 0
                   else
                       if (-e >= c) then
                           t = 1
                       else
                           t = -e / c
                           ! of region 3
                       endif
                   endif
                endif 
           else
               if (t < 0) then
                   ! region 5
                   t = 0
                   if (d >= 0) then
                       s = 0
                   else
                       if (-d >= a) then
                           s = 1
                       else
                           s = -d / a
                       endif
                   endif
               else
                   ! region 0
                   !write(*,*) 'region'
                   invDet = 1.0 / det
                   s = s * invDet
                   t = t * invDet
               endif
           endif
         else
           if (s < 0.0) then
               ! region 2
               tmp0 = b + d
               tmp1 = c + e
               if (tmp1 > tmp0) then  !minimum on edge s+t=1
                   numer = tmp1 - tmp0
                   denom = a - 2.0 * b + c
                   if (numer >= denom) then
                       s = 1.0
                       t = 0.0
                   else
                       s = numer / denom
                       t = 1 - s
                   endif
   
               else  ! minimum on edge s=0
                   s = 0.0
                   if (tmp1 <= 0.0) then
                       t = 1
                   else
                       if (e >= 0.0) then
                           t = 0.0
                       else
                           t = -e / c
                           ! of region 2
                       endif
                   endif
                endif
           else
               if (t < 0.0) then
                   ! region6
                   tmp0 = b + e
                   tmp1 = a + d
                   if (tmp1 > tmp0) then
                       numer = tmp1 - tmp0
                       denom = a - 2.0 * b + c
                       if (numer >= denom) then
                           t = 1.0
                           s = 0
                       else
                           t = numer / denom
                           s = 1 - t
                       endif
   
                   else
                       t = 0.0
                       if (tmp1 <= 0.0) then
                           s = 1
                       else
                           if (d >= 0.0) then
                               s = 0.0
                           else
                               s = -d / a
                           endif
                        endif
                   endif
               else
                   ! region 1
                   numer = c + e - b - d
                   if (numer <= 0) then
                       s = 0.0
                       t = 1.0
                   else
                       denom = a - 2.0 * b + c
                       if (numer >= denom) then
                           s = 1.0
                           t = 0.0
                       else
                           s = numer / denom
                           t = 1 - s
                       endif
                   endif
               endif
           endif
         endif
         sqrdistance = a*s**2+2*b*s*t+c*t**2+2.*d*s+2.*e*t+f
         !write(*,*) a,s,b,s,t,c,e,f,d
         !write(*,*) sqrdistance,f,'vert'
         if (sqrdistance <0) sqrdistance = 0
         
         PP0 = B0 + s * E0 + t * E1
         !scallar projection of distance between sphere and octahedron
         !into the normal vector of triangle facet
         u3x = (E1(2)*E0(3) - E1(3)*E0(2))
         u3y = (E1(3)*E0(1) - E1(1)*E0(3))
         u3z = (E1(1)*E0(2) - E1(2)*E0(1))
         !lenght 
         u3r2 = u3x**2 + u3y**2 + u3z**2
         dotur = tt(1)*u3x + tt(2)*u3y + tt(3)*u3z
         dscmin = dotur**2/u3r2
      
         !if (dscmin <0) then
         !  write(*,*) dscmin
         !  write(*,*) sqrdistance
         !endif
         !WRITE(*,*) dotur**2,u3r2,'DLAT'
           return
          END

          subroutine calcvertice(ri,ipart,VERT)
            include 'parameters.h'
            include 'coords.h'
            double precision VERT(6,3)
            double precision ri(3),u(3,3),dc
            integer ipart,i
         
            u(:,:) = axes(:,:,ipart) 
            dc = roverh11
      
           do i=1,3
            VERT(1,i) = ri(i)+dc*u(i,1) 
            VERT(2,i) = ri(i)+dc*u(i,2)
            VERT(3,i) = ri(i)+dc*u(i,3)
            VERT(4,i) = ri(i)-dc*u(i,1)
            VERT(5,i) = ri(i)-dc*u(i,2)
            VERT(6,i) = ri(i)-dc*u(i,3)
            
          end do
        
        
          END
    
          subroutine sortindex(NARRAY,ARRAY,sorted)
            INTEGER sorted(NARRAY),ii,jj
            real*8  ARRAY(NARRAY)
          
            do j=2,NARRAY
              r2 = ARRAY(j)
              nlj = sorted(j) 
              do ii=j-1,1,-1
                if(ARRAY(ii).le.r2) goto 10
                ARRAY(ii+1)=ARRAY(ii)
                sorted(ii+1) = sorted(ii)
              end do
              ii=0
       10        ARRAY(ii+1) = r2
              sorted(ii+1)=nlj
             end do
             return
            end
            FUNCTION TO CALCULATE MINIMUM IMAGE DISTANCES
                REAL*8 FUNCTION DMIN3(XX)
                IMPLICIT REAL*8(A-H,O-Z)
          
                DMIN3 = XX
                if(XX.gt.0.5) then
                    DMIN3= XX -1.0
                elseif(XX.lt.-0.5) then
                    DMIN3 = XX+1.0
                endif
          
                RETURN
                END

  
     end module overlapcodes