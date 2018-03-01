      program FindCOM

      integer stream,i,j,Nwater
      parameter (stream=50,Nwater=20)

      double precision point
      double precision com(3)

c     Initialize COM to 0
      do j=1,3
         com(j)=0.0D0
      enddo

c     Open Stream
      open(stream,FILE='points.20',STATUS='OLD')


c     Calculate COM on the fly (xyz Oxygen, then xyz H1, xyz H2)

      do i=1,Nwater
         do j=1,3
            read(stream,*)point
            com(j)=com(j)+dble(point)*15.9994D0
         enddo

         do j=1,3
            read(stream,*)point
            com(j)=com(j)+dble(point)*1.00794D0
         enddo

         do j=1,3
            read(stream,*)point
            com(j)=com(j)+dble(point)*1.00794D0
         enddo
      enddo
      
      do j=1,3
         com(j)=com(j)/(15.9994D0 + 2.0D0*1.00794D0)
      enddo

c     Output COM coordinates
      print*,"Water Cage Center is : "
      do j=1,3
         print*, com(j)
      enddo

      stop
      end
