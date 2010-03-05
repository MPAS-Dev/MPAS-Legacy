program scanner

integer iday, ndays
integer idata, ndata, nskip
character*80 :: tmp, rev

real, allocatable, dimension(:,:) :: h

open(12,file='summary',form='formatted',access='append')


!open(13,file='svn.info',form='formatted')
!do i=1,4
!read(13,*) tmp
!enddo
!read(13,*) tmp, idata
!write(8,*) idata
!close(8)
!read(8,'(a4)') rev
!rev = trim(rev)
!write(9,'(a4)') rev
!stop

read(5,*) ndays
read(5,*) ndata
read(5,*) nskip
allocate(h(ndata,ndays))
h=0

do i=1,nskip
read(10,*) 
enddo

do idays=1,ndays
do idata=1,ndata
read(10,*) h(idata,idays)
enddo
enddo

do idays=1,ndays
write(6,'(2e20.10)') minval(h(:,idays)), maxval(h(:,idays))
if(idays.eq.ndays) write(12,'(i8,2e20.10)') ndata, minval(h(:,idays)), maxval(h(:,idays))
enddo

end program scanner
