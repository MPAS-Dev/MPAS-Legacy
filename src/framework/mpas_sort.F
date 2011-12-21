module mpas_sort

   use mpas_kind_types

   interface quicksort
      module procedure mpas_quicksort_int
      module procedure mpas_quicksort_real
   end interface


   contains


   recursive subroutine mpas_mergesort(array, d1, n1, n2)
   
      implicit none
   
      ! Arguments
      integer, intent(in) :: n1, n2, d1
      integer, dimension(1:d1,n1:n2), intent(inout) :: array
   
      ! Local variables
      integer :: i, j, k
      integer :: rtemp
      integer, dimension(1:d1,1:n2-n1+1) :: temp
   
      if (n1 >= n2) return
   
      if (n2 - n1 == 1) then
        if (array(1,n1) > array(1,n2)) then
           do i=1,d1
              rtemp = array(i,n1)
              array(i,n1) = array(i,n2)
              array(i,n2) = rtemp
           end do
        end if
        return
      end if
   
      call mpas_mergesort(array(1:d1,n1:n1+(n2-n1+1)/2), d1, n1, n1+(n2-n1+1)/2)
      call mpas_mergesort(array(1:d1,n1+((n2-n1+1)/2)+1:n2), d1, n1+((n2-n1+1)/2)+1, n2)
   
      i = n1
      j = n1 + ((n2-n1+1)/2) + 1
      k = 1
      do while (i <= n1+(n2-n1+1)/2 .and. j <= n2)
        if (array(1,i) < array(1,j)) then
          temp(1:d1,k) = array(1:d1,i)
          k = k + 1
          i = i + 1
        else
          temp(1:d1,k) = array(1:d1,j)
          k = k + 1
          j = j + 1
        end if
      end do
   
      if (i <= n1+(n2-n1+1)/2) then
        do while (i <= n1+(n2-n1+1)/2)
          temp(1:d1,k) = array(1:d1,i)
          i = i + 1
          k = k + 1
        end do
      else
        do while (j <= n2)
          temp(1:d1,k) = array(1:d1,j)
          j = j + 1
          k = k + 1
        end do
      end if
   
      array(1:d1,n1:n2) = temp(1:d1,1:k-1)
   
   end subroutine mpas_mergesort


   subroutine mpas_quicksort_int(nArray, array)

      implicit none

      integer, intent(in) :: nArray
      integer, dimension(2,nArray), intent(inout) :: array

      integer :: i, j, top, l, r, pivot, s
      integer :: pivot_value
      integer, dimension(2) :: temp
      integer, dimension(1000) :: lstack, rstack

      if (nArray < 1) return

      top = 1
      lstack(top) = 1
      rstack(top) = nArray

      do while (top > 0)

         l = lstack(top)
         r = rstack(top)
         top = top - 1

         pivot = (l+r)/2

         pivot_value = array(1,pivot)
         temp(:) = array(:,pivot)
         array(:,pivot) = array(:,r)
         array(:,r) = temp(:)

         s = l
         do i=l,r-1
            if (array(1,i) <= pivot_value) then
               temp(:) = array(:,s)
               array(:,s) = array(:,i)
               array(:,i) = temp(:)
               s = s + 1
            end if
         end do

         temp(:) = array(:,s)
         array(:,s) = array(:,r)
         array(:,r) = temp(:)

         if (s-1 > l) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = l
            rstack(top) = s-1
         end if

         if (r > s+1) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = s+1
            rstack(top) = r
         end if
      end do

   end subroutine mpas_quicksort_int


   subroutine mpas_quicksort_real(nArray, array)

      implicit none

      integer, intent(in) :: nArray
      real (kind=RKIND), dimension(2,nArray), intent(inout) :: array

      integer :: i, j, top, l, r, pivot, s
      real (kind=RKIND) :: pivot_value
      real (kind=RKIND), dimension(2) :: temp
      integer, dimension(1000) :: lstack, rstack

      if (nArray < 1) return

      top = 1
      lstack(top) = 1
      rstack(top) = nArray

      do while (top > 0)

         l = lstack(top)
         r = rstack(top)
         top = top - 1

         pivot = (l+r)/2

         pivot_value = array(1,pivot)
         temp(:) = array(:,pivot)
         array(:,pivot) = array(:,r)
         array(:,r) = temp(:)

         s = l
         do i=l,r-1
            if (array(1,i) <= pivot_value) then
               temp(:) = array(:,s)
               array(:,s) = array(:,i)
               array(:,i) = temp(:)
               s = s + 1
            end if
         end do

         temp(:) = array(:,s)
         array(:,s) = array(:,r)
         array(:,r) = temp(:)

         if (s-1 > l) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = l
            rstack(top) = s-1
         end if

         if (r > s+1) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = s+1
            rstack(top) = r
         end if
      end do

   end subroutine mpas_quicksort_real


   integer function mpas_binary_search(array, d1, n1, n2, key)

      implicit none

      integer, intent(in) :: d1, n1, n2, key
      integer, dimension(d1,n1:n2), intent(in) :: array

      integer :: l, u, k

      mpas_binary_search = n2+1

      l = n1
      u = n2
      k = (l+u)/2
      do while (u >= l)
         if (array(1,k) == key) then
            mpas_binary_search = k
            exit   
         else if (array(1,k) < key) then
            l = k + 1
            k = (l+u)/2
         else   
            u = k - 1
            k = (l+u)/2
         end if 
      end do 

   end function mpas_binary_search

end module mpas_sort
