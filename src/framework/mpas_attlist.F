module mpas_attlist

   use mpas_kind_types

   ! Derived type for holding field attributes
   type att_list_type
      character (len=64) :: attName
      integer :: attType
      integer :: attValueInt
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=1024) :: attValueText
      type (att_list_type), pointer :: next => null()
   end type att_list_type

   interface mpas_add_att
      module procedure mpas_add_att_int0d
      module procedure mpas_add_att_int1d
      module procedure mpas_add_att_real0d
      module procedure mpas_add_att_real1d
      module procedure mpas_add_att_text
   end interface mpas_add_att

   interface mpas_get_att
      module procedure mpas_get_att_int0d
      module procedure mpas_get_att_int1d
      module procedure mpas_get_att_real0d
      module procedure mpas_get_att_real1d
      module procedure mpas_get_att_text
   end interface mpas_get_att


   !!!!! PRIVATE? !!!!!

   integer, parameter :: ATT_INT   = 1
   integer, parameter :: ATT_INTA  = 2
   integer, parameter :: ATT_REAL  = 3
   integer, parameter :: ATT_REALA = 4
   integer, parameter :: ATT_TEXT  = 5


contains


   subroutine mpas_add_att_int0d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      integer, intent(in) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_INT
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueInt = attValue

   end subroutine mpas_add_att_int0d


   subroutine mpas_add_att_int1d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      integer, dimension(:), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_INTA
      allocate(cursor % attValueIntA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueIntA(:) = attValue(:)

   end subroutine mpas_add_att_int1d


   subroutine mpas_add_att_real0d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_REAL
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueReal = attValue

   end subroutine mpas_add_att_real0d


   subroutine mpas_add_att_real1d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_REALA
      allocate(cursor % attValueRealA(size(attValue)))
      write(cursor % attName,'(a)') trim(attName)
      cursor % attValueRealA(:) = attValue(:)

   end subroutine mpas_add_att_real1d


   subroutine mpas_add_att_text(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      character (len=*), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      if (.not. associated(attList)) then
         allocate(attList)
         cursor => attList
      else
         cursor => attList
         do while (associated(cursor % next))
            cursor => cursor % next
         end do
         allocate(cursor % next)
         cursor => cursor % next
      end if
     
      cursor % attType = ATT_TEXT
      write(cursor % attName,'(a)') trim(attName)
      write(cursor % attValueText,'(a)') trim(attValue)

   end subroutine mpas_add_att_text


   subroutine mpas_get_att_int0d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      integer, intent(out) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_INT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueInt
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_int0d


   subroutine mpas_get_att_int1d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      integer, dimension(:), pointer :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_INTA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueIntA)))
               attValue(:) = cursor % attValueIntA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_int1d


   subroutine mpas_get_att_real0d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(out) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_REAL) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               attValue = cursor % attValueReal
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_real0d


   subroutine mpas_get_att_real1d(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), pointer :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_REALA) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               allocate(attValue(size(cursor % attValueRealA)))
               attValue(:) = cursor % attValueRealA(:)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_real1d


   subroutine mpas_get_att_text(attList, attName, attValue, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      character (len=*), intent(out) :: attValue
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            if (cursor % attType /= ATT_TEXT) then
               if (present(ierr)) ierr = 1        ! Wrong type
            else
               write(attValue,'(a)') trim(cursor % attValueText)
            end if
            return
         end if 
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_get_att_text


   subroutine mpas_remove_att(attList, attName, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      character (len=*), intent(in) :: attName
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor, cursor_prev

      if (present(ierr)) ierr = 0

      cursor => attList

      ! Item is at the head of the list
      if (trim(attName) == trim(cursor % attName)) then
         attList => cursor % next
         if (cursor % attType == ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         return
      end if

      cursor_prev => cursor
      cursor => cursor % next
      do while (associated(cursor))
         if (trim(attName) == trim(cursor % attName)) then
            cursor_prev % next => cursor % next

            if (cursor % attType == ATT_REALA) then
               deallocate(cursor % attValueRealA)
            else if (cursor % attType == ATT_INTA) then
               deallocate(cursor % attValueIntA)
            end if
            deallocate(cursor)
            
            return
         end if 

         cursor_prev => cursor
         cursor => cursor % next
      end do

      if (present(ierr)) ierr = 1    ! Not found

   end subroutine mpas_remove_att


   subroutine mpas_deallocate_attlist(attList, ierr)

      implicit none

      type (att_list_type), pointer :: attList
      integer, intent(out), optional :: ierr

      type (att_list_type), pointer :: cursor

      if (present(ierr)) ierr = 0

      cursor => attList
      do while (associated(cursor))
         attList => attList % next
         if (cursor % attType == ATT_REALA) then
            deallocate(cursor % attValueRealA)
         else if (cursor % attType == ATT_INTA) then
            deallocate(cursor % attValueIntA)
         end if
         deallocate(cursor)
         cursor => attList
      end do

   end subroutine mpas_deallocate_attlist

end module mpas_attlist
