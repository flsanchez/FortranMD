module tablas

  implicit none
  private

    public :: leer_tablas

  contains

    ! Esta subrutina lee el archivo en la ubicacion pathToTable (suponiendo
    ! mismo directorio que el ejecutable) y aloca en la memoria la variable
    ! tabla, con el tamaño del array dado por la primer linea de la tabla

    subroutine leer_tablas(tabla,pathToTable)

      integer :: j
      integer :: N
      real(16), dimension(:), allocatable, intent(inout) :: tabla
      character(len=*), intent(in) :: pathToTable

      open(unit = 100, file = pathToTable, status = 'old', action = 'read')
      read(100,*), N
      allocate(tabla(N))
      !do j = 1, N
      read(100,*) (tabla(j),j = 1,N)
      !end do
      close(unit = 100)

    end subroutine

end module

! program main
!
!   use tablas
!   double precision, dimension(:), allocatable :: tablas
!   character(len=10) :: pathToTable
!   integer :: N
!
!   pathToTable = 'tablas.txt'
!   call leer_tablas(tablas,pathToTable)
!   N = size(tablas)
!   do j = 1, N
!     write(*,*), tablas(j)
!   end do
!   write(*,*) 'N =',size(tablas)
!   deallocate(tablas)
!
! end program
