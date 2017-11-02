program main

  use tablas
  real(16), dimension(:), allocatable :: tablas
  character(len=10) :: pathToTable
  integer :: N

  pathToTable = 'tablas.txt'
  call leer_tablas(tablas,pathToTable)
  N = size(tablas)
  do j = 1, N
    write(*,*), tablas(j)
  end do
  write(*,*) 'N =',size(tablas)
  deallocate(tablas)

end program
