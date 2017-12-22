module grabar

implicit none
private

	public :: grabarXYZ
	public :: grabarVector

contains

  subroutine grabarXYZ(vector,unit,L)
    real(16), dimension(:,:), intent(in) :: vector
    integer(4), intent(in) :: unit
    real(16), intent(in) :: L
    integer(4) :: i !indice de particula
    integer(4) :: j
    integer(4) :: N
    integer(4) :: s = 9

    N = ubound(vector, dim=1)
    write(unit,"(I5.5)") N
    write(unit,*)
    do i = 1,N
      write(unit,"(A2,3F8.5)") "N ",vector(i,1:3)/L*s
    end do
  end subroutine

	subroutine grabarVector(vector,unit)
		integer :: i
		integer :: j
		integer, intent(in) :: unit
		real(16), dimension(:,:), intent(in) :: vector

		do i = 1,ubound(vector,dim=1)
			do j = 1,ubound(vector,dim=2)
				write(unit, "(f16.8A)", advance='no') vector(i,j),','
			enddo
			write(unit, "(A)", advance='no') ";"
		enddo
		write(unit,*)
	end subroutine

end module
