module integrador

implicit none
private

	!public :: TranformacionHa
!	public :: TranformacionHb
	public :: TranformacionHc
!vector  (q,y,x,p)  [12,N]
!definir y hacer la matrizhc
contains
  subroutine TranformacionHc(vector,matrizhc)
	real(16), dimension(:,:), intent(inout) :: vector
  real(16), dimension(12) :: nvector
  real(16), dimension(4,4), intent(in) :: matrizhc
  integer :: i
  integer :: k
  integer :: l
  nvector=0
  do i=1,2
    do l=1,3
      do k=1,4
        nvector(l)=nvector(l)+0.5*matrizhc(k,1)*vector(l+3*(k-1),i)
        nvector(l+3)=nvector(l+3)+0.5*matrizhc(k,2)*vector(l+3*(k-1),i)
        nvector(l+6)=nvector(l+6)+0.5*matrizhc(k,3)*vector(l+3*(k-1),i)
        nvector(l+9)=nvector(l+9)+0.5*matrizhc(k,4)*vector(l+3*(k-1),i)
      end do
    end do
    vector(:,i)=nvector
  end do
	end subroutine

end module
