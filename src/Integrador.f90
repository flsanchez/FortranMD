module integrador

implicit none
private

	!public :: TranformacionHa
!	public :: TranformacionHb
	public :: TranformacionHc
	public :: matrizhc
!vector  (q,y,x,p)  [12,N]

contains
  subroutine TranformacionHc(vector,matrizhc)
	real(16), dimension(:,:), intent(inout) :: vector
  real(16), dimension(12) :: nvector
  real(16), dimension(4,4), intent(in) :: matrizhc
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: m
  integer :: n

  nvector=0
  do i=1,2
    do l=1,3
      do k=1,4
        nvector(l+0)=nvector(l+0)+0.5*matrizhc(k,1)*vector(l+3*(k-1),i)
        nvector(l+3)=nvector(l+3)+0.5*matrizhc(k,2)*vector(l+3*(k-1),i)
        nvector(l+6)=nvector(l+6)+0.5*matrizhc(k,3)*vector(l+3*(k-1),i)
        nvector(l+9)=nvector(l+9)+0.5*matrizhc(k,4)*vector(l+3*(k-1),i)
      end do
    end do
!aca tenia que ordenar de nuevo  yo tenia (qpxy)--->(qyxp)
!entonces q-->q//p-->y...etc 
	do j=1,3	
    vector(j,i)=nvector(j)
    vector(j+6,i)=nvector(j+6)
	end do
	do m=1,3	
    vector(m+9,i)=nvector(m+3)
	end do
	do n=1,3	
    vector(m+3,i)=nvector(m+9)
	end do
  end do
  end subroutine

	real(16) function matrizhc(delta,omega)
		real(16), dimension(4,4)  :: matrizhc
		real(16) :: delta
		real(16) :: omega
		real(16) :: arg
		matrizhc=0
		arg=2*delta*omega
		matrizhc(:,1)=(1+cos(arg), -sin(arg)  , 1-cos(arg), sin(arg) )
		matrizhc(:,2)=(-sin(arg) ,  1-cos(arg),   sin(arg),1+cos(arg))
		matrizhc(:,3)=(1-cos(arg),  sin(arg)  , 1+cos(arg),-sin(arg) )
		matrizhc(:,4)=(sin(arg)  , 1+cos(arg) ,  -sin(arg),1-cos(arg))
	end function
  
end module
