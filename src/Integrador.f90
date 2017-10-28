module integrador

implicit none
private

	public :: TranformacionHa
  !public :: TranformacionHb
	public :: TranformacionHc
	public :: matrizhc
	!vector  (q,y,x,p)  [12,N]

contains
	! vector(i,j): i es numero de particula (1 <= i <= N),
	! j es q,y,x,p (1 <= j <= 12),
	! dim(vector(i,j)) = Nx12
	! el V0 = (3pi^2)^(2/3)/5*(7/5)^(5/2)*h_barra^2*rho^(2/3)

	subroutine TransformacionHa(vector,delta,V,LUT)
		real(16), dimension(:,:), intent(inout) :: vector
		real(16), dimension(:), intent(in) :: LUT
		real(16), intent(in) :: delta
		! real(16), intent(in) :: omega
		real(16), intent(in) :: V
		real(16), intent(out) :: qij2
		real(16), intent(out) :: yij2
		real(16), intent(out) :: qi
		real(16), intent(out) :: qj
		real(16), intent(out) :: yi
		real(16), intent(out) :: yj
		integer, intent(out) :: i ! indice de particula
		integer, intent(out) :: j ! indice de sumatoria
		integer, intent(out) :: k ! indice de componente vectorial
		integer, intent(out) :: N

		N = ubound(vector, dim=1)
		do i = 1, N
			! q e y quedan igual
			do k = 1, 6
				vector(i,k) = vector(i,k)
			enddo
			do j = 1, N
				qij2 = DistanciaCuad(vector,i,j,L,0)
				yij2 = DistanciaCuad(vector,i,j,L,1)
				do k = 7, 9
					yi = vector(i,k-3)
					yj = vector(j,k-3)
					vector(i,k) = vector(i,k) - delta/2*(yi-2*V*(yi-yj)*Valor_LUT(-qij2-yij2))
				enddo
				do k = 10, 12
					qi = vector(i,k-9)
					qj = vector(j,k-9)
					vector(i,k) = vector(i,k) + delta*V*(qi-qj)*Valor_LUT(-qij2-yij2))
				enddo
		enddo

	endsubroutine

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
