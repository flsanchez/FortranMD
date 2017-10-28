module integrador

implicit none
private

	public :: TranformacionHa
	public :: TranformacionHb
	public :: TranformacionHc

contains
	! vector(i,j): i es numero de particula (1 <= i <= N), (1 <= j <= 12)
	! es q,y,x,p, dim(vector(i,j)) = Nx12
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
