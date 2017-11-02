module observables

	use funciones

	implicit none

	private

		public :: EnergiaCinetica

	contains

		! Por default, estoy tomando como coordenadas (q,p) en lugar de (x,y) o alguna combinación de ambas

		real(16) function EnergiaCinetica(vector)
			real(16), dimension(:,:), intent(in) :: vector
			integer :: i
			EnergiaCinetica = 0
			do i=1,ubound(vector)
				EnergiaCinetica = EnergiaCinetica + vector(i,10)*vector(i,10) + vector(i,11)*vector(i,11) + vector(i,12)*vector(i,12)  ! Medio pajero, pero deberia funcionar
			end do
		end function

		real(16) function EnergiaPotencial(vector,L,LUT,V)
			real(16), dimension(:,:), intent(in) :: vector
			integer :: i
			integer :: j
			EnergiaPotencial = 0
			do i=2,ubound(vector)
				do j=1,(i-1)
					EnergiaPotencial = EnergiaPotencial + Valor_LUT(LUT,DistanciaCuad(vector,i,j,L,1)+DistanciaCuad(vector,i,j,L,4))
				end do
			end do
			EnergiaPotencial = V*EnergiaPotencial
		end function

end module