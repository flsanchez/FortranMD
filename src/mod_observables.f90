module observables

	use funciones

	implicit none

	private

		public :: EnergiaCinetica
		public :: EnergiaPotencial

	contains

		! Por default, estoy tomando como coordenadas (q,p) en lugar de (x,y) o alguna combinación de ambas

		real(16) function EnergiaCinetica(vector)
			real(16), dimension(:,:), intent(in) :: vector
			EnergiaCinetica = 0.5*sum(vector(:,10:12)**2)
		end function

		real(16) function EnergiaPotencial(vector,L,LUT,V)
			real(16), dimension(:,:), intent(in) :: vector
			real(16), dimension(:), intent(in) :: LUT
			real(16), intent(in) :: L
			real(16), intent(in) :: V
			integer :: i
			integer :: j
			EnergiaPotencial = 0
			do i=2,ubound(vector,1)
				do j=1,(i-1)
					EnergiaPotencial = EnergiaPotencial + Valor_LUT(LUT,(DistanciaCuad(vector,i,j,L,0)+DistanciaCuad(vector,i,j,L,3))/2)
				end do
			end do
			EnergiaPotencial = V*EnergiaPotencial
		end function

end module
