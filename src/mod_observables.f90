module observables

	use funciones

	implicit none

	private

		public :: EnergiaCinetica
		public :: EnergiaPotencial
		public :: HBoltzmannj
		public :: HBoltzmann

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
					EnergiaPotencial = EnergiaPotencial + Valor_LUT(LUT,(DistanciaCuad(vector,i,j,L,0)+DistanciaCuad(vector,i,j,L,3))*0.5)
					!EnergiaPotencial = EnergiaPotencial + exp(-(DistanciaCuad(vector,i,j,L,0)+DistanciaCuad(vector,i,j,L,3))/2)
				end do
			end do
			EnergiaPotencial = V*EnergiaPotencial
		end function

		real(16) function HBoltzmannj(vector,Nbins,j)
			! vector(i,j): i es numero de particula (1 <= i <= N),
			! j es q,y,x,p (1 <= j <= 12),
			! dim(vector(i,j)) = Nx12
			real(16), dimension(:,:), intent(in) :: vector
			integer, intent(in) :: j ! indice espacial 10,11,12 = x,y,z
			integer, intent(in) :: Nbins
			real(16), dimension(:), allocatable :: histo
			real(16) :: deltav
			integer :: N
			integer :: i
			integer :: idxHisto
			real(16) :: minP

			N = ubound(vector, dim=1)
			minP = minVal(vector(1:N,j))
			! write(*,*) minP, maxVal(vector(1:N,j)),(maxVal(vector(1:N,j))-minP)/(Nbins-1)
			deltav = (maxVal(vector(1:N,j))-minP)/(Nbins-1)
			allocate(histo(Nbins))
			histo = 0
			do i = 1,N
				idxHisto = floor((vector(i,j)-minP)/deltav)+1
				histo(idxHisto) = histo(idxHisto)+1
			end do
			histo = histo/N
			HBoltzmannj = sum(histo*log(histo)*deltav,mask=histo.gt.0)
			deallocate(histo)
		end function

		real(16) function HBoltzmann(vector,Nbins)
			real(16), dimension(:,:), intent(in) :: vector
			integer, intent(in) :: Nbins
			integer :: j

			HBoltzmann = 0
			do j = 10,12
				HBoltzmann = HBoltzmann + (1.0/3.0)*HBoltzmannj(vector,Nbins,j)
			end do
		end function

end module
