module observables

	use funciones

	implicit none

	private

		public :: EnergiaCinetica
		public :: EnergiaPotencial
		public :: HBoltzmannj
		public :: HBoltzmann
		public :: dist_vels
		public :: corr_T_virial
		public :: Reescalar_vel

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


		function dist_vels(vector,E,Nbins)
			real(16), dimension(:,:), intent(in) :: vector
			integer, intent(in) :: Nbins
			real(16), intent(in) :: E
			real(16), dimension(Nbins) :: dist_vels
			integer :: i
			integer :: j
			real(16) :: delta_p2
			delta_p2 = E/Nbins
			dist_vels(:) = 0
			do i = 1,ubound(vector,1)
				j = ceiling(sum(vector(i,10:12)**2)/delta_p2)
				dist_vels(j) = dist_vels(j)+1
			end do
		end function

		real(16) function corr_T_virial(vector,L,LUT,V)		! Correccion a la temperatura del gas ideal
			real(16), dimension(:,:), intent(in) :: vector
			real(16), dimension(:), intent(in) :: LUT
			real(16), intent(in) :: L
			real(16), intent(in) :: V
			real(16) :: prod_int
			integer(4) :: N
			integer :: i
			integer :: j
			N = ubound(vector, dim=1)
			corr_T_virial = 0
			do i=2,N
				do j=1,(i-1)
					prod_int = sum(vector(i,10:12)*(vector(i,10:12)-vector(j,10:12)))
					corr_T_virial = corr_T_virial + prod_int*Valor_LUT(LUT,(DistanciaCuad(vector,i,j,L,0)+DistanciaCuad(vector,i,j,L,3))*0.5)
				end do
			end do
			corr_T_virial = -V*corr_T_virial/(1.5*N)
		end function

		subroutine Reescalar_vel(vector,factor)
			real(16), dimension(:,:), intent(inout) :: vector
			real(16), intent(in) :: factor
			vector(:,10:12) = vector(:,10:12)*factor
			vector(:,4:6) = vector(:,4:6)*factor
		end subroutine


end module
