module funciones

implicit none
private

	public :: CdCP
	public :: DistanciaCuad
	public :: RestaFases
	public :: NormaCuadrado

contains

	! Esta subrutina transforma d según las CdC períodicas
	! Antes era la funcion Delta, pero le paso directo d para hacerla más directa
	! Aclaración: d es la diferencia entre 2 coordenadas espaciales (q o x) de 2 particulas
	subroutine CdCP(d,L)
		real, intent(in) :: L
		real, intent(inout) :: d
		if (d>0.5*L) then
			d = d-L
		else
			if (d<-0.5*L) then
				d = L+d
			end if
		end if
	end subroutine


	! Distancia en el espacio de posiciones tanto para q como para x
	! Creo que esto no sirve para nada, de última la borramos
	function DistanciaCuad(vector,i,j,L)
		real, dimension(:,:), intent(in) :: vector
		real, intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		real, dimension(2) :: DistanciaCuad
		real :: d
		integer :: k
		DistanciaCuad(:) = 0 ! La 1º coordenada es la distancia cuadrado en q y la 2º es la distancia cuadrado en x
		do k=1,3
			d = vector(k,i)-vector(k,j)
			call CdCP(d,L)
			DistanciaCuad(1) = DistanciaCuad(1)+d*d
			d = vector(k+3,i)-vector(k+3,j)
			call CdCP(d,L)
			DistanciaCuad(2) = DistanciaCuad(2)+d*d
		end do
	end function


	! Vector diferencia en el espacio de fases (qo=po=1)
	function RestaFases(vector,i,j,L)
		real, dimension(:,:), intent(in) :: vector
		real, intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		real, dimension(12) :: RestaFases
		integer :: k
		RestaFases(:) = 0
		do k=1,3   ! Vector q
			RestaFases(k) = vector(k,i)-vector(k,j)
			call CdCP(RestaFases(k),L)
		end do
		do k=4,6	! Vector y
			RestaFases(k) = vector(k,i)-vector(k,j)
		end do
		do k=7,9	! Vector x
			RestaFases(k) = vector(k,i)-vector(k,j)
			call CdCP(RestaFases(k),L)
		end do
		do k=10,12	! Vector p
			RestaFases(k) = vector(k,i)-vector(k,j)
		end do
	end function


!!! ------------------------- !!!

	! Calcula la norma cuadrado de un vector cualquiera
	real function NormaCuadrado(V)
		real, dimension(:), intent(in) :: V
		integer :: i
		NormaCuadrado = 0
		do i= lbound(V,1), ubound(V,1)
			NormaCuadrado = NormaCuadrado + V(i)*V(i)
		end do
	end function

	! Calcula el promedio de un vector cualquiera
	real function Esperanza(V)
		real, dimension(:), intent(in) :: V
		integer :: i
		Esperanza = 0
		do i= lbound(V,1), ubound(V,1)
			Esperanza = Esperanza + V(i)
		end do
		Esperanza = Esperanza/(ubound(V,1))
	end function

	real function Varianza(V)
		real, dimension(:), intent(in) :: V
		real :: mu
		mu = Esperanza(V)
		Varianza = NormaCuadrado(V)/ubound(V,1)-mu*mu
	end function


end module
