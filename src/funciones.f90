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
	! Aclaración: d es la diferencia entre 2 coordenadas espaciales de 2 particulas
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
	end sobroutine


	! Distancia en el espacio de posiciones tanto para q como para x
	! Creo que esto no sirve para nada, de última la borramos
	real function DistanciaCuad(vector,i,j,L) result (res)
		real, dimension(:), intent(in) :: vector
		real, intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		real, dimension(2) :: res
		real :: d
		integer :: k
		res(:) = 0 ! La 1º coordenada es la distancia cuadrado en q y la 2º es la distancia cuadrado en x
		do k=1,3
			d = vector(k,i)-vector(k,j)
			call CdCP(d,L)
			res(1) = res(1)+d*d
			d = vector(k+3,i)-vector(k+3,j)
			call CdCP(d,L)
			res(2) = res(2)+d*d
		end do
	end function


	! Vector diferencia en el espacio de fases (qo=po=1)
	real function RestaFases(vector,i,j,L) result (res)
		real, dimension(:), intent(in) :: vector
		real, intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		real, dimension(12) :: res
		real :: d
		integer :: k
		res(:) = 0
		do k=1,3   ! Vector q
			res(k) = vector(k,i)-vector(k,j)
			call CdCP(res(i),L)
		end do
		do k=4,6	! Vector p
			res(k) = vector(k,i)-vector(k,j)
		end do
		do k=7,9	! Vector x
			res(k) = vector(k,i)-vector(k,j)
			call CdCP(res(i),L)
		end do
		do k=10,12	! Vector y
			res(k) = vector(k,i)-vector(k,j)
		end do
	end function

	! Calcula la norma cuadrado de un vector cualquiera
	real function NormaCuadrado(V)
		real, dimension(:), intent(in) :: V
		integer :: i
		NormaCuadrado = 0
		do i= lbound(V,1), ubound(V,1)
			NormaCuadrado = NormaCuadrado + V(i)*V(i)
		end do
	end function


end module