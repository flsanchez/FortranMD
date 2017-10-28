module funciones

implicit none
private

	public :: CdCP
	public :: DistanciaCuad
	public :: RestaFases
	public :: NormaCuadrado
	public :: Valor_LUT
	public :: Esperanza
	public :: Varianza

contains

	! Esta subrutina transforma d según las CdC períodicas
	! Antes era la funcion Delta, pero le paso directo d para hacerla más directa
	! Aclaración: d es la diferencia entre 2 coordenadas espaciales (q o x) de 2 particulas
	subroutine CdCP(d,L)
		real(16), intent(in) :: L
		real(16), intent(inout) :: d
		if (d>0.5*L) then
			d = d-L
		else
			if (d<-0.5*L) then
				d = L+d
			end if
		end if
	end subroutine

	! El vector es Coordenada x Particula, vector(k,i) es la coordenada k de la i-esima particula


	! Distancia en el espacio de posiciones tanto para q como para x
	! Creo que esto no sirve para nada, de última la borramos
	real(16) function DistanciaCuad(vector,i,j,L,k)
		real(16), dimension(:,:), intent(in) :: vector
		real(16), intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer(4), intent(in) :: k
		real(16) :: d
		integer(4) :: m
		DistanciaCuad = 0 ! La 1º coordenada es la distancia cuadrado en q y la 2º es la distancia cuadrado en x
		if (k==1 .or. k==3) then   ! Es un impulso, no hay CdCP
			do m=1,3
				d = vector(k*3+m,i)-vector(k*3+m,j)
				DistanciaCuad = DistanciaCuad+d*d
			end do
		else  ! Es una posicion, tengo que aplicar CdCP
			do m=1,3
				d = vector(k*3+m,i)-vector(k*3+m,j)
				call CdCP(d,L)
				DistanciaCuad = DistanciaCuad+d*d
			end do
		end if
	end function


	! Vector diferencia en el espacio de fases (qo=po=1)
	function RestaFases(vector,i,j,L)
		real(16), dimension(:,:), intent(in) :: vector
		real(16), intent(in) :: L
		integer, intent(in) :: i
		integer, intent(in) :: j
		real(16), dimension(12) :: RestaFases
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

! Funcion para obtener el valor de la LUT, aun
	real(16) function Valor_LUT(LUT,s)
		real(16), intent(in) :: s
		real(16), dimension(:), intent(in) :: LUT
		integer :: i
		real(16) :: r_step
		if (s>=3.0) then
			Valor_LUT=0
		else
			r_step = 3.0/(ubound(LUT,1)-1)
			i = floor(s/r_step)+1
			Valor_LUT = (LUT(i+1)-LUT(i))*s/r_step+LUT(i)
		end if
	end function



!!! ------------------------- !!!

	! Calcula la norma cuadrado de un vector cualquiera
	real(16) function NormaCuadrado(V)
		real(16), dimension(:), intent(in) :: V
		integer :: i
		NormaCuadrado = 0
		do i= lbound(V,1), ubound(V,1)
			NormaCuadrado = NormaCuadrado + V(i)*V(i)
		end do
	end function

	! Calcula el promedio de un vector cualquiera
	real(16) function Esperanza(V)
		real(16), dimension(:), intent(in) :: V
		integer :: i
		Esperanza = 0
		do i= lbound(V,1), ubound(V,1)
			Esperanza = Esperanza + V(i)
		end do
		Esperanza = Esperanza/(ubound(V,1))
	end function

	real(16) function Varianza(V)
		real(16), dimension(:), intent(in) :: V
		real(16) :: mu
		mu = Esperanza(V)
		Varianza = NormaCuadrado(V)/ubound(V,1)-mu*mu
	end function


end module
