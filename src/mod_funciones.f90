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

	! vector(i,j): i es numero de particula (1 <= i <= N),
	! j es q,y,x,p (1 <= j <= 12),
	! dim(vector(i,j)) = Nx12

	! Esta subrutina transforma d según las CdC períodicas
	! Antes era la funcion Delta, pero le paso directo d para hacerla más directa
	! Aclaración: d es la diferencia entre 2 coordenadas espaciales (q o x) de 2 particulas
	subroutine CdCP(d,L)
		real(16), intent(in) :: L
		real(16), intent(inout) :: d
		if (d>0.5*L) then
			d = d-L
		else
			if (d>-0.5*L) then
				d = d
			else
				d = d+L
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
		integer(4) :: m = 0
		if (k==1 .or. k==3) then   ! Es un impulso, no hay CdCP
			m = 3*k
			DistanciaCuad = sum((vector(i,m:m+3)-vector(j,m:m+3))**2)
		else  ! Es una posicion, tengo que aplicar CdCP
			DistanciaCuad = 0
			do m=1,3
				d = vector(i,k*3+m)-vector(j,k*3+m)
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
		RestaFases = vector(i,:)-vector(j,:)
		do k=1,3   ! Vector q
			call CdCP(RestaFases(k),L)
		end do
		do k=7,9	! Vector x
			call CdCP(RestaFases(k),L)
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
		NormaCuadrado = sum(V**2)
	end function

	! Calcula el promedio de un vector cualquiera
	real(16) function Esperanza(V)
		real(16), dimension(:), intent(in) :: V
		Esperanza = sum(V)/ubound(V,1)
	end function

	real(16) function Varianza(V)
		real(16), dimension(:), intent(in) :: V
		real(16) :: mu
		mu = Esperanza(V)
		Varianza = NormaCuadrado(V)/ubound(V,1)-mu*mu
	end function


end module
