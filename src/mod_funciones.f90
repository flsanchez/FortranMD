module funciones

implicit none
private

	public :: CdCP
	public :: CCperiod
	public :: CCperiod2
	public :: DistanciaCuad
	public :: RestaFases
	public :: NormaCuadrado
	public :: Valor_LUT
	public :: Valor_LUT_der
	public :: Esperanza
	public :: Varianza
	public :: init_random_seed

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
			if (d<-0.5*L) then
				d = d+L
			end if
		end if
	end subroutine

	real(16) function CCperiod(d,L)
		real(16), intent(in) :: L
		real(16), intent(in) :: d
		CCperiod = d - L*floor(d/L)
	end function

	subroutine CCperiod2(vector,L)
		real(16), dimension(:,:), intent(inout) :: vector
		real(16), intent(in) :: L
		integer(4) :: i
		integer(4) :: j
		do i=1,ubound(vector,1)
			do j=1,3
				if(vector(i,j)<0 .and. vector(i,j+6)<0) then
					vector(i,j) = vector(i,j)+L
					vector(i,j+6) = vector(i,j+6)+L
				else
					if (vector(i,j)>L .and. vector(i,j+6)>L) then
						vector(i,j) = vector(i,j)-L
						vector(i,j+6) = vector(i,j+6)-L
					end if
				end if
			end do
		end do
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
			DistanciaCuad = sum((vector(i,m+1:m+3)-vector(j,m+1:m+3))**2)
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
		real(16) :: r_max = 8
		if (s>=r_max) then
			Valor_LUT=0
		else
			r_step = r_max/(ubound(LUT,1)-1)
			i = floor(s/r_step)
			Valor_LUT = (LUT(i+2)-LUT(i+1))*(s/r_step-i)+LUT(i+1)
		end if
	end function


	! real(16) function Valor_LUT_der(LUT,s)
	! 	real(16), intent(in) :: s
	! 	real(16), dimension(:), intent(in) :: LUT
	! 	integer :: i
	! 	real(16) :: r_step
	! 	real(16) :: r_max = 8
	! 	if (s>=r_max) then
	! 		Valor_LUT_der=0
	! 	else
	! 		r_step = r_max/(ubound(LUT,1)-1)
	! 		i = floor(s/r_step)
	! 		Valor_LUT_der = (LUT(i+2)-LUT(i+1))*(s/r_step-i)+LUT(i+1)+ 3.35575200841245E-04
	! 	end if
	! end function

	real(16) function Valor_LUT_der(LUT,s)
		real(16), intent(in) :: s
		real(16), dimension(:), intent(in) :: LUT
		integer :: i
		real(16) :: r_step
		real(16) :: r_max = 8
		real(16) :: r_spline = 6
		if (s>=r_max) then
			Valor_LUT_der=0
		else
			if(s>=r_spline) then
				Valor_LUT_der = 0.5*Valor_LUT(LUT,r_spline)*(r_max-s)
			else
				r_step = r_max/(ubound(LUT,1)-1)
				i = floor(s/r_step)
				Valor_LUT_der = (LUT(i+2)-LUT(i+1))*(s/r_step-i)+LUT(i+1)
			end if
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

	! 0 es leer la semilla, 1 es hacer una nueva random

	subroutine init_random_seed(rw)

				integer, intent(in) :: rw
	      INTEGER :: i, n, clock
	      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

	      CALL RANDOM_SEED(size = n)
	      ALLOCATE(seed(n))

				CALL SYSTEM_CLOCK(COUNT=clock)

				!write(*,*) "seed",clock

	      seed = clock + 37 * (/ (i - 1, i = 1, n) /)

				if(rw .eq. 1) then
					write(*,*) "GRABANDO SEED..."
					open(unit = 100, file="seed.txt")
					write (100, *) seed
					close(100)
				else
					write(*,*) "LEYENDO SEED..."
					open(unit = 100, file="seed.txt")
					read (100, *), seed
					close(100)
				endif

					CALL RANDOM_SEED(PUT = seed)

	      DEALLOCATE(seed)

	end

end module
