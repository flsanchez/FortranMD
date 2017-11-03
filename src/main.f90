program hola
  use observables
  use funciones
  use integrador
  use tablas

  real(16), parameter :: m_pi = 3.14159265359
  real(16), parameter :: h_barra = 25.0/47
  integer, parameter :: N = 100000    ! Cantidad de pasos
  ! Choque frontal (unidimensional) de 2 particulas
  real(16) :: L = 1000    ! Caja muy grande, donde las imagenes no interactuan
  real(16), dimension(2,12) :: vector   ! q,y,x,p
  real(16), dimension(:), allocatable :: LUT
  real(16) :: delta= 1E-4
  real(16) :: w = 1E2 ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V
  real(16), dimension(N) :: dq
  real(16), dimension(N) :: dp
  real(16), dimension(N) :: dx
  real(16), dimension(N) :: dy
  real(16), dimension(N) :: energia
  real(16), dimension(N) :: energiaC
  integer :: i
  !call leer_tablas(LUT,'tabla.txt')

  matriz = matrizhc(2*delta*w)

  V = (200)*h_barra**3 !(3*pi**2)**(2.0/3)/5*(1.4)**(2.5)*h_barra**2*(2/L**3)**(2.0/3)
  dp(1) = -1     ! Diferencia de impulso inicial
  dq(1) = 5.0    ! Las particulas arrancan a distancia 3, donde el potencial se anula
  dx(1) = dq(1)
  dy(1) = dp(1)
  vector(:,:) = 0           ! El choque es en una sola direccion
! Coordenadas q,p
  vector(1,1) = L/2
  vector(2,1) = L/2+dq(1)
  vector(1,10) = -dp(1)     ! Esta es la velocidad de la particula disparada, la otra esta quieta
! Coordenadas x,y
  vector(1,7) = L/2
  vector(2,7) = L/2+dq(1)
  vector(1,4) = -dp(1)

  energia(1) = EnergiaCinetica(vector) + EnergiaPotencial(vector,L,LUT,V)
  energiaC(1) = EnergiaCinetica(vector)
  do i=2,N
    call avanzar(vector,delta,V,L,LUT,matriz)
    dp(i) = vector(2,10)-vector(1,10)
    dq(i) = vector(2,1)-vector(1,1)
    dx(i) = vector(2,7)-vector(1,7)
    dy(i) = vector(2,4)-vector(1,4)
    energia(i) = EnergiaCinetica(vector) + EnergiaPotencial(vector,L,LUT,V)
    energiaC(i) = EnergiaCinetica(vector)
  end do

  open(unit = 100, file = "choque.txt")
  do i = 1,N
    write(100,*) dq(i),";",dp(i),";",dx(i),";",dy(i),";",energia(i),";",energiaC(i)
  enddo
  close(unit = 100)

end program
