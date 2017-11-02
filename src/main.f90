program hola
  use observables
  use funciones
  use integrador
  use tablas

  real(16), parameter :: pi = 3.14159265359
  real(16), parameter :: h_barra = 25/47
  integer, parameter :: N = 10000    ! Cantidad de pasos
  ! Choque frontal (unidimensional) de 2 particulas
  real(16) :: L = 10    ! Caja muy grande, donde las imagenes no interactuan
  real(16), dimension(2,24) :: vector
  real(16), dimension(:), allocatable :: LUT
  real(16) :: delta= 1E-3
  real(16) :: w =1E3  ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V
  real(16), dimension(N) :: dq
  real(16), dimension(N) :: dp
  integer :: i
  matriz = matrizhc(2*delta*w)
  V = 2E-2 !(3*pi**2)**(2.0/3)/5*(1.4)**(2.5)*h_barra**2*(2/L**3)**(2.0/3)
  dp(1) = 1   ! Diferencia de impulso inicial
  dq(1) = 3.0   ! Las particulas arrancan a distancia 3, donde el potencial se anula
  vector(:,:) = 0           ! El choque es en una sola direccion
  vector(1,1) = L/2
  vector(2,1) = L/2+dq(1)
  vector(1,10) =  dp(1)     ! Esta es la velocidad de la particula disparada, la otra esta quieta

  call leer_tablas(LUT,'tabla.txt')

  do i=2,N
    call avanzar(vector,delta,V,L,LUT,matriz)
    dp(i) = vector(1,10)-vector(2,10)
    dq(i) = vector(2,1)-vector(1,1)
  end do

  open(unit = 100, file = "choque.txt")
  write(100,*) (dq(i),i=1,N)
  write(100,*) (dp(i),i=1,N)
end program
