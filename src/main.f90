program menu
  integer(4) :: N_args
  character(len=10), dimension(:), allocatable :: args
  integer(4) :: i
  N_args = IARGC()
  allocate(args(N_args))
  ! Argumentos pasados por consola
  do i=1,N_args
    CALL getarg(i,args(i))
  end do
  if (N_args==0) then
    call Opcion1
  end if
  if (N_args==98347924) then
    call EspacioFases(args)
  end if
  if(args(1)=="i") then
    call main_inicializar(args)
  end if
  deallocate(args)
end program


!!! --------------------------------------------------- !!!
!!! ----------------> SUBRUTINAS <--------------------- !!!
!!! --------------------------------------------------- !!!

subroutine Opcion1()
  use observables
  use funciones
  use integrador
  use tablas
  real(16), parameter :: m_pi = 3.14159265359
  real(16), parameter :: h_barra = 25.0/47
  integer, parameter :: N = 50000    ! Cantidad de pasos
  ! Choque frontal (unidimensional) de 2 particulas
  real(16) :: L = 1000    ! Caja muy grande, donde las imagenes no interactuan
  real(16), dimension(2,12) :: vector   ! q,y,x,p
  real(16), dimension(:), allocatable :: LUT
  real(16) :: delta= 1E-4
  real(16) :: w = 1E3 ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V
  real(16), dimension(N) :: dq
  real(16), dimension(N) :: dp
  real(16), dimension(N) :: dx
  real(16), dimension(N) :: dy
  real(16), dimension(N) :: energia
  real(16), dimension(N) :: energiaC
  real(16), dimension(N) :: Hbolt
  integer :: i
  call leer_tablas(LUT,'tabla.txt')

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
    hbolt(i) = HBoltzmannj(vector,2,10)
  end do

  open(unit = 100, file = "choque.txt")
  do i = 1,N
    write(100,*) dq(i),";",dp(i),";",dx(i),";",dy(i),";",energia(i),";",energiaC(i),";",hbolt(i)
  enddo
  close(unit = 100)
end subroutine

subroutine EspacioFases(args)  ! Aun no esta terminada (de hecho, no tiene nada)
  use observables
  use funciones
  use integrador
  use tablas

  real(16), parameter :: m_pi = 3.14159265359
  real(16), parameter :: h_barra = 25.0/47
  character(len=10), dimension(3), intent(in) :: args
  character(len=13) :: archivo
  real(8) :: p_min
  real(8) :: p_max
  real(8) :: p_step
  integer(4) :: Cant_ps
  ! Variables del problema
  real(16), dimension(2,12) :: vector = 0   ! q,y,x,p
  real(16), dimension(:), allocatable :: LUT
  real(16), parameter :: L = 100
  real(16) :: delta= 1E-4
  real(16) :: w = 1E2 ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V = 10
  real(16), dimension(:), allocatable :: dq
  real(16), dimension(:), allocatable :: dp
  integer(4) :: i
  integer(8) :: j
  integer(8) :: k
  integer(8) :: N = 100000 ! Estimacion de la cantidad de pasos (de tanto correrlo)
  ! Asignacion argumentos
  read(args(1), *) p_min
  read(args(2), *) p_max
  read(args(3), *) Cant_ps
  call leer_tablas(LUT,'tabla.txt')
  p_step = (p_max-p_min)/(Cant_ps-1)
  write(*,*) p_min,p_max, Cant_ps, p_step

  matriz = matrizhc(2*delta*w)

  allocate(dq(N),dp(N))

  do i=0,Cant_ps-1
      dp(1) = -(p_min+p_step*i)
      write(archivo,"(A7,I2,A4)") "choque_", i,".txt"
      write(*,*) -dp(1)
      dq(1) = 3.0
      if(dp(1)>0) then  ! Si el impulso va para el otro lado,
        dq(1) = -dq(1)  ! pongo la particula para el otro lado
      end if
    ! Reseteo el vector
      vector(:,:) = 0
    ! Coordenadas q,p
      vector(1,1) = L/2
      vector(2,1) = L/2+dq(1)
      vector(1,10) = -dp(1)     ! Esta es la velocidad de la particula disparada, la otra esta quieta
    ! Coordenadas x,y
      vector(1,7) = L/2
      vector(2,7) = L/2+dq(1)
      vector(1,4) = -dp(1)

      j = 1
      do while((dq(j) .le. 3.0) .and. (dq(j) .ge. -3.0))
        j = j+1
        call avanzar(vector,delta,V,L,LUT,matriz)
        dq(j) = vector(2,7)-vector(1,7)
        dp(j) = vector(2,4)-vector(1,4)
      end do
    open(unit = 100, file = archivo)
    do k = 1,j
      write(100,*) dq(k),";",dp(k)
    enddo
    close(unit = 100)
  end do
end subroutine

subroutine main_inicializar(args)
  use observables
  use funciones
  use integrador
  use tablas
  use init

  character(len=10), dimension(4) :: args
  integer(4) :: N
  real(16) :: T
  real(16) :: rho
  real(16) :: L
  real(16), dimension(:,:), allocatable :: vector
  integer(4) :: i

  read(args(2), *) N
  read(args(3), *) T
  read(args(4), *) rho
  write(*,*) N,T,rho
  L = (N/rho)**(1.0/3)
  allocate(vector(N,12))

  call inicializar(vector,T,L)

  open(unit=100, file ="init.txt")
  do i=1,N
    !write(100,*) vector(i,1),vector(i,2),vector(i,3),vector(i,4),vector(i,5),vector(i,6)
    write(100,*) vector(i,1:6)
  end do
  close(100)
  open(unit=100, file ="gauss.txt")
  do i=1,100000
    !write(100,*) vector(i,1),vector(i,2),vector(i,3),vector(i,4),vector(i,5),vector(i,6)
    write(100,*) Rand_Gauss(T)
  end do
  close(100)

end subroutine
