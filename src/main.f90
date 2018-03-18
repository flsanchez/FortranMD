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
  if (args(1)=="o1") then
    call Opcion1
  end if
  if (args(1)=="f") then
    call EspacioFases(args)
  end if
  if(args(1)=="d") then
    call main_distribucion(args)
  end if
  if(args(1)=="t") then
    call main_temperatura(args)
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
  real(16) :: L = 10000    ! Caja muy grande, donde las imagenes no interactuan
  real(16), dimension(3,12) :: vector   ! q,y,x,p
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
  dq(1) = 2.0    ! Las particulas arrancan a distancia 4, donde el potencial se anula
  dx(1) = dq(1)
  dy(1) = dp(1)
  vector(:,:) = 0           ! El choque es en una sola direccion
  ! Coordenadas q,p
  vector(1,1) = L/2-dq(1)/2
  vector(2,1) = L/2+dq(1)/2
  vector(3,1) = L/2
  vector(3,2) = 1
  vector(1,10) = -dp(1)/2     ! Esta es la velocidad de la particula disparada, la otra esta quieta
  vector(2,10) = dp(1)/2     ! Esta es la velocidad de la particula disparada, la otra esta quieta
  ! Coordenadas x,y
  vector(1,7) = L/2-dq(1)/2
  vector(2,7) = L/2+dq(1)/2
  vector(3,7) = L/2
  vector(3,8) = 1
  vector(1,4) = -dp(1)/2
  vector(2,4) = dp(1)/2

  energia(1) = EnergiaCinetica(vector) + EnergiaPotencial(vector,L,LUT,V)
  energiaC(1) = EnergiaCinetica(vector)
  do i=2,N
    call avanzar(vector,delta,V,L,LUT,matriz)
    !dp(i) = vector(2,10)-vector(1,10)
    dp(i) = DistanciaCuad(vector,1,3,L,3)
    !dq(i) = vector(2,1)-vector(1,1)
    dq(i) = vector(3,2)
    !dx(i) = vector(2,7)-vector(1,7)
    dx(i) = DistanciaCuad(vector,1,3,L,0)
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


subroutine main_distribucion(args)
  use observables
  use funciones
  use integrador
  use tablas
  use init
  use grabar

  real(16), parameter :: m_pi = 3.14159265359
  real(16), parameter :: h_barra = 25.0/47
  real(16), parameter :: const_Vo = (h_barra**5)*(1.4**2.5)/5

  character(len=10), dimension(4) :: args
  character(len=40) :: archivo
  real(16), dimension(:,:), allocatable :: vector

  ! Variables el problema
  integer(4) :: N
  real(16) :: T
  real(16) :: rho
  real(16) :: L
  real(16) :: pmax

  ! Variables de muestreo
  integer(4) :: n_term = 2000
  integer(4) :: n_muestras = 500
  integer(4) :: n_desc = 100
  integer(4) :: pasos_term = 20
  real(16) :: factor = 0
  real(16) :: Tinit = 0.02
  real(16), dimension(:), allocatable :: distribucion

  ! Gilada
  real(16), dimension(:), allocatable :: LUT
  real(16) :: delta = 1E-2
  real(16) :: w = 1E2 ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V
  integer(4) :: i
  integer(4) :: j

  real(16) :: start
  real(16) :: stop

  read(args(2), *) N
  read(args(3), *) T
  read(args(4), *) L
  write(*,*) N,T,L
  write(*,*) "Muestras: ", n_muestras
  write(*,*) "Termalizacion:", n_term
  write(*,*) "Reescalamientos:", pasos_term
  write(*,*) "Descorrelacion:", n_desc
  write(archivo,"(A13,I3,A1,F6.4,A4)") 'distribucion_', N,'_', T,'.txt'
  write(*,*) archivo

  rho = N/(L**3)
!  V = ((6*m_pi*rho)**(2.0/3.0))*const_Vo  ! Sin spin
  V = ((3*m_pi*rho)**(2.0/3.0))*const_Vo  ! Con spin 1/2

  call init_random_seed()

  call leer_tablas(LUT,'tabla.txt')
  matriz = matrizhc(2*delta*w)

!  allocate(distribucion(Nbins))
  allocate(vector(N,12))
  call inicializar(vector,Tinit,L)
  factor = (T/2)**(1.0/pasos_term)

  call CPU_TIME(start)
  do i=1,pasos_term
    call avanzarN(vector,delta,V,L,LUT,matriz,n_term)   ! Tiempo de termalizacion
    vector(:,10:12) = factor*vector(:,10:12)
    vector(:,4:6) = factor*vector(:,4:6)
  end do
  call CPU_TIME(stop)
  write(*,*) "Termalizacion:", stop-start

  open(unit = 100, file=archivo)
! Metodo guardando los p
  call CPU_TIME(start)
  call avanzarN(vector,delta,V,L,LUT,matriz,n_term)   ! Tiempo de termalizacion
  write(100,*) EnergiaCinetica(vector)
  write(100,*) corr_T_virial(vector,L,LUT,V)
  do j=1,ubound(vector,1)
    !write(100,*) sum(vector(j,10:12)**2)
    write(100,*) vector(j,10)
    write(100,*) vector(j,11)
    write(100,*) vector(j,12)
  end do
  call CPU_TIME(stop)
  write(*,*) "Muestra 1: Tardo", stop-start
  factor = (stop-start)*(n_muestras-1)
  write(*,*) "Tiempo restante estimado:", factor
  do i = 2,n_muestras   ! Para cada muestra posterior tengo que descorrelacionar
    call avanzarN(vector,delta,V,L,LUT,matriz,n_desc)
    write(100,*) EnergiaCinetica(vector)
    write(100,*) corr_T_virial(vector,L,LUT,V)
    do j=1,ubound(vector,1)
      !write(100,*) sum(vector(j,10:12)**2)
      write(100,*) vector(j,10)
      write(100,*) vector(j,11)
      write(100,*) vector(j,12)
    end do
    write(*,*) "Muestra", i
  end do
  close(100)

end subroutine


! Para correr mientras se come Vitel Tone
subroutine main_temperatura(args)
  use observables
  use funciones
  use integrador
  use tablas
  use init
  use grabar

  real(16), parameter :: m_pi = 3.14159265359
  real(16), parameter :: h_barra = 25.0/47
  real(16), parameter :: const_Vo = (1.4**2.5)/5
  character(len=10), dimension(7) :: args
  character(len=30) :: name = "temperaturas_64.txt"

  real(16), dimension(:,:), allocatable :: vector
  real(16), dimension(:), allocatable :: T_viriales
  real(16), dimension(:), allocatable :: T_Boltzmann

  ! Variables el problema
  integer(4) :: N
  real(16) :: factor
  real(16) :: T_init
  real(16) :: T_final
  real(16) :: rho
  real(16) :: L = 8
  real(16) :: start
  real(16) :: stop

  ! Variables de muestreo
  integer(4) :: Nterm=100
  integer(4) :: Nmuestras=1000

  ! Gilada
  real(16), dimension(:), allocatable :: LUT
  real(16) :: delta = 1E-2
  real(16) :: w = 1E2 ! Parametro magico del integrador
  real(16), dimension(4,4) :: matriz
  real(16) :: V
  integer(4) :: i
  integer(4) :: j

  read(args(2), *) Nterm
  read(args(3), *) Nmuestras
  read(args(4), *) N
  read(args(5), *) T_init
  read(args(6), *) T_final
  read(args(7), *) N_temp
  write(*,*) N,T_init,T_final,N_temp

  rho = N/(L**3)
  V = 10*((3*m_pi*rho)**(2.0/3.0))*const_Vo

  call init_random_seed()

  call leer_tablas(LUT,'tabla.txt')
  matriz = matrizhc(2*delta*w)

  allocate(vector(N,12))
  allocate(T_Boltzmann(N_temp))
  allocate(T_viriales(N_temp))
  call inicializar(vector,T_init,L)
  T_Boltzmann(:) = 0
  T_viriales(:) = 0

  ! Termalizacion
  call CPU_TIME(start)
  call avanzarN(vector,delta,V,L,LUT,matriz,Nterm)
  call CPU_TIME(stop)
  write(*,*) "Segundo",stop-start,": Termalizacion"

  ! Primer temperatura
  do j=1,Nmuestras
    call avanzarN(vector,delta,V,L,LUT,matriz,Nterm/100)
    T_Boltzmann(1) = T_Boltzmann(1) + EnergiaCinetica(vector)/(1.5*N*Nmuestras)
    T_viriales(1) = T_viriales(1) + corr_T_virial(vector,L,LUT,V)/Nmuestras
  end do
  T_viriales(1) = T_viriales(1) + T_Boltzmann(1)
  call CPU_TIME(stop)
  write(*,*) "Segundo",stop-start,": T=", T_init

  do i = 2,N_temp
    factor = (T_init*(N_temp-1)-(i-1)*(T_init-T_final))/(T_init*(N_temp-1)-(i-2)*(T_init-T_final))
    call Reescalar_vel(vector,sqrt(factor))
    call avanzarN(vector,delta,V,L,LUT,matriz,Nterm/100)  ! Termalizo
    do j=1,Nmuestras
      call avanzarN(vector,delta,V,L,LUT,matriz,Nterm/100)
      T_Boltzmann(i) = T_Boltzmann(i) + EnergiaCinetica(vector)/(1.5*N*Nmuestras)
      T_viriales(i) = T_viriales(i) + corr_T_virial(vector,L,LUT,V)/Nmuestras
    end do
    T_viriales(i) = T_viriales(i) + T_Boltzmann(i)
    call CPU_TIME(stop)
    write(*,*) "Segundo",stop-start,": T=", T_init-(i-1)*(T_init-T_final)/(N_temp-1)
  end do

  open(unit = 100, file=name)
  do i = 1,N_temp
    write(100,*) T_init-(i-1)*(T_init-T_final)/(N_temp-1), ";", T_Boltzmann(i), ";",  T_viriales(i)
  end do
  close(100)

end subroutine
