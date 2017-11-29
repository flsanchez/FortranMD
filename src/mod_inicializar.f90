module init

implicit none
private

	public ::	inicializar
	public :: Rand_Gauss

contains

	subroutine inicializar(vector, T, L)
		real(16), dimension(:,:), intent(inout) :: vector
		real(16), intent(in) :: T
		real(16), intent(in) :: L
		integer(4) :: N
		integer(4) :: m
		real(16) :: a
		integer(4) :: x = 1
		integer(4) :: y = 1
		integer(4) :: z = 1
		integer(4) :: i
		N = ubound(vector, dim=1)
		m = N**(1.0/3)
		a = L/m

		do x=1,m
			do y=1,m
				do z=1,m
					i = x+(y-1)*m+m*m*(z-1);
					vector(i,1) = (a/2) + (x-1)*a;
					vector(i,2) = (a/2) + (y-1)*a;
					vector(i,3) = (a/2) + (z-1)*a;
					vector(i,7) = (a/2) + (x-1)*a;
					vector(i,8) = (a/2) + (y-1)*a;
					vector(i,9) = (a/2) + (z-1)*a;
				end do
			end do
		end do

		do i=1,N
			vector(i,4) = Rand_Gauss(T)
			vector(i,5) = Rand_Gauss(T)
			vector(i,6) = Rand_Gauss(T)
		end do
		vector(:,4) = vector(:,4) - sum(vector(:,4))/N
		vector(:,5) = vector(:,5) - sum(vector(:,5))/N
		vector(:,6) = vector(:,6) - sum(vector(:,6))/N
		vector(:,10:12) = vector(:,4:6)
	end subroutine



	real(16) function Rand_Gauss(T)
		real(16), intent(in) :: T
	  integer(4) ::  repe = 24
	  real(16) :: s
		integer :: i
		s = sqrt(12*T/repe);  ! s es la normalizacion EPS = 1
		Rand_Gauss = sum((/ (rand()-0.5, i=1,repe) /))*s
	end function

end module
