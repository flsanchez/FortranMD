module integrador

	use funciones

	implicit none

	private

		public :: TransformacionHa
		public :: TransformacionHb
		public :: TransformacionHc
		public :: matrizhc
		public :: avanzar
		!vector  (q,y,x,p)  [12,N]

	contains
		! vector(i,j): i es numero de particula (1 <= i <= N),
		! j es q,y,x,p (1 <= j <= 12),
		! dim(vector(i,j)) = Nx12
		! el V0 = (3pi^2)^(2/3)/5*(7/5)^(5/2)*h_barra^2*rho^(2/3)

		subroutine TransformacionHa(vector,delta,V,L,LUT)
			real(16), dimension(:,:), intent(inout) :: vector
			real(16), dimension(:), intent(in) :: LUT
			real(16), intent(in) :: delta
			real(16), intent(in) :: V
			real(16), intent(in) :: L
			real(16) :: qij2
			real(16) :: yij2
			real(16) :: qi
			real(16) :: qj
			real(16) :: yi
			real(16) :: yj
			integer :: i ! indice de particula
			integer :: j ! indice de sumatoria
			integer :: k ! indice de componente vectorial
			integer :: N

			N = ubound(vector, dim=1)
			do i = 1, N
				! q e y quedan igual
				do j = 1, N
					qij2 = DistanciaCuad(vector,i,j,L,0)
					yij2 = DistanciaCuad(vector,i,j,L,1)
					! para x
					do k = 7, 9
						yi = vector(i,k-3)
						yj = vector(j,k-3)
						vector(i,k) = vector(i,k) - delta/2*(yi-2*V*(yi-yj)*Valor_LUT(LUT,qij2+yij2))
					enddo
					! para p
					do k = 10, 12
						qi = vector(i,k-9)
						qj = vector(j,k-9)
						vector(i,k) = vector(i,k) + delta*V*(qi-qj)*Valor_LUT(LUT,qij2+yij2)
					enddo
				enddo
			enddo
		end subroutine

		subroutine TransformacionHb(vector,delta,V,L,LUT)
			real(16), dimension(:,:), intent(inout) :: vector
			real(16), dimension(:), intent(in) :: LUT
			real(16), intent(in) :: delta
			real(16), intent(in) :: V
			real(16), intent(in) :: L
			real(16) :: pij2
			real(16) :: xij2
			real(16) :: pi
			real(16) :: pj
			real(16) :: xi
			real(16) :: xj
			integer :: i ! indice de particula
			integer :: j ! indice de sumatoria
			integer :: k ! indice de componente vectorial
			integer :: N
			!q,y,x,p
			N = ubound(vector, dim=1)
			do i = 1, N
				! x y p quedan igual
				do j = 1, N
					pij2 = DistanciaCuad(vector,i,j,L,3)
					xij2 = DistanciaCuad(vector,i,j,L,2)
					! para q
					do k = 1, 3
						pi = vector(i,k+9)
						pj = vector(j,k+9)
						vector(i,k) = vector(i,k) + delta/2*(pi-2*V*(pi-pj)*Valor_LUT(LUT,pij2+xij2))
					enddo
					! para y
					do k = 4, 6
						xi = vector(i,k+3)
						xj = vector(j,k+3)
						vector(i,k) = vector(i,k) + delta*V*(xi-xj)*Valor_LUT(LUT,xij2+pij2)
					enddo
				enddo
			enddo
		end subroutine

	  subroutine TransformacionHc(vector,matrizhc)
		real(16), dimension(:,:), intent(inout) :: vector
	  real(16), dimension(12) :: nvector
	  real(16), dimension(4,4), intent(in) :: matrizhc
	  integer :: i
	  integer :: j
	  integer :: k
	  integer :: l
	  integer :: m
	  integer :: n

	  nvector=0
	  do i=1,ubound(vector, dim=1)
	    do l=1,3
	      do k=1,4
	        nvector(l+0)=nvector(l+0)+0.5*matrizhc(k,1)*vector(l+3*(k-1),i)
	        nvector(l+3)=nvector(l+3)+0.5*matrizhc(k,2)*vector(l+3*(k-1),i)
	        nvector(l+6)=nvector(l+6)+0.5*matrizhc(k,3)*vector(l+3*(k-1),i)
	        nvector(l+9)=nvector(l+9)+0.5*matrizhc(k,4)*vector(l+3*(k-1),i)
	      end do
	    end do
		!aca tenia que ordenar de nuevo  yo tenia (qpxy)--->(qyxp)
		!entonces q-->q//p-->y...etc
			do j=1,3
	    	vector(i,j)=nvector(j)
	    	vector(i,j+6)=nvector(j+6)
			end do
			do m=1,3
	    	vector(i,m+9)=nvector(m+3)
			end do
			do n=1,3
	    	vector(i,m+3)=nvector(m+9)
			end do
	  end do
	  end subroutine

	 function matrizhc(arg)
		real(16), dimension(4,4)  :: matrizhc
		real(16), intent(in) :: arg
		matrizhc=0
!		matrizhc(1,1)=1+cos(arg)
!		matrizhc(3,3)=1+cos(arg)
!		matrizhc(2,4)=1+cos(arg)
!		matrizhc(4,2)=1+cos(arg)

!		matrizhc(1,2)=-sin(arg)
!		matrizhc(2,1)=-sin(arg)
!		matrizhc(3,4)=-sin(arg)
!		matrizhc(4,3)=-sin(arg)

!		matrizhc(1,4)=sin(arg)
!		matrizhc(4,1)=sin(arg)
!		matrizhc(2,3)=sin(arg)
!		matrizhc(3,2)=sin(arg)

!		matrizhc(1,3)=1-cos(arg)
!		matrizhc(2,2)=1-cos(arg)
!		matrizhc(3,1)=1-cos(arg)
!		matrizhc(4,4)=1-cos(arg)

!creo que asi se solucionaba pero lo hice de la manera chota de arriba  igual habira que testear esto

		matrizhc(1,:)=(/ 1+cos(arg), -sin(arg)  , 1-cos(arg), sin(arg) /)
		matrizhc(2,:)=(/ -sin(arg) ,  1-cos(arg), sin(arg),1+cos(arg) /)
		matrizhc(3,:)=(/ 1-cos(arg),  sin(arg)  , 1+cos(arg),-sin(arg) /)
		matrizhc(4,:)=(/ sin(arg)  , 1+cos(arg) , -sin(arg),1-cos(arg) /)
	end function

	subroutine avanzar(vector,delta,V,L,LUT,matriz)
		real(16), dimension(:,:), intent(inout) :: vector
		real(16), dimension(:,:), intent(in) :: matriz
		real(16), dimension(:), intent(in) :: LUT
		real(16), intent(in) :: delta
		real(16), intent(in) :: V
		real(16), intent(in) :: L
		call TransformacionHa(vector,delta*0.5,V,L,LUT)
	  call TransformacionHb(vector,delta*0.5,V,L,LUT)
		call TransformacionHc(vector,matriz)
		call TransformacionHb(vector,delta*0.5,V,L,LUT)
		call TransformacionHa(vector,delta*0.5,V,L,LUT)
	end subroutine

end module
