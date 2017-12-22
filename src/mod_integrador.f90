module integrador

	use funciones

	implicit none

	private

		public :: TransformacionHa
		public :: TransformacionHb
		public :: TransformacionHc
		public :: matrizhc
		public :: avanzar
		!vector  (q,y,x,p)  [N,12]

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
			real(16) :: V_aux
			real(16) :: qij2
			real(16) :: yij2
			real(16) :: qij
			!real(16) :: qj
			real(16) :: yij
			real(16) :: inc
			!real(16) :: yj
			integer :: i ! indice de particula
			integer :: j ! indice de sumatoria
			integer :: k ! indice de componente vectorial
			integer :: N
			!q,y,x,p
			N = ubound(vector, dim=1)
			vector(:,7:9) = vector(:,7:9)+delta*vector(:,4:6)
			do i = 1, N
				! para x
				! q e y quedan igual
				do j = i+1, N
					qij2 = DistanciaCuad(vector,i,j,L,0)
					yij2 = DistanciaCuad(vector,i,j,L,1)
					V_aux = V*Valor_LUT_der(LUT,0.5*(qij2+yij2))
					! para x
					do k = 7, 9
						yij = vector(i,k-3)-vector(j,k-3)
						inc = delta*yij*V_aux
						!yj = vector(j,k-3)
						vector(i,k) = vector(i,k) - inc
						vector(j,k) = vector(j,k) + inc
					enddo
					! para p
					do k = 10, 12
						qij = vector(i,k-9)-vector(j,k-9)
						call CdCP(qij,L)
						inc = delta*qij*V_aux
						vector(i,k) = vector(i,k) + inc
						vector(j,k) = vector(j,k) - inc
				enddo
			enddo
		end do
		end subroutine

		subroutine TransformacionHb(vector,delta,V,L,LUT)
			real(16), dimension(:,:), intent(inout) :: vector
			real(16), dimension(:), intent(in) :: LUT
			real(16), intent(in) :: delta
			real(16), intent(in) :: V
			real(16), intent(in) :: L
			real(16) :: V_aux
			real(16) :: pij2
			real(16) :: xij2
			real(16) :: pij
			!real(16) :: pj
			real(16) :: xij
			real(16) :: inc
			!real(16) :: xj
			integer :: i ! indice de particula
			integer :: j ! indice de sumatoria
			integer :: k ! indice de componente vectorial
			integer :: N
			!q,y,x,p
			N = ubound(vector, dim=1)
			vector(:,1:3) = vector(:,1:3)+delta*vector(:,10:12)
			do i = 1, N
				! x y p quedan igual
				do j = i+1,N
					xij2 = DistanciaCuad(vector,i,j,L,2)
					pij2 = DistanciaCuad(vector,i,j,L,3)
					V_aux = V*Valor_LUT_der(LUT,0.5*(xij2+pij2))
					! para q
					do k = 1, 3
						pij = vector(i,k+9)-vector(j,k+9)
						!pj = vector(j,k+9)
						inc = delta*pij*V_aux
						vector(i,k) = vector(i,k) - inc
						vector(j,k) = vector(j,k) + inc
					enddo
					! para y
					do k = 4, 6
						xij = vector(i,k+3)-vector(j,k+3)
						call CdCP(xij,L)
						inc = delta*xij*V_aux
						vector(i,k) = vector(i,k) + inc
						vector(j,k) = vector(j,k) - inc
					enddo
				enddo
			enddo
		end subroutine

	subroutine TransformacionHc(vector,matrizhc,Lado)
		real(16), dimension(:,:), intent(inout) :: vector
	  real(16), dimension(12) :: nvector
	  real(16), dimension(4,4), intent(in) :: matrizhc
		real(16), intent(in) :: Lado
	  integer :: i
	  integer :: j
	  integer :: k
	  integer :: l

		nvector(:)=0
	  do i=1,ubound(vector, dim=1)
			nvector(:)=0
	    do l=1,3
	      do k=1,4
	        nvector(l+0)=nvector(l+0)+0.5*matrizhc(k,1)*vector(i,l+3*(k-1))
	        nvector(l+3)=nvector(l+3)+0.5*matrizhc(k,2)*vector(i,l+3*(k-1))
	        nvector(l+6)=nvector(l+6)+0.5*matrizhc(k,3)*vector(i,l+3*(k-1))
	        nvector(l+9)=nvector(l+9)+0.5*matrizhc(k,4)*vector(i,l+3*(k-1))
	      end do
	    end do
		!aca tenia que ordenar de nuevo  yo tenia (qpxy)--->(qyxp)
		!entonces q-->q//p-->y...etc
			do j=1,3
	    	vector(i,j)=nvector(j)
				vector(i,j+3)=nvector(j+9)
	    	vector(i,j+6)=nvector(j+6)
				vector(i,j+9)=nvector(j+3)
			end do
	  end do
	  end subroutine

	 function matrizhc(arg)
		real(16), dimension(4,4)  :: matrizhc
		real(16), intent(in) :: arg
		matrizhc=0

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
		integer(4) :: i
		integer(4) :: j
		call TransformacionHa(vector,delta*0.5,V,L,LUT)
	  call TransformacionHb(vector,delta*0.5,V,L,LUT)
		!call CCperiod2(vector,L)
		call TransformacionHc(vector,matriz,L)
		!call CCperiod2(vector,L)
		call TransformacionHb(vector,delta*0.5,V,L,LUT)
		call TransformacionHa(vector,delta*0.5,V,L,LUT)
		call CCperiod2(vector,L)
	end subroutine

end module
