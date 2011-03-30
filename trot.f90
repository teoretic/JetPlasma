module trot

    use names

  implicit none
  ! ===================================
  ! Вычисление матрицы преобразования T
  ! Toro, с. 573, раздел 16.7.3
  !       с. 103, раздел 3.2.1
  contains
  subroutine T_rot(n, T, Ti)
    implicit none

    real (kind=precis), dimension(1:2), intent(in):: n ! n = (cos(theta), sin(theta)) -- вектор нормали

    real (kind=precis), dimension(1:5,1:5), intent(out):: T ! матрица преобразования
    real (kind=precis), dimension(1:5,1:5), intent(out):: Ti ! обратная к ней

	T = 0.0

	T(1,1)=1.0
	T(4,4)=1.0
	T(5,5)=1.0
	T(2,2)=n(1)
	T(2,3)=n(2)
	T(3,2)=-n(2)
	T(3,3)=n(1)
    Ti = transpose(T)

  end subroutine T_rot

    function pres95(node) result (pressure)
  	real (kind=precis) :: pressure
  	real (kind=precis), dimension(2) :: node

  	real (kind=precis) :: rad, c1,c2

  	r_0 = sqrt(r_max**2+z_max**2)
  	c1 = p_0 - G*rho_init/r_0
  	c2 = c1 + 2.0*G*rho_init/gr_rad
  	rad = sqrt(sum(node**2))
 ! 	print *, c1, c2
!  	read *
  	if (rad<gr_rad) then
  			pressure = c2 - rho_init*G*rad/(gr_rad**2)
  		else
  			pressure = c1+G*rho_init/rad
  	end if

  	if (pressure<0.0) then
print *, 'ahhhaaa'
  		print *, pressure, r_max, z_max,r_0, node, rad, rho_init*G*rad/(gr_rad**2), G*rho_init/rad
  		read *
  	end if
  end function pres95

end module trot