module geometry

use names

implicit none

contains
  ! ===============================================
  ! Расчет правой части -- только для (r,z) геометрии
  ! Физические источники учитываются отдельно
  ! U -- вектор консервативных переменных
  ! SS -- правая часть для учета (r,z) геометрии
  ! если alpha = 0, то H должен возвращаться нулевым

  ! геометрия x,r - ось вращения х, координаты x,r,phi
  function H_GEOM(U) result (SS)
    implicit none

    real (kind=precis), dimension(1:5):: U
    real (kind=precis), dimension(1:5):: SS

    real (kind=precis) :: rho1, u1,w1,v1 ,k_e1, p1, E1 !I
 !   real (kind=precis), parameter:: gmm = 7.0/5.0 !gmm = 1.4 ! показатель адиабаты

    SS = 0.0

!I ->
    rho1 = U(1)
    u1 = U(2)/rho1
    v1 = U(3)/rho1
    w1 = U(4)/rho1
    E1 = U(5)

    k_e1 = u1*u1+v1*v1+w1*w1 !check
    p1 = (E1-0.5*rho1*k_e1)*(gmm-1.0)
! радиальная часть во второй компоненте скорости
    SS(1) = rho1*v1
    SS(2) = rho1*u1*v1
    SS(3) = rho1*(v1*v1-w1*w1)
    SS(4) = rho1*v1*w1
    SS(5) = (E1+p1)*v1

!I end
  end function H_GEOM

  function H_MAGN(U,B) result (SS)
    implicit none

    real (kind=precis), dimension(1:5):: U
	real (kind=precis), dimension(1:3):: B
    real (kind=precis), dimension(1:5):: SS

    real (kind=precis) :: rho1, u1,w1,v1 ,k_e1, p1, E1,B2 !I
 !   real (kind=precis), parameter:: gmm = 7.0/5.0 !gmm = 1.4 ! показатель адиабаты

    SS = 0.0

!I ->
    rho1 = U(1)
    u1 = U(2)/rho1
    v1 = U(3)/rho1
    w1 = U(4)/rho1
    E1 = U(5)
	B2 = B(1)*B(1) + B(2)*B(2) + B(3)*B(3)

    k_e1 = u1*u1+v1*v1+w1*w1 !check
    p1 = (E1-0.5*rho1*k_e1)*(gmm-1.0)
! радиальная часть во второй компоненте скорости
    SS(1) = 0.0
    SS(2) = B(1)*B(2)
    SS(3) = B(2)*B(2) - B(3)*B(3)
    SS(4) = B(2)*B(3)
    SS(5) = B(2)*(u1*B(1)+v1*B(2)+w1*B(3)) - v1*B2

	SS = -SS/(4.0*pi)

!I end
  end function H_MAGN

  function Right(n) result (SS)
    implicit none

	integer :: n
    real (kind=precis), dimension(1:5):: SS
	real (kind=precis) :: r, t

	r = elems_center(2,n)
	t = t_current

! радиальная часть во второй компоненте скорости
    SS(1) = 3.*(r**2)*t+t+r*t+t*(r**3)
    SS(2) = (r**2)*(pi*(r**2)+pi+6.0*pi*(r**2)*(t**2)+4.0*pi*(t**2)-1.)/pi
    SS(3) = r*(2.0*pi*(r**2)+6.0*pi+8.0*pi*(r**2)*(t**2)+4.0*pi*(t**2)+r**2)/(2.0*pi)
    SS(4) = r*(4.0*pi*(r**2)+4.0*pi+20.0*(r**2)*(t**2)*pi+12.0*pi*(t**2)-3.)/(4.0*pi)
    SS(5) = t*(-9.0*(r**4)*(t**2)+2.0*gmm-2.0*(r**2)-r**6-3*(r**4)+9.0*gmm*(r**4)*(t**2)+4*gmm*(r**2)*(t**2)+4.0*(r**6)*(t**2)*gmm-4.0*(r**2)*(t**2)-4.0*(r**6)*(t**2)+6.0*gmm*(r**2)+gmm*(r**6)+3.0*gmm*(r**4))/(gmm-1.0)

!I end
  end function Right

end module geometry