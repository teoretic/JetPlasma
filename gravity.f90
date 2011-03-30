! ==========================================================
! Содержит процедyры для гравитации
! ==========================================================

module gravity

  use names

  implicit none

  real (kind=precis) :: GR_x = 0.0
  real (kind=precis) :: GR_y = 0.0
  real (kind=precis) :: GR_z = 0.0 !-9.8

contains

pure function GRAV(Num, pU) result (GR)

  implicit none

  real (kind=precis), dimension(1:5), Intent(In):: pU
  integer, Intent(In) :: Num
  real (kind=precis), dimension(1:5):: GR

  real (kind=precis) :: dist, Mass, rho, zet, rad

  GR(:) = 0.0
  rho = pU(1)
  zet = elems_center(1,Num)
  rad = elems_center(2,Num)
  dist = sqrt(sum(elems_center(:,Num)**2))
  if (dist>gr_rad) then
    GR(2) = -G*rho*zet/(dist**3)
    GR(3) = -G*rho*rad/(dist**3)
  else
    GR(2) = -G*rho*zet/(gr_rad**3)
    GR(3) = -G*rho*rad/(gr_rad**3)
  end if
  GR(5) = dot_product(pU(2:3),GR(2:3))/pU(1)
end function GRAV

end module gravity
