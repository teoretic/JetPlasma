module initials

use names
use trot
use radiationTR
use solution

implicit none

contains

  ! ============================================
  ! процедура задает начальные условия
  ! вызывается один раз в начале программы перед
  ! циклом пересчета по времени
  ! -- только для внутренних ячеек !!!

  ! см. пример на странице 128 (Chapter 7) "Part 2_3"
  subroutine init_ic
!	call alfven
!	call axial95	!main
!	call rad_init
    call radjet2D
!    call simple

!	call stairs

!	call explosion
!	call test_visc
!	call dissip
!	call riemann_pb
!	call rotor_pb
!	call orszag
!	call blast
!	call cloud_shock_pb
  end subroutine init_ic

  subroutine simple
    integer:: I

    gmm = 1.4

    do I = 1,nElems
      U(1,I) = 1.0
!U(1,I) = 1.0
      U(2,I) = 0.0
      U(3,I) = 0.0
      U(4,I) = 0.0
      U(5,I) = 0.125 !1.7857 ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
    end do
  end subroutine simple

  function rad_test_dens(x) result (den)
    real (kind=precis) :: x, den

    real (kind=precis) :: rho0, rho1

    rho0 = 0.5
    rho1 = 1.0

    den = (rho1-rho0)*Exp(-600.*((x-0.2)**2))+rho0
  end function rad_test_dens

  function rad_test_pres(x) result (pres)
    real (kind=precis) :: x, pres

    pres = 2.*Pi*int0*(1.-exp(-kkap*x))/(3.0*lsp) + p0
  end function rad_test_pres

  subroutine rad_test
    integer :: I
    real (kind=precis) :: pres

    gmm = 3.0
    kkap = 1.0
    int0 = 0.05
!    int0 = 0.0
    p0 = 0.1
    vell = 0.2

    do I = 1,nElems
      U(1,I) = rad_test_dens(elems_center(2,I))
!U(1,I) = 1.0
      U(2,I) = 0.0
      U(3,I) = rad_test_dens(elems_center(2,I))*lsp
      U(4,I) = 0.0
      U(5,I) = U(3,I)*U(3,I)*0.5/U(1,I) + rad_test_pres(elems_center(2,I))/(gmm-1.0) !1.7857 ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))

      intensity_hat(i) = 2.*Pi*int0*exp(-kkap*elems_center(2,I))
      Poynting_hat(:,i) = (/0.0_precis, Pi*int0*exp(-kkap*elems_center(2,I)), 0.0_precis/)
      radstress(1,:,i) = (/ 2.*Pi*exp(-kkap*elems_center(2,I))/3., 0.0_precis, 0.0_precis /)
      radstress(2,:,i) = (/ 0.0_precis, 2.*Pi*exp(-kkap*elems_center(2,I))/3., 0.0_precis /)
      radstress(3,:,i) = (/ 0.0_precis, 0.0_precis, 2.*Pi*exp(-kkap*elems_center(2,I))/3. /)
      radstress(:,:,i) = int0*radstress(:,:,i)/lsp
      intensity(i) = intensity_hat(i)
      Poynting(:,i) = Poynting_hat(:,i)
    end do
  end subroutine rad_test

  subroutine rad_init
    integer:: I

    gmm = 1.4

    do I = 1,nElems
      U(1,I) = 0.50
!U(1,I) = 1.0
      U(2,I) = 0.0
      U(3,I) = 0.0
      U(4,I) = 0.0
      U(5,I) = 0.125 !1.7857 ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
    end do

    call radiation_solverTR

    do i=1,nElems
      intensity(i) = intensity_hat(i)
      Poynting(:,i) = Poynting_hat(:,i)
    end do
!read *
  end subroutine rad_init

  subroutine stairs
    integer:: I

    gmm = 1.4

    do I = 1,nElems
      U(1,I) = 0.50
!U(1,I) = 1.0
      U(2,I) = 0.0
      U(3,I) = 0.0
      U(4,I) = 0.0
      U(5,I) = 0.125 !1.7857 ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
    end do
  end subroutine stairs

  subroutine explosion
    integer:: I
    real (kind=precis) :: p, rho, r0,r
    real (kind=precis), dimension(1:2) :: centr = (/1.,1./)

    gmm = 1.4

    rho = 0.125
    p = 0.1
    do I = 1,nElems
      U(1,I) = 0.125
!U(1,I) = 1.0
      U(2,I) = 0.0
      U(3,I) = 0.0
      U(4,I) = 0.0
      U(5,I) = p/(gmm-1.0) !1.7857 ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
    end do

    rho = 1.0
    p = 1.
    r0 = 0.4
    do I = 1,nElems
      r = distance(elems_center(:,i), centr)
      if (r<r0) then
	U(1,I) = rho !rho
	U(5,I) = p/(gmm-1.0)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
      end if
    end do
  end subroutine explosion

  subroutine blast

    Character(Len=10) :: buf, buf1
	real (kind=precis) :: rho, p, vx, vy, vz, Bx, By, Bz
	real (kind=precis), dimension(1:2) :: centr = (/0.5,0.75/), center
    real (kind=precis) :: r0, r1, r
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

    gmm = 5./3.
	p = 0.1
	rho = 1.0
	Bx = sqrt(4.*pi)/sqrt(2.)
	By = sqrt(4.*pi)/sqrt(2.)
	Bz = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0
	t_final = 1.0

!!!!начальные условия в ячейках
	 do I = 1,nElems
		  U(1,I) = rho !rho
		  U(2,I) = rho*vx  !rho*u
		  U(3,I) = rho*vy  !rho*v
		  U(4,I) = rho*vz  !rho*w
		  U(5,I) = p/(gmm-1.0)+0.5*rho*(vx*vx+vy*vy+vz*vz)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = Bx !Bx
	      B(2,I) = By !By
	      B(3,I) = Bz !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	      V_nodes(1,I) = vx
		  V_nodes(2,I) = vy
		  B_nodes(1,I) = Bx !Bx
	      B_nodes(2,I) = By  !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	      B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
          B_edges(2,I) = Bz
		  V_edges(1,I) = vx*edge_n(1,I) + vy*edge_n(2,I) !Bn
          V_edges(2,I) = vz
	 end do

	p = 10.
	r0 = 0.1
	 do I = 1,nElems
		r = distance(elems_center(:,i), centr)
		if (r<r0) then
		  U(1,I) = rho !rho
		  U(2,I) = 0.0  !rho*u
		  U(3,I) = 0.0  !rho*v
		  U(5,I) = p/(gmm-1.0)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
		end if
	 end do
  end subroutine blast


  subroutine orszag

    Character(Len=10) :: buf, buf1
	real (kind=precis) :: rho, p, vx, vy, vz, Bx, By, Bz,x,y
	real (kind=precis), dimension(1:2) :: centr = (/0.5,0.5/), center
    real (kind=precis) :: r0, r1, r
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

    gmm = 5./3.
	p = 5./(12.*pi)
	rho = 25./(36.*pi)
	Bx = 5.0
	By = 0.0
	Bz = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0
	t_final = 0.5


!!!!начальные условия в ячейках
	 do I = 1,nElems
		x = elems_center(1,I)
		y = elems_center(2,I)
		  U(1,I) = rho !rho
		  U(2,I) = rho*(-sin(2.*pi*y))  !rho*u
		  U(3,I) = rho*sin(2.*pi*x)  !rho*v
		  U(4,I) = 0.0  !rho*w
		  U(5,I) = p/(gmm-1.0)+0.5*(sum(U(2:4,I)**2))/U(1,I)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = -sin(2.*pi*y) !Bx
	      B(2,I) = sin(4.*pi*x) !By
	      B(3,I) = Bz !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
		x = nodes(1,I)
		y = nodes(2,I)
	      V_nodes(1,I) = -sin(2.*pi*y)
		  V_nodes(2,I) = sin(2.*pi*x)
		  B_nodes(1,I) = -sin(2.*pi*y) !Bx
	      B_nodes(2,I) = sin(4.*pi*x)  !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
		x = 0.5*(nodes(1,edges(1,I))+nodes(1,edges(2,I)))
		y = 0.5*(nodes(2,edges(1,I))+nodes(2,edges(2,I)))
	      B_edges(1,I) = -sin(2.*pi*y)*edge_n(1,I) + sin(4.*pi*x)*edge_n(2,I) !Bn
          B_edges(2,I) = 0.0
		  V_edges(1,I) = -sin(2.*pi*y)*edge_n(1,I) + sin(2.*pi*x)*edge_n(2,I) !Bn
          V_edges(2,I) = 0.0
	 end do
  end subroutine orszag


  subroutine rotor_pb

    Character(Len=10) :: buf, buf1
	real (kind=precis) :: rho, p, vx, vy, vz, v0, Bx, By, Bz
	real (kind=precis), dimension(1:2) :: centr = (/0.5,0.5/), center
    real (kind=precis) :: r0, r1, r
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

!	t_final = 0.15
    t_final = 0.295
!    gmm = 1.4
	gmm = 5./3.
!	p = 1.0
	p = 0.5
	rho = 1.0
!	Bx = 5.0
	Bx = 2.5
	By = 0.0
	Bz = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0
!	v0 = 2.0
	v0 = 1.0


!!!!начальные условия в ячейках
	 do I = 1,nElems
		  U(1,I) = rho !rho
		  U(2,I) = rho*vx  !rho*u
		  U(3,I) = rho*vy  !rho*v
		  U(4,I) = rho*vz  !rho*w
		  U(5,I) = p/(gmm-1.0)+0.5*rho*(vx*vx+vy*vy+vz*vz)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = Bx !Bx
	      B(2,I) = By !By
	      B(3,I) = Bz !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	      V_nodes(1,I) = vx
		  V_nodes(2,I) = vy
		  B_nodes(1,I) = Bx !Bx
	      B_nodes(2,I) = By  !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	      B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
          B_edges(2,I) = Bz
		  V_edges(1,I) = vx*edge_n(1,I) + vy*edge_n(2,I) !Bn
          V_edges(2,I) = vz
	 end do

	rho = 10.0
	r0 = 0.1
	r1 = 0.115
	 do I = 1,nElems
		r = distance(elems_center(:,i), centr)
		if (r<r0) then
		  U(1,I) = rho !rho
		  U(2,I) = rho*(-v0/r0)*(elems_center(2,i)-0.5)  !rho*u
		  U(3,I) = rho*(v0/r0)*(elems_center(1,i)-0.5)  !rho*v
		  U(5,I) = p/(gmm-1.0)+0.5*sum(U(2:4,I)**2)/U(1,I)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
		else
		  if (r<r1) then
			U(1,I) = 1.0+9.0*ff(r) !rho
			U(2,I) = -U(1,I)*ff(r)*v0*(elems_center(2,i)-0.5)/r0  !rho*u
			U(3,I) = U(1,I)*ff(r)*v0*(elems_center(1,i)-0.5)/r0!rho*v
			U(5,I) = p/(gmm-1.0)+0.5*sum(U(2:4,I)**2)/U(1,I)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
		  end if
		end if
	 end do

!!!! начальные данные в точках
	 do I=1,nNodes
	   r = distance(nodes(:,I), centr)
		if (r<r0) then
		  V_nodes(1,I)= (-v0/r0)*(nodes(2,i)-0.5)  !rho*u
		  V_nodes(2,I)= (v0/r0)*(nodes(1,i)-0.5)  !rho*v
		else
		  if (r<r1) then
			V_nodes(1,I) = -ff(r)*v0*(nodes(2,i)-0.5)/r0  !rho*u
			V_nodes(2,I) = ff(r)*v0*(nodes(1,i)-0.5)/r0 !rho*v
		  end if
		end if
	 end do
!!!! начальные данные на гранях
	 do I=1,nEdges
	   center = 0.0
	   do J=1,2
	     center=center+nodes(:,edges(J,I))*0.5
	   end do
	   r = distance(center, centr)
		if (r<r0) then
		  vx = (-v0/r0)*(center(2)-0.5)  !rho*u
		  vy = (v0/r0)*(center(1)-0.5)  !rho*v
		  V_edges(1,I) = vx*edge_n(1,I) + vy*edge_n(2,I) !Bn
          V_edges(2,I) = 0.0
		else
		  if (r<r1) then
			vx = -ff(r)*v0*(center(2)-0.5)/r0  !rho*u
			vy = ff(r)*v0*(center(1)-0.5)/r0 !rho*v
			V_edges(1,I) = vx*edge_n(1,I) + vy*edge_n(2,I) !Bn
            V_edges(2,I) = 0.0
		  end if
		end if
	 end do
  end subroutine rotor_pb

  real (kind=precis) function ff(x)
    real (kind=precis) :: x
    real (kind=precis) :: r0=0.1, r1=0.115

    ff = (r1-x)/(r1-r0)
  end function ff

  subroutine riemann_pb
    Character(Len=10) :: buf, buf1
	real (kind=precis) :: rhoL, pL, uL, vL, wL, BxL, ByL, BzL
    real (kind=precis) :: rhoR, pR, uR, vR, wR, BxR, ByR, BzR
	real (kind=precis) :: bound = 0.5
	real (kind=precis), dimension(1:2) :: disc_center

	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J
	 real (kind=precis)   :: centr

!	t_final = 0.1
t_final = 10.
	open(unit=100, file='./data/input.dat', form='formatted')
	!read the input settings from input file  rho, p, u, v, w, By, Bz
    read (100,*) buf, buf1, gmm
    read (100,*) buf, buf1, rhoL
    read (100,*) buf, buf1, pL
    read (100,*) buf, buf1, uL
    read (100,*) buf, buf1, vL
    read (100,*) buf, buf1, wL
    read (100,*) buf, buf1, BxL
	read (100,*) buf, buf1, ByL
    read (100,*) buf, buf1, BzL
    read (100,*) buf, buf1, rhoR
    read (100,*) buf, buf1, pR
    read (100,*) buf, buf1, uR
    read (100,*) buf, buf1, vR
    read (100,*) buf, buf1, wR
	read (100,*) buf, buf1, BxR
    read (100,*) buf, buf1, ByR
    read (100,*) buf, buf1, BzR
    close (unit=100,status='Keep')

    BxL = sqrt(4.0*pi)*BxL
	BxR = sqrt(4.0*pi)*BxR
	ByL = sqrt(4.0*pi)*ByL
	ByR = sqrt(4.0*pi)*ByR
	BzL = sqrt(4.0*pi)*BzL
	BzR = sqrt(4.0*pi)*BzR

!   BxL = 0.0
!	BxR = 0.0
!	ByL = 0.0
!	ByR = 0.0
!	BzL = 0.0
!	BzR = 0.0

!	gmm = 5.0/3.0
!!!!начальные условия в ячейках
	 do I = 1,nElems
	   centr=0.
	   do J=1,3
	     centr=centr+nodes(2,elems(J,I))/3.0
	   end do
	   if (centr<bound) then
		  U(1,I) = rhoL !rho
!		  U(1,I) = 1.0
		  U(2,I) = rhoL*uL  !rho*u
		  U(3,I) = rhoL*vL  !rho*v
		  U(4,I) = rhoL*wL  !rho*w
		  U(5,I) = pL/(gmm-1.0)+0.5*rhoL*(uL*uL+vL*vL+wL*wL)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = BxL !Bx
	      B(2,I) = ByL  !By
	      B(3,I) = BzL  !Bz
		else
		  U(1,I) = rhoR !rho
!		  U(1,I) = 1.0
		  U(2,I) = rhoR*uR  !rho*u
		  U(3,I) = rhoR*vR  !rho*v
		  U(4,I) = rhoR*wR  !rho*w
		  U(5,I) = pR/(gmm-1.0)+0.5*rhoR*(uR*uR+vR*vR+wR*wR)!+(Bx*Bx+ByR*ByR+BzR*BzR)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = BxR !Bx
	      B(2,I) = ByR  !By
	      B(3,I) = BzR  !Bz
	   end if
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	   	if (nodes(2,I)<bound) then
	      V_nodes(1,I) = uL
		  V_nodes(2,I) = vL
		  B_nodes(1,I) = BxL !Bx
	      B_nodes(2,I) = ByL  !By
		else
		  V_nodes(1,I) = uR
		  V_nodes(2,I) = vR
	      B_nodes(1,I) = BxR !Bx
	      B_nodes(2,I) = ByR  !By
	   end if

	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	   centr=0.0
	   do J=1,2
	     centr=centr+nodes(2,edges(J,I))*0.5
	   end do
	   if (centr<bound) then
	      B_edges(1,I) = BxL*edge_n(1,I) + ByL*edge_n(2,I) !Bn
          B_edges(2,I) = BzL
		  V_edges(1,I) = uL*edge_n(1,I) + vL*edge_n(2,I) !Bn
          V_edges(2,I) = wL
		else
	      B_edges(1,I) = BxR*edge_n(1,I) + ByR*edge_n(2,I) !Bn
		  B_edges(2,I) = BzR
		  V_edges(1,I) = uR*edge_n(1,I) + vR*edge_n(2,I) !Bn
          V_edges(2,I) = wR
	   end if
	 end do

  end subroutine riemann_pb

  subroutine cloud_shock_pb
    Character(Len=10) :: buf, buf1
	real (kind=precis) :: rhoL, pL, uL, vL, wL, BxL, ByL, BzL
    real (kind=precis) :: rhoR, pR, uR, vR, wR, BxR, ByR, BzR
	real (kind=precis) :: bound = 0.5
	real (kind=precis), dimension(1:2) :: disc_center
	integer :: nx
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J
	 real (kind=precis)   :: centr

	open(unit=100, file='./data/input.dat', form='formatted')
	!read the input settings from input file  rho, p, u, v, w, By, Bz
    read (100,*) buf, buf1, gmm
    read (100,*) buf, buf1, rhoL
    read (100,*) buf, buf1, pL
    read (100,*) buf, buf1, uL
    read (100,*) buf, buf1, vL
    read (100,*) buf, buf1, wL
    read (100,*) buf, buf1, BxL
	read (100,*) buf, buf1, ByL
    read (100,*) buf, buf1, BzL
    read (100,*) buf, buf1, rhoR
    read (100,*) buf, buf1, pR
    read (100,*) buf, buf1, uR
    read (100,*) buf, buf1, vR
    read (100,*) buf, buf1, wR
	read (100,*) buf, buf1, BxR
    read (100,*) buf, buf1, ByR
    read (100,*) buf, buf1, BzR
    close (unit=100,status='Keep')

    BxL = sqrt(4.0*pi)*BxL
	BxR = sqrt(4.0*pi)*BxR
	ByL = sqrt(4.0*pi)*ByL
	ByR = sqrt(4.0*pi)*ByR
	BzL = sqrt(4.0*pi)*BzL
	BzR = sqrt(4.0*pi)*BzR

	gmm = 5.0/3.0
	t_final = 5.0

	nx = 1
    disc_center=(/0.8,0.5/)
	bound = 0.6
!!!!начальные условия в ячейках
	 do I = 1,nElems
	   centr=0.
	   do J=1,3
	     centr=centr+nodes(nx,elems(J,I))/3.0
	   end do
	   if (centr<bound) then
		  U(1,I) = rhoL !rho
!		  U(1,I) = 1.0
		  U(2,I) = rhoL*uL  !rho*u
		  U(3,I) = rhoL*vL  !rho*v
		  U(4,I) = rhoL*wL  !rho*w
		  U(5,I) = pL/(gmm-1.0)+0.5*rhoL*(uL*uL+vL*vL+wL*wL)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = BxL !Bx
	      B(2,I) = ByL  !By
	      B(3,I) = BzL  !Bz
		else
		  U(1,I) = rhoR !rho
!		  U(1,I) = 1.0
		  U(2,I) = rhoR*uR  !rho*u
		  U(3,I) = rhoR*vR  !rho*v
		  U(4,I) = rhoR*wR  !rho*w
		  U(5,I) = pR/(gmm-1.0)+0.5*rhoR*(uR*uR+vR*vR+wR*wR)!+(Bx*Bx+ByR*ByR+BzR*BzR)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = BxR !Bx
	      B(2,I) = ByR  !By
	      B(3,I) = BzR  !Bz
	   end if
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	   	if (nodes(nx,I)<bound) then
	      V_nodes(1,I) = uL
		  V_nodes(2,I) = vL
		  B_nodes(1,I) = BxL !Bx
	      B_nodes(2,I) = ByL  !By
		else
		  V_nodes(1,I) = uR
		  V_nodes(2,I) = vR
	      B_nodes(1,I) = BxR !Bx
	      B_nodes(2,I) = ByR  !By
	   end if

	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	   centr=0.0
	   do J=1,2
	     centr=centr+nodes(nx,edges(J,I))*0.5
	   end do
	   if (centr<bound) then
	      B_edges(1,I) = BxL*edge_n(1,I) + ByL*edge_n(2,I) !Bn
          B_edges(2,I) = BzL
		  V_edges(1,I) = uL*edge_n(1,I) + vL*edge_n(2,I) !Bn
          V_edges(2,I) = wL
		else
	      B_edges(1,I) = BxR*edge_n(1,I) + ByR*edge_n(2,I) !Bn
		  B_edges(2,I) = BzR
		  V_edges(1,I) = uR*edge_n(1,I) + vR*edge_n(2,I) !Bn
          V_edges(2,I) = wR
	   end if
	 end do


	do I=1,nElems
	  if (sum((elems_center(:,I)-disc_center)**2)<0.15**2) then
		  U(1,I) = 10.0 !rho
		  U(2,I) = 10.0*uR  !rho*u
		  U(3,I) = 10.0*vR  !rho*v
		  U(4,I) = 10.0*wR  !rho*w
		  U(5,I) = 1.0/(gmm-1.0)+0.5*10.0*(uR*uR+vR*vR+wR*wR)
	  end if
	end do

  end subroutine cloud_shock_pb


  subroutine dissip
    real (kind=precis) :: rho, p, ux, v, w, Bx, By, Bz
	real (kind=precis), dimension(2):: centr
	real (kind=precis) :: nx, ny, v0, B0, ksi, rho0, p0, phi,x,y, r
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

	r = 6.0
	gmm = 5.0/3.0
	nx = 1.0/sqrt(1+r**2)
	ny = r/sqrt(1+r**2)
	v0 = 0.0
	B0 = 1.0
	ksi = 0.2
	rho0 = 1.0
	p0 = 1.0
	t_final = 100.0

!!!!начальные условия в ячейках
	 do I = 1,nElems
		  x = elems_center(1,i)
		  y = elems_center(2,i)
		  phi = 2.0*pi*(nx*x+ny*y)/ny
		  U(1,I) = rho0 !rho
		  U(2,I) = rho0*(v0*nx-ksi*ny*cos(phi))  !rho*u
		  U(3,I) = rho0*(v0*ny+ksi*nx*cos(phi))  !rho*v
		  U(4,I) = rho0*ksi*sin(phi)  !rho*w
		  U(5,I) = p0/(gmm-1.0)+0.5*(sum(U(2:4,I)**2))/rho0 !+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = B0*nx+ksi*ny*sqrt(4.0*pi*rho0)*cos(phi) !Bx
	      B(2,I) = B0*ny-ksi*nx*sqrt(4.0*pi*rho0)*cos(phi)  !By
	      B(3,I) = -ksi*sqrt(4.0*pi*rho0)*sin(phi)  !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
		  x = nodes(1,i)
		  y = nodes(2,i)
		  phi = 2.0*pi*(nx*x+ny*y)/ny
	      V_nodes(1,I) = v0*nx-ksi*ny*cos(phi)
		  V_nodes(2,I) = v0*ny+ksi*nx*cos(phi)
		  B_nodes(1,I) = B0*nx+ksi*ny*sqrt(4.0*pi*rho0)*cos(phi) !Bx
	      B_nodes(2,I) = B0*ny-ksi*nx*sqrt(4.0*pi*rho0)*cos(phi) !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	   centr=0.0
	   do J=1,2
	     centr=centr+nodes(:,edges(J,I))*0.5
	   end do
		  x = centr(1)
		  y = centr(2)
		  phi = 2.0*pi*(nx*x+ny*y)/ny

		  ux = v0*nx-ksi*ny*cos(phi)
		  v = v0*ny+ksi*nx*cos(phi)
		  w = ksi*sin(phi)
		  Bx = B0*nx+ksi*ny*sqrt(4.0*pi*rho0)*cos(phi)
		  By = B0*ny-ksi*nx*sqrt(4.0*pi*rho0)*cos(phi)
		  Bz = -ksi*sqrt(4.0*pi*rho0)*sin(phi)
	      B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
          B_edges(2,I) = Bz
		  V_edges(1,I) = ux*edge_n(1,I) + v*edge_n(2,I) !Bn
          V_edges(2,I) = w

	 end do

  end subroutine dissip


  subroutine alfven
    real (kind=precis) :: rho, p, ux, v, w, Bx, By, Bz, alpha, Bpar, Bper, upar, uper
	real (kind=precis), dimension(2):: centr
	real (kind=precis) :: nx, ny, v0, B0, ksi, rho0, p0, phi,x,y, r
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

	t_final = 5.0
	alpha = pi/6.0
	r = 6.0
	gmm = 5.0/3.0

	rho0 = 1.0
	upar = 0.0
	Bpar = sqrt(4.0*pi) * 1.0
	p0 = 0.1

!!!!начальные условия в ячейках
	 do I = 1,nElems
		  x = elems_center(1,i)
		  y = elems_center(2,i)
		  ksi = x*cos(alpha)+y*sin(alpha)
		  uper = 0.1 * sin(2.0*pi*ksi)
		  Bper = sqrt(4.0*pi) * 0.1 * sin(2.0*pi*ksi)
		  U(1,I) = rho0 !rho
		  U(2,I) = rho0*(upar*cos(alpha)-uper*sin(alpha))  !rho*u
		  U(3,I) = rho0*(uper*cos(alpha)+upar*sin(alpha))  !rho*v
		  U(4,I) = rho0*0.1*cos(2.0*pi*ksi)  !rho*w
		  U(5,I) = p0/(gmm-1.0)+0.5*(sum(U(2:4,I)**2))/rho0 !+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = Bpar*cos(alpha)-Bper*sin(alpha) !Bx
	      B(2,I) = Bper*cos(alpha)+Bpar*sin(alpha)  !By
	      B(3,I) = sqrt(4.0*pi)*0.1*cos(2.0*pi*ksi)  !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
		  x = nodes(1,i)
		  y = nodes(2,i)
		  ksi = x*cos(alpha)+y*sin(alpha)
		  uper = 0.1 * sin(2.0*pi*ksi)
		  Bper = sqrt(4.0*pi) * 0.1 * sin(2.0*pi*ksi)

	      V_nodes(1,I) = upar*cos(alpha)-uper*sin(alpha)
		  V_nodes(2,I) = uper*cos(alpha)+upar*sin(alpha)
		  B_nodes(1,I) = Bpar*cos(alpha)-Bper*sin(alpha) !Bx
	      B_nodes(2,I) = Bper*cos(alpha)+Bpar*sin(alpha) !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	   centr=0.0
	   do J=1,2
	     centr=centr+nodes(:,edges(J,I))*0.5
	   end do
		  x = centr(1)
		  y = centr(2)
		  ksi = x*cos(alpha)+y*sin(alpha)
		  uper = 0.1 * sin(2.0*pi*ksi)
		  Bper = sqrt(4.0*pi) * 0.1 * sin(2.0*pi*ksi)

		  ux = upar*cos(alpha)-uper*sin(alpha)
		  v = uper*cos(alpha)+upar*sin(alpha)
		  w = 0.1*cos(2.0*pi*ksi)
		  Bx = Bpar*cos(alpha)-Bper*sin(alpha)
		  By = Bper*cos(alpha)+Bpar*sin(alpha)
		  Bz = sqrt(4.0*pi)*0.1*cos(2.0*pi*ksi)
	      B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
          B_edges(2,I) = Bz
		  V_edges(1,I) = ux*edge_n(1,I) + v*edge_n(2,I) !Bn
          V_edges(2,I) = w

	 end do

  end subroutine alfven

  subroutine axial

	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)
	 real (kind=precis) :: pres, lef, rig
	 integer:: I, J
	 real (kind=precis) :: cent
	 real (kind=precis) :: r_c1 = 0.6
	 real (kind=precis), dimension(2):: buuff

!!!!начальные условия в ячейках
	gmm = 1.4
	t_final = 50.0
	pres = 0.2
	 do I = 1,nElems
		  U(1,I) = 1.0 !rho
		  U(2,I) = 0.0  !rho*u
		  U(3,I) = 0.0  !rho*v
		  U(4,I) = 0.0  !rho*w
		  U(5,I) = pres/(gmm-1.0) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
!	      B(1,I) = 1.0  !Bx
!		B(1,I) = (0.9-elems_center(2,I))/(0.9-r_c1)  !Bx
B(1,I) = 0.0
		if (elems_center(2,I)<r_c1) B(1,I) = 1.0
!		if (elems_center(2,I)>0.9) B(1,I) = 0.0
	      B(2,I) = 0.0  !By
	      B(3,I) = 0.0  !Bz

		if (elems_center(2,I)>1.0) then
		  U(1,I) = 1.0 !rho
		  buuff = inflow(elems_center(1,I),1.0_precis)
		  U(2:3,I) = buuff  !rho*u
		  U(5,I) = pres/(gmm-1.0)+0.5*sum(buuff**2) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
		end if
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	      V_nodes(1,I) = 0.0
		  V_nodes(2,I) = 0.0
!		  B_nodes(1,I) = 1.0  !Bx
!		  B_nodes(1,I) = (0.9-nodes(2,I))/(0.9-r_c1)
B_nodes(1,I) = 0.0
if (nodes(2,I)<r_c1) B_nodes(1,I) = 1.0
!if (nodes(2,I)>0.9) B_nodes(1,I) = 0.0
	      B_nodes(2,I) = 0.0  !Br

		if (nodes(2,I)>1.0) then
		  buuff = inflow(nodes(1,I),1.0_precis)
	      V_nodes(:,I) = buuff
		end if

	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
!	      B_edges(1,I) = 1.0*edge_n(1,I)  !Bn
cent = 0.5*(nodes(2,edges(2,I))+nodes(2,edges(1,I)))
!B_edges(1,I) = ((0.9-cent)/(0.9-r_c1))*edge_n(1,I)  !Bn
!B_edges(1,I) = 0.0

lef = B_nodes(1, edges(1,I))
rig = B_nodes(1, edges(2,I))
B_edges(1,I) = 0.5*(lef+rig)*edge_n(1,I)

!if (cent<r_c1) B_edges(1,I) = 1.0*edge_n(1,I)  !Bn
!if (cent>0.9) B_edges(1,I) = 0.0  !Bn
          B_edges(2,I) = 0.0
		  V_edges(1,I) = 0.0 !Bn
          V_edges(2,I) = 0.0

		if (cent>1.0) then
		cent = 0.5*(nodes(1,edges(2,I))+nodes(1,edges(1,I)))
		  buuff = inflow(cent,1.0_precis)
	      V_edges(:,I) = dot_product(buuff,edge_n(1:2,I))
		end if

	 end do
  end subroutine axial

  subroutine axial1

	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)
	 real (kind=precis) :: pres, rad
	 integer:: I, J
!!!!начальные условия в ячейках
	t_final = 100.0
	pres = 0.2
	gmm = 5.0/3.0
	 do I = 1,nElems
		  rad = elems_center(2,I)
		  pres = rad**2+1.0
		  U(1,I) = rad**2+1.0 !rho
		  U(2,I) = 0.0  !rho*u
		  U(3,I) = 0.0  !rho*v
		  U(4,I) = 0.0  !rho*w
		  U(5,I) = pres/(gmm-1.0) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = rad**2  !Bx
	      B(2,I) = rad  !By
	      B(3,I) = rad  !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	      rad = nodes(2,I)
		  V_nodes(1,I) = 0.0
		  V_nodes(2,I) = 0.0
		  B_nodes(1,I) = rad**2  !Bx
	      B_nodes(2,I) = rad  !Br
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
		  rad = 0.5*(nodes(2,edges(1,I))+nodes(2,edges(2,I)))
	      B_edges(1,I) = (rad**2)*edge_n(1,I)+rad*edge_n(2,I)  !Bn
          B_edges(2,I) = rad
		  V_edges(1,I) = 0.0 !Bn
          V_edges(2,I) = 0.0


	 end do
  end subroutine axial1


  subroutine axial95

	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)
	 real (kind=precis) :: pres, lef, rig
	 integer:: I, J
	 real (kind=precis) :: cent
	 real (kind=precis) :: r_c1 = 0.6
	 real (kind=precis), dimension(2):: buuff

!!!!начальные условия в ячейках
	b_scale=sqrt(4.0*pi)
	gmm = 5.0/3.0
	t_final = 50.0
	pres = 0.01
	do I = 1,nElems
		if (elems_center(2,I)<r_c1) then
!				pres = 0.0025
				B(1,I) = 1.0*b_scale
!				U(1,I) = 0.5 !rho
			else
!				pres = 0.01
				B(1,I) = 0.0
!				U(1,I) = 1.0 !rho
		end if
		pres = pres95(elems_center(:,I))
		U(1,I) = 1.0 !rho
		U(2,I) = 0.0  !rho*u
		U(3,I) = 0.0  !rho*v
		U(4,I) = 0.0  !rho*w
		U(5,I) = pres/(gmm-1.0) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	    B(2,I) = 0.0  !By
	    B(3,I) = 0.0  !Bz

		if (elems_center(2,I)>1.5) then
		  U(1,I) = 1.0 !rho
!		  pres = 0.01
		  buuff = inflow(elems_center(1,I),1.5_precis)
		  U(2:3,I) = buuff  !rho*u
		  U(5,I) = pres/(gmm-1.0)+0.5*sum(buuff**2) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
		end if

	end do
!!!! начальные данные в точках

	 do I=1,nNodes
	    V_nodes(1,I) = 0.0
		V_nodes(2,I) = 0.0
		B_nodes(1,I) = 0.0
		if (nodes(2,I)<r_c1) B_nodes(1,I) = 1.0*b_scale
	    B_nodes(2,I) = 0.0  !Br

		if (nodes(2,I)>1.5) then
		  buuff = inflow(nodes(1,I),1.0_precis)
	      V_nodes(:,I) = buuff
		end if
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
		cent = 0.5*(nodes(2,edges(2,I))+nodes(2,edges(1,I)))
		lef = B_nodes(1, edges(1,I))
		rig = B_nodes(1, edges(2,I))
		B_edges(1,I) = 0.5*(lef+rig)*edge_n(1,I)
        B_edges(2,I) = 0.0
		V_edges(1,I) = 0.0 !Bn
        V_edges(2,I) = 0.0

		if (cent>1.5) then
		cent = 0.5*(nodes(1,edges(2,I))+nodes(1,edges(1,I)))
		  buuff = inflow(cent,1.5_precis)
	      V_edges(:,I) = dot_product(buuff,edge_n(1:2,I))
		end if
	 end do
  end subroutine axial95

  function inflow(zet,rr1) result (inf)
    real (kind=precis) :: zet, rad,rr,rr1
	real (kind=precis), dimension(2) :: inf

	rr = 1.5
	rad = sqrt(rr**2+zet**2)
	inf(1) = -2.2*zet/(rad**3)
	inf(2) = -2.2*rr/(rad**3)

!	inf = 2.0*inf
  end function inflow

  subroutine test_visc

	real (kind=precis) :: rho, p, vx, vy, vz, v0, Bx, By, Bz, B0
	real (kind=precis), dimension(1:2) :: centr = (/0.0,1.5/), center
    real (kind=precis) :: r0, r1,z1, r,dist, dist1
	! rho = 1, p = 1/gmm, gmm = 7.0/3.0
	! U(:,nElem) = (rho, rho*u, rho*v, E)

	 integer:: I, J

!	B0 = 4.0*pi
B0 = 10.0
	gmm = 5./3.
	p = 1.0
!	p = 0.01/(4.0*pi)
	rho = 1.0
    t_final = 2.0
!	Bx = 5.0
	Bx = 0.0
	By = 0.0
	Bz = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0

!!!!начальные условия в ячейках
	 do I = 1,nElems
		  U(1,I) = rho !rho
		  U(2,I) = rho*vx  !rho*u
		  U(3,I) = rho*vy  !rho*v
		  U(4,I) = rho*vz  !rho*w
		  U(5,I) = p/(gmm-1.0)+0.5*rho*(vx*vx+vy*vy+vz*vz)!+(Bx*Bx+ByL*ByL+BzL*BzL)/(8.0*pi) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
	      B(1,I) = Bx !Bx
	      B(2,I) = By !By
	      B(3,I) = Bz !Bz
	 end do
!!!! начальные данные в точках

	 do I=1,nNodes
	      V_nodes(1,I) = vx
		  V_nodes(2,I) = vy
		  B_nodes(1,I) = Bx !Bx
	      B_nodes(2,I) = By  !By
	 end do

!!!! начальные данные на гранях
	 do I=1,nEdges
	      B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
          B_edges(2,I) = Bz
		  V_edges(1,I) = vx*edge_n(1,I) + vy*edge_n(2,I) !Bn
          V_edges(2,I) = vz
	 end do

	r0 = 0.1
	 do I = 1,nElems
		z1 = elems_center(1,I)
		r1 = elems_center(2,I)
		dist = max(abs(z1-centr(1)) , abs(r1-centr(2)))
		z1 = elems_center(1,I) - centr(1)
		r1 = elems_center(2,I) - centr(2)
		dist1 = z1**2+r1**2
		if (dist1<r0) then
!	      B(1,I) = cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !Bx
!	      B(2,I) = -cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !By
	      B(1,I) = B0*r1*exp(-dist1/(r0**2))/sqrt(dist1) !Bx
	      B(2,I) = -B0*z1*exp(-dist1/(r0**2))/sqrt(dist1) !By
		end if
	 end do

!!!! начальные данные в точках
	 do I=1,nNodes
		z1 = nodes(1,I)
		r1 = nodes(2,I)
		dist = max(abs(z1-centr(1)) , abs(r1-centr(2)))
		z1 = nodes(1,I) - centr(1)
		r1 = nodes(2,I) - centr(2)
		dist1 = z1**2+r1**2
		if (dist1<r0) then
!	      B_nodes(1,I) = cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !Bx
!	      B_nodes(2,I) = -cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !By
	      B_nodes(1,I) = B0*r1*exp(-dist1/(r0**2))/sqrt(dist1) !Bx
	      B_nodes(2,I) = -B0*z1*exp(-dist1/(r0**2))/sqrt(dist1) !By
		end if
	 end do
!!!! начальные данные на гранях
	 do I=1,nEdges
	   center = 0.0
	   do J=1,2
	     center=center+nodes(:,edges(J,I))*0.5
	   end do
		z1 = center(1)
		r1 = center(2)
		dist = max(abs(z1-centr(1)) , abs(r1-centr(2)))
		z1 = center(1) - centr(1)
		r1 = center(2) - centr(2)
		dist1 = z1**2+r1**2
		if (dist1<r0) then
!	      Bx = cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !Bx
!	      By = -cos(pi*(z1-centr(1))/0.4)*cos(pi*(r1-centr(2))/0.4) !By
	      Bx = B0*r1*exp(-dist1/(r0**2))/sqrt(dist1) !Bx
	      By = -B0*z1*exp(-dist1/(r0**2))/sqrt(dist1) !By
		  B_edges(1,I) = Bx*edge_n(1,I) + By*edge_n(2,I) !Bn
		end if
	 end do
  end subroutine test_visc


  subroutine radjet2D

    real (kind=precis) :: pres, lef, rig
    integer:: I, J
    real (kind=precis) :: cent
    real (kind=precis), dimension(2):: buuff
    real (kind=precis) :: z_dim, r_dim
    logical :: cp_exists

!     z_dim = 2.5
!     r_dim = 2.5
!
!     do i=1,nNodes
!       nodes(1,i) = nodes(1,i)*z_dim/z_max
!       nodes(2,i) = nodes(2,i)*z_dim/z_max
!     end do

  t_final = 50.0
  tau_out = 0.05
!!!!начальные условия в ячейках
  b_scale = sqrt(4.0*pi)
  omega0 = 1.
  pres = 0.05
  p_0 = 0.01
  init_dens = 1.0
  r_disk = 0.6
  gr_rad = 0.1
  G = 0.5
  gmm = 5.0/3.0
  refp = (/1.5,1.5/)
  refsspeed = sqrt(gmm*pres/init_dens)
  refMach = 1.5
  lsp = 3.35E5
!*********************
!     absorb0 = 2.28E-8
!     int0 = 1.!1.
  cp_exists = read_control_point()
  if (.not. cp_exists) then
    do I = 1,nElems
      pres = pres95(elems_center(:,I))
      if (elems_center(2,I)<r_disk) then
	B(1,I) = b_scale*magnz(elems_center(2,I))
      else
	B(1,I) = 0.0
      end if
      U(1,I) = init_dens !rho
      U(2,I) = 0.0  !rho*u
      U(3,I) = 0.0  !rho*v
      U(4,I) = 0.0  !rho*w
      U(5,I) = pres/(gmm-1.0) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
      B(2,I) = 0.0  !By
      B(3,I) = 0.0  !Bz

      if (elems_center(2,I)>refp(2)) then
	U(1,I) = init_dens !rho

	buuff = inflow_radjet2D(elems_center(1,I),refp(2),refp,refMach,refsspeed)
	U(2:3,I) = buuff  !rho*u
	U(5,I) = pres/(gmm-1.0)+0.5*sum(buuff**2) !1.7857 ! ! E = rho*(0.5*(u**2+v**2)+p/(rho*(gmm-1)))
      end if
    end do

!!!! начальные данные в точках
    do I=1,nNodes
      V_nodes(1,I) = 0.0
      V_nodes(2,I) = 0.0
      B_nodes(1,I) = 0.0
      if (nodes(2,I)<r_disk) B_nodes(1,I) = b_scale*magnz(nodes(2,I))
	B_nodes(2,I) = 0.0  !Br
	if (nodes(2,I)>refp(2)) then
	  buuff = inflow_radjet2D(nodes(1,I),refp(2),refp,refMach,refsspeed)
	  V_nodes(:,I) = buuff
	end if
    end do

!!!! начальные данные на гранях
    do I=1,nEdges
      cent = 0.5*(nodes(2,edges(2,I))+nodes(2,edges(1,I)))
      lef = B_nodes(1, edges(1,I))
      rig = B_nodes(1, edges(2,I))
      B_edges(1,I) = 0.5*(lef+rig)*edge_n(1,I)
      B_edges(2,I) = 0.0
      V_edges(1,I) = 0.0 !Bn
      V_edges(2,I) = 0.0

      if (cent>refp(2)) then
	cent = 0.5*(nodes(1,edges(2,I))+nodes(1,edges(1,I)))
	buuff = inflow_radjet2D(cent,refp(2),refp,refMach,refsspeed)
	V_edges(:,I) = dot_product(buuff,edge_n(1:2,I))
      end if
    end do

  end if

  do i=1,nElems
    U_hat(:,i) = U(:,i)
    U_til(:,i) = U(:,i)
  end do
!!!! начальные данные для излучения
!    call radiation_solverTR

!     do i=1,nElems
!       intensity(i) = intensity_hat(i)
!       Poynting(:,i) = Poynting_hat(:,i)
!     end do

end subroutine radjet2D

  function inflow_radjet2D(x1,x2,refpoint,Mach0,sspeed0) result (inf)
    real (kind=precis) :: x1, rad,x2,Mach0,sspeed0
    real (kind=precis), dimension(2) :: inf, refpoint
    real (kind=precis) :: vel0

    vel0 = Mach0*sspeed0*sum(refpoint**2)
!print *, Mach0, sspeed0, vel0, sum(refpoint**2)
!read *
    rad = sqrt(x1**2+x2**2)
    inf(1) = -vel0*x1/(rad**3)
    inf(2) = -vel0*x2/(rad**3)
  end function inflow_radjet2D

  function magnz(radius) result (Hz0)
    real (kind=precis) :: radius, Hz0

    Hz0 = 1.1*(Pi/2 - ATan((radius - 0.55)/0.04))/Pi - 0.06
!    Hz0 = 1D0

!     if(radius<r_disk) then
!       Hz0 = 1D0
!     else
!       Hz0 = (r_disk-radius)*5D0
!       Hz0 = 0D0
!     end if

  end function

end module initials
