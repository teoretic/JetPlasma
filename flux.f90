! ===========================================================
! РЕШЕНИЕ ОДНОМЕРНОЙ ЗАДАЧИ О РАСПАДЕ РАЗРЫВА и вычисление потоков
!
! - расчет физических и численных потоков для решения 1D задачи
!   на границах ячеек
! - расчет собственных чисел и векторов
! - расчет правой части -- для (r,z) геометрии
! ===========================================================

! ==========================================================
! вектор консервативных переменных: U = (rho, rho*u, rho*v, E)
! ==========================================================

! ==========================================================
! Используется консервативая форма записи уравнений вида
! dU/dt + dF/dr + dG/dz + alpha/r*H = 0
! Toro, с. 28, раздел 1.6.3, (1.102)-(1.104)
! ==========================================================


module flux

  use names

  implicit none

  integer, parameter:: alpha = 0 ! 0 -- (x,y) геометрия, 1 -- (r,z) геометрия

contains

  ! ======================================================
  ! вычисление физического потока, собственных чисел и векторов
  ! входные переменные: UU
  ! выходные переменные: F, D, R, L
  ! UU -- вектор консервативных переменных
  ! F -- вектор физического потока
  ! D -- диагональная матрица собственных чисел
  ! R -- матрица правых собственных векторов
  ! L -- матрица левых собственных веторов, L = R^{-1}
  ! Toro, с. 103, раздел 3.2.1

pure subroutine FL(UU,F,D,R,L)
  implicit none

  ! входные переменные
  real (kind=precis), dimension(1:5), intent(in):: UU

  ! выходные переменные
  real (kind=precis), dimension(1:5), intent(out):: F
  real (kind=precis), dimension(1:5,1:5), intent(out):: D, R, L

  real (kind=precis):: rho, u, v, w, E, H2, a, p, k_e, mult
  real (kind=precis):: courant, verif

  ! вычисляем все простые переменные

  rho = UU(1)

  u = UU(2)/rho
  v = UU(3)/rho
  w = UU(4)/rho
  ! полная энергия (внутренняя плюс кинетическая)
  E = UU(5)
  k_e = u*u+v*v+w*w
  p = (E-0.5*rho*k_e)*(gmm-1.0)         ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
  H2 = (E+p)/rho								! энтальпия

  verif = p*gmm/rho
  if (verif.LT.0.0) then
    verif = 0.0
  end if

  ! ошибка вычисления SQRT из-за деления на нуль в предыдущих формулах  !!!!
  a = SQRT(verif)	! скорость звука


  F = (/ UU(2), UU(2)*u+p, UU(2)*v, UU(2)*w, u*(E+p) /)

  D = 0.0
  D(1,1) = u-a
  D(2,2) = u
  D(3,3) = u
  D(4,4) = u
  D(5,5) = u+a


  courant = D(5,5)*tau/0.2

  R(1,1) =1D0;    R(1,2) =0D0; R(1,3) =0D0; R(1,4) =1D0;     R(1,5) =1D0;
  R(2,1) =u-a;    R(2,2) =0D0; R(2,3) =0D0; R(2,4) =u;       R(2,5) =u+a;
  R(3,1) =v;      R(3,2) =1D0; R(3,3) =0D0; R(3,4) =v;       R(3,5) =v;
  R(4,1) =w;      R(4,2) =0D0; R(4,3) =1D0; R(4,4) =w;       R(4,5) =w;
  R(5,1) =H2-a*u; R(5,2) =v;   R(5,3) =w;   R(5,4) =0.5*k_e; R(5,5) =H2+a*u;

  mult = (gmm-1.0)/(2.0*a*a)

  L(1,1) =H2*mult+1.0/(a*2.0)*(u-a);  L(1,2) =-(u*mult+1.0/(a*2.0));   L(1,3) =-v*mult;    L(1,4) = -w*mult;       L(1,5) =1D0*mult;
  L(2,1) =-v;                         L(2,2) =0D0;                L(2,3) =1.0;        L(2,4) =0D0;            L(2,5) =0D0;
  L(3,1) =-w;                         L(3,2) =0D0;                L(3,3) =0D0;        L(3,4) = 1.0;      L(3,5) =0D0;
  L(4,1) =-2.0*H2*mult+2.0;           L(4,2) = 2.0*u*mult;        L(4,3) =2.0*v*mult; L(4,4) =2.0*w*mult;     L(4,5) =-2D0*mult;
  L(5,1) =H2*mult-1.0/(a*2.0)*(u+a);  L(5,2) =-u*mult+1.0/(a*2.0);L(5,3) =-v*mult;    L(5,4) =-w*mult;        L(5,5) =1D0*mult;

end subroutine FL

pure subroutine FL_hllc(UU,F,u,v,w,a,p,rho,E)
  implicit none

  ! входные переменные
  real (kind=precis), dimension(1:5), intent(in):: UU

  ! выходные переменные
  real (kind=precis), dimension(1:5), intent(out):: F
  real (kind=precis), intent(out):: u, v, w, a, p, rho, E

  real (kind=precis):: H2, k_e, mult
  real (kind=precis):: courant

  ! вычисляем все простые переменные

  rho = UU(1)

  u = UU(2)/rho
  v = UU(3)/rho
  w = UU(4)/rho
  E = UU(5)
  k_e = u*u+v*v+w*w
  p = (E-0.5*rho*k_e)*(gmm-1.0)    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
  H2 = (E+p)/rho  ! энтальпия

  if (p*gmm/rho<0.0) then
    p = 0.0
  end if

  a = SQRT(p*gmm/rho)						! скорость звука

  F(1)=UU(2)
  F(2)=UU(2)*u+p
  F(3)=UU(2)*v
  F(4)=UU(2)*w
  F(5)=u*(E+p)

end subroutine FL_hllc

  ! ===============================================
  ! Расчет правой части -- только для (r,z) геометрии
  ! Физические источники учитываются отдельно
  ! U -- вектор консервативных переменных
  ! H -- правая часть для учета (r,z) геометрии
  ! если alpha = 0, то H должен возвращаться нулевым


  ! ===============================================
  ! Расчет правой части -- для физики
  ! U -- вектор консервативных переменных
  ! H -- правая часть для учета (r,z) геометрии

  function H_PHYS(U) result (H)
    implicit none

    real (kind=precis), dimension(1:5):: U
    real (kind=precis), dimension(1:5):: H

    H = 0.0

  end function H_PHYS

! ======================================================
pure function flux_num_hllc(Ul,Ur) result(fr)
  implicit none

  real (kind=precis), dimension(1:5):: FR ! численный поток
  real (kind=precis), dimension(1:5), intent(in):: Ul, Ur
  real (kind=precis), dimension(1:5)::U_lr, Ur_et, Ul_et
  real (kind=precis), dimension(1:5):: F_l, F_r, F_lr
  real (kind=precis), dimension(1:5,1:5)::  R_lr, L_lr, D_lr

  real (kind=precis):: S_l, S_r, S_et, q_l, q_r

  real (kind=precis):: u_l,a_l,p_l,rho_l,E_l,v_l,w_l
  real (kind=precis):: u_r,a_r,p_r,rho_r,E_r,v_r,w_r
  real (kind=precis):: a_lr,p_lr,rho_lr,E_lr,v_lr,ulr, w_lr
  real (kind=precis):: u_et, p_et
  real (kind=precis):: numb1, numb2

  call FL_hllc(Ul,F_l,u_l,v_l,w_l,a_l,p_l,rho_l,E_l)
  call FL_hllc(Ur,F_r,u_r,v_r,w_r,a_r,p_r,rho_r,E_r)

  rho_lr = 0.5*(rho_l+rho_r)
  a_lr = 0.5*(a_l+a_r)
  p_et = 0.5*(p_l+p_r)+0.5*(u_l-u_r)*(rho_lr*a_lr)
  u_et = 0.5*(u_l+u_r)+0.5*(p_l-p_r)/(rho_lr*a_lr)

  if (p_et .LE. p_l) then
    q_l = 1.0
  else
    q_l = sqrt(1.0+(gmm+1.0)/(2.0*gmm)*(p_et/p_l-1.0))
  end if

  if (p_et .LE. p_r) then
    q_r = 1.0
  else
    q_r = sqrt(1.0+(gmm+1.0)/(2.0*gmm)*(p_et/p_r-1.0))
  end if

  S_l = u_l-a_l*q_l
  S_r = u_r+a_r*q_r
  S_et = u_et
  !S_et = (p_r-p_l+rho_l*u_l*(S_l-u_l)-rho_r*u_r*(S_r-u_r))/(rho_l*(S_l-u_l)-rho_r*(S_r-u_r))

  numb1 = rho_r*(S_r-u_r)/(S_r-S_et)
  numb2 = rho_l*(S_l-u_l)/(S_l-S_et)

  Ur_et(1) =numb1*1D0
  Ur_et(2) =numb1*S_et
  Ur_et(3) =numb1*v_r
  Ur_et(4) =numb1*w_r
  Ur_et(5) =numb1*(E_r/rho_r+(S_et-u_r)*(S_et+p_r/(rho_r*(S_r-u_r))))
  Ul_et(1) =numb2*1D0
  Ul_et(2) =numb2*S_et
  Ul_et(3) =numb2*v_l
  Ul_et(4) =numb2*w_l
  Ul_et(5) =numb2*(E_l/rho_l+(S_et-u_l)*(S_et+p_l/(rho_l*(S_l-u_l))))

  if (S_l .ge. 0D0) then
    FR = F_l
  end if

  if (S_l .le. 0D0 .and. S_et .ge. 0D0) then
    FR = F_l + S_l*(Ul_et-Ul)
  end if

  if (S_et .le. 0D0 .and. S_r .ge. 0D0) then
    FR = F_r + S_r*(Ur_et-Ur)
  end if

  if (S_r .le. 0D0) then
    FR= F_r
  end if

end function flux_num_hllc

  subroutine compute_fluxes(fmo,feo,fmoj,feoj,fmi,fei)
	real (kind=precis), intent(inout):: fmo,feo,fmoj,feoj,fmi,fei
	integer :: i
	real (kind=precis), dimension(1:2) :: norr1
	real (kind=precis) :: l_b, dist

	fmo = 0.0
	feo = 0.0
	fmoj = 0.0
	feoj = 0.0
	fmi = 0.0
	fei = 0.0

	do nElem=1,nElems_b
		j=bc_elems(nElem)
		norr1 = (/1.0,0.0/)
		if (dot_product(n(:,j,bc_elem_edge(nElem)),norr1)>0.1) then
			l_b = edge_vol(bc_edges(nElem))
			dist = sqrt(sum(elems_center(:,j)**2))
			if (U(2,j)>0.0) then
				fmo = fmo + 2.0*pi*U(2,j)*elems_center(2,j)*l_b
				feo = feo + 2.0*pi*U(2,j)*elems_center(2,j)*l_b*(U(5,j)-U(1,j)*G/dist)/U(1,j)
				if (B(1,j)>0.1*b_scale) then
					fmoj = fmoj + 2.0*pi*U(2,j)*elems_center(2,j)*l_b
					feoj = feoj + 2.0*pi*U(2,j)*elems_center(2,j)*l_b*(U(5,j)-U(1,j)*G/dist)/U(1,j)
				end if
			end if
		end if

		norr1 = (/0.0,1.0/)
		if (dot_product(n(:,j,bc_elem_edge(nElem)),norr1)>0.1) then
			l_b = edge_vol(bc_edges(nElem))
			dist = sqrt(sum(elems_center(:,j)**2))
			if (U(3,j)<0.0) then
				fmi = fmi - 2.0*pi*U(3,j)*elems_center(2,j)*l_b
				fei = fei - 2.0*pi*U(3,j)*elems_center(2,j)*l_b*(U(5,j)-U(1,j)*G/dist)/U(1,j)
			end if
		end if
	end do

  end subroutine compute_fluxes

  subroutine compute_fluxes1(jet_rad,fmo,feo,fmoj,feoj,fmi,fei)
	real (kind=precis), intent(inout):: fmo,feo,fmoj,feoj,fmi,fei,jet_rad
	integer :: i,j, dots, i1, i2, imin,n_jet_full
	real (kind=precis), dimension(1:2) :: norr1
	real (kind=precis) :: l_b, dist, buf,minim, jet_bound, delta_r, dens, veloc, rad, energ
	real (kind=precis), dimension(1:20000) :: rads, density, velocity, magn_field, energy, dists

	fmo = 0.0
	feo = 0.0
	fmoj = 0.0
	feoj = 0.0
	fmi = 0.0
	fei = 0.0

	i=0
	do nElem=1,nElems_b
		j=bc_elems(nElem)
		if (bc_type(nElem)==2) then
			i = i+1
			rads(i) = elems_center(2, j)
			density(i) = U(1,j)
			velocity(i) = U(2,j)/U(1,j)
			magn_field(i) = B(1,j)
			energy(i) = U(5,j)
			dists(i) = sqrt(sum(elems_center(:,j)**2))
		end if
	end do
	dots = i
	do i=1,dots-1
		imin=i
		minim = rads(i)
		do j=i+1, dots
			if (minim>rads(j)) then
				imin = j
				minim = rads(j)
			end if
		end do
		if (imin>i) then
			buf = rads(i)
			rads(i) = rads(imin)
			rads(imin) = buf

			buf = density(i)
			density(i) = density(imin)
			density(imin) = buf

			buf = velocity(i)
			velocity(i) = velocity(imin)
			velocity(imin) = buf

			buf = magn_field(i)
			magn_field(i) = magn_field(imin)
			magn_field(imin) = buf

			buf = energy(i)
			energy(i) = energy(imin)
			energy(imin) = buf

			buf = dists(i)
			dists(i) = dists(imin)
			dists(imin) = buf
		end if
	end do

	i=strt_jet_num
	jet_bound = 0.15*b_scale
	do while ((magn_field(i)-jet_bound)*(magn_field(i+1)-jet_bound)>1E-7)
		i = i+1
	end do
	n_jet_full = i
	delta_r = rads(n_jet_full+1)-rads(n_jet_full)
	if (abs(magn_field(n_jet_full+1)-magn_field(n_jet_full))>1E-7) then
			jet_rad = rads(n_jet_full) + delta_r*(jet_bound-magn_field(n_jet_full))/ (magn_field(n_jet_full+1)-magn_field(n_jet_full))
		else
			jet_rad = rads(n_jet_full)
	end if

	fmoj = 0.0
	feoj = 0.0
	print *, n_jet_full, jet_rad
	do i=1,n_jet_full-1
		dens = 0.5*(density(i)+density(i+1))
		veloc = 0.5*(velocity(i)+velocity(i+1))
		rad = 0.5*(rads(i)+rads(i+1))
		energ = 0.5*(energy(i)+energy(i+1))
		delta_r = rads(i+1)-rads(i)
!		print *, 'delta =', delta_r
		dist = 0.5*(dists(i)+dists(i+1))
		fmoj = fmoj + 2.0*pi*dens*veloc*rad*delta_r
!		print *, 'dfmoj =', 2.0*pi*dens*veloc*rad*delta_r
		feoj = feoj + 2.0*pi*dens*veloc*rad*delta_r*(energ-dens*G/dist)/dens
	end do

	i = n_jet_full
	delta_r = jet_rad-rads(i)
!	print *, '************'
!	print *, 'delta_r', delta_r
	dens = density(i)+0.5*delta_r*(density(i+1)-density(i))/(rads(n_jet_full+1)-rads(n_jet_full))
	veloc = velocity(i)+0.5*delta_r*(velocity(i+1)-velocity(i))/(rads(n_jet_full+1)-rads(n_jet_full))
	energ = energy(i)+0.5*delta_r*(energy(i+1)-energy(i))/(rads(n_jet_full+1)-rads(n_jet_full))
	dist = dists(i)+0.5*delta_r*(dists(i+1)-dists(i))/(rads(n_jet_full+1)-rads(n_jet_full))
	rad = 0.5*(rads(i)+jet_rad)
	fmoj = fmoj + 2.0*pi*dens*veloc*rad*delta_r
!	print *, 'dfmoj =', 2.0*pi*dens*veloc*rad*delta_r
	feoj = feoj + 2.0*pi*dens*veloc*rad*delta_r*(energ-dens*G/dist)/dens

	do nElem=1,nElems_b
		j=bc_elems(nElem)
		norr1 = (/1.0,0.0/)
		if (dot_product(n(:,j,bc_elem_edge(nElem)),norr1)>0.1) then
			l_b = edge_vol(bc_edges(nElem))
			dist = sqrt(sum(elems_center(:,j)**2))
			if (U(2,j)>0.0) then
				fmo = fmo + 2.0*pi*U(2,j)*elems_center(2,j)*l_b
				feo = feo + 2.0*pi*U(2,j)*elems_center(2,j)*l_b*(U(5,j)-U(1,j)*G/dist)/U(1,j)
			end if
		end if

		norr1 = (/0.0,1.0/)
		if (dot_product(n(:,j,bc_elem_edge(nElem)),norr1)>0.1) then
			l_b = edge_vol(bc_edges(nElem))
			dist = sqrt(sum(elems_center(:,j)**2))
			if (U(3,j)<0.0) then
				fmi = fmi - 2.0*pi*U(3,j)*elems_center(2,j)*l_b
				fei = fei - 2.0*pi*U(3,j)*elems_center(2,j)*l_b*(U(5,j)-U(1,j)*G/dist)/U(1,j)
			end if
		end if
	end do

  end subroutine compute_fluxes1

  subroutine compute_flux(jet_rad, fmoj, big_jet_rad, fmo, max_vel)
    real (kind=precis), intent(inout):: fmo, fmoj, jet_rad, big_jet_rad, max_vel
    integer :: i,j,imin,n_jet_full,n_injet_full
    real (kind=precis), dimension(1:2) :: norr1
    real (kind=precis) :: l_b, dist, jet_bound, delta_r, dens, veloc, rad, energ, mag_b
    real (kind=precis), dimension(1:nElems_b) :: rads, density, velocity, magn_field, energy, dists

    n_injet_full = 0
    n_jet_full = 0
    mag_b = 1E8
    jet_bound = 0.05*b_scale
    do nElem=1,outb_N
      j=outb_nums(nElem)
      rads(nElem) = elems_center(2,j)
      density(nElem) = U(1,j)
      velocity(nElem) = U(2,j)/U(1,j)
      magn_field(nElem) = B(1,j)
      energy(nElem) = U(5,j)
      dists(nElem) = sqrt(sum(elems_center(:,j)**2))
      if (B(3,j)<mag_b) then
        n_injet_full = nElem
        mag_b = B(3,j)
        jet_rad = rads(n_injet_full)
      end if
      if (nElem>strt_jet_num-1 .and. n_jet_full<1) then
        if ((magn_field(nElem-1)-jet_bound)*(magn_field(nElem)-jet_bound)>1E-7) then
          n_jet_full = nElem
          big_jet_rad = rads(n_jet_full)
        end if
      end if
    end do

    fmoj = 0.0
    do i=1,n_injet_full
      dens = 0.5*(density(i)+density(i+1))
      veloc = 0.5*(velocity(i)+velocity(i+1))
      rad = 0.5*(rads(i)+rads(i+1))
      energ = 0.5*(energy(i)+energy(i+1))
      delta_r = rads(i+1)-rads(i)
      dist = 0.5*(dists(i)+dists(i+1))
      fmoj = fmoj + 2.0*pi*dens*veloc*rad*delta_r
    end do

    fmo = 0.0
    do i=1,n_jet_full
      dens = 0.5*(density(i)+density(i+1))
      veloc = 0.5*(velocity(i)+velocity(i+1))
      rad = 0.5*(rads(i)+rads(i+1))
      energ = 0.5*(energy(i)+energy(i+1))
      delta_r = rads(i+1)-rads(i)
      dist = 0.5*(dists(i)+dists(i+1))
      fmo = fmo + 2.0*pi*dens*veloc*rad*delta_r
    end do

    veloc = 0.0
    do nElem=1,nElems
      l_b = sqrt(sum(U(2:4,nElem)**2))/U(1,nElem)
      if (veloc<l_b) veloc = l_b
    end do
    max_vel = veloc
  end subroutine compute_flux


  subroutine strt_jets

	integer :: i,j, dots, i1, i2, imin
	real (kind=precis) :: buf, minim, jet_rad, jet_bound, rad
	real (kind=precis), dimension(1:20000) :: rads, magn_field

	i=0
	do nElem=1,nElems_b
		j=bc_elems(nElem)
		if (bc_type(nElem)==2) then
			i = i+1
			rads(i) = elems_center(2, j)
			magn_field(i) = B(1,j)
		end if
	end do
	dots = i
	do i=1,dots-1
		imin=i
		minim = rads(i)
		do j=i+1, dots
			if (minim>rads(j)) then
				imin = j
				minim = rads(j)
			end if
		end do
		if (imin>i) then
			buf = rads(i)
			rads(i) = rads(imin)
			rads(imin) = buf

			buf = magn_field(i)
			magn_field(i) = magn_field(imin)
			magn_field(imin) = buf
		end if
	end do

	i=1
	jet_bound = 0.15*b_scale
	do while ((magn_field(i)-jet_bound)*(magn_field(i+1)-jet_bound)>1E-7)
		i = i+1
	end do
	strt_jet_num = i
  end subroutine strt_jets


  subroutine err_compute

	real (kind=precis), dimension(1:2) :: norr1

	norr1 = (/1.0,0.0/)
	do nElem=1,nElems_b
		j=bc_elems(nElem)
		if (dot_product(n(:,j,bc_elem_edge(nElem)),norr1)>0.1 .and. U(2,j)>1E-5) errr = 0.0
	end do

  end subroutine err_compute

end module flux