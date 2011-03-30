module pressure_correct

	use names
	use trot

implicit none

contains


  function press_correct(NumEl) result (pressure)

	integer :: NumEl, NearN
	real (kind=precis) :: vx,vy,vz, rho, p, k_e, rho1, E, A
	real (kind=precis) :: pressure
	real (kind=precis), dimension(1:5) :: Ur, Ul, Up, Unear, Up_up, U_cons, flow_edge, flow_elem

	U_cons = U(:,NumEl)
	rho = U(1,NumEl)
	rho1 = rho**(1.0-gmm)

	vx = U_cons(2)/rho
	vy = U_cons(3)/rho
	vz = U_cons(4)/rho
	E = U_cons(5)
	k_e = vx*vx+vy*vy+vz*vz
	p = (E-0.5*rho*k_e)*(gmm-1.0)		    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!

	Up(1:4) = U(1:4, NumEl)
	Up(5) = p * rho1
    flow_elem = 0.0
	do nEdge=1,3
	    call T_rot(n(:,NumEl,nEdge), T, Ti)

		NearN = elems_elems(nEdge, NumEl)
		rho = U(1,NearN)
		rho1 = rho**(1.0-gmm)

		vx = U(2,NearN)/rho
		vy = U(3,NearN)/rho
		vz = U(4,NearN)/rho
		E = U(5,NearN)
		k_e = vx*vx+vy*vy+vz*vz
		p = (E-0.5*rho*k_e)*(gmm-1.0)		    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!

		Unear(1:4) = U(1:4, NearN)
		Unear(5) = p * rho1
	    do i=1,5
	      Ul(i) = Up(i)
		  Ur(i) = Unear(i)
	    end do

        Ul = matmul(T,Ul)
        Ur = matmul(T,Ur)
  	    flow_edge = flux_num_hll_pres(Ul,Ur)
        flow_edge = matmul(Ti,flow_edge)
        flow_edge = flow_edge*h(NumEl,nEdge)
        flow_elem = flow_elem + flow_edge
    end do
	A = -tau/elems_vol(NumEl)
	Up_up = Up + flow_elem*A
	rho = Up_up(1)
	rho1 = rho**(1.0-gmm)
	pressure = Up_up(5)/rho1
  end function press_correct

  function flux_num_hll_pres(Ul,Ur) result (FR)
      implicit none

    real (kind=precis), dimension(1:5):: FR ! численный поток
    real (kind=precis), dimension(1:5):: Ul, Ur, U_lr
	real (kind=precis), dimension(1:5):: F_l, F_r, F_lr
	real (kind=precis), dimension(1:5,1:5):: D_r, R_r, L_r
	real (kind=precis), dimension(1:5,1:5):: D_l, R_l, L_l
	real (kind=precis), dimension(1:5,1:5):: D_lr, R_lr, L_lr

	integer:: I
	integer:: J

	call FL_pres(Ul,F_l,D_l)
	call FL_pres(Ur,F_r,D_r)

    U_lr = 0.5*(Ur + Ul)
	call FL_pres(U_lr,F_lr,D_lr)

	if (D_lr(1,1).ge.0.0) then
		FR= F_l
	end if

	if (D_lr(1,1).le.0.0 .and. D_lr(5,5).ge.0.0) then
		FR = (D_lr(5,5)*F_l - D_lr(1,1)*F_r + D_lr(5,5)*D_lr(1,1)*(Ur-Ul))/(D_lr(5,5)-D_lr(1,1))
	end if

	if (D_lr(5,5).le.0.0) then
		FR= F_r
	end if

  end function flux_num_hll_pres

subroutine FL_pres(UU,F,D)
    implicit none

    ! входные переменные
    real (kind=precis), dimension(1:5), intent(in):: UU
	! UU = (rho, rho u, rho v, rho w, p rho^(1-gmm))

    ! выходные переменные
    real (kind=precis), dimension(1:5), intent(out):: F
    real (kind=precis), dimension(1:5,1:5), intent(out)::D
	real (kind=precis):: u, v, w, a, p, rho, E

    real (kind=precis):: H2, k_e, mult, rho1
	real (kind=precis):: courant

    ! вычисляем все простые переменные

	 rho = UU(1)
	 rho1 = rho**(1.0-gmm)

     u = UU(2)/rho
     v = UU(3)/rho
     w = UU(4)/rho									! полная энергия (внутренняя плюс кинетическая)
	 p = UU(5)/rho1
!	 E = UU(5)
	 k_e = u*u+v*v+w*w
	 E = p/(gmm-1.0)+0.5*k_e
 !    p = (E-0.5*rho*k_e)*(gmm-1.0)		    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
	 H2 = (E+p)/rho								! энтальпия


	if (p*gmm/rho.LT. 0.0) then
	  print*, ntime
	  print *, nElem
	  print *, elems_center(:,nElem)
	  print*, p
	  p = 0.0
	end if
	 ! ошибка вычисления SQRT из-за деления на нуль в предыдущих формулах  !!!!
	 a = SQRT(p*gmm/rho)						! скорость звука


    F = (/ UU(2), UU(2)*u+p, UU(2)*v, UU(2)*w, u*p*rho1 /)

    D = 0.0
	D(1,1) = u-a
	D(2,2) = u
	D(3,3) = u
	D(4,4) = u
	D(5,5) = u+a

!	R = RESHAPE((/1.0, u-a, v, w, H2-a*u, 0.0, 0.0, 1.0, 0.0, v, &
!		 0.0, 0.0, 0.0, 1.0, w,  1.0, u, v, w, 0.5*k_e, 1.0, u+a, v, w, H2+a*u/), SHAPE= (/5, 5/))

!	mult = (gmm-1.0)/(2.0*a*a)

!	L = mult*RESHAPE((/H2+1.0/(a*2.0*mult)*(u-a), -v/mult, -w/mult, -2.0*H2+2.0/mult, H2-1.0/(a*2.0*mult)*(u+a), &
 !               -(u+1.0/(a*2.0*mult)), 0.0, 0.0, 2.0*u, -u+1.0/(a*2.0*mult), &
!				-v, 1.0/mult, 0.0, 2.0*v, -v, -w, 0.0, 1.0/mult, 2.0*w, -w, &
!				  1.0, 0.0, 0.0, -2.0, 1.0/), SHAPE= (/5, 5/))

  end subroutine FL_pres
! ======================================================

  subroutine correct_press_1
  	integer :: i
	real (kind=precis):: u1, v1, w1, p, rho, E, k_e, eps1

    ! вычисляем все простые переменные
     eps1 = 1E-5

    do i=1,nElems
		rho = U(1,i)
	    u1 = U(2,i)/rho
	    v1 = U(3,i)/rho
	    w1 = U(4,i)/rho									! полная энергия (внутренняя плюс кинетическая)
		E = U(5,i)
		k_e = u1*u1+v1*v1+w1*w1
   		p = (E-0.5*rho*k_e)*(gmm-1.0)		    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
		p = max(p,eps1)
		U(5,i) = p/(gmm-1.0)+0.5*rho*k_e
  	end do

  end subroutine correct_press_1

  function clean_fl(upb) result(fl_clean)
  	integer :: fl_clean
  	integer :: i,j
  	real (kind=precis) :: eps_clean
  	real (kind=precis) :: upb

  	eps_clean = 1E-6
  	fl_clean = 0
  	do i=1,nElems
  		if (elems_center(2,i)>upb) then
  			do j=1,3
  				if (abs(B_hat(j,i))>eps_clean) fl_clean = 1
  			end do
  		end if
  	end do

  end function clean_fl

  subroutine clean1
	integer :: i,j, ed
	real (kind=precis) :: magmod
	real (kind=precis) :: epsil

	epsil = 1E-4
	do i=1,nElems
		magmod = sqrt(sum(B_hat(:,i)**2))
		if (magmod<epsil) then
			U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
			B_hat(:,i) = 0.0
			do j=1,3
				ed = elem_edges(j,i)
				B_edges_hat(:,ed) = 0.0
			end do
		end if
	end do

  end subroutine clean1

subroutine total_clean
  integer :: i,j, ed
  real (kind=precis) :: magmod
  real (kind=precis) :: epsil, rtc, ztc

  if (t_current>5.0) then
!    if (t_current - cleanprev > 0.5) then
      cleanprev = t_current
      do i=1,nElems
        ztc = elems_center(1,i)
        rtc = elems_center(2,i)
        if (B_hat(1,i)/b_scale<0.0001 ) then
          U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
          B_hat(:,i) = 0.0
          do j=1,3
            ed = elem_edges(j,i)
            B_edges_hat(:,ed) = 0.0
          end do
        end if
      end do
!    end if
    do i=1,nElems
      ztc = elems_center(1,i)
      rtc = elems_center(2,i)
      if (rtc>2.0) then
        U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
        B_hat(:,i) = 0.0
        do j=1,3
          ed = elem_edges(j,i)
          B_edges_hat(:,ed) = 0.0
        end do
      end if
    end do
    if (t_current<5.0+2.*tau) then
      do i=1,nElems
        ztc = elems_center(1,i)
        rtc = elems_center(2,i)
        if ((ztc<1.15 .and. rtc>0.7) .or. (ztc>1.15 .and. rtc>(ztc-1.15)*(1.35-0.65)/(2.5-1.15)+0.65 )) then
          U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
          B_hat(:,i) = 0.0
          do j=1,3
            ed = elem_edges(j,i)
            B_edges_hat(:,ed) = 0.0
          end do
        end if
        if (rtc<0.7 .and. ztc<1.15 .and. B_hat(1,i)/b_scale<0.06) then
          U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
          B_hat(:,i) = 0.0
          do j=1,3
            ed = elem_edges(j,i)
            B_edges_hat(:,ed) = 0.0
          end do
        end if
      end do
    end if
  end if

end subroutine total_clean

  subroutine total_clean2
	integer :: i,j, ed
	real (kind=precis) :: magmod
	real (kind=precis) :: epsil, rtc, ztc
	real (kind=precis) :: jet_rad

	jet_rad = r_max
	do i=1,nElems_b
		j = bc_elems(i)
		if (bc_type(i)==2) then
			if(abs(B_hat(1,j))/b_scale<0.005 .and. elems_center(2,j)<jet_rad) jet_rad = elems_center(2,j)
		end if
	end do

	do i=1,nElems
		ztc = elems_center(1,i)
		rtc = elems_center(2,i)
		if (rtc>jet_rad .and. ztc>z_max-2.0*hmin) then
			U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
			B_hat(:,i) = 0.0
			do j=1,3
				ed = elem_edges(j,i)
				B_edges_hat(:,ed) = 0.0
			end do
		end if
	end do

  end subroutine total_clean2

  subroutine clean
  	integer :: i
  	real (kind=precis) :: upbn,upb1

  	upbn = 0.75*z_max
   	upb1 = 0.9*z_max
  	if (clean_fl(upb1)>0.5) then
		do i=1,nElems
			if (elems_center(2,i)>upbn) then
				U_til(5,i) = U_til(5,i) +  sum(B_hat(:,i)**2)/(8.0*pi)
				B(:,i) = 0.0
				B_hat(:,i) = 0.0
			end if
		end do
		do i=1,nEdges
			if (nodes(2,edges(1,i))>upbn) then
				B_edges(:,i) = 0.0
				B_edges_hat(:,i) = 0.0
			end if
		end do
  	end if
  end subroutine clean
end module pressure_correct