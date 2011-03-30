! =========================================================
! В модуле содержится информация, относящаяяся к временной сетке
!   - шаг по времени
!   - момент времени начала расчета
!   - момент времени окончания расчета
!   - номер временного слоя
!   и т.д.
! =========================================================

module time
  use names
  use trot
  use solution

  implicit none

contains

  ! ===================================
  ! Инициализация временных параметров
  subroutine init_time
    implicit none

    t_start = 0.0
    t_final = 4.0
    t_current = t_start
    tau = 1.0E-3
    ntime = 0
    tau_out = 0.1
  end subroutine init_time

  ! ==============================================
  ! Проверка условия окончания расчета по времени
  function is_finished_time() result (flag)
    implicit none

    logical:: flag

    flag = .FALSE.

    if (t_current.GE.t_final) flag = .TRUE.

! !     print *, "is_finished_time"
! !     print *, tau, t_current, t_final

  end function is_finished_time

  subroutine time_step_compute

!    tau = min(CFL*hmin/SpecRadius_gas(),1E-3)
! !     if (flag_magnetic>0.5) then
! !       tau = CFL*hmin/SpecRadius()
! !     else
      tau = CFL*hmin/SpecRadius_gas()
!      print *, CFL,hmin,SpecRadius_gas()
! !     endif

  end subroutine time_step_compute

  function SpecRadius() result (SpecRadius1)
    real (kind=precis) :: SpecRadius1
    integer :: nElem, nEdge
	real (kind=precis) :: rho, v1, v2, v3, E, p, k_e, asound,afast, B2, a1,a2
	real (kind=precis), dimension (1:5) :: UU
	real (kind=precis), dimension (1:3) :: BB
	real (kind=precis), dimension (1:3,1:3) :: Tb
	real (kind=precis), dimension (1:2) :: nT

	SpecRadius1 = 0.0
	do nElem=1,nElems
	  if (elems_center(1,nElem)<0.9*z_max .and. elems_center(2,nElem)<0.9*r_max) then
	  do nEdge=1,3
	    UU(:) = U(:,nElem)
		BB(:) = B(:,nElem)
		nT = n(:,nElem,nEdge)
	    call T_rot(n(:,nElem,nEdge), T, Ti)
		Tb(1,:) = (/ nT(1), nT(2), 0._precis/)
		Tb(2,:) = (/-nT(2), nT(1), 0._precis/)
		Tb(3,:) = (/ 0._precis , 0._precis , 1._precis/)
		UU = matmul(T,UU)
		BB = matmul(Tb,BB)
		B2 = sum(BB**2)
	    rho = UU(1)
		v1 = UU(2)/rho
		v2 = UU(3)/rho
		v3 = UU(4)/rho
		E = UU(5)
		k_e = v1*v1+v2*v2+v3*v3 !CHECK!
		p = (E-0.5*rho*k_e)*(gmm-1.0)
		if (p*gmm/rho<0.0) then
		  print *, '***************'
		  print *, 'Sound speed < 0'
		  print *, '***************'
		  print *, 'ntime = ',ntime
		  print *, 'tau = ',tau
		  print *, 'nElem = ',nELem
		  print *, 'Elems_center = ',elems_center(:,nElem)
		  print *, '***************'
		  print *, 'rho = ', U(1,nElem)
		  print *, 'u = ', U(2,nElem)/U(1,nElem)
		  print *, 'v = ', U(3,nElem)/U(1,nElem)
		  print *, 'w = ', U(4,nElem)/U(1,nElem)
		  print *, 'E = ', U(5,nElem)
		  print *, '***************'
		  print *, 'p = ', p
		  print *, '***************'
		  print *, 'Bx = ', B(1,nElem)
		  print *, 'By = ', B(2,nElem)
		  print *, 'Bz = ', B(3,nElem)
		  print *, '***************'
		  print *, 'Bn1 = ', B_nodes(:,elems(1,nElem))
		  print *, 'Bn2 = ', B_nodes(:,elems(2,nElem))
		  print *, 'Bn3 = ', B_nodes(:,elems(3,nElem))
		  print *, '***************'
!		  call save_solution
		  read *
!		  p = -p
		end if
		asound = p*gmm/rho
		a1 = asound + B2/(4.*pi*rho) + abs(BB(1))*sqrt(asound)/sqrt(pi*rho)
		a2 = asound + B2/(4.*pi*rho) - abs(BB(1))*sqrt(asound)/sqrt(pi*rho)
		afast = 0.5 * (sqrt(a1)+sqrt(a2))
		if (abs(v1)+abs(afast)>SpecRadius1) then
			SpecRadius1 = abs(v1)+abs(afast)
			nmsr = nElem
		end if
	  end do
	  end if
	end do
  end function


  function SpecRadius_gas() result (SpecRadius1)
    real (kind=precis) :: SpecRadius1
    integer :: nElem, nEdge
    real (kind=precis) :: rho, v1, v2, v3, E, p, k_e, asound
    real (kind=precis), dimension (1:5) :: UU
    real (kind=precis), dimension (1:3,1:3) :: Tb
    real (kind=precis), dimension (1:2) :: nT

    SpecRadius1 = 0.0
    do nElem=1,nElems
      UU(:) = U(:,nElem)
      rho = UU(1)
      v1 = UU(2)/rho
      v2 = UU(3)/rho
      v3 = UU(4)/rho
      E = UU(5)
      k_e = v1*v1+v2*v2+v3*v3 !CHECK!
      p = (E-0.5*rho*k_e)*(gmm-1.0)
      if (p*gmm/rho<0.0) then
        print *, '***************'
        print *, 'Sound speed < 0'
        print *, '***************'
	print *, 'ntime = ',ntime
	print *, 'tau = ',tau
        print *, 'nElem = ',nELem
        print *, 'Elems_center = ',elems_center(:,nElem)
        print *, '***************'
        print *, 'rho = ', U(1,nElem)
        print *, 'u = ', U(2,nElem)/U(1,nElem)
        print *, 'v = ', U(3,nElem)/U(1,nElem)
        print *, 'w = ', U(4,nElem)/U(1,nElem)
        print *, 'E = ', U(5,nElem)
        print *, '***************'
        print *, 'p = ', p
        print *, '***************'
	print *, 'Solution print'
	call save_vtk_solution
        print *, '***************'
!		  call save_solution
        read *
!		  p = -p
      end if
      if (SpecRadius1<sqrt(k_e)+p*gmm/rho) then
	SpecRadius1 = sqrt(k_e)+p*gmm/rho
        nmsr = nElem
      end if
    end do
  end function

end module time