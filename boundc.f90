! ==========================================================
! Содержит процедyру и данные, необходимые для задания
! граничных и начальных условий
! ==========================================================

module boundc
  use names
  use trot
  use limiters
  use flux
  use initials
  use gravity
  use geometry
!  use boundc_general
!  use boundc_axis
!  use boundc_periodic
!  use boundc_radjet2d

  implicit none

contains

  ! ==============================================
  ! восстановление газовых переменных на границе
  ! и V в точках на верхнем слое
  ! ex update_U_hat
subroutine update_U_til
  implicit none

  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis), dimension(1:2):: xy1 = 0.0, xy2 = 0.0, e1 = 0.0, n1 = 0.0,infl
  real (kind=precis):: V(1:5) = 0.0
  real (kind=precis) :: V2(1:2, 1:5) = 0.0
  integer :: k
  integer:: nElem,nNode
  integer:: NG1 = 0, NG2 = 0
  real (kind=precis):: h1 = 0.0, nor = 0.0
  integer :: i,j, j1,nrb1
  real (kind=precis), dimension(1:2) :: buf
  real (kind=precis), dimension(1:5) :: U_l, U_r, U_l1, VWind,VWall
  real (kind=precis)    :: uL, vL, uR, vR, uLR,vLR, rhoL,rhoR
  real (kind=precis), dimension(1:2) :: x01, x02, xc,norm1,norm2
  real (kind=precis) :: rho1, pres, rad, B_poloidal_mod, B_switch, sw_wei, wind_vel,rho2
  real (kind=precis), dimension(2) :: wind
  real (kind=precis), dimension(1:3):: wind_dir
    ! граничное условие U_til
    ! цикл по граничным ячейкам
  do nElem = 1,nElems_b
    i = bc_elems(nElem)
    select case (bc_type(nElem))

      case (5)
          ! экваториальная симметрия
                !!!!!!!!!! гран условие - все различие, остальное одинаково
        V = (/U(1,i), 0.0_precis, U(3,i), U(4,i), U(5,i) - 0.5*(U(2,i)**2)/U(1,i) /)
        flow_elem = BoundFlow(nElem, V) !?????????????????????????????????????
        A = -tau/elems_vol(i)

        U_til(:,i) = U(:,i) + flow_elem*A
        U_til(:,i) = U_til(:,i) - tau*H_Geom(U_til(:, i))/elems_center(2,i)
        U_til(:,i) = U_til(:,i) + tau*GRAV(i,U_til(:,i))

        U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
        U_til(2,nElems + nElem) = -U_til(2,bc_elems(nElem))
        U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
        U_til(4,nElems + nElem) = U_til(4,bc_elems(nElem))
        U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))

      case (1)
                  ! диск
        B_poloidal_mod = sqrt(B(1,i)**2 + B(2,i)**2)
        B_switch = Abs(B(2,i)/B_poloidal_mod)
        ! вес для диска: 1 - ветра нет, 0 - только ветер
        sw_wei = max(1.08*(Pi/2D0 - ATan((B_switch - 0.5)/0.05))/Pi - 0.05,0D0)

        rho1 = U_til(1,i)
        pres = (gmm-1.0)*(U_til(5,i)-0.5*sum(U_til(2:4,i)**2)/rho1)
        VWall(1) = rho1
        VWall(2) = - U_til(2,i)
        VWall(3) = 0D0
        VWall(4) = rho1*omeg(elems_center(2,i))*elems_center(2,i)
        VWall(5) = pres/(gmm-1.0) + 0.5*sum(VWall(2:4)**2)/rho1

!        rho1 = 1D0
!         rho2 = rho1+(0.01-rho1)*Exp(-10D0*(rad**2))
!         rho1 = max(rho1,rho2)
        wind_dir = B(:,i)/norm2_3(B(:,i))
        wind_vel = sqrt(sum(U_til(2:4,i)**2)/(rho1**2))
!        wind_vel = sqrt(sum(B(:,i)**2)/(4D0*pi))
        VWind(1) = rho1
        VWind(2:4) = rho1*wind_vel*wind_dir
        VWind(4) = VWind(4) + rho1*omeg(elems_center(2,i))*elems_center(2,i)
        VWind(5) = pres/(gmm-1.0) + 0.5*sum(VWind(2:4)**2)/rho1

        U_til(:,nElems + nElem) = sw_wei * VWall + (1D0-sw_wei) * VWind

!         if (-B(2,i)/B_poloidal_mod>0.5 .and. U(2,i)>1E-5) then
!           rho1 = U_til(1,bc_elems(nElem))
!           U_til(1,nElems + nElem) = rho1
!           U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
!           U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
!           U_til(4,nElems + nElem) = rho1*omeg(elems_center(2,i))*elems_center(2,i)
!           pres = (gmm-1.0)*(U_til(5,bc_elems(nElem))-0.5*sum(U_til(2:4,bc_elems(nElem))**2)/rho1)
!           U_til(5,nElems + nElem) = pres/(gmm-1.0)+0.5*sum(U_til(2:4,nElems + nElem)**2)/rho1
!         else
!                 !!!!!!!!!! гран условие - все различие, остальное одинаково
!           pres = (gmm-1.0)*(U(5,i)-0.5*sum(U(2:4,i)**2)/U(1,i))
!           V = (/ U(1,i), 0.0_precis, 0.0_precis, U(1,i)*omeg(elems_center(2,i))*elems_center(2,i), 0.0_precis /)
!           V(5) = pres/(gmm-1.0)
!           flow_elem = BoundFlow(nElem, V) !?????????????????????????????????????
!           A = -tau/elems_vol(i)
!           U_til(:,bc_elems(nElem)) = U(:,i) + flow_elem*A
!           U_til(:,nElems + nElem) = V
!         end if

      case (2)
                         ! симметрия
        U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
        U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
        U_til(3,nElems + nElem) = -U_til(3,bc_elems(nElem))
        U_til(4,nElems + nElem) = -U_til(4,bc_elems(nElem))
        U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))

      case (3)

           ! вытекание
        j = elems_elems(Q(bc_elem_edge(nElem)),i)
        j1 = elems_elems(Q(1+bc_elem_edge(nElem)),i)
        U_r = U_til(:,i)
        U_l = U_til(:,nElems+nElem)
                !                OutFlow(Num,    pU,     Ud,    Un1,     Un2,  bU, nrb)
        call OutFlow(nElem, U_r, U(:,i), U(:,j), U(:,j1), U_l, nrb1)
        if (nrb1 == 0) U_til(:,i) = U_r
        U_til(:,nElems+nElem) = U_l

      case (4)
                         ! аккреция
        rho1 = init_dens
!                       pres = 0.01
        pres = pres95(elems_center(:,i))
        infl = inflow_radjet2D(elems_center(1,i),refp(2),refp,refMach,refsspeed)
        U_l(1) = rho1
        U_l(2) = rho1*infl(1)
        U_l(3) = rho1*infl(2)
        U_l(4) = 0.0_precis
        U_l(5) = pres/(gmm-1.0)+0.5*sum(infl**2)/rho1

        U_til(1:5,nElems + nElem) = U_l

      case default
        print *, 'ERROR: update_u_til: incorrect BC type'
        read *
     end select
  end do

  call compute_upstr
    ! Скорость в точках
  do i=1,nNodes
    nor = 0.0
    V_nodes(1:2,i) = 0.0
    do j=1,node_elems_count(i)
      buf = (/ U_til(2,node_elem(j,i)), U_til(3,node_elem(j,i))/)
      buf = buf/U_til(1,node_elem(j,i))
      h1 = 1.0
      V_nodes(1:2,i) = V_nodes(1:2,i) + buf*h1
      nor = nor + h1
    end do
    if (nor<epsil) then
      nor = 0.0
      V_nodes(1:2,i) = 0.0
      do j=1,node_elems_count(i)
        buf = (/ U_til(2,node_elem(j,i)), U_til(3,node_elem(j,i))/)
        buf = buf/U_til(1,node_elem(j,i))
        h1 = elems_vol(node_elem(j,i))
        V_nodes(1:2,i) = V_nodes(1:2,i) + buf*h1
        nor = nor + h1
      end do
    end if
    V_nodes(1:2,i) = V_nodes(1:2,i)/nor
  end do

  ! ==============================================
  ! вычисление скорости на ребрах на
  ! промежуточном временном слое
  do i=1,nEdges
    V_edges(2,i) = 0.5*(U_til(4,edge_elems(1,i))/U_til(1,edge_elems(1,i))+U_til(4,edge_elems(2,i))/U_til(1,edge_elems(2,i)))
    uL = U_til(2,edge_elems(1,i))/U_til(1,edge_elems(1,i))
    vL = U_til(3,edge_elems(1,i))/U_til(1,edge_elems(1,i))
    uR = U_til(2,edge_elems(2,i))/U_til(1,edge_elems(2,i))
    vR = U_til(3,edge_elems(2,i))/U_til(1,edge_elems(2,i))
    rhoL = U_til(1,edge_elems(1,i))
    rhoR = U_til(1,edge_elems(2,i))
    uLR = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
    vLR = (sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
    V_edges(1,i) = uLR*edge_n(1,i) + vLR*edge_n(2,i)
  end do

  do nElem = 1,nElems_b
    i = bc_elems(nElem)
!         rad = elems_center(2,bc_elems(nElem))
    select case (bc_type(nElem))
      case (5)
        V_edges(1,bc_edges(nElem)) = 0.0
        do k=1,2
          rad = nodes(2,bc_nodes(k,nElem))
          V_nodes(1,bc_nodes(k,nElem)) = 0.0
        end do

      case (1)
          ! диск


!         B_poloidal_mod = sqrt(B(1,i)**2 + B(2,i)**2)
!         if (-B(2,i)/B_poloidal_mod>0.5 .and. U(2,i)>1E-5) then
          V_edges(1,bc_edges(nElem)) = edge_n(1,bc_edges(nElem))*U_til(2,nElems+nElem)+edge_n(2,bc_edges(nElem))*U_til(3,nElems+nElem)
          V_edges(1,bc_edges(nElem)) = V_edges(1,bc_edges(nElem))/U_til(1,nElems+nElem)
          V_edges(2,bc_edges(nElem)) = U_til(4,nElems+nElem)/U_til(1,nElems+nElem)
!         else
!           V_edges(1,bc_edges(nElem)) = 0.0
!         end if
!         V_edges(2,bc_edges(nElem)) = omeg(elems_center(2,i))*elems_center(2,i)
        do k=1,2
          V_nodes(1,bc_nodes(k,nElem)) = 0.0
          V_nodes(2,bc_nodes(k,nElem)) = 0.0
        end do

      case (3)
           ! вытекание
        V_edges(1,bc_edges(nElem)) = edge_n(1,bc_edges(nElem))*U_til(2,nElems+nElem)+edge_n(2,bc_edges(nElem))*U_til(3,nElems+nElem)
        V_edges(2,bc_edges(nElem)) = U_til(4,nElems+nElem)
        V_edges(:,bc_edges(nElem)) = V_edges(:,bc_edges(nElem))/U_til(1,nElems+nElem)

      case (2)
             ! симметрия
        V_edges(1,bc_edges(nElem)) = 0.0
        V_edges(2,bc_edges(nElem)) = 0.0
        V_nodes(2,bc_nodes(1,nElem)) = 0.0
        V_nodes(2,bc_nodes(2,nElem)) = 0.0
      case (4)
                         ! аккреция
        rho1 = sqrt(sum(B(:,i)**2))
        if(rho1>epsil) then
          V_edges(1,bc_edges(nElem)) = edge_n(1,bc_edges(nElem))*U_til(2,bc_elems(nElem))+edge_n(2,bc_edges(nElem))*U_til(3,bc_elems(nElem))
          V_edges(2,bc_edges(nElem)) = U_til(4,bc_elems(nElem))
          V_edges(:,bc_edges(nElem)) = V_edges(:,bc_edges(nElem))/U_til(1,bc_elems(nElem))
        else
          infl = inflow_radjet2D(elems_center(1,i),refp(2),refp,refMach,refsspeed)
          V_edges(1,bc_edges(nElem)) = -infl(2)
          V_edges(2,bc_edges(nElem)) = 0.0
          infl = inflow_radjet2D(nodes(1,bc_nodes(1,nElem)),refp(2),refp,refMach,refsspeed)
          V_nodes(:,bc_nodes(1,nElem)) = infl
          infl = inflow_radjet2D(nodes(1,bc_nodes(2,nElem)),refp(2),refp,refMach,refsspeed)
          V_nodes(:,bc_nodes(2,nElem)) = infl
        end if

      case default
        print *, 'ERROR: update_bc: incorrect BC type'
        read *
    end select
  end  do
end subroutine update_U_til

subroutine update_B_hat

  implicit none
  integer :: i,j
  real (kind=precis), dimension(1:2) :: tau_r,tau_l
  real (kind=precis)    :: A, A1, A2, bet1, bet2, dist1
  real (kind=precis)    :: BzL,BzR, rad, rhoL, rhoR, nor, rho1
  real (kind=precis), dimension(1:2) :: x01, x02, xc, x03
  real (kind=precis):: V(1:5) = 0.0
  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis) :: V2(1:2, 1:5) = 0.0
  integer :: k

  ! Bx,y на границе
  do nElem = 1,nElems_b
    i = bc_elems(nElem)
    select case (bc_type(nElem))
      case (1,3,5)
      case (2)
                  ! симметрия
        B_edges_hat(1,bc_edges(nElem)) = 0.0
        B_edges_hat(2,bc_edges(nElem)) = 0.0
      case (4)
                  ! аккреция
        rho1 = sqrt(sum(B(:,i)**2))
        if(rho1<epsil) then
          B_edges_hat(1,bc_edges(nElem)) = 0.0
          B_edges_hat(2,bc_edges(nElem)) = 0.0
        end if
      case default
        print *, 'ERROR: update_b_hat: incorrect BC type'
        read *
    end select
  end do

  ! Bx,y в ячейках
  do i=1,nElems
    B_hat(1:2,i) = 0.0
    do j=1,3
      A = dot_product(n(:,i,j),edge_n(:,elem_edges(j,i)))*B_edges_hat(1,elem_edges(j,i))
      A = A*h(i,j)
      tau_r = A*(e(:,i,Q(Q(j)))-e(:,i,Q(j)))
      B_hat(1:2,i) = B_hat(1:2,i) + tau_r
    end do
    B_hat(1:2,i) = B_hat(1:2,i)/(6.0*elems_vol(i))
  end do

  ! Bx,y на границе
  do nElem = 1,nElems_b
    i = bc_elems(nElem)
    rad = elems_center(2,bc_elems(nElem))
    select case (bc_type(nElem))
      case (5)
        ! экваториальная симметрия
        B_hat(1,nElems + nElem) = B_hat(1,bc_elems(nElem))
        B_hat(2,nElems + nElem) = B_hat(2,bc_elems(nElem))
        B_hat(3,nElems + nElem) = B_hat(3,bc_elems(nElem))
      case (1)
        ! диск
        B_hat(1,nElems + nElem) = b_scale*magnz(nodes(2,I))
        B_hat(2,nElems + nElem) = B_hat(2,bc_elems(nElem))
        B_hat(3,nElems + nElem) = B_hat(3,bc_elems(nElem))
      case (2)
        ! симметрия
        B_hat(1,nElems + nElem) = B_hat(1,bc_elems(nElem))
        B_hat(2,nElems + nElem) = -B_hat(2,bc_elems(nElem))
        B_hat(3,nElems + nElem) = -B_hat(3,bc_elems(nElem))
      case (3)
        ! вытекание
        B_hat(1,nElems + nElem) = B_hat(1,bc_elems(nElem))
        B_hat(2,nElems + nElem) = B_hat(2,bc_elems(nElem))
        B_hat(3,nElems + nElem) = B_hat(3,bc_elems(nElem))
      case (4)
        ! аккреция
        rho1 = sqrt(sum(B(:,i)**2))
        if(rho1>epsil) then
          B_hat(1:3,nElems + nElem) = B_hat(1:3,bc_elems(nElem))
        else
          B_hat(1,nElems + nElem) = 0.0
          B_hat(2,nElems + nElem) = 0.0
          B_hat(3,nElems + nElem) = 0.0
        end if
      case default
        print *, 'ERROR: update_b_hat: incorrect BC type'
        read *
     end select
  end do

  call limit_B_hat
  do i=1,nElems_b
    B_hat_angles(:,:,nElems+i) = 0.0
  end do
    ! Bx,y в точках
    ! противопоток, проба 2
  call compute_upstr

  do i=1,nNodes
    do j=1,3
      B_nodes(j,i) = 0.0
    end do

    A1 = 0.0
    xc = nodes(:,i)
    do j=1,node_elems_count(i)
      if (node_elem(j,i)<nElems+1) then
        x01 = elems_center(:,node_elem(j,i))
        x01 = xc - x01
        bet1 = upstr_nodes(j,i)
        A1 = A1 + bet1
        do k=1,2
          B_nodes(k,i) = B_nodes(k,i) + (B_hat(k,node_elem(j,i))+dot_product(B_hat_angles(:,k,node_elem(j,i)),x01))*bet1
        end do
      end if
    end do

    if (A1>1E-7) then
      do k=1,2
        B_nodes(k,i) = B_nodes(k,i)/A1
      end do
    else
      do j=1,2
        B_nodes(j,i) = 0.0
      end do
      A = 0.0
      xc = nodes(:,i)
      do j=1,node_elems_count(i)
        if (node_elem(j,i)<nElems+1) then
          x01 = elems_center(:,node_elem(j,i))
          rhoL = 1.0
          A = A + rhoL
          do k=1,2
            B_nodes(k,i) = B_nodes(k,i) + rhoL*(B_hat(k,node_elem(j,i)) + dot_product(B_hat_angles(:,k,node_elem(j,i)),xc-x01))
          end do
        endif
      end do
      B_nodes(1:2,i) = B_nodes(1:2,i)/A
    end if

    B_nodes(3,i) = 0.0
    A = 0.0
    xc = nodes(:,i)
    do j=1,node_elems_count(i)
      if (node_elem(j,i)<nElems+1) then
        x01 = elems_center(:,node_elem(j,i))
        rhoL = 1.0
        A = A + rhoL
        B_nodes(3,i) = B_nodes(3,i) + rhoL*B_hat(3,node_elem(j,i))
      endif
    end do
    B_nodes(3,i) = B_nodes(3,i)/A
  end do

! Bz противопоточно

  do i=1,nEdges
    x01 = elems_center(:,edge_elems(1,i))
    x02 = elems_center(:,edge_elems(2,i))
    xc  = 0.5*(nodes(:,edges(1,i))+nodes(:,edges(2,i)))
    BzL = B_hat(3,edge_elems(1,i)) + dot_product(B_hat_angles(:,3,edge_elems(1,i)),xc-x01)
    BzR = B_hat(3,edge_elems(2,i)) + dot_product(B_hat_angles(:,3,edge_elems(2,i)),xc-x02)
    bet1 = up_koef_z(elem_edgeN(i,1) ,edge_elems(1,i))
    bet2 = up_koef_z(elem_edgeN(i,2) ,edge_elems(2,i))
    A1 = bet1+bet2
    if (A1>1E-7 .and. max(edge_elems(1,i),edge_elems(2,i))<nElems+1 ) then
      B_edges_hat(2,i) = (BzL*bet1+BzR*bet2)/A1
    else
      rhoL = U_til(1,edge_elems(1,i))
      rhoR = U_til(1,edge_elems(2,i))
      B_edges_hat(2,i) = (BzL*sqrt(rhoR)+sqrt(rhoL)*BzR)/(sqrt(rhoL)+sqrt(rhoR))
    end if
  end do

  do nElem = 1,nElems_b
    i = bc_elems(nElem)
    select case (bc_type(nElem))
      case (1)
        B_nodes(1,bc_nodes(1,nElem)) = b_scale*magnz(nodes(2,I))
        B_nodes(1,bc_nodes(2,nElem)) = b_scale*magnz(nodes(2,I))
      case (2)
                  ! симметрия
        B_nodes(2,bc_nodes(1,nElem)) = 0.0
        B_nodes(2,bc_nodes(2,nElem)) = 0.0
        B_edges_hat(2,bc_edges(nElem)) = 0.0
      case (4)
                ! аккреция
        rho1 = sqrt(sum(B(:,i)**2))
        if(rho1<epsil) then
          B_edges_hat(1,bc_edges(nElem)) = 0.0
          B_edges_hat(2,bc_edges(nElem)) = 0.0
          B_nodes(:,bc_nodes(1,nElem)) = 0.0
          B_nodes(:,bc_nodes(2,nElem)) = 0.0
        end if
      case (3,5)

      case default
           print *, 'ERROR: update_b_hat: incorrect BC type'
           read *
     end select
  end do


end subroutine update_B_hat

subroutine bound_B
  implicit none

  real (kind=precis):: V(1:5) = 0.0
  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis), dimension(1:2):: xy1 = 0.0, xy2 = 0.0
  real (kind=precis) :: V2(1:2, 1:5) = 0.0
  integer :: k

  integer:: nElem,nNode
  integer:: NG1 = 0, NG2 = 0
  real (kind=precis):: rad, rho1

    ! цикл по граничным ячейкам
  do nElem = 1,nElems_b
    i = bc_elems(nElem)
    select case (bc_type(nElem))
      case (5)
        B(1,nElems + nElem) = B(1,bc_elems(nElem)) !Bx
        B(2,nElems + nElem) = B(2,bc_elems(nElem)) !By
        B(3,nElems + nElem) = B(3,bc_elems(nElem)) !Bz
      case (1)
                ! диск
        B(1,nElems + nElem) = b_scale*magnz(nodes(2,I)) !Bx
        B(2,nElems + nElem) = B(2,bc_elems(nElem)) !By
        B(3,nElems + nElem) = B(3,bc_elems(nElem)) !Bz
        B_nodes(1,bc_nodes(1,nElem)) = b_scale*magnz(nodes(2,I))
        B_nodes(1,bc_nodes(2,nElem)) = b_scale*magnz(nodes(2,I))

      case (2)
             ! симметрия
        B(1,nElems + nElem) = B(1,bc_elems(nElem))
        B(2,nElems + nElem) = -B(2,bc_elems(nElem))
        B(3,nElems + nElem) = -B(3,bc_elems(nElem))

        B_edges(1,bc_edges(nElem)) = 0.0
        B_edges(2,bc_edges(nElem)) = 0.0
        B_nodes(2,bc_nodes(1,nElem)) = 0.0
        B_nodes(2,bc_nodes(2,nElem)) = 0.0

      case (3)
         ! вытекание
        B(1,nElems + nElem) = B(1,bc_elems(nElem))
        B(2,nElems + nElem) = B(2,bc_elems(nElem))
        B(3,nElems + nElem) = B(3,bc_elems(nElem))
        B_edges(2,bc_edges(nElem)) = B(3,bc_elems(nElem))

      case (4)

        rho1 = sqrt(sum(B(:,i)**2))
        if(rho1>epsil) then
          B(:,nElems + nElem) = B(:,bc_elems(nElem))
          B_edges(2,bc_edges(nElem)) = B(3,bc_elems(nElem))
        else
          B(:,nElems + nElem) = 0.0
          B_edges(:,bc_edges(nElem)) = 0.0
          B_nodes(1,bc_nodes(1,nElem)) = 0.0
          B_nodes(2,bc_nodes(1,nElem)) = 0.0
          B_nodes(1,bc_nodes(2,nElem)) = 0.0
          B_nodes(2,bc_nodes(2,nElem)) = 0.0
        end if

      case default
        print *, 'ERROR: update_bc: incorrect BC type'
        read *
    end select
  end  do
end subroutine bound_B

subroutine update_bc
  implicit none

  real (kind=precis) :: pres, rho1,rad,rho2
  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis), dimension(1:2):: xy1 = 0.0, xy2 = 0.0, e = 0.0, n1 = 0.0, norm1, norm2, infl
  real (kind=precis):: V(1:5) = 0.0, U_r(1:5), U_l(1:5)
  integer:: nElem, i,j,j1,nrb1
  integer:: NG1 = 0, NG2 = 0
  real (kind=precis):: h = 0.0, signum
  real (kind=precis) :: B_poloidal_mod
  real (kind=precis) :: V2(1:2, 1:5) = 0.0
  integer :: k
  real (kind=precis) :: dens
  real (kind=precis), dimension(2) :: wind
  real (kind=precis) :: B_switch, sw_wei, wind_vel
  real (kind=precis), dimension(1:3):: wind_dir
  real (kind=precis), dimension(1:5):: VWall, VWind
    ! цикл по граничным ячейкам
  do nElem = 1,nElems_b

    i = bc_elems(nElem)
                ! U(:,nElems + nElem) = (rho, rho*u, rho*v, rho*w, E)
    rad = elems_center(2,bc_elems(nElem))

    select case (bc_type(nElem))

      case (5)
          ! экваториальная симметрия
          U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
          U_til(2,nElems + nElem) = -U_til(2,bc_elems(nElem))
          U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
          U_til(4,nElems + nElem) = U_til(4,bc_elems(nElem))
          U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))

      case (1)
           !       диск

        B_poloidal_mod = sqrt(B(1,i)**2 + B(2,i)**2)
        B_switch = Abs(B(2,i)/B_poloidal_mod)
        ! вес для диска: 1 - ветра нет, 0 - только ветер
        sw_wei = max(1.08*(Pi/2D0 - ATan((B_switch - 0.5)/0.05))/Pi - 0.05,0D0)

        rho1 = U_til(1,i)
        pres = (gmm-1.0)*(U_til(5,i)-0.5*sum(U_til(2:4,i)**2)/rho1)
        VWall(1) = rho1
        VWall(2) = - U_til(2,i)
        VWall(3) = 0D0
        VWall(4) = rho1*omeg(elems_center(2,i))*elems_center(2,i)
        VWall(5) = pres/(gmm-1.0) + 0.5*sum(VWall(2:4)**2)/rho1

!        rho2 = rho1+(0.01-rho1)*Exp(-10D0*(rad**2))
!        rho1 = max(rho1,rho2)
        wind_dir = B(:,i)/norm2_3(B(:,i))
        wind_vel = sqrt(sum(U_til(2:4,i)**2)/(rho1**2))
!        wind_vel = sqrt(sum(B(:,i)**2)/(4D0*pi))
        VWind(1) = rho1
        VWind(2:4) = rho1*wind_vel*wind_dir
        VWind(4) = VWind(4) + rho1*omeg(elems_center(2,i))*elems_center(2,i)
        VWind(5) = pres/(gmm-1.0) + 0.5*sum(VWind(2:4)**2)/rho1

        U_til(:,nElems + nElem) = sw_wei * VWall + (1D0-sw_wei) * VWind

!           B_poloidal_mod = sqrt(B(1,i)**2 + B(2,i)**2)
!           if (-B(2,i)/B_poloidal_mod>0.5 .and. U_til(2,i)>1E-5) then
!             rho1 = U_til(1,bc_elems(nElem))
!             U_til(1,nElems + nElem) = rho1
!             U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
!             U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
!             U_til(4,nElems + nElem) = rho1*omeg(elems_center(2,i))*elems_center(2,i)
!             pres = (gmm-1.0)*(U_til(5,bc_elems(nElem))-0.5*sum(U_til(2:4,bc_elems(nElem))**2)/rho1)
!             U_til(5,nElems + nElem) = pres/(gmm-1.0)+0.5*sum(U_til(2:4,nElems + nElem)**2)/rho1
!           else
!             rho1 = U_til(1,bc_elems(nElem))
!             U_til(1,nElems + nElem) = rho1
!             U_til(2,nElems + nElem) = -U_til(2,bc_elems(nElem))
!             U_til(3,nElems + nElem) = 0.0
!             U_til(3,bc_elems(nElem)) = 0.0
!             U_til(4,nElems + nElem) = rho1*omeg(elems_center(2,i))*elems_center(2,i)
!             pres = (gmm-1.0)*(U_til(5,bc_elems(nElem))-0.5*sum(U_til(2:4,bc_elems(nElem))**2)/rho1)
!             U_til(5,nElems + nElem) = pres/(gmm-1.0)+0.5*sum(U_til(2:4,nElems + nElem)**2)/rho1
!           end if

! !         U_til(1,nElems + nElem) = 1.0
! !         U_til(2,nElems + nElem) = 0.1
! !         U_til(3,nElems + nElem) = -0.5 * (rad**2)
! !         U_til(4,nElems + nElem) = 0D0
! !         U_til(5,nElems + nElem) = 0.5*sum(U_til(2:4,nElems + nElem)**2)+0.125

      case (2)
                         ! симметрия
        U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
        U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
        U_til(3,nElems + nElem) = -U_til(3,bc_elems(nElem))
        U_til(4,nElems + nElem) = U_til(4,bc_elems(nElem))
        U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))

      case (3)
!         U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
!         U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
!         U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
!         U_til(4,nElems + nElem) = U_til(4,bc_elems(nElem))
!         U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))
         ! вытекание
        j = Q(bc_elem_edge(nElem))
        j1 = Q(1+bc_elem_edge(nElem))
        U_r = U_til(:,i)
        U_l = U_til(:,nElems+nElem)
        call OutFlow(nElem, U_r, U_til(:,i), U_til(:,j), U_til(:,j1), U_l, nrb1)
        if (nrb1 == 0) then
          U_til(:,nElems+nElem) = U_til(:,i)
        else
          U_til(:,nElems+nElem) = U_l
        end if

      case (4)
                         ! аккреция
          rho1 = init_dens
!                       pres = 0.01
          pres = pres95(elems_center(:,i))
          infl = inflow_radjet2D(elems_center(1,i),refp(2),refp,refMach,refsspeed)
          U_l(1) = rho1
          U_l(2) = rho1*infl(1)
          U_l(3) = rho1*infl(2)
          U_l(4) = 0.0
          U_l(5) = pres/(gmm-1.0)+0.5*sum(infl**2)/rho1

          U_til(1:5,nElems + nElem) = U_l
!         U_til(1,nElems + nElem) = U_til(1,bc_elems(nElem))
!         U_til(2,nElems + nElem) = U_til(2,bc_elems(nElem))
!         U_til(3,nElems + nElem) = U_til(3,bc_elems(nElem))
!         U_til(4,nElems + nElem) = U_til(4,bc_elems(nElem))
!         U_til(5,nElems + nElem) = U_til(5,bc_elems(nElem))

        case default
          print *, 'ERROR: update_bc: incorrect BC type', nElem, bc_type(nElem)
          read *
      end select
    end  do
end subroutine update_bc

subroutine OutFlow(Num, pU, Ud, Un1, Un2, bU, nrb)

  implicit none

  real (kind=precis), dimension(1:5), Intent(In):: Un1, Un2, Ud
  real (kind=precis), dimension(1:5), Intent(InOut):: pU, bU
  real (kind=precis), dimension(1:5) :: Flow, U_l, U_r, V, pU1
  real (kind=precis), dimension(1:5,1:2) :: Ubb
  integer, Intent(In) :: Num
  integer, Intent(InOut) :: nrb
  integer :: i,j
  real (kind=precis), dimension(1:2):: norm1, norm2
  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis):: vx, vy, vz, c_a, pres, rho, En, k_e
  real (kind=precis) :: a_g, rho_g, E_g, u_g, v_g, w_g, ke_g, p_g, M_g, vv_g
  real (kind=precis) :: a_0, rho_0, E_0, u_0, v_0, w_0, ke_0, p_0, M_0, vv_0

    !*****************************
  nrb = 0
  i = bc_elems(Num)
  norm1 = n(:,i,bc_elem_edge(Num))
  call T_rot(norm1, T, Ti)

  pU1 = matmul(T,pU)

  rho_g = pU1(1)
  u_g = pU1(2)/rho_g
  v_g = pU1(3)/rho_g
  w_g = pU1(4)/rho_g
  E_g = pU1(5)
  ke_g = 0.5*(u_g**2+v_g**2+w_g**2)
  vv_g = sqrt(2.0*ke_g)
  p_g = (E_g - ke_g * rho_g)*(gmm-1.0)
  a_g = sqrt(gmm*p_g/rho_g)
  M_g = u_g/a_g
  if (M_g<1.0 .and. u_g>1E-2 .and. t_current>1.0) then
    a_0 = (gmm-1.0)*(u_g+2.0*a_g/(gmm-1.0))/(gmm+1.0)
    u_0 = a_g
    rho_0 = rho_g*((a_0/a_g)**(2.0/(gmm-1.0)))
    v_0 = v_g
    w_0 = w_g
    ke_0 = 0.5*(u_0**2+v_0**2+w_0**2)
    p_0 = rho_0*(a_0**2)/gmm
    E_0 = p_0/(gmm-1.0)+ke_0*rho_0
    bU(1) = rho_0
    bU(2) = rho_0*u_0
    bU(3) = rho_0*v_0
    bU(4) = rho_0*w_0
    bU(5) = E_0
    bU = matmul(Ti,bU)
    nrb = 1
  else
    Flow(:) = 0.0
    Ubb(:,1) = Un1
    Ubb(:,2) = Un2

    do nEdge = 1,2
      j = Q(nEdge-1+bc_elem_edge(Num))
      norm2 = n(:,i,j)-2.0*norm1*dot_product(n(:,i,j),norm1)
      call T_rot(norm2, T, Ti)
      U_l = Ud
      U_r = Ubb(:,nEdge)
      U_l = matmul(T,U_l)
      U_r = matmul(T,U_r)
      flow_edge = flux_num_hllc(U_l,U_r)
      flow_edge = matmul(Ti,flow_edge)
      flow_edge = flow_edge*h(i,j)
      Flow = Flow + flow_edge

      call T_rot(n(:,i,j), T, Ti)
      U_l = Ud
      U_r = Ubb(:,nEdge)
      U_l = matmul(T,U_l)
      U_r = matmul(T,U_r)
      flow_edge = flux_num_hllc(U_l,U_r)
      flow_edge = matmul(Ti,flow_edge)
      flow_edge = flow_edge*h(i,j)
      Flow = Flow + flow_edge
    end do
    A = -0.5*tau/elems_vol(i)
    pU = Ud + Flow*A
    if (flag_axial>0.5) pU = pU - tau*H_GEOM(pU)/elems_center(2,bc_elems(Num))
    pU = pU + tau*GRAV(bc_elems(Num),pU)
    bU = pU
  end if
end subroutine OutFlow

function omeg(radius) result (Omm)
  real (kind=precis) :: radius, Omm

  Omm = omega0*(1.0-(radius/r_disk)**2)
!       Omm = 0.0
  if (radius>r_disk) Omm = 0.0
end function omeg

function BoundFlow(Num, pU) result (Flow)

  implicit none

  real (kind=precis), dimension(1:5):: Flow, pU, bU, U_l, U_r
  integer :: Num
  real (kind=precis), dimension(1:5,1:5):: T = 0.0, Ti = 0.0
  real (kind=precis):: vx, vy, vz, c_a, pres, rho, En, k_e

    !*****************************
  i = bc_elems(Num)
  Flow(:) = 0.0

  do nEdge = 1,3
    call T_rot(n(:,i,nEdge), T, Ti)
    if (bc_elem_edge(Num)==nEdge) then
      bU = matmul(T,pU)
      rho = bU(1)
      vx = bU(2)/rho
      vy = bU(3)/rho
      vz = bU(4)/rho              ! полная энергия (внутренняя плюс кинетическая)
      En = bU(5)
      k_e = vx**2+vy**2+vz**2
      pres = (En-0.5*rho*k_e)*(gmm-1.0)
      flow_edge = (/bU(2), bU(2)*vx+pres, bU(2)*vy, bU(2)*vz, vx*(En+pres)/)
      flow_edge = matmul(Ti,flow_edge)
      flow_edge = flow_edge * h(i, nEdge)
    else
      U_l = U(:,i)
      U_r = U(:,elems_elems(nEdge,i))
      U_l = matmul(T,U_l)
      U_r = matmul(T,U_r)
      flow_edge = flux_num_hllc(U_l,U_r)
      flow_edge = matmul(Ti,flow_edge)
      flow_edge = flow_edge*h(i, nEdge)
    end if
    Flow = Flow + flow_edge
  end do

end function BoundFlow






subroutine info_read_create
  implicit none
  logical:: EXT1
INQUIRE (FILE='./data/information.dat', EXIST=EXT1)
  IF (.NOT. EXT1) THEN
         open (unit=1,file='./data/information.dat',form='formatted')

			write(1,*)  ntime
			write(1,*)  t_current
			write(1,*)  tau
			write(1,*)  t_start
			write(1,*)  t_final
			do I=1,nElems
				write(1,*)  U(1,I)
			end do
			do I=1,nElems
				write(1,*)  U(2,I)
			end do
			do I=1,nElems
				write(1,*)  U(3,I)
			end do
			do I=1,nElems
				write(1,*)  U(4,I)
			end do
			do I=1,nElems
				write(1,*)  U(5,I)
			end do
	!		write(1,*)  U(6,1:nElems)
	!		write(1,*)  U(7,1:nElems)
	!		write(1,*)  U(8,1:nElems)

		 close (unit=1,status='Keep')
		 write (*,*) 'File ./data/information.dat has been created'
  else
		open (unit=1,file='./data/information.dat',form='formatted')

			read(1,*)  ntime
			read(1,*)  t_current
			read(1,*)  tau
			read(1,*)  t_start
			read(1,*)  t_final
			do I=1,nElems
				read(1,*)  U(1,I)
			end do
			do I=1,nElems
				read(1,*)  U(2,I)
			end do
			do I=1,nElems
				read(1,*)  U(3,I)
			end do
			do I=1,nElems
				read(1,*)  U(4,I)
			end do
			do I=1,nElems
				read(1,*)  U(5,I)
			end do
	!		read(1,*)  U(6,1:nElems)
	!		read(1,*)  U(7,1:nElems)
	!		read(1,*)  U(8,1:nElems)

		 close (unit=1,status='Keep')
		 write (*,*) 'File ./data/information.dat has been read'

  END IF
  end subroutine info_read_create


 ! ==============================================

  subroutine info_write
  implicit none

	open (unit=1,file='./data/information.dat',form='formatted')

			write(1,*)  ntime
			write(1,*)  t_current
			write(1,*)  tau
			write(1,*)  t_start
			write(1,*)  t_final
			do I=1,nElems
				write(1,*)  U(1,I)
			end do
			do I=1,nElems
				write(1,*)  U(2,I)
			end do
			do I=1,nElems
				write(1,*)  U(3,I)
			end do
			do I=1,nElems
				write(1,*)  U(4,I)
			end do
			do I=1,nElems
				write(1,*)  U(5,I)
			end do
	!		write(1,*)  U(6,1:nElems)
	!		write(1,*)  U(7,1:nElems)
	!		write(1,*)  U(8,1:nElems)

	close (unit=1,status='Keep')
	!write (*,*) 'File ./data/information.dat has been written'

  end subroutine info_write

! ==============================================

end module boundc