module solvers

  use geometry
  use gravity
  use names
  use limiters
  use flux
  use trot
  use multies
  use boundc
  use solution
  use pressure_correct

  implicit none

contains

subroutine gas_hllc

  implicit none

  integer :: nEdge,elem1,elem2
  real (kind=precis), dimension(1:5) :: flow_edge
  real (kind=precis), dimension(1:2) :: edge_normal

  real (kind=precis) ::vx,vy
  real (kind=precis), dimension(1:2) :: x01, x02, xc,norm1,norm2
  real (kind=precis), dimension(1:5) :: buf1


  call update_bc

  !$OMP PARALLEL DEFAULT(none) Shared(U_til,nEdges,gas_edge_flux,edge_n,edge_vol,elem_edgeN,edge_elems,n,elems_elems) PRIVATE(nEdge,elem1,elem2,flow_edge,Ul,Ur,edge_normal,vx,vy)
  !$OMP DO SCHEDULE(auto)
  do nEdge = 1,nEdges
    ! вычисление потока через ребро
    ! цикл по ребрам
    elem1=edge_elems(1,nEdge)
    elem2=edge_elems(2,nEdge)

    edge_normal(:)=edge_n(:,nEdge)
    Ul(:) = U_til(:,elem1)
    Ur(:) = U_til(:,elem2)

    !поворот
    vx=Ul(2)
    vy=Ul(3)
    Ul(2)= edge_normal(1)*vx + edge_normal(2)*vy
    Ul(3)=-edge_normal(2)*vx + edge_normal(1)*vy

    vx=Ur(2)
    vy=Ur(3)
    Ur(2)= edge_normal(1)*vx + edge_normal(2)*vy
    Ur(3)=-edge_normal(2)*vx + edge_normal(1)*vy

    ! ==================================
    ! вычисляем поток из распадника
    flow_edge=flux_num_hllc(Ul,Ur)

    vx=flow_edge(2)
    vy=flow_edge(3)
    flow_edge(2)= edge_normal(1)*vx - edge_normal(2)*vy
    flow_edge(3)= edge_normal(2)*vx + edge_normal(1)*vy

    flow_edge = flow_edge*edge_vol(nEdge)
    gas_edge_flux(nEdge,:)= flow_edge(:)
  end do
  !$OMP END DO nowait
  !$omp end parallel

  !$OMP PARALLEL DEFAULT(none) Shared(U_hat,elems_vol,elems_center,U_til,edge_normSign,elem_edges,gas_edge_flux,nElems,tau) PRIVATE(nEdge,flow_elem,flow_edge,A,buf1)
  !$OMP DO SCHEDULE(auto)
  do nElem = 1,nElems
    flow_elem = 0D0
    do nEdge = 1,3
      ! ================================================================
      ! добавляем поток через данное ребро в полный поток через границу ячейки
      flow_elem = flow_elem + edge_normSign(nEdge,nElem)*gas_edge_flux(elem_edges(nEdge,nElem),:)
    end do
    ! ====================================================
    ! вычисляем новое значение в ячейке
    ! A -- коэфиициент, зависит от шага сетки и площади ячейки
    A = -tau/elems_vol(nElem)

    U_hat(:,nElem) = U_til(:,nElem) + flow_elem(:)*A

    ! осевая симметрия
    U_hat(:,nElem) = U_hat(:,nElem) - tau*H_Geom(U_hat(:, nElem))/elems_center(2,nElem)

    ! гравитация
    U_hat(:,nElem) = U_hat(:,nElem) + tau*GRAV(nElem,U_hat(:, nElem))

  end do
 !$OMP END DO nowait
 !$omp end parallel
end subroutine gas_hllc

subroutine magnetic
  implicit none
  integer :: i,j
  real (kind=precis),dimension(1:3) :: B_lr, B_l, B_r
  real (kind=precis), dimension(2):: buf, tau_r
  real (kind=precis) :: dist
  real (kind=precis) :: r1, r2, r_mid

  call update_U_til

  !========================================================
  ! Изменение магнитного поля по закону Фарадея
  call bound_B

  ! Пересчет B
  ! = Bn

  !$OMP PARALLEL DEFAULT(none) Shared(B_edges_hat,nodes,edges,V_nodes,B_nodes,tau,edge_vol,nEdges,B_edges) PRIVATE(r1,r2,r_mid,circ,nEdge,A)
  !$OMP DO SCHEDULE(auto)
  do nEdge=1,nEdges
    r1 = nodes(2,edges(1,nEdge))
    r2 = nodes(2,edges(2,nEdge))
    r_mid = 0.5*(r1+r2)
    circ = r2*(V_nodes(1,edges(2,nEdge))*B_nodes(2,edges(2,nEdge)) - B_nodes(1,edges(2,nEdge))* V_nodes(2,edges(2,nEdge)))
    circ = circ - r1*((V_nodes(1,edges(1,nEdge))*B_nodes(2,edges(1,nEdge)) - B_nodes(1,edges(1,nEdge))* V_nodes(2,edges(1,nEdge))))
    if (r_mid>1E-7) then
      A = tau/(edge_vol(nEdge)*r_mid)
      B_edges_hat(1,nEdge) = B_edges(1,nEdge) + A*circ
    else
      B_edges_hat(1,nEdge) = 0.0
    end if
  end do
  !$OMP END DO nowait
  !$omp end parallel

  ! = B_phi
  !$OMP PARALLEL DEFAULT(none) Shared(h,edge_normSign,V_edges,B_edges,elems_vol,tau,B,B_hat,elem_edges,nElems) PRIVATE(nElem,circ,nEdge,A)
  !$OMP DO SCHEDULE(auto)
  do nElem=1,nElems
    circ = 0D0
    do nEdge=1,3
      A = h(nElem,nEdge)*edge_normSign(nEdge,nElem)
      circ = circ - A*( V_edges(1,elem_edges(nEdge,nElem))*B_edges(2,elem_edges(nEdge,nElem))- V_edges(2,elem_edges(nEdge,nElem))*B_edges(1,elem_edges(nEdge,nElem)))
    end do
    A = tau/elems_vol(nElem)
    B_hat(3,nElem) = B(3,nElem) + A*circ
  end do
  !$OMP END DO nowait
  !$omp end parallel

  call update_B_hat

end subroutine magnetic

subroutine rotB_incell
  real (kind=precis) :: A1, A, rhoL, rhoR, bet1, bet2, circ
  real (kind=precis), dimension(2)::B_L, B_R, tang, xc, x01,x02, tau_r
  integer :: i
  real (kind=precis) :: r1, r2, r_mid

  ! rotBn на гранях
  !$OMP PARALLEL DEFAULT(none) Shared(nodes,edges,rotB_edges,B_nodes,nEdges,edge_vol) PRIVATE(nEdge,r1,r2,r_mid)
  !$OMP DO SCHEDULE(auto)
  do nEdge=1,nEdges
    r1 = nodes(2,edges(1,nEdge))
    r2 = nodes(2,edges(2,nEdge))
    r_mid = 0.5*(r1+r2)
    if (r_mid>0.0) then
      rotB_edges(nEdge) = (r2*B_nodes(3,edges(2,nEdge))-r1*B_nodes(3,edges(1,nEdge)))/(edge_vol(nEdge)*r_mid)
    else
      rotB_edges(nEdge) = 0.0
    end if
  end do
  !$OMP END DO nowait
  !$omp end parallel

  ! rotBxy в ячейках
  !$OMP PARALLEL DEFAULT(none) Shared(nElems,rotB,n,edge_n,rotB_edges,elem_edges,e,elems_vol,h,edge_normSign) PRIVATE(i,A,tau_r)
  !$OMP DO SCHEDULE(auto)
  do i=1,nElems
    rotB(1:2,i) = 0.0
    do j=1,3
      A = edge_normSign(j,i)*rotB_edges(elem_edges(j,i))*h(i,j)
      ! dot_product(n(:,i,j),edge_n(:,elem_edges(j,i)))*
      tau_r = A*(e(:,i,Q(Q(j)))-e(:,i,Q(j)))
      rotB(1:2,i) = rotB(1:2,i) + tau_r
    end do
    rotB(1:2,i) = rotB(1:2,i)/(6.0*elems_vol(i))
  end do
  !$OMP END DO nowait
  !$omp end parallel

  ! Bt iUW на гранях
  !$OMP PARALLEL DEFAULT(none) Shared(nEdges,elems_center,edge_elems,nodes,edges,B_hat,B_hat_angles,B_edges_hat,edge_vol) PRIVATE(i,j,B_L,B_R,x01,x02,xc,tang,rhoL,rhoR)
  !$OMP DO SCHEDULE(auto)
  do i=1,nEdges
    x01 = elems_center(:,edge_elems(1,i))
    x02 = elems_center(:,edge_elems(2,i))
    xc  = 0.5*(nodes(:,edges(1,i))+nodes(:,edges(2,i)))
    do j=1,2
      B_L(j) = B_hat(j,edge_elems(1,i)) + dot_product(B_hat_angles(:,j,edge_elems(1,i)),xc-x01)
      B_R(j) = B_hat(j,edge_elems(2,i)) + dot_product(B_hat_angles(:,j,edge_elems(2,i)),xc-x02)
    end do

    tang = nodes(:,edges(2,i))-nodes(:,edges(1,i))
    tang = tang/edge_vol(i)

    rhoL = 1.0
    rhoR = 1.0
    B_edges_hat(3,i) = dot_product(tang,B_L*sqrt(rhoR)+sqrt(rhoL)*B_R)/(sqrt(rhoL)+sqrt(rhoR))
  end do
  !$OMP END DO nowait
  !$omp end parallel

  ! rotBz в ячейках
  !$OMP PARALLEL DEFAULT(none) Shared(nElems,edge_normSign,elem_edges,edge_vol,elems_vol,rotB,B_edges_hat) PRIVATE(nElem,circ,nEdge,A)
  !$OMP DO SCHEDULE(auto)
  do nElem=1,nElems
    circ = 0.0
    do nEdge=1,3
      A = edge_normSign(nEdge,nElem)!dot_product(n(:,nElem, nEdge), edge_n(:,elem_edges(nEdge,nElem)))
      circ = circ + B_edges_hat(3,elem_edges(nEdge,nElem))*edge_vol(elem_edges(nEdge,nElem))*A
    end do
    A = 1.0/elems_vol(nElem)
    rotB(3,nElem) = A*circ
  end do
  !$OMP END DO nowait
  !$omp end parallel
end subroutine rotB_incell

  ! ======================================================
  ! Восполнение правой части газовых уравнений
subroutine magnetic_kinematics
  implicit none

  real (kind=precis), dimension(1:3) :: f_lorenz,normal
  real (kind=precis) :: pres1, energ1,flow,B_l,B_r,B_lr
  real (kind=precis), dimension(1:2) :: x01, x02, xc
  integer :: i,j
  !=====================================================
  ! поправка в гидродинамическую часть от магнитных сил

  call rotB_incell

  normal(3) = 0.0

  do nElem = 1,nElems
    f_lorenz = vect_mult(rotB(:,nElem),B_hat(:,nElem))/mag_fact
    RS(1) = 0.0
    RS(2:4) = f_lorenz
      ! костыль на фитую компоненту - разобраться в будущем
    RS(4) = 0.0
    x01 = elems_center(:,nElem)
    flow = 0.0
    do nEdge=1,3
      x02 = elems_center(:,elems_elems(nEdge,nElem))
      xc  = elem_edges_center(:,nEdge,nElem)
      B_l = B_hat(3,nElem) + dot_product(B_hat_angles(:,3,nElem),xc-x01)
      B_r = B_hat(3,elems_elems(nEdge, nElem)) + dot_product(B_hat_angles(:,3,elems_elems(nEdge, nElem)),xc-x02)
      B_lr = 0.5*(B_l+B_r)
      A = h(nElem,nEdge)*dot_product(n(:,nElem,nEdge),edge_n(:,elem_edges(nEdge,nElem)))
      flow = A*B_edges_hat(1,elem_edges(nEdge,nElem))*B_lr/(mag_fact*elems_vol(nElem))
      RS(4) = RS(4) + flow
    end do

    RS(4) = RS(4) - B_hat(2,nElem)*B_hat(3,nElem)/(elems_center(2,nElem)*mag_fact)

    RS(5) = dot_product(RS(2:4), U_til(2:4,nElem))/U_til(1,nElem)

    U_hat(:,nElem) = U_til(:,nElem)+tau*RS

    energ1 = U_hat(5,nElem)
    pres1 = (energ1-0.5*(sum(U_hat(2:4,nElem)**2))/U_hat(1,nElem))*(gmm-1.0)
    if (pres1<0.05*energ1) then
      press = press_correct(nElem)
      U_hat(5,nElem) = press/(gmm-1.0)+0.5*(sum(U_hat(2:4,nElem)**2))/U_hat(1,nElem)
    end if
  end do
  !=====================================================

!	call compute_bc_rs
  !=====================================================
end subroutine magnetic_kinematics

end module solvers