! ========================================
! Описатели сетки и работа с ними.
! Работа с геометрией.
! ========================================

! =====================================================================
! 1. Узлы КЭ обходятся против часовой стрелки
! 2. Ребра КЭ обходятся против часовой стрелки -- см. комментарии ниже
! 3. Внешняя часть расчетной области соотвествует КЭ с номером 0 (ноль) --
!    всюду, где это может потребоваться
! 4. Файлы с данными для этой программы -- в форматном (текстовом) виде --
!    это нужно для переносимости программы между системами и компиляторами
! 5. Для (r,z) геометрии ось r направлена горизонтально вправо,
!    ось z -- вертикально вверх
!    Для (x,y) направление осей обычное.
!    x = r, y = z
! ======================================================================

! =================================================================================
! Раскомментированны массивы, которые _уже_ используются в главной программе
! на данный момент.
! Для расчета восполнения могут понадобиться и другие.
! =================================================================================

module mesh
  use names
  use multies

  implicit none

contains

subroutine scalemesh(scale_x, scale_y)
  integer :: i
  real (kind=precis),Intent(In) :: scale_x, scale_y

  z_max = 0D0
  r_max = 0D0
  do i=1,nNodes
    nodes(1,i) = scale_x*nodes(1,i)
    nodes(2,i) = scale_y*nodes(2,i)
    if (nodes(1,i)>z_max) z_max = nodes(1,i)
    if (nodes(2,i)>r_max) r_max = nodes(2,i)
  end do
end subroutine scalemesh
  ! ================================
  ! Инициали зация всех описателей сетки
subroutine init_mesh()

  ! used:
  !	nodes, elems, elems_elems, elems_vol, bc_elems, bc_nodes,
  !	bc_type, n, h, is_bc_elem

  implicit none

  integer :: node, elem,edge
  real (kind=precis), dimension(1:2) :: ed1,ed2
  real (kind=precis) :: r1,r2,z1,z2,radc,delz, bbb

  print *
  print *, 'init_mesh()'
  print *, '==========='

  open(101, file='./data/mesh/nodes.dat')
  open(102, file='./data/mesh/elems.dat')
  open(103, file='./data/mesh/bc_elems.dat')

  read(101,*) nNodes
  print *, '  > NNODES   =', NNODES
  read(102,*) nElems
  print *, '  > NELEMS   =', NELEMS
  read(103,*) nElems_b
  print *, '  > NELEMS_B =', NELEMS_B

  allocate(nodes(1:2,1:nNodes), stat = err)
  print *, "Nodes allocate stat = ", err
  allocate(is_bc_node(1:nNodes), stat = err)
  print *, "Nodes allocate stat = ", err
  allocate(elems(1:3,1:nElems+nElems_b), stat = err)
  print *, "elems allocate stat = ", err
  allocate(elems_elems(1:3,1:nElems+nElems_b), stat = err)
  print *, "elems_elems allocate stat = ", err
  allocate(elems_vol(1:nElems+nElems_b), stat = err)
  print *, "elems_vol allocate stat = ", err
  allocate(elems_3dvol(1:nElems+nElems_b), stat = err)
  print *, "elems_3dvol allocate stat = ", err
  allocate(node_edges(1:10,1:nNodes), stat = err)
  print *, "node_edges allocate stat = ", err
  allocate(node_edges_count(1:nNodes), stat = err)
  print *, "node_edges_count allocate stat = ", err
  is_bc_node(:) = 0
  node_edges(:,:)=0
  node_edges_angles(:,:)=0.0
  node_edges_count(:)=0

  allocate(bc_elems(1:nElems_b), stat = err)
  print *, "bc_elems allocate stat = ", err
  allocate(bc_nodes(1:2,1:nElems_b), stat = err)
  print *, "bc_nodes allocate stat = ", err
  allocate(bc_type(1:nElems_b), stat = err)
  print *, "bc_type allocate stat = ", err

  allocate(elem_edges(1:3,1:nElems+nElems_b), stat = err)
  print *, "elem_edges allocate stat = ", err
  allocate(elems_center(1:2,1:nElems+nElems_b), stat = err)
  print *, "elems_center allocate stat = ", err
  allocate(elem_edges_center(1:2,1:3,1:nElems), stat = err)
  print *, "elem_edges_center allocate stat = ", err
  allocate(elems_vol1(1:nElems), stat = err)
  print *, "elems_vol1 allocate stat = ", err
  elem_edges(:,:) = 0
  elems_center(:,:) = 0.0
  elem_edges_center(:,:,:) = 0.0
  elems_vol1(:)= 0.0

  allocate(is_bc_elem(1:nElems+nElems_b), stat = err)
  print *, "is_bc_elem allocate stat = ", err
  allocate(near_ort_node(1:nElems+nElems_b), stat = err)
  print *, "near_ort_node allocate stat = ", err
  is_bc_elem(:) = 0
  near_ort_node(:) = 0

  allocate(bc_elem_edge(1:nElems_b), stat = err)
  print *, "bc_elem_edge allocate stat = ", err
  allocate(not_bc_node(1:nElems_b), stat = err)
  print *, "not_bc_node allocate stat = ", err
  allocate(bc_edges(1:nElems_b), stat = err)
  print *, "bc_edges allocate stat = ", err
  allocate(bc_edges_type(1:nElems_b), stat = err)
  print *, "bc_edges_type allocate stat = ", err
  allocate(bc_node(1:2*nElems_b), stat = err)
  print *, "bc_node allocate stat = ", err
  allocate(bc_node_type(1:2*nElems_b), stat = err)
  print *, "bc_node_type allocate stat = ", err
  allocate(period(1:nElems_b), stat = err)
  print *, "period allocate stat = ", err

  bc_elem_edge(:)= 0
  not_bc_node(:)= 0
  bc_edges(:)= 0
  bc_edges_type(:)= 0
  bc_node(:)= 0
  bc_node_type(:)=0
  period(:) = 0

  ! =======================
  ! Reading NODES.DAT (101)
  do node = 1,nNodes
    read (101,*) nodes(1,node), nodes(2,node)
    if (nodes(1,node)>z_max) z_max = nodes(1,node)
    if (nodes(2,node)>r_max) r_max = nodes(2,node)
  enddo

  ! =======================
  ! Reading ELEMS.DAT (102)
  do elem = 1,nElems
    read (102,*) elems(1,elem), elems(2,elem), elems(3,elem), &
               elems_elems(1,elem), elems_elems(2,elem), elems_elems(3,elem), &
               elems_vol(elem)
  enddo

  ! =======================
  ! Reading BC_ELEMS.DAT
  do elem = 1,nElems_b
    read (103,*) bc_nodes(1,elem), bc_nodes(2,elem), bc_elems(elem), bc_type(elem)
  enddo

  print *, bc_elems(nElems_b)
  print *, bc_nodes(:,nElems_b)
  print *, bc_type(nElems_b)

! Close files
  close(101)
  close(102)
  close(103)

!!!!вычисление площадей элементов
! ?????????????????? объемов для цилиндрии?
!!!!

  allocate(h(nElems,3), stat = err)
  print *, 'h allocate stat=', err
  allocate(n(2,nElems,3), stat = err)
  print *, 'n allocate stat=', err
  allocate(e(2,nElems,3), stat = err)
  print *, 'e allocate stat=', err

  do nElem=1, nElems
    do nEdge=1,3
      ! ===========================
      ! локальные номера вершин ребра
      NL1 = nEdge  ! это так в силу согласованности
      NL2 = Q(NL1) ! нумерации узлов и ребер
      ! ============================
      ! глобальные номера вершин ребра
      NG1 = elems(NL1,nElem)
      NG2 = elems(NL2,nElem)
      ! ===============================
      ! координаты узлов, образующих ребро
      xy1 = nodes(1:2,NG1)
      xy2 = nodes(1:2,NG2)
      ! ===============================================================
      ! вектор, соответствующий ребру, ориентированный против часовой стрелки,
      ! от вершины с меньшим номером к вершине с большим номером
      e(:,nElem,nEdge) = xy2 - xy1
      ! внешняя нормаль нормаль, длина которой равна длине ребра
      n(:,nElem,nEdge) = (/ e(2,nElem,nEdge), -e(1,nElem,nEdge) /)
      ! длина ребра
      h(nElem,nEdge) = sqrt(sum(n(:,nElem,nEdge)**2.0))
      ! единичная нормаль
      n(:,nElem,nEdge) = n(:,nElem,nEdge)/h(nElem,nEdge)
    end do
  end do

  do elem=1,nElems_b
    elems_vol(elem+nElems) = elems_vol(bc_elems(elem))
  end do


  do elem=1,nElems
    elems_3dvol(elem) = 0D0
    elems_vol1(elem) = 0.0
    do edge=1,3
      r1 = nodes(2,elems(edge,elem))
      z1 = nodes(1,elems(edge,elem))
      r2 = nodes(2,elems(Q(edge),elem))
      z2 = nodes(1,elems(Q(edge),elem))
      delz = z2-z1
      radc = r1*r1+r1*r2+r2*r2
      elems_vol1(elem) = elems_vol1(elem) - pi*delz*radc/3D0
      elems_3dvol(elem) = elems_3dvol(elem) - pi*delz*radc/3D0
    end do
  end do

  do elem=1,nElems_b
    is_bc_elem(bc_elems(elem)) = 1
  end do

  print *
  print *, 'init_mesh() DONE'
  print *, '================'

end subroutine init_mesh

subroutine init_mesh2
  implicit none

! used static
!	edges, edge_vol, edge_n, edge_elems, elem_edges
!	node_edges, bc_node, bc_edges...

  integer:: node, edge, elem
  integer:: i,j,k,nEdge_global
  real (kind=precis), dimension(1:2) :: buf1, buf2
  real (kind=precis)   :: sina, cosa, buf, lam, mju, nju, aside, bside
  integer:: sort1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *
  print *, 'init_mesh2()'
  print *, '==========='

  open(1, file='./data/mesh/edges.dat')

  read(1,*) nEdges
  print *, '  > nEdges   =', nEdges

  allocate(edges(2,nEdges), stat = err)
  print *, 'edges allocate stat=', err
  allocate(edge_elems(2,nEdges), stat = err)
  print *, 'edge_elems allocate stat=', err
  allocate(edge_n(2,nEdges), stat = err)
  print *, 'edge_n allocate stat=', err
  allocate(edge_vol(nEdges), stat = err)
  print *, 'edge_vol allocate stat=', err
  allocate(elem_edgeN(nEdges, 2), stat = err)
  print *, 'elem_edgeN allocate stat=', err
  allocate(edge_normSign(1:3,1:nElems), stat = err)
  print *, "edge_fluxSigns allocate stat = ", err
  allocate(edge_3dsurf(nEdges), stat = err)
  print *, 'edge_3dsurf allocate stat=', err


  do edge = 1,nEdges
    read (1,*) edges(1,edge), edges(2,edge), edge_elems(1,edge), edge_elems(2,edge), elem_edgeN(edge,1), elem_edgeN(edge,2), edge_vol(edge)
    edge_n(:,edge) = n(:,edge_elems(1,edge),elem_edgeN(edge,1))
    elem_edges(elem_edgeN(edge,1),edge_elems(1,edge)) = edge
    elem_edges(elem_edgeN(edge,2),edge_elems(2,edge)) = edge
  enddo
  close(1)

  do nElem = 1,nElems
    do nEdge = 1,3
      nEdge_global=elem_edges(nEdge,nElem)
      edge_normSign(nEdge,nElem) = dot_product(edge_n(:,nEdge_global),n(:,nElem,nEdge))
    enddo
  enddo

  allocate(node_elems_count(nNodes), stat = err)
  print *, 'node_elems_count allocate stat=', err
  allocate(node_elem(10,nNodes), stat = err)
  print *, 'node_elem allocate stat=', err
  allocate(elem_nodeN(10,nNodes), stat = err)
  print *, 'elem_nodeN allocate stat=', err

  node_elems_count(:) = 0
  node_edges_count(:) = 0

  do elem=1,nElems
    do edge=1,3
      node = elems(edge,elem)
      node_elems_count(node) = node_elems_count(node) + 1
      node_elem(node_elems_count(node),node) = elem
      elem_nodeN(node_elems_count(node),node) = edge
    end do
  end do

  k = 1
  do edge=1,nEdges
    buf1 = nodes(:,edges(1,edge))
    buf2 = nodes(:,edges(2,edge))
    edge_3dsurf(edge) = pi*distance(buf1,buf2)*(buf1(2)+buf2(2))
!     print *, buf1
!     print *, buf2
!     print *, edge_3dsurf(edge)
!     read *
    do j=1,2
      node = edges(j,edge)
      node_edges_count(node) = node_edges_count(node) + 1
      node_edges(node_edges_count(node),node) = edge
      if (node_edges_count(node)>k) k = node_edges_count(node)
    end do
  end do
  print *, 'node_edges_max = ', k

  bc_node(1) = bc_nodes(1,1)
  is_bc_node(bc_nodes(1,1)) = 1
  bc_node(2) = bc_nodes(2,1)
  is_bc_node(bc_nodes(2,1)) = 2
  i = 3
  do elem = 2, nElems_b
    do j = 1,2
      k = 1
      do while (bc_node(k)/=bc_nodes(j,elem) .and. k<i)
        k = k + 1
      end do
      if (k==i) then
        bc_node(i) = bc_nodes(j,elem)
        is_bc_node(bc_nodes(j,elem)) = i
        i = i + 1
      end if
    end do
  end do
  nNodes_b = i - 1

  print *, '  > nNodes_b =', nNodes_b
!!!! неверно!
  do elem=1,nNodes_b
    bc_node_type(elem) = 2
  end do

  ! ребра на границе
  do edge=1,nEdges
    if (edge_elems(2,edge)>nElems) then
      elem = edge_elems(2,edge)-nElems
      bc_edges(elem) = edge
      bc_edges_type(elem) = bc_type(elem)
      bc_elem_edge(elem) = elem_edgeN(edge,1)

      node = edges(1,edge)
      node_elems_count(node) = node_elems_count(node) + 1
      node_elem(node_elems_count(node),node) = edge_elems(2,edge)
      elem_nodeN(node_elems_count(node), node) = 1

      node = edges(2,edge)
      node_elems_count(node) = node_elems_count(node) + 1
      node_elem(node_elems_count(node),node) = edge_elems(2,edge)
      elem_nodeN(node_elems_count(node), node) = 2
    end if
  end do

  print *, '  > nNodes_b >'

  allocate(bc_node_neib(1:nNodes_b,1:2,1:3), stat = err)
  print *, "bc_node_neib allocate stat = ", err

  do i=1,nNodes_b
    bc_node_neib(i,:,:) = 0
  end do

  do elem=1,nElems_b
    k = 1
    do while (bc_nodes(1,elem)==elems(k,bc_elems(elem)) .or. bc_nodes(2,elem)==elems(k,bc_elems(elem)))
      k = k + 1
    end do
    not_bc_node(elem) = k

    do j=1,2
      k=1
      if (bc_node_neib(is_bc_node(bc_nodes(j,elem)),k,1)/=0) k=2

      bc_node_neib(is_bc_node(bc_nodes(j,elem)),k,1) = bc_elems(elem)
      bc_node_neib(is_bc_node(bc_nodes(j,elem)),k,2) = Q(not_bc_node(elem))
      bc_node_neib(is_bc_node(bc_nodes(j,elem)),k,3) = elems(Q(not_bc_node(elem)),bc_elems(elem))
      if (elems(Q(not_bc_node(elem)),bc_elems(elem))==bc_nodes(j,elem)) then
	bc_node_neib(is_bc_node(bc_nodes(j,elem)),k,3) = elems(Q(Q(not_bc_node(elem))),bc_elems(elem))
      end if
    end do
!      bc_node_neib(bc_nodes(2,elem),:,:) =
  end do

  print *
  print *, 'init_mesh2() DONE'
  print *, '================'
end subroutine init_mesh2

subroutine init_mesh3
  ! вычисление центров ячеек, центров ребер
  integer:: node, edge, elem
  real (kind=precis), dimension(2):: norm, try,try1
  real (kind=precis) :: sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! used
!	elems_center
!   elem_edges_center
  print *
  print *, 'init_mesh3()'
  print *, '==========='

  do elem=1,nElems
    elems_center(:,elem) = (/0.0,0.0/)
    do node=1,3
      elems_center(:,elem) = elems_center(:,elem) + (1.0/3.0)*nodes(:,elems(node,elem))
    end do
    do edge=1,3
      elem_edges_center(:,edge,elem) = 0.5*(nodes(:,elems(edge,elem))+nodes(:,elems(Q(edge),elem)))
    end do
  end do
  do elem=1,nElems_b
    try = elems_center(:,bc_elems(elem)) - nodes(:,bc_nodes(1,elem))
    norm = nodes(:,bc_nodes(1,elem)) - nodes(:,bc_nodes(2,elem))
    norm = (/norm(2),-norm(1)/)
    sh = sqrt(norm(1)**2+norm(2)**2)
    norm =  norm/sh
    elems_center(:,nElems+elem) = elems_center(:,bc_elems(elem)) - 2.0*norm*dot_product(norm,try)
  end do

  print *
  print *, 'init_mesh3() DONE'
  print *, '================'
end subroutine init_mesh3

subroutine init_mesh4
  ! near_ort_node
  integer:: node, edge, elem,nn, i41, i42
  real (kind=precis), dimension(2):: norm, try
  real (kind=precis), dimension(3):: try1
  real (kind=precis) :: sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *
  print *, 'init_mesh4()'
  print *, '==========='

  do elem=1,nElems_b
    norm = nodes(:,bc_nodes(1,elem))-nodes(:,bc_nodes(2,elem))
    try = (/-norm(2),norm(1)/)
    norm = try
    do node=1,2
      nn = bc_nodes(node,elem)
      do i41=1,node_elems_count(nn)-2
        do i42=1,3
          if (elems(i42,node_elem(i41,nn))/=nn) then
            try = nodes(:,nn)-nodes(:,elems(i42,node_elem(i41,nn)))
            if (abs(vect_mult2d(try,norm))<1E-7) near_ort_node(nn) = elems(i42,node_elem(i41,nn))
          end if
        end do
      end do
    end do
  end do

  print *
  print *, 'init_mesh4() DONE'
  print *, '================'
end subroutine init_mesh4

subroutine init_outbound
  integer:: node, edge, elem,nn, i41, i42
  real (kind=precis), dimension(2):: norm, try
  real (kind=precis), dimension(3):: try1
  real (kind=precis) :: sh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *
  print *, 'init_outbound()'
  print *, '==========='

  allocate(outb_nums(nElems_b), stat = err)
  print *, 'outb_nums allocate stat=', err
  outb_nums = 0
  outb_N = 0
  do elem=1,nElems_b
    i41=bc_elems(elem)
    if (bc_type(elem)==3) then
      outb_N = outb_N + 1
      outb_nums(outb_N) = i41
    end if
  end do

  do i41=1,outb_N
    do i42=outb_N,i41+1,-1
      if(elems_center(2, outb_nums(i41))>elems_center(2, outb_nums(i42))) then
        nn = outb_nums(i41)
        outb_nums(i41) = outb_nums(i42)
        outb_nums(i42) = nn
      end if
    end do
  end do

!   print *, outb_N
!   do i41=1,outb_N
!     print *, outb_nums(i41), elems_center(2, outb_nums(i41))
!   end do
!   read *

  print *
  print *, 'init_outbound() DONE'
  print *, '================'
end subroutine init_outbound

subroutine init_mesh6
! четвертый тип ГУ
  real (kind=precis):: zet,rad
  real (kind=precis), dimension(2) :: norm1,norm2
  integer :: i

  print *
  print *, 'init_mesh6()'
  print *, '==========='

  norm1 = (/0.0,1.0/)
  do nElem=1,nElems_b
    i = bc_elems(nElem)
    norm2 = n(:,i,bc_elem_edge(nElem))
    if (dot_product(norm1,norm2)>1E-5) bc_type(nElem) = 4
  end do
  print *
  print *, 'init_mesh6() DONE'
  print *, '================'
end subroutine init_mesh6

subroutine init_mesh7
! четвертый тип ГУ
  real (kind=precis):: zet,rad
  real (kind=precis), dimension(2) :: norm1,norm2, norm3, norm4, norm5
  integer :: i

  print *
  print *, 'init_mesh7()'
  print *, '==========='

  norm1 = (/0.0,1.0/)
  norm2 = (/0.0,-1.0/)
  norm3 = (/1.0,0.0/)
  norm4 = (/-1.0,0.0/)
  do nElem=1,nElems_b
    i = bc_elems(nElem)
!	  bc_type(nElem) = 2
!	  if (elems_center(2,i)>1.0) bc_type(nElem) = 4
    bc_type(nElem) = 5
    norm5 = n(:,i,bc_elem_edge(nElem))
    if (dot_product(norm1,norm5)>1.0-1E-5) bc_type(nElem) = 4
    if (dot_product(norm2,norm5)>1.0-1E-5) bc_type(nElem) = 2
    if (dot_product(norm3,norm5)>1.0-1E-5) bc_type(nElem) = 3
    if (dot_product(norm4,norm5)>1.0-1E-5) bc_type(nElem) = 5
    if (dot_product(norm4,norm5)>1.0-1E-5 .and. elems_center(2,i)<r_disk) bc_type(nElem) = 1
  end do
  print *
  print *, 'init_mesh7() DONE'
  print *, '================'
end subroutine init_mesh7

subroutine init_mesh_test
! четвертый тип ГУ
  real (kind=precis):: zet,rad
  real (kind=precis), dimension(2) :: norm2, norm5
  integer :: i

    print *
    print *, 'init_mesh8()'
    print *, '==========='

	norm2 = (/0.0,-1.0/)
	do nElem=1,nElems_b
	  norm5 = n(:,i,bc_elem_edge(nElem))
	  if (dot_product(norm2,norm5)>1.0-1E-5) then
	  		bc_type(nElem) = 3
	  	else
	  		bc_type(nElem) = 2
	  end if
	end do
    print *
    print *, 'init_mesh8() DONE'
    print *, '================'
  end subroutine init_mesh_test

subroutine hmin_compute
  integer :: nEdge
  real (kind=precis) :: sqrt3

  print *
  print *, 'hmin_compute()'
  print *, '==========='

  hmin = 1.0
  sqrt3 = sqrt(3.0)
  do nEdge=1,nEdges
    if (edge_vol(nEdge)/sqrt3 < hmin) then
      hmin = edge_vol(nEdge)/sqrt3
    endif
  end do

  print *, 'hmin = ', hmin
  print *
  print *, 'hmin_compute() DONE'
  print *, '================'
end subroutine hmin_compute

end module mesh
