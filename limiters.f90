module limiters

  use names

implicit none

contains

  real (kind=precis) function minmod(a,b,c)
    real (kind=precis) :: a,b,c

	minmod = 0.0
	if (min(a,b,c)*max(a,b,c)>0.0) then
	  minmod = sign(1.0, a)*min(abs(a),abs(b),abs(c))
	end if
  end function minmod

  real (kind=precis) function limiter(a)
    real (kind=precis),dimension(3) :: a
	limiter = minmod(a(1),a(2),a(3))
!   limiter = superbee(a,b,c)
  end function limiter

  real (kind=precis) function det(A)
    real (kind=precis), dimension(2,2)::A

	det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
  end function det

  subroutine limit_B
	real (kind=precis), dimension(3) :: alf, bet
	real (kind=precis) :: det1
	real (kind=precis), dimension(2) :: x0,xi,xii, rs
	real (kind=precis),dimension(2,2) :: A, A1
    integer :: i,j,nElem

	do nElem = 1,nElems
      x0 = elems_center(:,nElem)
      do j=1,3
	    do i=1,3
	      xi = elems_center(:,elems_elems(i,nElem))
	      xii = elems_center(:,elems_elems(Q(i),nElem))
		  A(1,:) = xi - x0
		  A(2,:) = xii - x0
		  det1 = det(A)
		  rs(1) = B(j,elems_elems(i,nElem)) - B(j,nElem)
		  rs(2) = B(j,elems_elems(Q(i),nElem)) - B(j,nElem)
		  A1 = A
		  A1(:,1) = rs
          alf(i) = det(A1)/det1
		  A1 = A
		  A1(:,2) = rs
		  bet(i) = det(A1)/det1
	    end do
        B_angles(1,j,nElem) = limiter(alf)
        B_angles(2,j,nElem) = limiter(bet)
	  end do
	end do
	do nElem = 1, nElems_b
	  B_angles(:,:,nElem+nElems) = 0.0
	end do
  end subroutine limit_B

  subroutine limit_U
	real (kind=precis), dimension(3) :: alf, bet
	real (kind=precis) :: det1
	real (kind=precis), dimension(2) :: x0,xi,xii, rs
	real (kind=precis),dimension(2,2) :: A, A1
    integer :: i,j,nElem

	do nElem = 1,nElems
      x0 = elems_center(:,nElem)
      do j=1,5
	    do i=1,3
	      xi = elems_center(:,elems_elems(i,nElem))
	      xii = elems_center(:,elems_elems(Q(i),nElem))
		  A(1,:) = xi - x0
		  A(2,:) = xii - x0
		  det1 = det(A)
		  rs(1) = U(j,elems_elems(i,nElem)) - U(j,nElem)
		  rs(2) = U(j,elems_elems(Q(i),nElem)) - U(j,nElem)
		  A1 = A
		  A1(:,1) = rs
          alf(i) = det(A1)/det1
		  A1 = A
		  A1(:,2) = rs
		  bet(i) = det(A1)/det1
	    end do
        U_angles(1,j,nElem) = limiter(alf)
        U_angles(2,j,nElem) = limiter(bet)
	  end do
	end do
	do nElem = 1, nElems_b
	  U_angles(:,:,nElem+nElems) = 0.0
	end do
  end subroutine limit_U

subroutine limit_B_hat
  real (kind=precis), dimension(3) :: alf, bet
  real (kind=precis) :: det1
  real (kind=precis), dimension(2) :: x0,xi,xii, rs
  real (kind=precis),dimension(2,2) :: A, A1
  integer :: i,j,nElem

  do nElem = 1,nElems
    x0 = elems_center(:,nElem)
    do j=1,3
      do i=1,3
        xi = elems_center(:,elems_elems(i,nElem))
        xii = elems_center(:,elems_elems(Q(i),nElem))
        A(1,:) = xi - x0
        A(2,:) = xii - x0
        det1 = det(A)
        rs(1) = B_hat(j,elems_elems(i,nElem)) - B_hat(j,nElem)
        rs(2) = B_hat(j,elems_elems(Q(i),nElem)) - B_hat(j,nElem)
        A1 = A
        A1(:,1) = rs
        alf(i) = det(A1)/det1
        A1 = A
        A1(:,2) = rs
        bet(i) = det(A1)/det1
      end do
      B_hat_angles(1,j,nElem) = limiter(alf)
      B_hat_angles(2,j,nElem) = limiter(bet)
!        B_hat_angles(1,j,nElem) = 0.0
!        B_hat_angles(2,j,nElem) = 0.0

    end do
  end do
  do nElem = 1, nElems_b
    B_hat_angles(:,:,nElem+nElems) = 0.0
  end do
end subroutine limit_B_hat

  subroutine limit_U_hat
	real (kind=precis), dimension(3) :: alf, bet
	real (kind=precis) :: det1
	real (kind=precis), dimension(2) :: x0,xi,xii, rs
	real (kind=precis),dimension(2,2) :: A, A1
    integer :: i,j,nElem

	do nElem = 1,nElems
      x0 = elems_center(:,nElem)
      do j=1,5
	    do i=1,3
	      xi = elems_center(:,elems_elems(i,nElem))
	      xii = elems_center(:,elems_elems(Q(i),nElem))
		  A(1,:) = xi - x0
		  A(2,:) = xii - x0
		  det1 = det(A)
		  rs(1) = U_hat(j,elems_elems(i,nElem)) - U_hat(j,nElem)
		  rs(2) = U_hat(j,elems_elems(Q(i),nElem)) - U_hat(j,nElem)
		  A1 = A
		  A1(:,1) = rs
          alf(i) = det(A1)/det1
		  A1 = A
		  A1(:,2) = rs
		  bet(i) = det(A1)/det1
	    end do
        U_hat_angles(1,j,nElem) = limiter(alf)
        U_hat_angles(2,j,nElem) = limiter(bet)
	  end do
	end do
	do nElem = 1, nElems_b
	  U_hat_angles(:,:,nElem+nElems) = 0.0
	end do
  end subroutine limit_U_hat


  subroutine compute_upkoef
	real (kind=precis) :: buf1, buf2, rho
	integer :: nNode, i, nnum, nbcel, nNode1
	real (kind=precis), dimension(2) :: nor_i, lam, nor1, nor2

! up_koef
	do nElem=1, nElems
	  do nNode=1,3
		nor_i = -n(:,nElem, Q(nNode))

		!!!
!		nor_i = nodes(:,elems(nNode,nElem))-elems_center(:,nElem)
!		nor_i = nor_i/(sqrt(sum(nor_i**2)))
		!!!

		rho = U_til(1,nElem)
!		lam = (/0.0, U_til(3,nElem) /)
		lam = (/U_til(2,nElem), U_til(3,nElem) /)
		lam = lam/rho
		up_koef(1,nNode,nElem) = max(0.0,dot_product(nor_i, lam))
!		lam = (/ U_til(2,nElem), 0.0/)
		lam = (/ U_til(2,nElem), U_til(3,nElem)/)
		lam = lam/rho
		up_koef(2,nNode,nElem) = max(0.0,dot_product(nor_i, lam))
	  end do
	  if (sum(up_koef(1,:,nElem))>1E-7) up_koef(1,:,nElem) = up_koef(1,:,nElem)/sum(up_koef(1,:,nElem))
	  if (sum(up_koef(2,:,nElem))>1E-7) up_koef(2,:,nElem) = up_koef(2,:,nElem)/sum(up_koef(2,:,nElem))
	end do

	do i=1, nElems_b
	  nElem = nElems+i
	  nbcel = bc_elems(i)
	  up_koef(:,:,nElem) = 0.0

	  nor_i = n(:,nbcel,bc_elem_edge(i))

	  nnum = Q(bc_elem_edge(i))
	  nor1 = n(:,nbcel,nnum)-2.0*nor_i*dot_product(n(:,nbcel,nnum),nor_i)
	  nor1 = -nor1

	  nnum = Q(Q(bc_elem_edge(i)))
	  nor2 = n(:,nbcel,nnum)-2.0*nor_i*dot_product(n(:,nbcel,nnum),nor_i)
	  nor2 = -nor2

	  rho = U_til(1,nElem)
	  lam = (/U_til(2,nElem), U_til(3,nElem) /)
	  lam = lam/rho

	  up_koef(1,1,nElem) = max(0.0,dot_product(nor1, lam))
	  up_koef(1,2,nElem) = max(0.0,dot_product(nor2, lam))
	  lam = (/ U_til(2,nElem), U_til(3,nElem)/)
	  lam = lam/rho
	  up_koef(2,1,nElem) = max(0.0,dot_product(nor1, lam))
	  up_koef(2,2,nElem) = max(0.0,dot_product(nor2, lam))

	  if (sum(up_koef(1,:,nElem))>1E-7) up_koef(1,:,nElem) = up_koef(1,:,nElem)/sum(up_koef(1,:,nElem))
	  if (sum(up_koef(2,:,nElem))>1E-7) up_koef(2,:,nElem) = up_koef(2,:,nElem)/sum(up_koef(2,:,nElem))

	  if (bc_type(i)==10) then
	  	  nElem = bc_elems(i)
	  	  nnum = bc_elem_edge(i)
  		  do nNode1=nnum,nnum+1
  		  	nNode = nNode1
  		  	if (nNode1>3) nNode = 1
			nor_i = -n(:,nElem, Q(nNode))
			rho = U_til(1,nElem)
			lam = (/U_til(2,nElem), U_til(3,nElem) /)
			lam = lam/rho
			up_koef(1,nNode,nElem) = max(0.5*sqrt(sum(lam**2)), abs(dot_product(nor_i, lam)))
			lam = (/ U_til(2,nElem), U_til(3,nElem)/)
			lam = lam/rho
			up_koef(2,nNode,nElem) = max(0.5*sqrt(sum(lam**2)), abs(dot_product(nor_i, lam)))
		  end do

  		    nNode = Q(Q(nnum))
			nor_i = -n(:,nElem, Q(nNode))
			rho = U_til(1,nElem)
			lam = (/U_til(2,nElem), U_til(3,nElem) /)
			lam = lam/rho
			up_koef(1,nNode,nElem) = max(0.0, dot_product(nor_i, lam))
			lam = (/ U_til(2,nElem), U_til(3,nElem)/)
			lam = lam/rho
			up_koef(2,nNode,nElem) = max(0.0, dot_product(nor_i, lam))

	  	  if (sum(up_koef(1,:,nElem))>1E-7) up_koef(1,:,nElem) = up_koef(1,:,nElem)/sum(up_koef(1,:,nElem))
	  	  if (sum(up_koef(2,:,nElem))>1E-7) up_koef(2,:,nElem) = up_koef(2,:,nElem)/sum(up_koef(2,:,nElem))

		  nElem = nElems+i
		  up_koef(1,:,nElem) = 0.0
		  up_koef(2,:,nElem) = 0.0

	  end if
	end do

! up_koef_z
	do nElem=1, nElems
	  do nEdge=1,3
		nor_i = n(:,nElem, nEdge)
		rho = U_til(1,nElem)
		lam = (/U_til(2,nElem), U_til(3,nElem) /)
		lam = lam/rho
		up_koef_z(nEdge,nElem) = max(0.0,dot_product(nor_i, lam))
	  end do
	  if (sum(up_koef_z(:,nElem))>1E-7) up_koef_z(:,nElem) = up_koef_z(:,nElem)/sum(up_koef_z(:,nElem))
	end do

  end subroutine compute_upkoef

subroutine compute_upstr
  !upstr_nodes(:,:) !upstr_nodes(nElem_of_node,nNode)
  integer :: nNode,numEl,node_local
  real (kind=precis) :: buf1, buf2, rho
  integer :: i, nnum, nbcel, nNode1
  real (kind=precis), dimension(2) :: nor_i, lam, nor1, nor2


  do nNode=1,nNodes
    do nElem=1,node_elems_count(nNode)
      numEl = node_elem(nElem,nNode)
      node_local = elem_nodeN(nElem,nNode)
      if (numEl>nElems) then
        nor_i = -n(:,bc_elems(numEl-nElems), bc_elem_edge(numEl-nElems))
      else
        nor_i = -n(:,numEl, Q(node_local))
      end if
		!!!
!			nor_i = nodes(:,elems(nNode,nElem))-elems_center(:,nElem)
!			nor_i = nor_i/(sqrt(sum(nor_i**2)))
		!!!
      rho = U_til(1,numEl)
      lam = (/U_til(2,numEl), U_til(3,numEl) /)
      lam = lam/rho
      upstr_nodes(nElem,nNode) = max(0.0,dot_product(nor_i, lam))
    end do
    if (sum(upstr_nodes(:,nNode))>1E-7) upstr_nodes(:,nNode) = upstr_nodes(:,nNode)/sum(upstr_nodes(:,nNode))
  end do

! up_koef_z
  do nElem=1, nElems
    do nEdge=1,3
      nor_i = n(:,nElem, nEdge)
      rho = U_til(1,nElem)
      lam = (/U_til(2,nElem), U_til(3,nElem) /)
      lam = lam/rho
      up_koef_z(nEdge,nElem) = max(0.0,dot_product(nor_i, lam))
    end do
    if (sum(up_koef_z(:,nElem))>1E-7) up_koef_z(:,nElem) = up_koef_z(:,nElem)/sum(up_koef_z(:,nElem))
  end do
end subroutine compute_upstr

end module limiters