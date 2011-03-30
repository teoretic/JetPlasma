module radiationTR

use names
use multies
use for77test
use mesh

implicit none

logical :: useFiles = .true.
logical :: debug = .true.
logical :: useBoundary = .false.
logical :: markers = .true.
integer :: testtask = 3

integer :: nDir

integer, parameter :: nAngular = 24! четное!!!
integer, dimension(1:nAngular/2) :: nAnglgr_phi
real (kind=precis), dimension(-nAngular/2:nAngular/2) :: gamma_ints
real (kind=precis), dimension(1:nAngular/2) :: delta_gamma
real (kind=precis), dimension(1:nAngular/2,0:2*nAngular) :: phi_ints
real (kind=precis), dimension(1:nAngular/2) :: delta_phis

integer :: nDirects = nAngular*(2+nAngular)
real(kind=precis), allocatable :: directnodes(:,:)
integer :: nDirectsUp

integer , allocatable :: directelems(:,:)
integer , allocatable :: directelems_buf(:,:)

integer, allocatable :: NS(:,:)
real (kind=precis), allocatable :: directns(:,:)
real (kind=precis), allocatable :: directnsND(:,:)

integer :: ns_iter
real (kind=precis) :: bodyAngle
real (kind=precis), allocatable :: bodyAngleTR(:)
integer , allocatable :: minGmmN(:)

real :: mid_gmm, mid_gmmp, mid_phi

real (kind=precis), allocatable :: intensityCell(:)
real (kind=precis), allocatable :: intensityNodes(:)
real (kind=precis), allocatable :: scatsource(:,:)
integer :: niterDG

real (kind=precis), allocatable :: waysLen(:,:), waysNodeDirPhi(:), waysSrcNode(:,:)
integer, allocatable :: waysElems(:,:), waysElemsCount(:)
integer :: waysNDirs

real (kind=precis), allocatable :: intensityTR(:,:)
integer, allocatable :: is_source_node(:,:) ! is_source_node(1:nNodes_b, 1:nDirects)

real (kind=precis) :: hmm = 50_precis
real (kind=precis) :: epsilTR = 1D-7

integer :: nDirectNodes

  contains

!!!!!!!!! задаются группы по направлениям и частотам ===============================================================
subroutine InitGroupsSn

  integer :: node

  allocate(directns(1:nDirects,1:3), stat = err)
  allocate(NS(1:nDirects,2), stat = err)

  do i=1,nAngular/2
    nAnglgr_phi(i) = 2*(nAngular - 2*i + 2)
  end do

  do i=1,nAngular/2
    delta_phis(i)=2.*Pi/nAnglgr_phi(i)
    do j=0,nAnglgr_phi(i)
      phi_ints(i,j) = delta_phis(i)*j !+0.25*delta_phis(i)*(1.+(-1.)**j)
    end do
  end do

  t92_nonz = 3*nAngular/2-2
  t92_neq = nAngular/2
  t92_nvar = nAngular/2

  t92_nac = 4*t92_nonz
  t92_nrend = t92_nac
  t92_mode = (/0,1,0,0/)
  t92_mxp = 10
  t92_ipd = 10
  t92_maxit = 16
  t92_stab = 0.5
  t92_bar = .15E-4
  t92_stpbar = 2.
  t92_epsin = 4E-7

  allocate(t92_row(t92_neq+1), stat = err)
  allocate(t92_h(9,-1:t92_neq), stat = err)
  allocate(t92_vt(t92_neq), stat = err)
  allocate(t92_c0(t92_nonz), stat = err)
  allocate(t92_c(0:t92_nac), stat = err)
  allocate(t92_r(0:t92_nrend), stat = err)
  allocate(t92_a0(t92_nonz), stat = err)
  allocate(t92_b(t92_neq), stat = err)
  allocate(t92_x(t92_nvar), stat = err)
  allocate(t92_d(t92_nvar), stat = err)
  allocate(t92_diag(t92_nvar), stat = err)
  allocate(t92_rr(t92_neq), stat = err)
  allocate(t92_v(t92_neq), stat = err)
  allocate(t92_cna0(t92_neq), stat = err)
  allocate(t92_barm(t92_neq), stat = err)
  allocate(t92_a(t92_nac), stat = err)

  do i=1,t92_nac
    t92_a(i)=0.
    t92_c(i)=0
  end do
  do i=1,t92_nrend
    t92_r(i)=0
  end do
  do i=1,t92_neq
    t92_rr(i)=0.
  end do
  do i=1,t92_nvar
    t92_x(i)=0.
    t92_b(i)=0.
  end do

  do i=1,nAngular/2-1
    t92_a0((i-1)*2+1)=delta_phis(i)
    t92_a0((i-1)*2+2)=-delta_phis(i+1)
    t92_c0((i-1)*2+1)=i
    t92_c0((i-1)*2+2)=i+1
    t92_b(i) = 0.0
    t92_row(i)=1+2*(i-1)
  end do

  do i=1,nAngular/2
    t92_a0(2*(nAngular/2-1)+i)=1.
    t92_c0(2*(nAngular/2-1)+i)=i
  end do
  t92_b(t92_neq) = 1.
  t92_row(t92_neq)=1+2*(nAngular/2-1)
  t92_row(t92_neq+1)=t92_nonz+1

  call t92(t92_mode, t92_a0,t92_c0,t92_row, t92_nvar,t92_neq,t92_nonz, t92_b, t92_x,t92_d,t92_rr,t92_v,t92_vt,t92_cna0,t92_barm,t92_mxp,t92_ipd,t92_stab,t92_bar,t92_stpbar,t92_epsin,t92_epsout,t92_maxit,t92_iter,t92_a,t92_nac,t92_c,t92_r,t92_nrend,t92_h(1,-1),t92_diag,t92_kstep, t92_merr)

  gamma_ints(0) = 0.
  do i=1,nAngular/2
    delta_gamma(i) = t92_x(i)
    gamma_ints(i) = gamma_ints(i-1)+t92_x(i)
    gamma_ints(-i) = -gamma_ints(i)
  end do

  bodyAngle = delta_gamma(1)*delta_phis(1)

  ns_iter=0
  do i=1,nAngular/2
    do j=1,nAnglgr_phi(i)
      ns_iter = ns_iter + 1
      NS(ns_iter,1) = i
      NS(ns_iter,2) = j
      NS(ns_iter+nAngular*(2+nAngular)/2,1) = -i
      NS(ns_iter+nAngular*(2+nAngular)/2,2) = j
    end do
  end do

  print *, "nDirects =",nDirects
  do i=1,nDirects
    mid_gmm = 0.5*(signum(NS(i,1))*gamma_ints(Abs(NS(i,1))-1) + signum(NS(i,1))*gamma_ints(Abs(NS(i,1))))
    mid_gmmp = sqrt(1.-mid_gmm**2)
    mid_phi = 0.5*(phi_ints(Abs(NS(i,1)),NS(i,2)-1) + phi_ints(Abs(NS(i,1)),NS(i,2)))
    directns(i,:) = (/ mid_gmmp*Cos(mid_phi), mid_gmmp*Sin(mid_phi), mid_gmm/)
  end do

  deallocate(t92_row,t92_h,t92_vt,t92_c0,t92_c,t92_r,t92_a0,t92_b,t92_x,t92_d,t92_diag,t92_rr,t92_v,t92_cna0,t92_barm,t92_a)

end subroutine InitGroupsSn

subroutine InitGroupsTri

  CHARACTER(LEN = max_filename_length):: ffname

  real :: foo
  integer :: node, elem, nDir
  real(kind=precis), dimension(1:2) :: bufcent
  real(kind=precis), dimension(1:3,1:2) :: xp
  real(kind=precis), dimension(1:3) :: buf_vec, xp0, xdir
  real(kind=precis), dimension(1:3,1:3) :: drot
  real(kind=precis) :: bbb, buf_real, cang, sang
  integer :: itr, jtr, ktr, buf_int,itrac

  integer(kind=2), allocatable :: buf_nums(:)
  real(kind=precis) :: diskDivVol
  logical :: flag

  open(11, file='./data/mesh/sphere_elems.dat')
  open(12, file='./data/mesh/sphere_nodes.dat')

  read(11,*) nDirects
  print *, '  > nDirects   =', nDirects

  read(12,*) nDirectNodes
  print *, '  > nDirectNodes =', nDirectNodes

  allocate(directnodes(1:2,1:nDirectNodes+1500), stat = err)
  print *, "Allocate stat = ", err
  allocate(directelems(1:3,1:nDirects+1500), stat = err)
  print *, "Allocate stat = ", err
  allocate(directelems_buf(1:3,1:nDirects+1500), stat = err)
  print *, "Allocate stat = ", err

  allocate(directns(1:nDirects+1500,1:3), stat = err)
  print *, "Allocate stat = ", err
  allocate(bodyAngleTR(1:nDirects+1500), stat = err)
  print *, "Allocate stat = ", err
  allocate(directnsND(1:nDirectNodes+1500,1:3), stat = err)
  print *, "Allocate stat = ", err

  bodyAngle = -1.0
  do nDir = 1, nDirects
    read (11,*) directelems_buf(1,nDir), directelems_buf(2,nDir), directelems_buf(3,nDir), foo, foo, foo, bodyAngleTR(nDir)
!    if (bodyAngleTR(nDir)>bodyAngle) bodyAngle = bodyAngleTR(nDir)
  enddo

  do node = 1,nDirectNodes
    read (12,*) directnodes(1,node), directnodes(2,node)
    mid_gmm = Cos(directnodes(2,node))
    mid_gmmp = sqrt(1.-mid_gmm**2)
    mid_phi = directnodes(1,node)
    directnsND(node,:) = (/ mid_gmmp*Cos(mid_phi), mid_gmm, -mid_gmmp*Sin(mid_phi)/)
!    directnsND(node,:) = (/ mid_gmmp*Cos(mid_phi), mid_gmmp*Sin(mid_phi), mid_gmm/)
!    directnsND(node,:) = (/ mid_gmm, mid_gmmp*Sin(mid_phi), -mid_gmmp*Cos(mid_phi)/)
  enddo

  close(11)
  close(12)

  buf_int = 0
  do nDir = 1,nDirects
    buf_vec(1) = distance(directnsND(directelems_buf(1,nDir),:),directnsND(directelems_buf(2,nDir),:))
    buf_vec(2) = distance(directnsND(directelems_buf(3,nDir),:),directnsND(directelems_buf(2,nDir),:))
    buf_vec(3) = distance(directnsND(directelems_buf(1,nDir),:),directnsND(directelems_buf(3,nDir),:))
    if (buf_vec(1)>epsilTR .and. buf_vec(2)>epsilTR .and. buf_vec(3)>epsilTR) then
      buf_int = buf_int + 1
      directelems(:,buf_int) = directelems_buf(:,nDir)
    end if
  end do
  nDirects = buf_int

  do nDir=1,nDirects
    buf_vec = (/0D0, 0D0, 0D0/)
    do node=1,3
      buf_vec = buf_vec + (1.0/3.0)*directnsND(directelems(node,nDir),:)
    end do
    directns(nDir,:) = buf_vec
  end do

  buf_real = 0D0
  buf_vec = (/0D0, 0D0, 0D0/)
  bodyAngle = 0D0
  do nDir=1,nDirects
    bodyAngleTR(nDir) = bodyAng(buf_vec, directnsND(directelems(1,nDir),:),directnsND(directelems(2,nDir),:),directnsND(directelems(3,nDir),:))
    buf_real = buf_real + bodyAngleTR(nDir)
    if (bodyAngleTR(nDir)>bodyAngle) bodyAngle = bodyAngleTR(nDir)
  end do

  print *, buf_real

  do itr=1,nDirects
    do jtr=nDirects,itr+1,-1
      if(directns(itr,3)>directns(jtr,3)) then
        buf_vec = directns(itr,:)
        directns(itr,:) = directns(jtr,:)
        directns(jtr,:) = buf_vec
        do ktr=1,3
          buf_int = directelems(ktr,itr)
          directelems(ktr,itr) = directelems(ktr,jtr)
          directelems(ktr,jtr) = buf_int
        end do
        buf_real = bodyAngleTR(itr)
        bodyAngleTR(itr) = bodyAngleTR(jtr)
        bodyAngleTR(jtr) = buf_real
      end if
    end do
  end do

  ffname = "./data/anglegrid.vtk"

  open (unit=1,file=ffname)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', nDirectNodes, ' float'

  do i=1,nDirectNodes
    write (1,*) directnodes(1,i), directnodes(2,i), 0
  end do

  write (1,"(A,I7,A,I7)") 'CELLS ', nDirects, ' ', nDirects+3*nDirects

  do i=1,nDirects
    write (1,*) 3, directelems(1,i)-1, directelems(2,i)-1, directelems(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', nDirects
  do i=1,nDirects
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'CELL_DATA ', nDirects

  write (1,"(A)") 'SCALARS foo float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nDirects
    write(1,"(I7)") i
  end do

  close (unit=1,status='Keep')
  write (*,*) 'AngleGrid file has been created'
  write (*,*) '************************'

!read *
end subroutine InitGroupsTri

subroutine InitGroupsTriRaz

  CHARACTER(LEN = max_filename_length):: ffname

  real(kind=precis), allocatable :: directnodes(:,:)
  integer (kind=precis), allocatable :: directelems(:,:), direct_bc_nodes(:,:)
  real :: foo
  integer :: nDirectNodes, node, elem, nDir, nDirects_b
  real(kind=precis), dimension(1:2) :: bufcent
  real(kind=precis) :: bbb

  open(1, file='./data/mesh/sphere_elems.dat')

  read(1,*) nDirects
  print *, '  > nDirects   =', nDirects

  allocate(directns(1:nDirects,1:3), stat = err)
  allocate(bodyAngleTR(1:nDirects), stat = err)
  allocate(directelems(1:3,1:nDirects), stat = err)
!  bbb = 0.0_precis
  do nDir = 1, nDirects
    read (1,*) directelems(1,nDir), directelems(2,nDir), directelems(3,nDir), foo, foo, foo, foo
    bbb = bbb + bodyAngleTR(nDir)
  enddo

!  print *, bbb
!  read *

  close(1)

  open(1, file='./data/mesh/sphere_nodes.dat')

  read(1,*) nDirectNodes
  print *, '  > nDirectNodes =', nDirectNodes
  allocate(directnodes(1:2,1:nDirectNodes), stat = err)

  do node = 1,nDirectNodes
    read (1,*) directnodes(1,node), directnodes(2,node)
  enddo

!bodyAngleTR(nDir)

  close(1)

  open(1, file='./data/mesh/sphere_bc_elems.dat')

  read(1,*) nDirects_b
  print *, '  > nDirects_b =', nDirects_b

  allocate(direct_bc_nodes(1:2,nDirects_b), stat = err)

  do nDir = 1,nDirects_b
    read (1,*) direct_bc_nodes(1,nDir), direct_bc_nodes(2,nDir), foo, foo
    do node=1,2
!      print *, "ooooooo"
!      print *, nDir, direct_bc_nodes(node,nDir)
!      print *, directnodes(:,direct_bc_nodes(node,nDir))

      if (abs(directnodes(1,direct_bc_nodes(node,nDir)))>1E-9) then
        bbb = directnodes(2,direct_bc_nodes(node,nDir))/directnodes(1,direct_bc_nodes(node,nDir))
        directnodes(1,direct_bc_nodes(node,nDir)) = signumr(directnodes(1,direct_bc_nodes(node,nDir)))*pi/sqrt(1.+(bbb*Pi)**2)
        directnodes(2,direct_bc_nodes(node,nDir)) = directnodes(1,direct_bc_nodes(node,nDir)) * bbb

      end if

!      print *, directnodes(:,direct_bc_nodes(node,nDir))
!      read *
    end do
  enddo

  do nDir=1,nDirectNodes
        bbb = sqrt(1.-directnodes(2,nDir)**2)
        if (abs(bbb)>1E-9) directnodes(1,nDir) = directnodes(1,nDir)/bbb
  end do

  close(1)

bbb = 0.0
  do nDir=1,nDirects
    bufcent = (/0.0_precis,0.0_precis/)
    do node=1,3
      bufcent = bufcent + (1.0/3.0)*directnodes(:,directelems(node,nDir))
    end do
    mid_gmm = bufcent(2)
    mid_gmmp = sqrt(1.-mid_gmm**2)
    mid_phi = bufcent(1)
    directns(nDir,:) = (/ mid_gmmp*Cos(mid_phi), mid_gmmp*Sin(mid_phi), mid_gmm/)
    bodyAngleTR(nDir) = TriangleVol(directnodes(:,directelems(2,nDir))-directnodes(:,directelems(1,nDir)), directnodes(:,directelems(3,nDir))-directnodes(:,directelems(1,nDir)))
    bbb = bodyAngleTR(nDir)+bbb
!    print *, bodyAngleTR(nDir)
!    read *
  end do
!print *, bbb, 4.*pi
!read *
  ffname = "./data/anglegrid.vtk"

  open (unit=1,file=ffname)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', nDirectNodes, ' float'

  do i=1,nDirectNodes
    write (1,*) directnodes(1,i), directnodes(2,i), 0
  end do

  write (1,"(A,I7,A,I7)") 'CELLS ', nDirects, ' ', nDirects+3*nDirects

  do i=1,nDirects
    write (1,*) 3, directelems(1,i)-1, directelems(2,i)-1, directelems(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', nDirects
  do i=1,nDirects
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'CELL_DATA ', nDirects

  write (1,"(A)") 'SCALARS foo float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nDirects
    write(1,"(I7)") i
  end do

  close (unit=1,status='Keep')
  write (*,*) 'AngleGrid file has been created'
  write (*,*) '************************'

!  read *

end subroutine InitGroupsTriRaz


subroutine IntensAllocate

  allocate(is_source_node(1:nNodes_b, 1:nDirects), stat=err)
  print *, 'is_source_node allocate stat=', err
  print *, 'nNodes_b = ', nNodes_b

  call source_node_axial_comp

  allocate(intensityCell(nElems), stat = err)
  print *, 'full intensity allocate stat=', err

  allocate(intens_p(nElems), stat = err)
  print *, 'intens_p allocate stat=', err
  allocate(scatsource(nElems,nDirects), stat = err)
  print *, 'full scatsource with directions allocate stat=', err
  allocate(absorp_koef(nElems), stat = err)
  print *, 'absorp_koef allocate stat=', err
  allocate(disper_koef(nElems), stat = err)
  print *, 'disper_koef allocate stat=', err

end subroutine IntensAllocate

subroutine source_node_comp
  real (kind=precis), dimension(1:2,1:2) :: bufn
  real (kind=precis), dimension(1:2) :: central = (/0.5, 0.5/)
  integer :: node

  do node=1,nNodes_b

    j = 0
    do nElem=1,node_elems_count(bc_node(node))
      if (node_elem(nElem,bc_node(node))<nElems+1) then
        do i=1,3
          if (elems(i,node_elem(nElem,bc_node(node)))/=bc_node(node) .and. is_bc_node(elems(i,node_elem(nElem,bc_node(node))))>0) then
            j = j+1
            bufn(j,:) = nodes(:,elems(i,node_elem(nElem,bc_node(node))))
            bufn(j,:) = nodes(:,bc_node(node)) - bufn(j,:)
            bufn(j,:) = (/ bufn(j,2), -bufn(j,1) /)
            bufn(j,:) = - signumr(dot_product(nodes(:,bc_node(node)) - central, bufn(j,:)))*bufn(j,:)
          endif
        end do
      endif
    end do

    do nDir=1,nDirects
      is_source_node(node,nDir) = 0
      if (dot_product(bufn(1,:), directns(nDir,1:2))>0. .and. dot_product(bufn(2,:), directns(nDir,1:2))>0.) then
        is_source_node(node,nDir) = 1
      end if
      if (dot_product(bufn(1,:), directns(nDir,1:2))*dot_product(bufn(2,:), directns(nDir,1:2))<1E-8) then
        is_source_node(node,nDir) = -1
      end if
      if (nDir==595 .and. bc_node(node) == 1 .and. debug) then
        print *, ":::::::"
        print *, nodes(:,bc_node(node))
        print *, bufn(1,:)
        print *, bufn(2,:)
        print *, directns(nDir,:)
        print *, dot_product(bufn(1,:), directns(nDir,1:2)), dot_product(bufn(2,:), directns(nDir,1:2))
        print *, is_source_node(node,nDir)
        read *
      end if

    end do
  end do

end subroutine source_node_comp

subroutine source_node_axial_comp
  real (kind=precis), dimension(1:2,1:2) :: bufn
  real (kind=precis), dimension(1:2) :: snode, central = (/0.5, 0.5/)
  real (kind=precis), dimension(1:3) :: dirt
  integer :: node

  do node=1,nNodes_b
    snode = nodes(:,bc_node(node))

    if (snode(2)<epsilTR .and. (snode(1)>epsilTR .and. snode(1)<z_max-epsilTR)) then
      do nDir=1,nDirects
        is_source_node(node,nDir) = 0
      end do
    else
      do nDir = 1, nDirects
        is_source_node(node,nDir) = 0
        dirt = directns(nDir,:)
        if (snode(1)<epsilTR .and. dirt(3)>epsilTR) is_source_node(node,nDir) = 1
        if (snode(1)>z_max-epsilTR .and. dirt(3)<-epsilTR) is_source_node(node,nDir) = 1
        if (snode(2)>r_max-epsilTR .and. dirt(1)<-epsilTR) is_source_node(node,nDir) = 1

        if ((snode(2)>r_max-epsilTR .and. snode(1)<epsilTR) .or. (snode(2)>r_max-epsilTR .and. snode(1)>z_max-epsilTR)) then
          is_source_node(node,nDir) = 0
          if (snode(1)>z_max-epsilTR .and. (dirt(1)<-epsilTR .and. dirt(3)<-epsilTR)) is_source_node(node,nDir) = 1
          if (snode(1)<epsilTR .and. (dirt(1)<-epsilTR .and. dirt(3)>epsilTR)) is_source_node(node,nDir) = 1
          if (snode(1)>z_max-epsilTR .and. (dirt(1)*dirt(3)<-epsilTR)) is_source_node(node,nDir) = -1
          if (snode(1)<epsilTR .and. (dirt(1)*dirt(3)>epsilTR)) is_source_node(node,nDir) = -1
        end if
      end do
    end if
  end do
end subroutine source_node_axial_comp

function absorbtion(Uc) result (abcof)
  real (kind=precis), dimension(1:5) :: Uc
  real (kind=precis) :: abcof
  real (kind=precis) :: densy, temper, pressur

  densy = Uc(1)
  pressur = (gmm-1.)*(Uc(5)-0.5*sum(Uc(2:4)**2)/Uc(1))
  temper = pressur/densy
  abcof = (densy**2)/(temper**3.5)
end function absorbtion

subroutine tracing(directt,noden,elemn)

  real (kind=precis), dimension(1:2), Intent(In) :: directt
  real (kind=precis), dimension(1:2) :: direct
  integer, Intent(In) :: noden, elemn
  real (kind=precis), dimension(1:2,1:2) :: trans, trans1, trans2
  real (kind=precis), dimension(1:2) :: centb, normd, bas1, bas2, basn1,basn2,nod, nodp, nodn
  integer :: nElem_next, nElem_prev, nribp, stepi, dotn1,dotn2, min_dist_n,itrac,jtrac,nElemtrac
  real (kind=precis) :: collecter, bufr, bufi, min_dist
  logical :: cond1, cond2, cond3, flag_br

  direct = directt/norm2(directt)

if (elemn==251 .and. noden == 14 .and. debug)  print *, "{{{{{{{{"
if (elemn==251 .and. noden == 14 .and. debug)  print *, nodes(:,noden)
if (elemn==251 .and. noden == 14 .and. debug)  print *, direct
if (elemn==251 .and. noden == 14 .and. debug) print *, is_source_node(is_bc_node(noden),elemn)
if (elemn==251 .and. noden == 14 .and. debug)  read*

  normd = (/direct(2), -direct(1)/)
  itrac=1
  nElemtrac = node_elem(itrac,noden)
  dotn1 = elems(Q(elem_nodeN(itrac,noden)),nElemtrac)
  dotn2 = elems(Q(Q(elem_nodeN(itrac,noden))),nElemtrac)
  cond1 = dot_product(nodes(:,dotn1)-nodes(:,noden),normd)*dot_product(nodes(:,dotn2)-nodes(:,noden),normd)>1E-14
  cond2 = (dot_product(nodes(:,dotn1)-nodes(:,noden),direct)>0. .and. dot_product(nodes(:,dotn2)-nodes(:,noden),direct)>0.)
  do while ((cond1 .or. cond2 .or. nElemtrac>nElems) .and.  itrac<=node_elems_count(noden))
    itrac = itrac + 1
    nElemtrac = node_elem(itrac,noden)
    dotn1 = elems(Q(elem_nodeN(itrac,noden)),nElemtrac)
    dotn2 = elems(Q(Q(elem_nodeN(itrac,noden))),nElemtrac)

    cond1 = dot_product(nodes(:,dotn1)-nodes(:,noden),normd)*dot_product(nodes(:,dotn2)-nodes(:,noden),normd)>1E-14
    cond2 = (dot_product(nodes(:,dotn1)-nodes(:,noden),direct)>0. .and. dot_product(nodes(:,dotn2)-nodes(:,noden),direct)>0.)

if (elemn==251 .and. noden == 14 .and. debug)    print *, itrac
  end do
  nElemtrac = node_elem(itrac,noden)
  nElem_next = elems_elems(Q(elem_nodeN(itrac,noden)),nElemtrac)
  nod = nodes(:,noden)
  basn1 = nodes(:,elems(Q(elem_nodeN(itrac,noden)),nElemtrac))
  basn2 = nodes(:,elems(Q(Q(elem_nodeN(itrac,noden))),nElemtrac))

  min_dist = 10000_precis
  flag_br = .true.

  do while(nElemtrac<nElems+1 .and.  flag_br)
    bas1 = basn2 - basn1
    bas1 = bas1/norm2(bas1)
    bas2 = direct

    trans(:,1) = bas1
    trans(:,2) = bas2

    trans1(:,1) = (/ trans(2,2), -trans(2,1)/)
    trans1(:,2) = (/ -trans(1,2), trans(1,1)/)
    trans1 = trans1/(trans(1,1)*trans(2,2)-trans(2,1)*trans(1,2))

    nodn = matmul(trans1,nod - basn1)
    nodn = (/nodn(1), 0.0_precis /)
    nodp = matmul(trans, nodn) + basn1

    waysElemsCount(elemn) = waysElemsCount(elemn) + 1
    waysElems(elemn,waysElemsCount(elemn)) = nElemtrac
    waysLen(elemn,waysElemsCount(elemn)) = distance(nodp,nod)

if (elemn==251 .and. noden == 14 .and. debug)    print *, "*******"
if (elemn==251 .and. noden == 14 .and. debug)    print *, waysElems(elemn,waysElemsCount(elemn)), waysLen(elemn,waysElemsCount(elemn))
if (elemn==251 .and. noden == 14 .and. debug)    print *, nodp
if (elemn==251 .and. noden == 14 .and. debug)    print *, "--"
if (elemn==251 .and. noden == 14 .and. debug) then
do i=1,3
      print *, nodes(:,elems(itrac,waysElems(elemn,waysElemsCount(elemn))))
    end do
end if


    nElem_prev = nElemtrac
    nElemtrac = nElem_next

    min_dist = 10000_precis

    do stepi=1,3
      if(distance(nodp,nodes(:,elems(stepi,nElemtrac)))<min_dist) then
        min_dist = distance(nodp,nodes(:,elems(stepi,nElemtrac)))
        min_dist_n = elems(stepi,nElemtrac)
      end if
    end do

    flag_br = .true.
    if (min_dist<1E-8 .and. is_bc_node(min_dist_n)>0.5) flag_br = .false.

    if (nElemtrac<nElems+1 .and. flag_br) then
      nod = nodp
      jtrac=1
      do while (elems_elems(jtrac,nElemtrac)/=nElem_prev)
        jtrac = jtrac + 1
      end do

if (elemn==251 .and. noden == 14 .and. debug)      print *, nodes(:,elems(Q(jtrac),nElemtrac))
if (elemn==251 .and. noden == 14 .and. debug)      print *, nodes(:,elems(Q(jtrac),nElemtrac))-nod
if (elemn==251 .and. noden == 14 .and. debug)      print *, dot_product(nodes(:,elems(Q(jtrac),nElemtrac))-nod,normd), signumr(dot_product(nodes(:,elems(Q(jtrac),nElemtrac))-nod,normd))
if (elemn==251 .and. noden == 14 .and. debug)      print *, "-------"
if (elemn==251 .and. noden == 14 .and. debug)      print *, nodes(:,elems(Q(Q(jtrac)),nElemtrac))
if (elemn==251 .and. noden == 14 .and. debug)         print *, nodes(:,elems(Q(Q(jtrac)),nElemtrac))-nod
if (elemn==251 .and. noden == 14 .and. debug) print *, dot_product(nodes(:,elems(Q(Q(jtrac)),nElemtrac))-nod,normd), signumr(dot_product(nodes(:,elems(Q(Q(jtrac)),nElemtrac))-nod,normd))

      if (signumr(dot_product(nodes(:,elems(Q(jtrac),nElemtrac))-nod,normd))*signumr(dot_product(nodes(:,elems(Q(Q(jtrac)),nElemtrac))-nod,normd))>-1E-14) then
        nribp = Q(Q(jtrac))
      else
        nribp = Q(jtrac)
      endif

!      if (abs(dot_product(nodes(:,elems(Q(j),nElem))-nod,normd))<5E-8 .or. abs(dot_product(nodes(:,elems(Q(Q(j)),nElem))-nod,normd))<5E-8) then

!      end if
!      print *, nribp
      basn1 = nodes(:,elems(nribp,nElemtrac))
      basn2 = nodes(:,elems(Q(nribp),nElemtrac))

!      print *, basn1
!      print *, basn2
      nElem_next = elems_elems(nribp,nElemtrac)
    end if
if (elemn==251 .and. noden == 14 .and. debug)    read *

  end do
  waysSrcNode(elemn,:) = nodp

end subroutine tracing
!!!!!!!! первое приближение

subroutine meshway

  real (kind=precis), dimension(1:2) :: direc, centb, normd, nod
  integer :: node, elem, stepi
  real (kind=precis) :: collecter, bufr, bufi,  int_curr
  integer :: add_direct

  real (kind=precis) :: min_dist

  real (kind=precis), dimension(1:2) :: min_dist_vect

  integer :: lentest = 1

  print *, "meshway start"
  print *, "==============="

  do elem=1,nElems
    intensityCell(elem) = 0.0_precis
  end do

  allocate(waysLen(2*nElems_b,4*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(waysElems(2*nElems_b,4*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(waysSrcNode(2*nElems_b,2), stat = err)
  print *, 'allocate stat=', err
  allocate(waysElemsCount(2*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(waysNodeDirPhi(2*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(intensityTR(nNodes,2*nElems_b), stat = err)
  print *, 'allocate stat=', err

  print *, "hmm = ", hmm*hmin
  do node = 1,nNodes
    if (mod(node,1000)==0) print *, "meshway node=",node

    waysNDirs = nElems_b
    min_dist = 100.0
    do elem = 1,nElems_b
      centb = 0.5*(nodes(:,bc_nodes(1,elem))+nodes(:,bc_nodes(2,elem)))
      direc = nodes(:,node) - centb
      direc = direc/norm2(direc)
      nod = edge_n(:,bc_edges(elem))

      waysElemsCount(elem) = 0
      waysNodeDirPhi(elem) = 0
      waysElemsCount(elem+nElems_b) = -1
      waysNodeDirPhi(elem+nElems_b) = 0

      if (is_bc_node(node)<0.5 .or. abs(dot_product(direc,nod))>1E-7) then
        waysNodeDirPhi(elem) = signumr(direc(2))*Acos(direc(1)) + (1.0-signumr(direc(2)))*Pi
        call tracing(direc,node,elem)
        if (lentest<waysElemsCount(elem)) lentest=waysElemsCount(elem)
        if (distance(nodes(:,node),centb)<min_dist) then
          min_dist = distance(nodes(:,node),centb)
          min_dist_vect = nodes(:, node) - centb
        end if

        if (is_bc_node(node)>0.5) then
          waysNDirs = 2*nElems_b
          waysNodeDirPhi(nElems_b + elem) = signumr(-direc(2))*Acos(-direc(1)) + (1.0-signumr(-direc(2)))*Pi
          waysElemsCount(nElems_b + elem) = 0
          waysSrcNode(nElems_b + elem,:) = nodes(:,node)
        end if
      endif
    end do

    if (min_dist<hmm*hmin .and. is_bc_node(node)<0.5) then
      do add_direct = 1, nElems_b
        centb = 0.5*(nodes(:,bc_nodes(1,add_direct))+nodes(:,bc_nodes(2,add_direct)))
        direc = - nodes(:,node) + centb
        direc = direc/norm2(direc)
        if (dot_product(direc, min_dist_vect)>0.0) then
          waysNodeDirPhi(nElems_b+add_direct) = signumr(direc(2))*Acos(direc(1)) + (1.0-signumr(direc(2)))*Pi
          call tracing(direc,node,nElems_b+add_direct)
          if (lentest<waysElemsCount(nElems_b+add_direct)) lentest=waysElemsCount(nElems_b+add_direct)
        else
          waysNodeDirPhi(nElems_b+add_direct) = 0.0
          waysElemsCount(nElems_b+add_direct) = -1
        end if
      end do
      waysNDirs = 2*nElems_b
    end if

    do elem=1,waysNDirs
      int_curr = bintens(waysSrcNode(elem,:),waysNodeDirPhi(elem), 0.0_precis)
      do stepi=waysElemsCount(elem),1,-1
        int_curr=intty_step(int_curr, absorp_koef(waysElems(elem,stepi)), absorp_koef(waysElems(elem,stepi))*intens_p(waysElems(elem,stepi)), waysLen(elem,stepi))
      end do
      intensityTR(node,elem) = int_curr
      if (waysElemsCount(elem)==-1) intensityTR(node,elem) = 0.0_precis
    end do

    do i=1,waysNDirs
      do j=waysNDirs,i+1,-1
        if(waysNodeDirPhi(i)>waysNodeDirPhi(j)) then
          bufi = waysNodeDirPhi(i)
          waysNodeDirPhi(i) = waysNodeDirPhi(j)
          waysNodeDirPhi(j) = bufi
          bufr = intensityTR(node,i)
          intensityTR(node,i) = intensityTR(node,j)
          intensityTR(node,j) = bufr
        end if
      end do
    end do

    collecter = 0.0_precis

    do nDir=1,waysNDirs-1
      if (abs(intensityTR(node,nDir)+intensityTR(node,nDir+1))>1E-4 .and. abs(waysNodeDirPhi(nDir+1)-waysNodeDirPhi(nDir))>1E-4) then
        collecter = collecter + (intensityTR(node,nDir)+intensityTR(node,nDir+1)) * (waysNodeDirPhi(nDir+1)-waysNodeDirPhi(nDir))
      end if
    end do
    collecter = collecter + (intensityTR(node,1)+intensityTR(node,waysNDirs))*(waysNodeDirPhi(1)-waysNodeDirPhi(waysNDirs)+2.*Pi)

    do i=1,node_elems_count(node)
      if(node_elem(i,node)<nElems+1) then
        intensityCell(node_elem(i,node)) = intensityCell(node_elem(i,node)) + collecter*0.3333333333333333333333333
      end if
    end do
  end do

  print *, lentest, "/", nElems_b

  print *, "meshway finish"
  print *, "==============="

end subroutine meshway


subroutine meshway_angr

  real (kind=precis), dimension(1:2) :: direc, centb, normd, nod
  integer :: node, elem, stepi
  real (kind=precis) :: collecter, bufr, bufi,  int_curr, gmmd, sgmmd!, bang

  integer :: lentest = 1
  integer :: flag

  print *, "meshway_angr start"
  print *, "==============="

  allocate(intensityTR(nNodes,nDirects), stat = err)
  print *, 'allocate stat=', err

!print *, nDirects



  allocate(waysLen(nDirects,4*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(waysElems(nDirects,4*nElems_b), stat = err)
  print *, 'allocate stat=', err
  allocate(waysSrcNode(nDirects,2), stat = err)
  print *, 'allocate stat=', err
  allocate(waysElemsCount(nDirects), stat = err)
  print *, 'allocate stat=', err
  allocate(waysNodeDirPhi(nDirects), stat = err)
  print *, 'allocate stat=', err

  print *, "meshway_calc start"

  do elem=1,nElems
    intensityCell(elem) = 0.0_precis
  end do

  do node = 1,nNodes
    if (mod(node,1000)==0) print *, "meshway_angr node=",node
!print *, "1hello"
    collecter = 0.0_precis
!    bang =  0.0_precis
!print *, "2hello"
    do nDir = 1,nDirects
!print *,nDir, "1hello"
      direc = directns(nDir,1:2)
      gmmd = directns(nDir,3)
      sgmmd = sqrt(1.0_precis - gmmd**2)

if (debug)      print *, direc

!print *,nDir, "2hello"
      waysElemsCount(nDir) = 0
      waysNodeDirPhi(nDir) = signumr(direc(2))*Acos(direc(1)/sgmmd) + (1.0-signumr(direc(2)))*Pi

if (debug)      print *, node, nDir, " hi1"

      flag = 0
!print *,nDir, "3hello"
      if (is_bc_node(node)<0.5) then
if (debug)      print *, node, nDir, " hi11"

        call tracing(direc,node,nDir)
        if (lentest<waysElemsCount(nDir)) lentest=waysElemsCount(nDir)
      else
        if (is_source_node(is_bc_node(node),nDir)>0.5 .or. is_source_node(is_bc_node(node),nDir)<-0.5) then
if (debug)      print *, node, nDir, " hi12"
          if (is_source_node(is_bc_node(node),nDir)<-0.5) flag = 1
          waysElemsCount(nDir) = 0
          waysSrcNode(nDir,:) = nodes(:,node)
        else
if (debug)      print *, node, nDir, " hi13"
!      print *, direc
!      print *, nodes(:, node)

          call tracing(direc,node,nDir)
          if (lentest<waysElemsCount(nDir)) lentest=waysElemsCount(nDir)
        end if
      endif

if (debug)      print *, node, nDir, " hi2"

      int_curr = bintens(waysSrcNode(nDir,:),waysNodeDirPhi(nDir), gmmd)
      do stepi=waysElemsCount(nDir),1,-1
          int_curr=intty_step(int_curr, absorp_koef(waysElems(nDir,stepi))+disper_koef(waysElems(nDir,stepi)), absorp_koef(waysElems(nDir,stepi))*intens_p(waysElems(nDir,stepi))+scatsource(waysElems(nDir,stepi),nDir), waysLen(nDir,stepi)/sqrt(1.0_precis-directns(nDir,3)**2))
      end do
      intensityTR(node,nDir) = int_curr
      if (flag>0.5) intensityTR(node,nDir) = 0.0

if (debug)      print *, node, nDir, " hi3"
      if (useFiles) then
        collecter = collecter + intensityTR(node,nDir)*bodyAngleTR(nDir)
!       bang = bang + bodyAngleTR(nDir)
      else
!        collecter = collecter + intensityTR(node,nDir)*bodyAngleTR(nDir)
        collecter = collecter + intensityTR(node,nDir)*bodyAngle
!       bang = bang + bodyAngle
      end if
    end do

!    print *, "*****", bang
!    read *

    do i=1,node_elems_count(node)
      if(node_elem(i,node)<nElems+1) then
        intensityCell(node_elem(i,node)) = intensityCell(node_elem(i,node)) + collecter*0.3333333333333333333333333
      end if
    end do
  end do

  print *, lentest, "/", nElems_b

  print *, "meshway_angr finish"
  print *, "==============="

end subroutine meshway_angr


pure function bintens(nod, directPhi, directGam) result(binten)
  real(kind=precis), dimension(1:2), Intent(In) :: nod
  real(kind=precis) :: binten, phhh
  real(kind=precis), Intent(In) :: directPhi, directGam
  real(kind=precis), dimension(2) :: centr,centr1

  centr=(/0.1,0.0/)
  centr1=(/0.5,0.0/)
!  phhh=13.*Pi/18.0
  phhh = Acos(1./sqrt(5.))
  binten = 0.0_precis !.and. abs(directPhi-phhh)<0.05

  select case (testtask)
    case (1)
      if (distance(centr, nod)<0.03 .and. Cos(directPhi-phhh)>0.9999 .and. abs(directGam)<0.1) binten = 0.1
    case (2)
      if (distance(centr1, nod)<0.0101) binten = 0.05_precis
    case (3)
      if (distance(centr1, nod)<0.3) binten = 1.0_precis! - (distance(centr1, nod)/0.2)**2
  end select
end function

pure function intty_step (i0, kap, sour, step) result(intty)
  real(kind=precis), Intent(In) :: i0, kap, sour, step
  real(kind=precis) :: intty

  intty = i0*Exp(-kap*step)
  if (abs(kap)>1E-7) intty = intty + (1.0_precis-Exp(-kap*step))*sour/kap
end function intty_step

subroutine radiation_solverTR
  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)
  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour

  if (useFiles) then
    call InitGroupsTri
!    call InitGroupsTriRaz

!    call sphereTesselate
    print *, 'InitGroupsTri done'
  else
    call InitGroupsSn
!    call sphereTesselate
    print *, 'InitGroups done'
  end if

  call IntensAllocate
  print *, 'IntensAllocate done'

  call rad_param_calc_TR
  print *, 'rad_param_calc_TR done'

  if (useBoundary) then
    call meshway
  else
!    call meshway_angr_pp
    call meshway_axial_pp
  end if
  print *, 'meshway done'

  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)

  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour - ddtime

  print *, 'Computation time (msec): ', ddtime

  call save_vtk_solution_radTR("-TRSN-"//i2c(nDirects))

  read *
end subroutine radiation_solverTR

subroutine rad_param_calc_TR
  real (kind=precis), dimension(2) :: centr=(/0.5,0.5/)
  real (kind=precis), dimension(2) :: centr1=(/0.5,0.0/)
  real (kind=precis), dimension(2) :: cnt, sec
  real (kind=precis) :: r0= 0.1
  integer :: nDir1
  integer :: seed=100
  real :: scat, scatsum

  do nElem=1,nElems
    absorp_koef(nElem) = 0.0_precis !absorb0*absorbtion(U_til(:,nElem))
!    if (distance(centr,elems_center(:,nElem))<0.1) absorp_koef(nElem) = 1000.
    intens_p(nElem) = 0.0_precis
    do nDir=1,nDirects
      scatsource(nElem, nDir) = 0.0_precis
    end do
    disper_koef(nElem) = 0.0_precis
  end do

end subroutine rad_param_calc_TR

subroutine save_vtk_solution_radTR(param)

  CHARACTER(LEN=*), INTENT(in):: param


  real (kind=precis) :: collecter, coll1
  integer :: ccc

  call filename_se_sol_vtk_radTR(niterDG,current_sol_name,param)

  open (unit=1,file=current_sol_name)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', nNodes, ' float'

  do i=1,nNodes
    write (1,*) nodes(1,i), nodes(2,i), 0
  end do

  write (1,"(A,I7,A,I7)") 'CELLS ', nElems, ' ', nElems+3*nElems

  do i=1,nElems
    write (1,*) 3, elems(1,i)-1, elems(2,i)-1, elems(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', nElems
  do i=1,nElems
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'CELL_DATA ', nElems

  write (1,"(A)") 'SCALARS intensity float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nElems
    write(1,*) intensityCell(i)
  end do

  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'

end subroutine save_vtk_solution_radTR

subroutine sphereTesselate

  real(kind=precis), allocatable :: directnodes(:,:)
  real(kind=precis), allocatable :: directelems(:,:)
  real :: foo
  integer :: nDirectNodes, node, elem, nDir
  real(kind=precis), dimension(1:3) :: bufcent
  real(kind=precis) :: bbb

  open(1, file='./data/mesh/sphere_elems.dat')

  read(1,*) nDirects
  print *, '  > nDirects   =', nDirects

  allocate(directns(1:nDirects,1:3), stat = err)
  allocate(bodyAngleTR(1:nDirects), stat = err)
  allocate(directelems(1:3,1:nDirects), stat = err)
!  bbb = 0.0_precis
  do nDir = 1, nDirects
    read (1,*) directelems(1,nDir), directelems(2,nDir), directelems(3,nDir)
  enddo

  close(1)

  open(1, file='./data/mesh/sphere_nodes.dat')

  read(1,*) nDirectNodes
  print *, '  > nDirectNodes =', nDirectNodes
  allocate(directnodes(1:3,1:nDirectNodes), stat = err)

  do node = 1,nDirectNodes
    read (1,*) directnodes(1,node), directnodes(2,node), directnodes(3,node)
  enddo

  close(1)

  bbb = 0.0_precis
  do nDir=1,nDirects
    bufcent = (/0.0_precis,0.0_precis,0.0_precis/)
    do node=1,3
      bufcent = bufcent + (1.0/3.0)*directnodes(:,directelems(node,nDir))
    end do
    directns(nDir,:) = bufcent/norm2_3(bufcent)
    bodyAngleTR(nDir) = TriangleVol3(directnodes(:,directelems(1,nDir)),directnodes(:,directelems(2,nDir)),directnodes(:,directelems(3,nDir)))
    bbb = bbb + bodyAngleTR(nDir)
  end do

  print *, "Sum body angle = ", bbb
!  read *

end subroutine sphereTesselate

subroutine meshway_angr_pp

  real (kind=precis), dimension(1:2) :: direc
  integer :: node, elem
  real (kind=precis) :: collecter, bufr, bufi,  int_curr, gmmd, sgmmd!, bang

  integer :: lentest = 1
  integer :: flag, stepin

  integer :: ktrac

  real (kind=precis), dimension(4000) :: waysLenp
  real (kind=precis) :: waysNodeDirPhip
  real (kind=precis), dimension(2) :: waysSrcNodep
  integer, dimension(4000) :: waysElemsp
  integer :: waysElemsCountp
  integer :: waysNDirsp

  integer :: OMP_GET_THREAD_NUM

  print *, "meshway_angr_pp start"
  print *, "==============="

  allocate(intensityTR(nNodes,nDirects), stat = err)
  print *, 'allocate stat=', err

  do elem=1,nElems
    intensityCell(elem) = 0.0_precis
  end do

  print *, nDirects

!$omp parallel

  ktrac =0
!$omp do SCHEDULE (GUIDED ,100) private(nElem,node,collecter,stepin,nDir,direc,gmmd,sgmmd,flag, bufr,bufi,int_curr,i,j,waysLenp,waysNodeDirPhip,waysSrcNodep,waysElemsp,waysElemsCountp, waysNDirsp)

  do node = 1,nNodes
    ktrac = ktrac + 1
    if (mod(ktrac,1000)==0) then
      print *, "meshway_angr ktrac=",ktrac," node=",node
  !    read *
    end if

    collecter = 0.0_precis

    do nDir = 1,nDirects

      direc = directns(nDir,1:2)
      gmmd = directns(nDir,3)
      sgmmd = sqrt(1.0_precis - gmmd**2)

      waysElemsCountp = 0
      waysNodeDirPhip = signumr(direc(2))*Acos(direc(1)/sgmmd) + (1.0-signumr(direc(2)))*Pi

      flag = 0

      if (is_bc_node(node)<0.5) then
        call tracing_pp(node,direc/norm2(direc), waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
        if (lentest<waysElemsCountp) lentest=waysElemsCountp
      else
       if (is_source_node(is_bc_node(node),nDir)>0.5 .or. is_source_node(is_bc_node(node),nDir)<-0.5) then
          if (is_source_node(is_bc_node(node),nDir)<-0.5) flag = 1
          waysElemsCountp = 0
          waysSrcNodep(:) = nodes(:,node)
       else
          call tracing_pp(node,direc/norm2(direc), waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
          if (lentest<waysElemsCountp) lentest=waysElemsCountp
       end if
      endif

      int_curr = bintens(waysSrcNodep(:),waysNodeDirPhip, gmmd)
      do stepin=waysElemsCountp,1,-1
          int_curr=intty_step(int_curr, absorp_koef(waysElemsp(stepin))+disper_koef(waysElemsp(stepin)), absorp_koef(waysElemsp(stepin))*intens_p(waysElemsp(stepin))+scatsource(waysElemsp(stepin),nDir), waysLenp(stepin)/sqrt(1.0_precis-directns(nDir,3)**2))
      end do
      intensityTR(node,nDir) = int_curr
      if (flag>0.5) intensityTR(node,nDir) = 0.0

      if (useFiles) then
        collecter = collecter + intensityTR(node,nDir)*bodyAngleTR(nDir)
      else
        collecter = collecter + intensityTR(node,nDir)*bodyAngle
      end if
    end do

    do stepin=1,node_elems_count(node)
      if(node_elem(stepin,node)<nElems+1) then
        intensityCell(node_elem(stepin,node)) = intensityCell(node_elem(stepin,node)) + collecter*0.3333333333333333333333333
      end if
    end do
  end do
!$omp end do
!$omp end parallel

  print *, lentest, "/", nElems_b

  print *, "meshway_angr finish"
  print *, "==============="

end subroutine meshway_angr_pp

pure subroutine tracing_pp(nnode,direct, waysSrcNodepp, waysElemsCountpp, waysElemspp, waysLenpp)
      logical :: cond1, cond2, flag_br
      integer :: nElem_next, nElem_prev, nribp, dotn1, dotn2, min_dist_n, itrac, jtrac, nElemtrac, bbf, stepi
      real (kind=precis), dimension(1:2) :: normd, bas1, bas2, basn1, basn2, nod, nodp, nodn
      real (kind=precis), dimension(1:2,1:2) :: trans, trans1
      real (kind=precis) :: min_dist

      real (kind=precis), dimension(1:2), Intent(In) :: direct
      real (kind=precis), dimension(2), Intent(InOut) :: waysSrcNodepp
      integer, Intent(InOut) :: waysElemsCountpp
      integer, dimension(4000), Intent(InOut) :: waysElemspp
      real (kind=precis), dimension(4000), Intent(InOut) :: waysLenpp
      integer, Intent(In) :: nnode

        normd = (/direct(2), -direct(1)/)
        itrac=1
        nElemtrac = node_elem(itrac,nnode)
        dotn1 = elems(Q(elem_nodeN(itrac,nnode)),nElemtrac)
        dotn2 = elems(Q(Q(elem_nodeN(itrac,nnode))),nElemtrac)
        cond1 = dot_product(nodes(:,dotn1)-nodes(:,nnode),normd)*dot_product(nodes(:,dotn2)-nodes(:,nnode),normd)>1E-14
        cond2 = (dot_product(nodes(:,dotn1)-nodes(:,nnode),direct)>0. .and. dot_product(nodes(:,dotn2)-nodes(:,nnode),direct)>0.)

        do while ((cond1 .or. cond2 .or. nElemtrac>nElems) .and.  itrac<=node_elems_count(nnode))
          itrac = itrac + 1
          nElemtrac = node_elem(itrac,nnode)
          dotn1 = elems(Q(elem_nodeN(itrac,nnode)),nElemtrac)
          dotn2 = elems(Q(Q(elem_nodeN(itrac,nnode))),nElemtrac)

          cond1 = dot_product(nodes(:,dotn1)-nodes(:,nnode),normd)*dot_product(nodes(:,dotn2)-nodes(:,nnode),normd)>1E-14
          cond2 = (dot_product(nodes(:,dotn1)-nodes(:,nnode),direct)>0. .and. dot_product(nodes(:,dotn2)-nodes(:,nnode),direct)>0.)
        end do

        nElemtrac = node_elem(itrac,nnode)
        nElem_next = elems_elems(Q(elem_nodeN(itrac,nnode)),nElemtrac)
        nod = nodes(:,nnode)
        basn1 = nodes(:,elems(Q(elem_nodeN(itrac,nnode)),nElemtrac))
        basn2 = nodes(:,elems(Q(Q(elem_nodeN(itrac,nnode))),nElemtrac))

        min_dist = 10000_precis
        flag_br = .true.

        do while(nElemtrac<nElems+1 .and.  flag_br)
          bas1 = basn2 - basn1
          bas1 = bas1/norm2(bas1)
          bas2 = direct

          trans(:,1) = bas1
          trans(:,2) = bas2

          trans1(:,1) = (/ trans(2,2), -trans(2,1)/)
          trans1(:,2) = (/ -trans(1,2), trans(1,1)/)
          trans1 = trans1/(trans(1,1)*trans(2,2)-trans(2,1)*trans(1,2))

          nodn = matmul(trans1,nod - basn1)
          nodn = (/nodn(1), 0.0_precis /)
          nodp = matmul(trans, nodn) + basn1

          waysElemsCountpp = waysElemsCountpp + 1
          waysElemspp(waysElemsCountpp) = nElemtrac
          waysLenpp(waysElemsCountpp) = distance(nodp,nod)

          nElem_prev = nElemtrac
          nElemtrac = nElem_next

          min_dist = 10000.0_precis

          do stepi=1,3
            if(distance(nodp,nodes(:,elems(stepi,nElemtrac)))<min_dist) then
              min_dist = distance(nodp,nodes(:,elems(stepi,nElemtrac)))
              min_dist_n = elems(stepi,nElemtrac)
            end if
          end do

          flag_br = .true.
          if (min_dist<1E-8 .and. is_bc_node(min_dist_n)>0.5) flag_br = .false.

          if (nElemtrac<nElems+1 .and. flag_br) then
            nod = nodp
            jtrac=1
            do while (elems_elems(jtrac,nElemtrac)/=nElem_prev)
              jtrac = jtrac + 1
            end do

            if (signumr(dot_product(nodes(:,elems(Q(jtrac),nElemtrac))-nod,normd))*signumr(dot_product(nodes(:,elems(Q(Q(jtrac)),nElemtrac))-nod,normd))>-1E-14) then
              nribp = Q(Q(jtrac))
            else
              nribp = Q(jtrac)
            endif

            basn1 = nodes(:,elems(nribp,nElemtrac))
            basn2 = nodes(:,elems(Q(nribp),nElemtrac))

            nElem_next = elems_elems(nribp,nElemtrac)
          end if
        end do
        waysSrcNodepp(:) = nodp
end subroutine tracing_pp

!!!!!!!!!!! axials

pure subroutine hyperb_line_intersect(x0,dirt,x1,x2,inter_flag,inter_points2D,inter_points3D,par_t,vortin)

  real (kind=precis), dimension(1:2,1:2), Intent(InOut) :: inter_points2D
  real (kind=precis), dimension(1:2,1:3), Intent(InOut) :: inter_points3D
  real (kind=precis), dimension(1:2), Intent(InOut) :: par_t
  integer, dimension(1:2), Intent(InOut) :: vortin
  real (kind=precis), dimension(1:2), Intent(In) :: x1,x2
  real (kind=precis), dimension(1:3), Intent(In) :: x0, dirt
  integer, Intent(InOut) :: inter_flag

  real (kind=precis), dimension(1:2,1:2) :: vertices
  real (kind=precis), dimension(1:2) :: par_t1, par_tt,ribdir,int_p
  real (kind=precis), dimension(1:3) :: buff
  real (kind=precis) :: buf1,buf2,buf3,buf4, buf5, zz, rr
  integer :: itt,jtt,ktt

  vertices(1,:)=x1
  vertices(2,:)=x2

  ribdir = x2-x1

  par_t1 = 0.0
  par_tt = 0.0

  if (Abs(ribdir(1))<1E-9) then
    zz = x2(1)
    par_tt(1) = (x0(3) - zz)/dirt(3)
    buff = x0 - par_tt(1) * dirt
    rr = Sqrt(Sum(buff(1:2)**2))
    par_t1(1) = (rr - x1(2))/ribdir(2)
    par_t1(2) = -1.0
    par_tt(2) = -1.0
  else
    buf2 = -4*((dirt(1)**2 + dirt(2)**2)*ribdir(1)**2 - dirt(3)**2*ribdir(2)**2)*(ribdir(1)**2*(x0(1)**2 + x0(2)**2 - x1(2)**2) -ribdir(2)**2*(x0(3) - x1(1))**2 + 2*ribdir(1)*ribdir(2)*x1(2)*(-x0(3) + x1(1))) + 4*(dirt(1)*ribdir(1)**2*x0(1) + dirt(2)*ribdir(1)**2*x0(2) - dirt(3)*ribdir(2)*(ribdir(2)*x0(3) + ribdir(1)*x1(2) - ribdir(2)*x1(1)))**2
    if (buf2>1E-30) then
      buf2 = Sqrt(buf2)
      buf5 = (dirt(1)**2 + dirt(2)**2)*(ribdir(1)**2) -(dirt(3)**2)*(ribdir(2)**2)
      buf1 = buf5*ribdir(1)
      buf3 = (dirt(3)**2)*ribdir(1)*ribdir(2)*x1(2) - (ribdir(1)**2)*(dirt(1)*dirt(3)*x0(1) + (dirt(1)**2)*(-x0(3) + x1(1)) + dirt(2)*(dirt(3)*x0(2) + dirt(2)*(-x0(3) + x1(1))))

      buf4 = (ribdir(1)**2)*(dirt(1)*x0(1) + dirt(2)*x0(2)) - dirt(3)*ribdir(1)*ribdir(2)*x1(2) + dirt(3)*(ribdir(2)**2)*(-x0(3) + x1(1))

      par_t1(1) = (buf3 + dirt(3)*buf2/2.)/buf1

      par_tt(1) = (buf4 - buf2/2.)/buf5

      par_t1(2) = (buf3 - dirt(3)*buf2/2.)/buf1

      par_tt(2) = (buf4 + buf2/2.)/buf5
    end if
  end if

  par_t = 0.0
  inter_flag = 0
  vortin = 0
  do itt=1,2
    if (par_t1(itt)>-epsilTR .and. par_t1(itt)<1D0+epsilTR .and.  par_tt(itt)>epsilTR) then
      inter_flag = inter_flag + 1
      par_t(inter_flag) = par_tt(itt)
      inter_points2D(inter_flag,:) = x1 + par_t1(itt) * ribdir
      inter_points3D(inter_flag,:) = x0 - par_tt(itt) * dirt
      int_p = x1 + par_t1(itt) * ribdir
      do jtt=1,2
        if (distance(int_p,vertices(jtt,:))<1D-5) then
          vortin(inter_flag) = jtt
        end if
      end do
    end if
  end do

end subroutine hyperb_line_intersect

pure subroutine tracing_axial_par(nnode,direct,waysSrcNodepp, waysElemsCountpp, waysElemspp, waysLenpp)

  real (kind=precis), dimension(1:3), Intent(In) :: direct
  real (kind=precis), dimension(2), Intent(InOut) :: waysSrcNodepp
  integer, Intent(InOut) :: waysElemsCountpp
  integer, dimension(10000), Intent(InOut) :: waysElemspp
  real (kind=precis), dimension(10000), Intent(InOut) :: waysLenpp
  integer, Intent(In) :: nnode

  integer :: itrac, jtrac, nElemtrac, nElemtrac1, nEdgetrac, nElem_next, nElem_prev, inter_flg, vnode

  real (kind=precis), dimension(1:3) :: x0, xb
  real (kind=precis), dimension(1:2,1:2):: ints_points2D
  real (kind=precis), dimension(1:2,1:3):: ints_points3D
  real (kind=precis), dimension(1:2) :: ints_par, probe2D
  real (kind=precis), dimension(1:2) :: next_int2D
  real (kind=precis), dimension(1:3) :: next_int3D, prev_int3D, probe3D, debv
  real (kind=precis), dimension(1:2) :: x1,x2, x02D
  real (kind=precis) :: direct_par, dirp, stepint
  integer, dimension(1:2) :: vortinter, vertex

  x02D = nodes(:,nnode)
  x0 = (/x02D(2), 0._precis, x02D(1)/)
  direct_par = 0.0
  waysElemsCountpp = 0

  dirp = 1E10

  stepint = 0.1
  probe2D = (/ z_max*2., r_max*2./)
  nElemtrac=0
  do while (probe2D(1)>z_max .or. probe2D(2)>r_max .or. nElemtrac==0)
    probe3D = x0 - stepint*direct
    probe2D = (/probe3D(3),Sqrt(Sum(probe3D(1:2)**2))/)
    do itrac=1, node_elems_count(nnode)
      nElemtrac1 = node_elem(itrac,nnode)
      if ( pintr(probe2D, nodes(:,elems(1,nElemtrac1)), nodes(:,elems(2,nElemtrac1)), nodes(:,elems(3,nElemtrac1))) ) nElemtrac = nElemtrac1
    end do
    stepint = 0.5*stepint
  end do
  waysElemsCountpp = waysElemsCountpp + 1
  waysElemspp(waysElemsCountpp) = nElemtrac
  waysLenpp(waysElemsCountpp) = norm2_3(probe3D-prev_int3D)
  next_int3D = probe3D
  dirp = 2.0*stepint

  do while(nElemtrac<nElems+1)
    direct_par = dirp
    prev_int3D = next_int3D
    dirp = 1E10
    do nEdgetrac=1,3
      x1 = nodes(:,elems(nEdgetrac,nElemtrac))
      x2 = nodes(:,elems(Q(nEdgetrac),nElemtrac))
      call hyperb_line_intersect(x0,direct,x1,x2,inter_flg,ints_points2D,ints_points3D,ints_par,vortinter)
      do jtrac=1,inter_flg
         if (ints_par(jtrac)<dirp .and. ints_par(jtrac)>direct_par+epsilTR*direct_par) then
          dirp = ints_par(jtrac)
          next_int2D = ints_points2D(jtrac,:)
          next_int3D = ints_points3D(jtrac,:)
          if (vortinter(jtrac)>0.5) then
            vertex = (/ elems(nEdgetrac,nElemtrac), elems(Q(nEdgetrac),nElemtrac) /)
            vnode = vertex(vortinter(jtrac))
            if (is_bc_node(vnode)>0.5) then
              nElem_next=nElems+1
            else
              stepint = 0.1
!               print *, "aaaaaaa"
              probe2D = (/ z_max*2., r_max*2./)
              do while (probe2D(1)>z_max .or. probe2D(2)>r_max .or. nElem_next == 0)
                probe3D = x0 - (stepint+dirp)*direct
                probe2D = (/probe3D(3),Sqrt(Sum(probe3D(1:2)**2))/)
!               print *, "probe2D", probe2D
                nElem_next = 0
                do itrac=1, node_elems_count(vnode)
                  nElemtrac1 = node_elem(itrac,vnode)
                  if ( pintr(probe2D, nodes(:,elems(1,nElemtrac1)), nodes(:,elems(2,nElemtrac1)), nodes(:,elems(3,nElemtrac1))) ) then
                    nElem_next = nElemtrac1
                  end if
                end do
                stepint = 0.5*stepint
              end do
            endif
            dirp = dirp + 2.0*stepint
          else
            nElem_next = elems_elems(nEdgetrac,nElemtrac)
          end if
        end if
      end do
    end do

    if (dirp > 1000D0) then
      stepint = 0.00001
!      dirp = direct_par + stepint
      probe3D = x0 - (direct_par + stepint)*direct
      probe2D = (/probe3D(3),Sqrt(Sum(probe3D(1:2)**2))/)
!       if (probe2D(1)<-epsilTR .or. probe2D(1)>z_max .or. probe2D(2)<-epsilTR .or. probe2D(2)>r_max) then
!         nElem_next = nElems+1
!         dirp = direct_par + stepint
!         next_int3D = probe3D
!       else
!        nElem_next = nElems+1
        do jtrac = 1, nElems
          if (pintr(probe2D, nodes(:,elems(1,jtrac)), nodes(:,elems(2,jtrac)), nodes(:,elems(3,jtrac)))) then
            nElem_next = jtrac
            dirp = direct_par + stepint
            next_int3D = probe3D
            exit
          end if
        end do
!       end if
    end if

    waysElemsCountpp = waysElemsCountpp + 1
    waysElemspp(waysElemsCountpp) = nElemtrac
    waysLenpp(waysElemsCountpp) = norm2_3(next_int3D-prev_int3D)

!debv = (/ 0.169657553918323   ,    0.449330292956042   ,   -0.877108090391651  /)
! debv = (/ 0.183095907194069  ,     0.982828652789592    ,   2.288510486271478E-002 /)
!  debv = (/ 0.296235120733546   ,    0.110296025454408     ,  0.948725218391997 /)
!  if (nnode==153 .and. distance(debv,direct)<1E-6) then
!      print *, "****************"
!      print *, nElemtrac, nDir, nElem_next, nElems
!      do itrac=1,3
!        print *, "x={",nodes(2,elems(itrac,nElemtrac)),",0,",nodes(1,elems(itrac,nElemtrac)),"};"
!      enddo
!      print *, dirp
!      print *, "x={",next_int2D(2),",0,",next_int2D(1),"};"
!      print *, vortinter
!      print *, "probe2D=", probe2D
!      print *, waysLenpp(waysElemsCountpp)
!
!       read *
!  end if

    nElemtrac = nElem_next
  end do

  waysSrcNodepp = next_int2D

end subroutine tracing_axial_par

subroutine meshway_axial_pp

  real (kind=precis), dimension(1:2) :: direc
  integer :: node, elem, buf_int
  real (kind=precis) :: collecter, bufr, bufi,  int_curr, gmmd, sgmmd!, bang

  integer :: lentest = 1
  integer :: flag, stepin
  integer :: ktrac

  real (kind=precis), dimension(10000) :: waysLenp
  real (kind=precis), dimension(2) :: waysSrcNodep
  integer, dimension(10000) :: waysElemsp
  integer :: waysElemsCountp
  integer :: waysNDirsp

  integer :: OMP_GET_THREAD_NUM

  print *, "meshway_axial_pp start"
  print *, "==============="

  allocate(intensityTR(nNodes,nDirects), stat = err)
  print *, 'allocate stat=', err

  do elem=1,nElems
    intensityCell(elem) = 0.0_precis
  end do

  print *, nDirects

!$omp parallel

  ktrac =0
!$omp do SCHEDULE (GUIDED ,100) private(nElem,node,collecter,stepin,nDir,direc,gmmd,sgmmd,flag, bufr,bufi,int_curr,i,j,waysLenp,waysSrcNodep,waysElemsp,waysElemsCountp, waysNDirsp)

  do node = 1,nNodes
    ktrac = ktrac + 1
    if (mod(ktrac,200)==0) then
      print *, "meshway_axial ktrac=",ktrac," node=",node
    end if
    collecter = 0.0_precis

    do nDir = 1,nDirects

      direc = directns(nDir,1:2)
      gmmd = directns(nDir,3)
      sgmmd = sqrt(1.0_precis - gmmd**2)

      waysElemsCountp = 0
      flag = 0

      if (is_bc_node(node)<0.5) then
        call tracing_axial_par(node,directns(nDir,:), waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
        if (lentest<waysElemsCountp) lentest=waysElemsCountp
      else
        if (is_source_node(is_bc_node(node),nDir)>0.5 .or. is_source_node(is_bc_node(node),nDir)<-0.5) then
          if (is_source_node(is_bc_node(node),nDir)<-0.5) flag = 1
          waysElemsCountp = 0
          waysSrcNodep(:) = nodes(:,node)
        else
          call tracing_axial_par(node,directns(nDir,:), waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
          if (lentest<waysElemsCountp) lentest=waysElemsCountp
            end if
      endif
      int_curr = bintens_axial(waysSrcNodep(:),directns(nDir,:))
      do stepin=waysElemsCountp,1,-1
          int_curr=intty_step(int_curr, absorp_koef(waysElemsp(stepin))+disper_koef(waysElemsp(stepin)), absorp_koef(waysElemsp(stepin))*intens_p(waysElemsp(stepin))+scatsource(waysElemsp(stepin),nDir), waysLenp(stepin))
      end do
      intensityTR(node,nDir) = int_curr
      if (flag>0.5) intensityTR(node,nDir) = 0.0
      if (useFiles) then
        collecter = collecter + intensityTR(node,nDir)*bodyAngleTR(nDir)
      else
        collecter = collecter + intensityTR(node,nDir)*bodyAngle
      end if
    end do
    do stepin=1,node_elems_count(node)
      if(node_elem(stepin,node)<nElems+1) then
        intensityCell(node_elem(stepin,node)) = intensityCell(node_elem(stepin,node)) + collecter*0.3333333333333333333333333
      end if
    end do
  end do
!$omp end do
!$omp end parallel
  print *, lentest, "/", nElems_b

  print *, "meshway_axial finish"
  print *, "==============="

end subroutine meshway_axial_pp

pure function bintens_axial(nod, direct) result(binten)
  real(kind=precis), dimension(1:2), Intent(In) :: nod
  real(kind=precis), dimension(1:3), Intent(In) :: direct
  real(kind=precis) :: binten, phhh
  real(kind=precis), dimension(2) :: centr,centr1

  binten = 0.0_precis

  if (nod(2)<r_disk .and. nod(1)<0.05) binten = 1.0_precis! - (distance(centr1, nod)/0.2)**2
end function

pure function bodyAng(x0,x1,x2,x3) result(ang)
  real(kind=precis), dimension(1:3), Intent(In) :: x0,x1,x2,x3
  real(kind=precis) :: ang

  real(kind=precis), dimension(1:3) :: n1,n2,n3
  real(kind=precis) :: alp, bet, gam

  n1 = vect_mult(x1-x0,x2-x0)
  n1 = n1 / norm2_3(n1)
  n2 = vect_mult(x2-x0,x3-x0)
  n2 = n2 / norm2_3(n2)
  n3 = vect_mult(x3-x0,x1-x0)
  n3 = n3 / norm2_3(n3)

  alp = Acos(-dot_product(n1,n2))
  bet = Acos(-dot_product(n2,n3))
  gam = Acos(-dot_product(n3,n1))

  ang = alp + bet + gam - pi

end function

end module radiationTR