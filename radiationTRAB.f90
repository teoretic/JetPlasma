module radiationTRAB

use radiationTR

use iso_c_binding


implicit none


logical :: notDiskOnly = .true.

!! добавочные направления в аккреционной задаче
integer :: diskNNodes, diskNElems
integer , allocatable :: disk_elems(:,:)
real(kind=precis), allocatable :: disk_nodes(:,:)
real(kind=precis), allocatable :: disk_srcnodes(:,:)
real (kind=precis), allocatable :: bodyAngleTRAB(:)
real (kind=precis), allocatable :: directnsAB(:,:)

integer, allocatable :: addDirsDivRange(:,:)
integer, allocatable :: addDirsDiv(:)
integer :: addDirsDivN
integer :: addDirsMax = 100
integer :: addKoef = 20

integer, parameter :: fullDirs = 10000
integer, parameter :: fullNodes = 10000
integer :: dirNAB1, dirNDnAB1
real (kind=precis), dimension(1:fullDirs) :: bodyAngAB1
real (kind=precis), dimension(1:fullDirs,1:3) :: dirsAB1
real (kind=precis), dimension(1:fullNodes,1:3) :: directnsNDAB1
integer, dimension(1:3,1:fullDirs) :: directelemsAB1
logical, dimension(1:fullDirs) :: toDiskAB1
logical, allocatable :: diskAddDiv(:)

integer, dimension(1:fullNodes) :: nDirectsAB

real (kind=precis) :: minDirMax = 2D-1

contains

subroutine InitGroupsTriAB

  CHARACTER(LEN = max_filename_length):: ffname

  real :: foo
  integer :: node, elem, nDir
  real(kind=precis), dimension(1:2) :: bufcent
  real(kind=precis), dimension(1:3,1:2) :: xp
  real(kind=precis), dimension(1:3) :: buf_vec, xp0, xdir
  real(kind=precis), dimension(1:3,1:3) :: drot
  real(kind=precis) :: bbb, buf_real, cang, sang
  integer :: itr, jtr, ktr, buf_int,itrac

  integer, allocatable :: buf_nums(:)
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
  allocate(directelems_buf(1:3,1:nDirects), stat = err)
  print *, "Allocate stat = ", err

  allocate(directns(1:nDirects+1500,1:3), stat = err)
  print *, "Allocate stat = ", err
  allocate(bodyAngleTR(1:nDirects+1500), stat = err)
  print *, "Allocate stat = ", err
  allocate(directnsND(1:nDirectNodes+1500,1:3), stat = err)
  print *, "Allocate stat = ", err

  do nDir = 1, nDirects
    read (11,*) directelems_buf(1,nDir), directelems_buf(2,nDir), directelems_buf(3,nDir), foo, foo, foo, foo
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

!$omp parallel
!$omp do SCHEDULE (GUIDED ,100) private(nDir, buf_vec, node)
  do nDir=1,nDirects
    buf_vec = (/0D0, 0D0, 0D0/)
    do node=1,3
      buf_vec = buf_vec + (1.0/3.0)*directnsND(directelems(node,nDir),:)
    end do
    directns(nDir,:) = buf_vec/norm2_3(buf_vec)
  end do
!$omp end do
!$omp end parallel

  buf_real = 0D0
  buf_vec = (/0D0, 0D0, 0D0/)
  bodyAngle = 0D0

  do nDir=1,nDirects
    bodyAngleTR(nDir) = bodyAng(buf_vec, directnsND(directelems(1,nDir),:),directnsND(directelems(2,nDir),:),directnsND(directelems(3,nDir),:))
    buf_real = buf_real + bodyAngleTR(nDir)
  end do
  bodyAngle=maxval(bodyAngleTR(1:nDirects))

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

    !===================
    ! разбиение аккреционного диска
  open(101, file='./data/mesh/circ_elems.dat')
  open(102, file='./data/mesh/circ_nodes.dat')

  read(101,*) diskNElems
  print *, '  > diskNElems   =', diskNElems

  read(102,*) diskNNodes
  print *, '  > diskNNodes =', diskNNodes

  allocate(disk_nodes(1:3,1:diskNNodes), stat = err)
  print *, "Allocate stat = ", err
  allocate(disk_elems(1:3,1:diskNElems), stat = err)
  print *, "Allocate stat = ", err
  allocate(disk_srcnodes(1:diskNElems,1:3), stat = err)
  print *, "Allocate stat = ", err

  do node = 1,diskNNodes
    read (102,*) disk_nodes(1,node), disk_nodes(2,node)
    disk_nodes(3,node) = 0D0
  enddo

  do nDir = 1, diskNElems
    read (101,*) disk_elems(1,nDir), disk_elems(2,nDir), disk_elems(3,nDir), foo, foo, foo, foo
    disk_srcnodes(nDir,3) = 0D0
    disk_srcnodes(nDir,1:2) = 0D0
    do node=1,3
      disk_srcnodes(nDir,1:2) = disk_srcnodes(nDir,1:2)  + disk_nodes(1:2,disk_elems(node,nDir))/3D0
    end do
  enddo

  close(101)
  close(102)

  allocate(diskAddDiv(1:nNodes), stat = err)
  print *, "Allocate stat = ", err

  do node = 1, nNodes
    diskAddDiv(node) = .false.
    xp0 = (/nodes(2,node), 0D0, nodes(1,node)/)
    if (xp0(3)>epsilTR) then
      buf_real = 0D0
      do nDir = 1, diskNElems
        do itr=1,3
          drot(itr,:) = disk_nodes(:,disk_elems(itr,nDir))
        end do
        bbb = bodyAng(xp0, drot(1,:), drot(2,:), drot(3,:))
        if (bbb>buf_real) buf_real = bbb
      end do
      if (buf_real>bodyAngle) diskAddDiv(node) = .true.
    end if
  end do

  write (*,*) 'disk_division file start'
  write (*,*) '************************'
  ffname = "./data/disk_division.vtk"
  open (unit=1,file=ffname)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', diskNNodes, ' float'

  do i=1,diskNNodes
    write (1,*) disk_nodes(1,i), disk_nodes(2,i), 0D0
  end do
  write (1,"(A,I7,A,I7)") 'CELLS ', diskNElems, ' ', diskNElems+3*diskNElems

  do i=1,diskNElems
    write (1,*) 3, disk_elems(1,i)-1, disk_elems(2,i)-1, disk_elems(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', diskNElems
  do i=1,diskNElems
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'CELL_DATA ', diskNElems
  write (1,"(A)") 'SCALARS Numb float'
  write (1,"(A)") 'LOOKUP_TABLE default'

  do i=1,diskNElems
    write(1,*) i
  end do

  close (unit=1,status='Keep')
  write (*,*) 'disk_division file has been created'
  write (*,*) '************************'
    !read *

    !===================
    ! расчет количества дополнительных направлений
  diskDivVol = pi*(r_disk**2)/addDirsMax
  allocate(buf_nums(1:1000000), stat = err)
  allocate(addDirsDivRange(1:nNodes,2), stat = err)
  addDirsDivN = 0
  bufcent = (/ 0D0, 0D0 /)

  do node=1,nNodes
    addDirsDivRange(node,1) = addDirsDivN + 1
    addDirsDivRange(node,2) = addDirsDivN
    xp0 = (/nodes(2,node), 0D0, nodes(1,node)/)
    nDir = nDirects
    if (xp0(3)>epsilTR) then
      do while (maxval(directnsND(directelems(1:3,nDir),3))>-0.1)!epsilTR)
        if(minval(directnsND(directelems(1:3,nDir),3))>epsilTR) then
          do itr=1,3
            buf_real = xp0(3)/directnsND(directelems(itr,nDir),3)
            xp(itr,:) = xp0(1:2) - buf_real * directnsND(directelems(itr,nDir),1:2)
          end do
          if (triInCirc(xp(1,:),xp(2,:),xp(3,:),bufcent,r_disk)) then!) then
            addDirsDivN = addDirsDivN +1
            buf_nums(addDirsDivN) = nDir
          end if
        else
           flag = .false.
           do itr=1,3
             if (directnsND(directelems(itr,nDir),3)>epsilTR) then
               buf_real = xp0(3)/directnsND(directelems(itr,nDir),3)
               xp(itr,:) = xp0(1:2) - buf_real * directnsND(directelems(itr,nDir),1:2)
               if (distance(xp(itr,:),bufcent)<r_disk-epsilTR) flag = .true.
             end if
           end do
           if (flag) then
            addDirsDivN = addDirsDivN +1
            buf_nums(addDirsDivN) = nDir
          end if
        end if
        nDir = nDir - 1
      end do
      addDirsDivRange(node,2) = addDirsDivN
    end if
    if (addDirsDivRange(node,1) > addDirsDivRange(node,2)) addDirsDivRange(node,:) = 0
  end do

  print *, "addDirsDivN =", addDirsDivN

  allocate(addDirsDiv(1:addDirsDivN), stat = err)
  print *, "addDirsDiv allocate stat=", err
  do itr = 1, addDirsDivN
    addDirsDiv(itr) = buf_nums(itr)
  end do

  deallocate(buf_nums)

!  call saveDirectsSphere(153)
!  read *

  do node=1,nNodes
    call formDirsAB(dirNAB1,dirNDnAB1,dirsAB1,bodyAngAB1,directelemsAB1,directnsNDAB1,toDiskAB1,node)
    nDirectsAB(node) = dirNAB1
  end do

  nDirectsUp = maxval(nDirectsAB(1:nNodes))

  ffname = "./data/addDirsDivN.vtk"

  open (unit=1,file=ffname)

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

  write (1,"(A,I7)") 'POINT_DATA ', nNodes
  write (1,"(A)") 'SCALARS addDirs float'
  write (1,"(A)") 'LOOKUP_TABLE default'

  do i=1,nNodes
    write(1,*) addDirsDivRange(i,2)-addDirsDivRange(i,1)
  end do

  write (1,"(A)") 'SCALARS nDirectsAB int'
  write (1,"(A)") 'LOOKUP_TABLE default'

  do i=1,nNodes
    write(1,*) nDirectsAB(i)
  end do

  close (unit=1,status='Keep')
  write (*,*) 'addDirsDivN file has been created'
  write (*,*) '************************'
end subroutine InitGroupsTriAB

subroutine IntensAllocateAB

  allocate(intensityNodes(nNodes), stat = err)
  print *, 'Nodal intensity allocate stat=', err
  allocate(Poynting(1:3,nNodes), stat = err)
  print *, 'Poynting vector allocate stat=', err
  allocate(radstress(1:3,1:3,nNodes), stat = err)
  print *, 'Radstress tensor allocate stat=', err
  allocate(intens_p(nElems), stat = err)
  print *, 'intens_p allocate stat=', err
  allocate(absorp_koef(nElems), stat = err)
  print *, 'absorp_koef allocate stat=', err

end subroutine IntensAllocateAB

function absorbtionAB(Uc) result (abcof)
  real (kind=precis), dimension(1:5) :: Uc
  real (kind=precis) :: abcof
  real (kind=precis) :: densy, temper, pressur

  densy = Uc(1)
  pressur = (gmm-1.)*(Uc(5)-0.5*sum(Uc(2:4)**2)/Uc(1))
  temper = pressur/densy
  abcof = (densy**2)/(temper**3.5)
end function absorbtionAB

subroutine radiation_solverTRAB
  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)
  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour

  call InitGroupsTriAB
  print *, 'InitGroupsTriAB done'

  call IntensAllocateAB
  print *, 'IntensAllocate done'

  call rad_param_calc_TRAB
  print *, 'rad_param_calc_TR done'

  call meshway_axial_AB

  print *, 'meshway done'

  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)

  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour - ddtime

  print *, 'Computation time (msec): ', ddtime

  call save_vtk_solution_radTRAB("-TRSN-"//i2c(nDirects))

!  read *
end subroutine radiation_solverTRAB

subroutine radiation_kinetics_TRAB
  integer :: node, elems

  real (kind=precis), allocatable :: divW(:)

  allocate(divW(nElems), stat = err)
  print *, 'absorp_koef allocate stat=', err

  do nElem=1,nElems
    divW(nElem) = 0D0
    do nEdge = 1,3

    end do
  end do
end subroutine radiation_kinetics_TRAB

subroutine rad_param_calc_TRAB
  real (kind=precis), dimension(2) :: centr=(/0.5,0.5/)
  real (kind=precis), dimension(2) :: centr1=(/0.5,0.0/)
  real (kind=precis), dimension(2) :: cnt, sec
  real (kind=precis) :: r0= 0.1
  integer :: nDir1
  integer :: seed=100
  real :: scat, scatsum

  do nElem=1,nElems
    absorp_koef(nElem) = 0.001_precis !absorb0*absorbtion(U_til(:,nElem))
!    if (distance(centr,elems_center(:,nElem))<0.1) absorp_koef(nElem) = 1000.
    intens_p(nElem) = 0.0_precis
  end do

end subroutine rad_param_calc_TRAB

subroutine save_vtk_solution_radTRAB(param)

  CHARACTER(LEN=*), INTENT(in):: param


  real (kind=precis) :: collecter, coll1
  integer :: ccc, node

  call filename_se_sol_vtk_radTR(niterDG,current_sol_name,param)

  open (unit=1,file=current_sol_name)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', nNodes, ' float'

  do i=1,nNodes
    write (1,*) nodes(2,i), 0, nodes(1,i)
  end do

  write (1,"(A,I7,A,I7)") 'CELLS ', nElems, ' ', nElems+3*nElems

  do i=1,nElems
    write (1,*) 3, elems(1,i)-1, elems(2,i)-1, elems(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', nElems
  do i=1,nElems
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'POINT_DATA ', nNodes

  write (1,"(A)") 'SCALARS Nodal_intensity float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nNodes
    write(1,*) intensityNodes(i)
  end do

  write (1,"(A)") 'SCALARS P_x float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nNodes
    write(1,*) Poynting(1,i)
  end do

  write (1,"(A)") 'SCALARS P_y float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nNodes
    write(1,*) Poynting(2,i)
  end do

  write (1,"(A)") 'SCALARS P_z float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do i=1,nNodes
    write(1,*) Poynting(3,i)
  end do

  write (1,"(A)") 'VECTORS Poynting float'
  do node=1,nNodes
    write(1,*) Poynting(:,node)
  end do

  write (1,"(A)") 'TENSORS radstress float'
  do node=1,nNodes
    do i=1,3
      write(1,*) radstress(i,:,node)
    end do
  end do

  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'

end subroutine save_vtk_solution_radTRAB

pure function bintens_axialAB(nod, direct) result(binten)
  real(kind=precis), dimension(1:2), Intent(In) :: nod
  real(kind=precis), dimension(1:3), Intent(In) :: direct
  real(kind=precis) :: binten, phhh
  real(kind=precis), dimension(2) :: centr,centr1

  binten = 0.0_precis
  centr1 = 0D0

  if (nod(2)<r_disk .and. nod(1)<0.05) binten = 1._precis - (distance(centr1, nod)/r_disk)**2
end function

subroutine meshway_axial_AB

  real (kind=precis), dimension(1:2) :: direc
  integer :: node, elem, buf_int
  real (kind=precis) :: collecter, bufr, bufi,  int_curr!, bang
  real (kind=precis), dimension(1:3) :: poynt_collect, dirn
  real (kind=precis), dimension(1:3,1:3) :: radstress_collect

  integer :: lentest = 1
  integer :: flag, stepin
  integer :: ktrac

  real (kind=precis), dimension(10000) :: waysLenp
  real (kind=precis), dimension(2) :: waysSrcNodep
  integer, dimension(10000) :: waysElemsp
  integer :: waysElemsCountp
  integer :: waysNDirsp
  integer :: is_src_node

  integer :: OMP_GET_THREAD_NUM

   integer :: control
   control = 0

  print *, "meshway_axial_AB start"
  print *, "==============="

  allocate(intensityTR(nNodes,fullDirs), stat = err)
  print *, 'allocate stat=', err

  do node=1,nNodes
    intensityNodes(node) = 0.0_precis
    Poynting(:,node) = 0.0_precis
  end do

  print *, nDirects
!$omp parallel
  ktrac =0
!$omp do SCHEDULE (GUIDED,100) private(nElem,node,collecter,stepin,nDir,direc,flag, bufr,bufi,int_curr,i,j,waysLenp,waysSrcNodep,waysElemsp,waysElemsCountp, waysNDirsp,dirNAB1, dirNDnAB1, dirsAB1, bodyAngAB1, directelemsAB1, directnsNDAB1, toDiskAB1,poynt_collect,dirn,radstress_collect)

  do node = 1,nNodes
    ktrac = ktrac + 1

    if (mod(ktrac,200)==0) then
      print *, "meshway_axial ktrac=",ktrac," node=",node
    end if

    call formDirsAB(dirNAB1,dirNDnAB1,dirsAB1,bodyAngAB1,directelemsAB1,directnsNDAB1,toDiskAB1,node)

    collecter = 0.0_precis
    poynt_collect = 0.0_precis
    radstress_collect = 0.0_precis

    do nDir = 1,dirNAB1
! if (markers .and. node>control)   print *, node, nDir

      dirn = dirsAB1(nDir,:)
      waysElemsCountp = 0
      waysSrcNodep = (/z_max,r_max/)
      flag = 0
      if (is_bc_node(node)<0.5) then
! if (markers .and. node>control)  print *, '1a', nDir, node
! if (markers .and. node>control)  print *, nodes(:, node)
! if (markers .and. node>control)  print *, dirsAB1(nDir,:)

        if (toDiskAB1(nDir) .or. notDiskOnly) then
          call tracing_axial_par(node, dirn, waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
          if (lentest<waysElemsCountp) lentest=waysElemsCountp
        endif
! if (markers .and. node>control)  print *, 'b', nDir, node

      else
! if (markers .and. node>control)  print *, '2a', nDir, node

        is_src_node = bound_source_node(node, dirn)
        if (is_src_node>0.5 .or. is_src_node<-0.5) then
          if (is_src_node<-0.5) flag = 1
          waysElemsCountp = 0
          waysSrcNodep(:) = nodes(:,node)
! if (markers .and. node>control)  print *, 'b', nDir, node

        else
! if (markers .and. node>control)  print *, '3a', nDir, node, is_bc_node(node)
! if (markers .and. node>control)  print *, nodes(:, node)
! if (markers .and. node>control)  print *, dirsAB1(nDir,:)

          if (toDiskAB1(nDir) .or. notDiskOnly) then
            call tracing_axial_par(node, dirn, waysSrcNodep, waysElemsCountp, waysElemsp, waysLenp)
            if (lentest<waysElemsCountp) lentest=waysElemsCountp
          end if
! if (markers .and. node>control)  print *, 'b', nDir

        end if
      endif
! if (markers .and. node>control)  print *, '4a', nDir
      int_curr = bintens_axialAB(waysSrcNodep, dirn)
! if (markers .and. node>control)  print *, '5a', nDir,waysElemsCountp
      do stepin=waysElemsCountp,1,-1
          int_curr=intty_step(int_curr, absorp_koef(waysElemsp(stepin)), absorp_koef(waysElemsp(stepin))*intens_p(waysElemsp(stepin)), waysLenp(stepin))
      end do
! if (markers .and. node>control)  print *, '6a', nDir

      intensityTR(node,nDir) = int_curr
      if (flag>0.5) intensityTR(node,nDir) = 0.0
! if (markers .and. node>control)  print *, '7a', nDir

      collecter = collecter + intensityTR(node,nDir)*bodyAngAB1(nDir)
      poynt_collect = poynt_collect + intensityTR(node,nDir)*bodyAngAB1(nDir)*dirn
      dirn(2) = -dirn(2)
      collecter = collecter + intensityTR(node,nDir)*bodyAngAB1(nDir)
      poynt_collect = poynt_collect + intensityTR(node,nDir)*bodyAngAB1(nDir)*dirn
! if (markers .and. node>control)  print *, '8a', nDir

      do i=1,3
        do j=1,3
          radstress_collect(i,j) = radstress_collect(i,j) + intensityTR(node,nDir)*bodyAngAB1(nDir)*dirn(i)*dirn(j)
        end do
      end do
! if (markers .and. node>control)  print *, '9a', nDir

    end do
    intensityNodes(node) = collecter
    Poynting(:,node) = poynt_collect

!!! поделить на безразмерную скорость света!
    radstress(:,:,node) = radstress_collect
  end do
!$omp end do
!$omp end parallel
  print *, lentest, "/", nElems_b

  print *, "meshway_axial_AB finish"
  print *, "==============="

end subroutine meshway_axial_AB

pure function bound_source_node(noden, dir) result (bsn)
  integer, Intent(In) :: noden
  real(kind=precis), dimension(1:3), Intent(In) :: dir

  integer :: bsn

  real (kind=precis), dimension(1:2) :: snode

  snode = nodes(:,noden)

  bsn = 0
  if (snode(2)<epsilTR .and. (snode(1)>epsilTR .and. snode(1)<z_max-epsilTR)) then
    bsn = 0
  else
    if (snode(1)<epsilTR .and. dir(3)>epsilTR) bsn = 1
    if (snode(1)>z_max-epsilTR .and. dir(3)<-epsilTR) bsn = 1
    if (snode(2)>r_max-epsilTR .and. dir(1)<-epsilTR) bsn = 1

    if ((snode(2)>r_max-epsilTR .and. snode(1)<epsilTR) .or. (snode(2)>r_max-epsilTR .and. snode(1)>z_max-epsilTR)) then
          bsn = 0
      if (snode(1)>z_max-epsilTR .and. (dir(1)<-epsilTR .and. dir(3)<-epsilTR)) bsn = 1
      if (snode(1)<epsilTR .and. (dir(1)<-epsilTR .and. dir(3)>epsilTR)) bsn = 1
      if (snode(1)>z_max-epsilTR .and. (dir(1)*dir(3)<-epsilTR)) bsn = -1
      if (snode(1)<epsilTR .and. (dir(1)*dir(3)>epsilTR)) bsn = -1
    end if
  end if

end function

pure subroutine formDirsAB(dirNAB,dirNDnAB,dirsAB,bodyAngAB,directelemsAB,directnsNDAB,toDiskAB,nnode)
  integer, Intent(InOut) :: dirNAB, dirNDnAB
  real (kind=precis), dimension(1:fullDirs,1:3), Intent(InOut) :: dirsAB
  real (kind=precis), dimension(1:fullDirs), Intent(InOut) :: bodyAngAB
  real (kind=precis), dimension(1:fullNodes,1:3), Intent(InOut) :: directnsNDAB
  integer, dimension(1:3,1:fullDirs), Intent(InOut) :: directelemsAB
  logical, dimension(1:fullDirs), Intent(InOut) :: toDiskAB

  integer, Intent(In) :: nnode

  integer :: nDir,itrac,jtrac,ktrac,tst,strtdsk
  real (kind=precis), dimension(1:2) :: bufcent
  real(kind=precis) :: buf_real
  real(kind=precis), dimension(1:3) :: buf_vec, buf_cent, x0

  integer,dimension(1:fullDirs) :: elemsToDivide, newElemsToDivide
  integer :: nElemsToDivide, newNElemsToDivide, resElems
  real(kind=precis), dimension(1:3,1:3) :: addInCirc, xp, xpp

  logical :: flag

  xp = 0D0
  addInCirc = 0D0
  buf_cent = (/0D0,0D0,0D0/)
  bufcent =  (/0D0,0D0/)

  x0 = (/ nodes(2,nnode), 0D0, nodes(1,nnode) /)
  toDiskAB = .false.

! копируем уже известные треугольники
  dirNAB = nDirects
  do nDir=1,dirNAB
    dirsAB(nDir,1:3) = directns(nDir,1:3)
    bodyAngAB(nDir) = bodyAngleTR(nDir)
    directelemsAB(:,nDir) = directelems(:,nDir)
  end do

  dirNDnAB = nDirectNodes
  do itrac=1,dirNDnAB
    directnsNDAB(itrac,:) = directnsND(itrac,:)
  end do

  nElemsToDivide = addDirsDivRange(nnode,2) - addDirsDivRange(nnode,1) + 1
  if (addDirsDivRange(nnode,2)==0) nElemsToDivide = 0
  do nDir = addDirsDivRange(nnode,1), addDirsDivRange(nnode,2)
    elemsToDivide(nDir-addDirsDivRange(nnode,1)+1) = addDirsDiv(nDir)
  end do

  do while (nElemsToDivide>0)
    newNElemsToDivide = 0
    do jtrac=1,nElemsToDivide
      nDir = elemsToDivide(jtrac)
      if (minval(directnsNDAB(directelemsAB(1:3,nDir),3))>epsilTR) then
! если треугольник над экватором
        xp = 0D0
        do itrac=1,3
          buf_real = x0(3)/directnsNDAB(directelemsAB(itrac,nDir),3)
          xp(itrac,1:2) = x0(1:2) - buf_real * directnsNDAB(directelemsAB(itrac,nDir),1:2)
        end do
        if (.not.fullInCirc(xp(1,1:2),xp(2,1:2),xp(3,1:2),bufcent,r_disk)) then
          do itrac=1,3
! разбиваем треугольник на диске
            addInCirc(itrac,:) = 0.5*(xp(itrac,:)+xp(Q(itrac),:))
! добавляем точки-направления
            dirNDnAB = dirNDnAB + 1
            buf_vec = x0 - addInCirc(itrac,:)
            directnsNDAB(dirNDnAB,:) = buf_vec/norm2_3(buf_vec)
          end do
! добавляем угловые треугольники
          do itrac=1,3
! добавляем элемент
            dirNAB = dirNAB + 1
            directelemsAB(:,dirNAB) = (/ directelemsAB(itrac,nDir), dirNDnAB-3+itrac, dirNDnAB-3+Q(Q(itrac)) /)
! добавляем направление и телесный угол
            dirsAB(dirNAB,:) = (directnsNDAB(directelemsAB(1,dirNAB),:)+directnsNDAB(directelemsAB(2,dirNAB),:)+directnsNDAB(directelemsAB(3,dirNAB),:))/3D0
            dirsAB(dirNAB,:) = dirsAB(dirNAB,:)/norm2_3(dirsAB(dirNAB,:))
            bodyAngAB(dirNAB) = bodyAng(x0,xp(itrac,:), addInCirc(itrac,:), addInCirc(Q(Q(itrac)),:))
! определяем, разбивать ли его при следующем заходе
            if (triInCirc(xp(itrac,1:2),addInCirc(Q(itrac),1:2),addInCirc(itrac,1:2),bufcent,r_disk)  .and. bodyAngAB(dirNAB) > minDirMax*bodyAngle) then
            newNElemsToDivide = newNElemsToDivide + 1
              newElemsToDivide(newNElemsToDivide) = dirNAB
            end if
          end do
! обновляем центральный треугольник
! добавляем элемент
          directelemsAB(:,nDir) = (/ dirNDnAB-2, dirNDnAB-1, dirNDnAB /)
! добавляем направление и телесный угол
          dirsAB(nDir,:) = (directnsNDAB(dirNDnAB-2,:)+directnsNDAB(dirNDnAB-1,:)+directnsNDAB(dirNDnAB,:))/3D0
          dirsAB(nDir,:) = dirsAB(nDir,:)/norm2_3(dirsAB(nDir,:))
          bodyAngAB(nDir) = bodyAngAB(nDir) - (bodyAngAB(dirNAB-2)+bodyAngAB(dirNAB-1)+bodyAngAB(dirNAB))
! определяем, разбивать ли его при следующем заходе
          if (triInCirc(addInCirc(1,1:2),addInCirc(2,1:2),addInCirc(3,1:2),bufcent,r_disk) .and. bodyAngAB(nDir) > minDirMax*bodyAngle) then
          newNElemsToDivide = newNElemsToDivide + 1
            newElemsToDivide(newNElemsToDivide) = nDir
          end if
        end if
      else
!
! если треугольник у экватора
        xp = 0D0
        do itrac=1,3
          xp(itrac,:) = directnsNDAB(directelemsAB(itrac,nDir),:)
        end do
        do itrac=1,3
! разбиваем треугольник на диске
          addInCirc(itrac,:) = 0.5*(xp(itrac,:)+xp(Q(itrac),:))
! добавляем точки-направления
          dirNDnAB = dirNDnAB + 1
          buf_vec = addInCirc(itrac,:)
          directnsNDAB(dirNDnAB,:) = buf_vec/norm2_3(buf_vec)
        end do
! добавляем угловые треугольники
        do itrac=1,3
! добавляем элемент
          dirNAB = dirNAB + 1
          directelemsAB(:,dirNAB) = (/ directelemsAB(itrac,nDir), dirNDnAB-3+itrac, dirNDnAB-3+Q(Q(itrac)) /)
! добавляем направление и телесный угол
          dirsAB(dirNAB,:) = (directnsNDAB(directelemsAB(1,dirNAB),:)+directnsNDAB(directelemsAB(2,dirNAB),:)+directnsNDAB(directelemsAB(3,dirNAB),:))/3D0
          dirsAB(dirNAB,:) = dirsAB(dirNAB,:)/norm2_3(dirsAB(dirNAB,:))
          bodyAngAB(dirNAB) = bodyAng(buf_cent,xp(itrac,:), addInCirc(itrac,:), addInCirc(Q(Q(itrac)),:))
! определяем, разбивать ли его при следующем заходе
          if (bodyAngAB(dirNAB) > minDirMax*bodyAngle) then
            if (minval(directnsNDAB(directelemsAB(1:3,dirNAB),3))>epsilTR) then
              newNElemsToDivide = newNElemsToDivide + 1
              newElemsToDivide(newNElemsToDivide) = dirNAB
            else
              flag = .false.
              do ktrac=1,3
                xpp = 0D0
                if (directnsNDAB(directelemsAB(ktrac,dirNAB),3)>epsilTR) then
                  buf_real = x0(3)/directnsNDAB(directelemsAB(ktrac,dirNAB),3)
                  xpp(ktrac,1:2) = x0(1:2) - buf_real * directnsNDAB(directelemsAB(ktrac,dirNAB),1:2)
                  if (distance(xpp(ktrac,1:2),bufcent)<r_disk-epsilTR) flag = .true.
                end if
              end do
              if (flag) then
                newNElemsToDivide = newNElemsToDivide + 1
                newElemsToDivide(newNElemsToDivide) = dirNAB
              end if
            end if
          end if
        end do
! обновляем центральный треугольник
! добавляем элемент
        directelemsAB(:,nDir) = (/ dirNDnAB-2, dirNDnAB-1, dirNDnAB /)
! добавляем направление и телесный угол
        dirsAB(nDir,:) = (directnsNDAB(dirNDnAB-2,:)+directnsNDAB(dirNDnAB-1,:)+directnsNDAB(dirNDnAB,:))/3D0
        dirsAB(nDir,:) = dirsAB(nDir,:)/norm2_3(dirsAB(nDir,:))
        bodyAngAB(nDir) = bodyAngAB(nDir) - (bodyAngAB(dirNAB-2)+bodyAngAB(dirNAB-1)+bodyAngAB(dirNAB))
 ! определяем, разбивать ли его при следующем заходе
        if (bodyAngAB(nDir) > minDirMax*bodyAngle) then
          if (minval(directnsNDAB(directelemsAB(1:3,nDir),3))>epsilTR) then
            newNElemsToDivide = newNElemsToDivide + 1
            newElemsToDivide(newNElemsToDivide) = nDir
         else
           flag = .false.
           do ktrac=1,3
             xpp = 0D0
             if (directnsNDAB(directelemsAB(ktrac,nDir),3)>epsilTR) then
               buf_real = x0(3)/directnsNDAB(directelemsAB(ktrac,nDir),3)
               xpp(ktrac,1:2) = x0(1:2) - buf_real * directnsNDAB(directelemsAB(ktrac,nDir),1:2)
               if (distance(xpp(ktrac,1:2),bufcent)<r_disk-epsilTR) flag = .true.
             end if
           end do
           if (flag) then
             newNElemsToDivide = newNElemsToDivide + 1
             newElemsToDivide(newNElemsToDivide) = nDir
           end if
          end if
        end if
      endif
    end do
    nElemsToDivide = newNElemsToDivide
    do itrac = 1, newNElemsToDivide
      elemsToDivide(itrac) = newElemsToDivide(itrac)
    end do
  end do

  if (x0(3)>epsilTR) then
    resElems = 0
    do nDir=1,dirNAB
      flag = .false.
      if (dirsAB(nDir,3)>epsilTR) then
        buf_real = x0(3)/dirsAB(nDir,3)
        buf_vec(:) =  x0(:) - buf_real * dirsAB(nDir,:)
        if (distance(buf_vec(1:2),bufcent)>r_disk) flag = .true.
      else
        flag = .true.
      end if
      if (flag) then
        resElems = resElems + 1
        dirsAB(resElems,1:3) = dirsAB(nDir,1:3)
        bodyAngAB(resElems) = bodyAngAB(nDir)
        directelemsAB(:,resElems) = directelemsAB(:,nDir)
      end if
    end do
    dirNAB = resElems

    do nDir=1,diskNNodes
      dirNDnAB = dirNDnAB + 1
      buf_vec = x0 - disk_nodes(:,nDir)
      directnsNDAB(dirNDnAB,:) = buf_vec/norm2_3(buf_vec)
    end do

    nElemsToDivide = 0
    do nDir = 1, diskNElems
      dirNAB = dirNAB +1
      buf_vec =  x0 - disk_srcnodes(nDir,1:3)
      dirsAB(dirNAB,1:3) = buf_vec/norm2_3(buf_vec)
      bodyAngAB(dirNAB) = bodyAng(x0,disk_nodes(:,disk_elems(1,nDir)),disk_nodes(:,disk_elems(2,nDir)),disk_nodes(:,disk_elems(3,nDir)))
      directelemsAB(:,dirNAB) = disk_elems(:,nDir)+dirNDnAB-diskNNodes
      toDiskAB(dirNAB) = .true.
      if (diskAddDiv(nnode) .and. bodyAngAB(dirNAB)>bodyAngle) then
        nElemsToDivide = nElemsToDivide + 1
        elemsToDivide(nElemsToDivide) = dirNAB
      end if
    end do

    do while (nElemsToDivide>0)
      newNElemsToDivide = 0
      do jtrac=1,nElemsToDivide
        nDir = elemsToDivide(jtrac)
        xp = 0D0
        do itrac=1,3
          buf_real = x0(3)/directnsNDAB(directelemsAB(itrac,nDir),3)
          xp(itrac,1:2) = x0(1:2) - buf_real * directnsNDAB(directelemsAB(itrac,nDir),1:2)
        end do
        do itrac=1,3
! разбиваем треугольник на диске
          addInCirc(itrac,:) = 0.5*(xp(itrac,:)+xp(Q(itrac),:))
! добавляем точки-направления
          dirNDnAB = dirNDnAB + 1
          buf_vec = x0 - addInCirc(itrac,:)
          directnsNDAB(dirNDnAB,:) = buf_vec/norm2_3(buf_vec)
        end do
! добавляем угловые треугольники
        do itrac=1,3
! добавляем элемент
          dirNAB = dirNAB + 1
!           print *, "dirNAB",dirNAB
          directelemsAB(:,dirNAB) = (/ directelemsAB(itrac,nDir), dirNDnAB-3+itrac, dirNDnAB-3+Q(Q(itrac)) /)
! добавляем направление и телесный угол
          dirsAB(dirNAB,:) = (directnsNDAB(directelemsAB(1,dirNAB),:)+directnsNDAB(directelemsAB(2,dirNAB),:)+directnsNDAB(directelemsAB(3,dirNAB),:))/3D0
          dirsAB(dirNAB,:) = dirsAB(dirNAB,:)/norm2_3(dirsAB(dirNAB,:))
          bodyAngAB(dirNAB) = bodyAng(x0,xp(itrac,:), addInCirc(itrac,:), addInCirc(Q(Q(itrac)),:))
          toDiskAB(dirNAB) = .true.
! определяем, разбивать ли его при следующем заходе
          if (bodyAngAB(dirNAB)>bodyAngle) then
            newNElemsToDivide = newNElemsToDivide + 1
            newElemsToDivide(newNElemsToDivide) = dirNAB
          end if
        end do
! обновляем центральный треугольник
! добавляем элемент
        directelemsAB(:,nDir) = (/ dirNDnAB-2, dirNDnAB-1, dirNDnAB /)
! добавляем направление и телесный угол
        dirsAB(nDir,:) = (directnsNDAB(dirNDnAB-2,:)+directnsNDAB(dirNDnAB-1,:)+directnsNDAB(dirNDnAB,:))/3D0
        dirsAB(nDir,:) = dirsAB(nDir,:)/norm2_3(dirsAB(nDir,:))
        bodyAngAB(nDir) = bodyAngAB(nDir) - (bodyAngAB(dirNAB-2)+bodyAngAB(dirNAB-1)+bodyAngAB(dirNAB))
! определяем, разбивать ли его при следующем заходе
        if (bodyAngAB(nDir) > bodyAngle) then
          newNElemsToDivide = newNElemsToDivide + 1
          newElemsToDivide(newNElemsToDivide) = nDir
        end if
      end do
      nElemsToDivide = newNElemsToDivide
      do itrac = 1, newNElemsToDivide
        elemsToDivide(itrac) = newElemsToDivide(itrac)
      end do
    end do
  end if

end subroutine

subroutine saveDirectsSphere(nnode)

  integer, Intent(In) :: nnode
  CHARACTER(LEN = max_filename_length):: ffname
  real(kind=precis), dimension(1:3) :: buf_vec, xp0
  integer :: node

  !call formDirectsAB(dirNAB1,dirNDnAB1,dirsAB1,bodyAngAB1,directelemsAB1,directnsNDAB1,node)
  call formDirsAB(dirNAB1,dirNDnAB1,dirsAB1,bodyAngAB1,directelemsAB1,directnsNDAB1,toDiskAB1,nnode)

  print *, diskAddDiv(nnode)
  print *, dirNAB1, dirNDnAB1, fullDirs
  write (*,*) 'direcsphere file start'
  write (*,*) '************************'
  ffname = "./data/direcsphere.vtk"
  xp0 = (/nodes(2,nnode), 0D0, nodes(1,nnode)/)
  print *, "xp0=", xp0
  open (unit=1,file=ffname)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A)") 'Unstructured Grid Example'
  write(1,"(A)") 'ASCII'
  write(1,*) ''
  write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1,"(A,I7,A)") 'POINTS ', dirNDnAB1, ' float'

  do i=1,dirNDnAB1
    buf_vec = xp0
! !buf_vec = (/0D0,0D0,0D0/)
!        if (directnsNDAB1(i,3)>0.05) then
!              buf_real = xp0(3)/directnsNDAB1(i,3)
!              buf_vec(:) = xp0(:) - buf_real * directnsNDAB1(i,:)
!  end if
!        write (1,*) buf_vec(1), buf_vec(2), buf_vec(3)
    write (1,*) directnsNDAB1(i,1), directnsNDAB1(i,2), directnsNDAB1(i,3)
  end do
  write (1,"(A,I7,A,I7)") 'CELLS ', dirNAB1, ' ', dirNAB1+3*dirNAB1

  do i=1,dirNAB1
    write (1,*) 3, directelemsAB1(1,i)-1, directelemsAB1(2,i)-1, directelemsAB1(3,i)-1
  end do

  write (1,"(A,I7)") 'CELL_TYPES ', dirNAB1
  do i=1,dirNAB1
    write (1,*) 5
  end do

  write (1,"(A,I7)") 'CELL_DATA ', dirNAB1
  write (1,"(A)") 'SCALARS bAng float'
  write (1,"(A)") 'LOOKUP_TABLE default'

  do i=1,dirNAB1
    write(1,*) bodyAngAB1(i)
  end do

  write (1,"(A)") 'VECTORS dirs float'

  do i=1,dirNAB1
    write(1,*) dirsAB1(i,:)
  end do


  close (unit=1,status='Keep')
  write (*,*) 'direcsphere file has been created'
  write (*,*) '************************'
    read *
end subroutine

pure function triInCirc(x1,x2,x3,x0,radius) result (res)
  real(kind=precis), dimension(1:2), intent(In) :: x1,x2,x3,x0
  real(kind=precis), Intent(In) :: radius
  logical :: res, cond1,cond2,cond3,cond4

  cond1 = distance(x1,x0)<radius .or. distance(x2,x0)<radius .or. distance(x3,x0)<radius

  cond2 = abs(vect_mult2d(x2-x0,x1-x0))/distance(x1,x2)<radius .and. abs(vect_mult2d(x2-x0,x1-x0))>epsilTR .and. dot_product(x1-x0,x2-x0)>epsilTR .and. dot_product(x1-x2,x0-x2)>epsilTR .and. dot_product(x0-x1,x2-x1)>epsilTR

  cond3 = abs(vect_mult2d(x3-x0,x1-x0))/distance(x1,x3)<radius .and. abs(vect_mult2d(x3-x0,x1-x0))>epsilTR .and. dot_product(x1-x0,x3-x0)>epsilTR .and. dot_product(x1-x3,x0-x3)>epsilTR .and. dot_product(x0-x1,x3-x1)>epsilTR

  cond4 = abs(vect_mult2d(x2-x0,x3-x0))/distance(x3,x2)<radius .and. abs(vect_mult2d(x2-x0,x3-x0))>epsilTR .and. dot_product(x3-x0,x2-x0)>epsilTR .and. dot_product(x3-x2,x0-x2)>epsilTR .and. dot_product(x0-x3,x2-x3)>epsilTR

  res = cond1 .or. cond2 .or. cond3 .or. cond4
end function

pure function fullInCirc(x1,x2,x3,x0,radius) result (res)
  real(kind=precis), dimension(1:2), intent(In) :: x1,x2,x3,x0
  real(kind=precis), Intent(In) :: radius
  logical :: res

  res = distance(x1,x0)<radius .and. distance(x2,x0)<radius .and. distance(x3,x0)<radius
end function

end module radiationTRAB
