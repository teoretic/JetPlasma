! ==========================================================================
! Модуль для массивов, описывающих решение

! ==========================================================================



module solution
  use names
!  use initials

  implicit none

contains

subroutine save_solution
    implicit none

    integer::I,K, nElem

    real (kind=precis):: rho, ux, v, w, E, p, k_e, sounds, Mach




	 call filename_se_sol(ntime,current_sol_name)

!	if (MOD (ntime, 10000).eq.0) then

	 open (unit=1,file=current_sol_name,form='formatted')

	 write (1,*) 'TITLE = "U graph t = ',t_current,'"'
	 if (flag_axial==1) then
			write (1,*) 'VARIABLES = "Z", "R", "rho", "Vz", "Vr","Vphi", "E", "Bz","Br","Bphi","p","M", "B2/2", "omega"'
		else
			write (1,*) 'VARIABLES = "X", "Y", "rho", "u", "v","w", "E", "Bx","By","Bz","p","M", "B2/2"'
	 end if
	 write (1,*) 'ZONE N=',nNodes, ',E=',nElems
	 write (1,*) 'DATAPACKING = BLOCK'
     write (1,*) 'VARLOCATION=([3,4,5,6,7,8,9,10,11,12,13,14]=CELLCENTERED),ZONETYPE=FETRIANGLE'



    do K=1,nNodes
      write(1,*)  nodes(1,K)
    end do
    do K=1,nNodes
      write(1,*)  nodes(2,K)
    end do

    do K=1,5
      do I=1,nElems
        if (K>1 .and. K<5) then
	      write(1,*)  U(K,I)/U(1,I)
	    else
          write(1,*)  U(K,I)
	    end if
      end do
    end do

    do K=1,3
      do I=1,nElems
	write(1,*) B(K,I)/b_scale
      end do
    end do

    do I=1,nElems
   ! вычисляем все простые переменные
      rho = U(1,I)
      ux = U(2,I)/rho
      v = U(3,I)/rho
      w = U(4,I)/rho									! полная энергия (внутренняя плюс кинетическая)
      E = U(5,I)
      k_e = ux*ux+v*v+w*w
      p = (E-0.5*rho*k_e)*(gmm-1.0)

      write(1,*) p
    end do

    do I=1,nElems
   ! вычисляем все простые переменные
      rho = U(1,I)
      ux = U(2,I)/rho
      v = U(3,I)/rho
      w = U(4,I)/rho									! полная энергия (внутренняя плюс кинетическая)
      E = U(5,I)
      k_e = ux*ux+v*v+w*w
      p = (E-0.5*rho*k_e)*(gmm-1.0)
      sounds = sqrt(gmm*p/rho)
      Mach = sqrt(ux*ux+v*v+w*w)/sounds
      write(1,*) Mach
    end do

    do I=1,nElems
      write(1,*) sum(B(:,I)**2)/(2.0*(b_scale**2))
    end do

    if (flag_axial==1) then
      do I=1,nElems
	write(1,*) U(4,I)/(U(1,I)*elems_center(2,I))
      end do
    end if

    write(1,*) ''

    do nElem=1,nElems
      write (1,*) elems(1,nElem),elems(2,nElem),elems(3,nElem)
    end do

    close (unit=1,status='Keep')
    write (*,*) 'File ',current_sol_name,' has been created'
	write (*,*) 't=',t_current
    write (*,*) '************************'
end subroutine save_solution

subroutine save_vtk_solution
  implicit none

  integer::I,K, j, nElem
  real (kind=precis):: rho, ux, v, w, Ek, p, k_e, sounds, Mach, radius, magn, vorticity

  call filename_se_sol_vtk(ntime,current_sol_name)

  open (unit=1,file=current_sol_name)

  write(1,"(A)") '# vtk DataFile Version 3.0'
  write(1,"(A,f6.3)") 'RadPlasma output, t= ',t_current
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

  write (1,"(A,I7)") 'CELL_DATA ', nElems

  write (1,"(A)") 'SCALARS density float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    write(1,*)  U(1,I)
  end do

  write (1,"(A)") 'VECTORS velocity float'
  do i=1,nElems
    write(1,*) U(3,i)/U(1,i), 0D0, U(2,i)/U(1,i)
  end do

  write (1,"(A)") 'SCALARS v_phi float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    write(1,*)  U(4,i)/U(1,i)
  end do

  write (1,"(A)") 'SCALARS omega float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    radius = sum(nodes(2, elems(1:3,i)))/3D0
    write(1,*) U(4,I)/(U(1,I)*radius)
  end do

  write (1,"(A)") 'SCALARS energy float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    write(1,*)  U(5,I)
  end do

  write (1,"(A)") 'SCALARS pressure float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    rho = U(1,I)
    ux = U(2,I)/rho
    v = U(3,I)/rho
    w = U(4,I)/rho
    Ek = U(5,I)
    k_e = ux*ux+v*v+w*w
    p = (Ek-0.5*rho*k_e)*(gmm-1.0)

    write(1,*) p
  end do

  write (1,"(A)") 'SCALARS Mach float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    rho = U(1,I)
    ux = U(2,I)/rho
    v = U(3,I)/rho
    w = U(4,I)/rho
    Ek = U(5,I)
    k_e = ux*ux+v*v+w*w
    p = (Ek-0.5*rho*k_e)*(gmm-1.0)
    sounds = sqrt(gmm*p/rho)
    Mach = sqrt(ux*ux+v*v+w*w)/sounds
    write(1,*) Mach
  end do

  write (1,"(A)") 'SCALARS Alfven_Mach float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    rho = U(1,I)
    magn = sqrt(sum(B(:,i)**2))
    ux = U(2,I)/rho
    v = U(3,I)/rho
    w = U(4,I)/rho
    sounds = magn/sqrt(4D0*pi*rho)
    Mach = 0D0
    if (sounds>1E-1) Mach = sqrt(ux*ux+v*v+w*w)/sounds
    write(1,*) Mach
  end do

  write (1,"(A)") 'VECTORS B float'
  do i=1,nElems
    write(1,*) B(2,i), 0D0, B(1,i)
  end do

  write (1,"(A)") 'SCALARS Bz float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    write(1,*) B(1,i)
  end do

  write (1,"(A)") 'SCALARS B_phi float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    write(1,*) B(3,i)
  end do

  write (1,"(A)") 'SCALARS Vorticity float'
  write (1,"(A)") 'LOOKUP_TABLE default'
  do I=1,nElems
    vorticity = 0D0
    do k=1,3
      j=elems_elems(k,i)
!   r=x          phi      z
!    2            3       1
! U(3,i)/U(1,i), 0D0, U(2,i)/U(1,i)
      v = 0D0
      if (U(1,j)>1D-7) then
        v = e(1,i,k)*(U(2,i)/U(1,i)+U(2,j)/U(1,j))+e(2,i,k)*(U(3,i)/U(1,i)+U(3,j)/U(1,j))
      else
        v = e(1,i,k)*(U(2,i)/U(1,i)+U(2,i)/U(1,i))+e(2,i,k)*(U(3,i)/U(1,i)+U(3,i)/U(1,i))
      end if
      v = v*0.5
      vorticity = vorticity + v
    end do
    vorticity = vorticity/elems_vol(i)
    write(1,*) vorticity
  end do


!   do nElem=1,nElems
!     circ = 0.0
!     do nEdge=1,3
!       A = edge_normSign(nEdge,nElem)
!       circ = circ + B_edges_hat(3,elem_edges(nEdge,nElem))*edge_vol(elem_edges(nEdge,nElem))*A
!     end do
!     A = 1.0/elems_vol(nElem)
!     rotB(3,nElem) = A*circ
!   end do

! !   if (flag_radiative>0.5) then
! !     write (1,"(A)") 'SCALARS absorp_koef float'
! !     write (1,"(A)") 'LOOKUP_TABLE default'
! !     do i=1,nElems
! !       write(1,*) absorp_koef(i)
! !     end do
! !
! !     write (1,"(A)") 'SCALARS intensity float'
! !     write (1,"(A)") 'LOOKUP_TABLE default'
! !     do i=1,nElems
! !       write(1,*) intensity(i)
! !     end do
! !
! !     write (1,"(A)") 'VECTORS poynting float'
! !     do i=1,nElems
! !       write(1,*) Poynting(1,i), Poynting(2,i), 0.! Poynting(3,i)
! !     end do
! !
! !     if (flag_debug>0.5) then
! !       write (1,"(A)") 'VECTORS radimpulse float'
! !       do i=1,nElems
! ! 	write(1,*) radimpulse(1,i), radimpulse(2,i), 0.!radimpulse(3,i)
! !       end do
! !       write (1,"(A)") 'SCALARS radenerg float'
! !       write (1,"(A)") 'LOOKUP_TABLE default'
! !       do i=1,nElems
! ! 	write(1,*) radenerg(i)
! !       end do
! !     end if
! !
! !   end if

  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'
end subroutine save_vtk_solution

subroutine save_control_point
  implicit none

  integer::I,K, nElem
  real (kind=precis):: rho, ux, v, w, E, p, k_e, sounds, Mach

  current_sol_name = root // '/cp-data' // extention

  open (unit=1,file=current_sol_name)

  write (1,*) t_current
  write (1,*) tau
  write (1,*) ntime

  do K=1,5
    write(1,*) U(K,:)
  end do

  do K=1,3
    write(1,*) B(K,:)
  end do

  do I=1,nEdges
    write (1,*) B_edges(:,I)
  end do

  do i=1,nNodes
    write (1,*) B_nodes(:,i)
  end do

  do I=1,nEdges
    write (1,*) V_edges(:,I)
  end do

  do i=1,nNodes
    write (1,*) V_nodes(:,i)
  end do

  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'
end subroutine save_control_point

function read_control_point result(success)
  implicit none

  integer::I,K, nElem
  real (kind=precis):: rho, ux, v, w, E, p, k_e, sounds, Mach
  logical :: success

  current_sol_name = root // '/cp-data' // extention

  inquire(file=current_sol_name, EXIST=success)

  if (success) then
    open (unit=1,file=current_sol_name)

    read (1,*) t_current
    read (1,*) tau
    read (1,*) ntime

    do K=1,5
      read (1,*) U(K,:)
    end do

    do K=1,3
      read (1,*) B(K,:)
    end do

    do I=1,nEdges
      read (1,*) B_edges(:,I)
    end do

    do i=1,nNodes
      read (1,*) B_nodes(:,i)
    end do

    do I=1,nEdges
      read (1,*) V_edges(:,I)
    end do

    do i=1,nNodes
      read (1,*) V_nodes(:,i)
    end do

    close (unit=1,status='Keep')
    write (*,*) 'File ',current_sol_name,' has been read'
    write (*,*) 't=',t_current
    write (*,*) '************************'
  end if
end function read_control_point

subroutine save_reservation
    implicit none

    integer::I,K, nElem

    real (kind=precis):: rho, ux, v, w, E, p, k_e, sounds, Mach




!	 call filename_se_sol(ntime,current_sol_name)

!	if (MOD (ntime, 10000).eq.0) then
	 current_sol_name = root // '/reserve' // extention

	 open (unit=1,file=current_sol_name,form='formatted')

	 write (1,*) 'TITLE = "U graph t = ',t_current,'"'
	 if (flag_axial==1) then
			write (1,*) 'VARIABLES = "Z", "R", "rho", "Vz", "Vr","Vphi", "E", "Bz","Br","Bphi","p","M", "B2/8pi", "omega"'
		else
			write (1,*) 'VARIABLES = "X", "Y", "rho", "u", "v","w", "E", "Bx","By","Bz","p","M", "B2/8pi"'
	 end if
	 write (1,*) 'ZONE N=',nNodes, ',E=',nElems
	 write (1,*) 'DATAPACKING = BLOCK'
     write (1,*) 'VARLOCATION=([3,4,5,6,7,8,9,10,11,12,13,14]=CELLCENTERED),ZONETYPE=FETRIANGLE'

	do K=1,nNodes
	  write(1,*)  nodes(1,K)
	end do
	do K=1,nNodes
      write(1,*)  nodes(2,K)
	end do
    do K=1,5
      do I=1,nElems
        if (K>1 .and. K<5) then
	      write(1,*)  U(K,I)/U(1,I)
	    else
          write(1,*)  U(K,I)
	    end if
      end do
    end do

	do K=1,3
      do I=1,nElems
	    write(1,*) B(K,I)/b_sc
	  end do
	end do

    do I=1,nElems
   ! вычисляем все простые переменные
	 rho = U(1,I)
     ux = U(2,I)/rho
     v = U(3,I)/rho
     w = U(4,I)/rho									! полная энергия (внутренняя плюс кинетическая)
	 E = U(5,I)
	 k_e = ux*ux+v*v+w*w
     p = (E-0.5*rho*k_e)*(gmm-1.0)

      write(1,*) p
    end do

    do I=1,nElems
   ! вычисляем все простые переменные
	 rho = U(1,I)
     ux = U(2,I)/rho
     v = U(3,I)/rho
     w = U(4,I)/rho									! полная энергия (внутренняя плюс кинетическая)
	 E = U(5,I)
	 k_e = ux*ux+v*v+w*w
     p = (E-0.5*rho*k_e)*(gmm-1.0)
	 sounds = sqrt(gmm*p/rho)
	 Mach = sqrt(ux*ux+v*v+w*w)/sounds
      write(1,*) Mach
    end do

    do I=1,nElems
      write(1,*) sum(B(:,I)**2)/2.0
    end do

    do I=1,nElems
      write(1,*) U(4,I)/(U(1,I)*elems_center(2,I))
    end do

    write(1,*) ''

    do nElem=1,nElems
      write (1,*) elems(1,nElem),elems(2,nElem),elems(3,nElem)
    end do

    close (unit=1,status='Keep')
    write (*,*) 'File ',current_sol_name,' has been created'
	write (*,*) 't=',t_current
    write (*,*) '************************'
end subroutine save_reservation

subroutine save_bufer
    implicit none

    integer::I,K, nElem

    real (kind=precis):: rho, ux, v, w, E, p, k_e, sounds, Mach


!	 call filename_se_sol(ntime,current_sol_name)

!	if (MOD (ntime, 10000).eq.0) then
	 current_sol_name = root // '/save' // extention

	 open (unit=1,file=current_sol_name,form='formatted')

	write (1,*) t_current
	write (1,*) tau
	write (1,*) ntime


    do K=1,5
		do I=1,nElems
			write(1,*)  U(K,I)
		end do
    end do

	do K=1,3
		do I=1,nElems
			write(1,*) B(K,I)
		end do
	end do

	do I=1,nEdges
		write(1,*) B_edges(1,I)
	end do

    close (unit=1,status='Keep')
    write (*,*) 'File ',current_sol_name,' has been created'
	write (*,*) 't=',t_current
    write (*,*) '************************'
end subroutine save_bufer


subroutine read_bufer
    implicit none

    integer::I,K, nElem


	current_sol_name = root // '/save' // extention

	open (1,file=current_sol_name)

	read (1,*) t_current
	read (1,*) tau
	read (1,*) ntime


    do K=1,5
		do I=1,nElems
			read(1,*)  U(K,I)
		end do
    end do

	do K=1,3
		do I=1,nElems
			read(1,*) B(K,I)
		end do
	end do

	do I=1,nEdges
		read(1,*) B_edges(1,I)
	end do

	close(1)

    write (*,*) 'Read file ',current_sol_name,' has been created'
	write (*,*) 't=',t_current
    write (*,*) '************************'
end subroutine read_bufer

subroutine elems_ce

	open (unit=1,file="./data/elems_ce.dat",form='formatted')
	write(1,*) nElems
	do i=1,nElems
		write(1,*) elems_center(:,i)
	end do

end subroutine elems_ce


subroutine save_bper

    real (kind=precis) :: rho, p, ux, v, w, Bx, By, Bz, alpha, Bpar, Bper, upar, uper

	integer :: k

	 open (unit=100,file='./data/bper.dat',form='formatted')

	 write (100,*) 'TITLE = "U graph"'
	 write (100,*) 'VARIABLES = "X", "Y", "bper", "Bz", "uper", "w"'
	 write (100,*) 'ZONE N=',nNodes, ',E=',nElems
	 write (100,*) 'DATAPACKING = BLOCK'
     write (100,*) 'VARLOCATION = ([3,4,5,6]=CELLCENTERED), &
	              ZONETYPE=FETRIANGLE '

	do K=1,nNodes
	  write(100,*)  nodes(1,K)
	end do
	do K=1,nNodes
      write(100,*)  nodes(2,K)
	end do

	alpha = pi/6.0


	do K=1,nElems
	  Bper = -B(1,k)*sin(alpha)+B(2,k)*cos(alpha)
	  write(100,*) Bper
	end do

	do K=1,nElems
	  write(100,*) B(3,K)
	end do

	do K=1,nElems
	  uper = -U(2,k)*sin(alpha)+U(3,k)*cos(alpha)
	  write(100,*) uper/U(1,K)
	end do

	do K=1,nElems
	  write(100,*) U(4,K)/U(1,K)
	end do

    write(100,*) ''

    do nElem=1,nElems
      write (100,*) elems(1,nElem),elems(2,nElem),elems(3,nElem)
    end do

    close (unit=100,status='Keep')
end subroutine save_bper

end module solution
