! ========================================
! Главная программа
! ========================================
! ====================================================================================
! - Базовая схема для газа -- см. Toro, с. 573, раздел 16.7.3
! ====================================================================================

program plasma

  use alloc
  use trot
  use time
  use mesh
  use solution
  use flux
  use boundc
  use names
  use multies
  use solvers
  use geometry
  use radiationTRAB

  implicit none

  logical :: fexist

  call mmap_test

  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)
  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour


!  integer :: flag1 = 1
!  integer :: n_dis = 0
!  real :: maxxx
  inquire(file='./data/fluxes.dat', EXIST=fexist)

  if (fexist) then
    open (unit=301,file='./data/fluxes.dat',form='formatted',access='APPEND')
  else
    open (unit=301,file='./data/fluxes.dat',form='formatted')
    write (301,*) 'VARIABLES = "t", "jet_rad", "Fmass_out_jet", "bigjet_rad", "Fmass_out", "maxVel"'
    write (301,*) 'ZONE I=',500, "F=POINT"
  end if

!t_current, jet_rr, f_mass_out_jet, bigjet_rr, f_mass_out, maxivel
!t_current, jet_rr, f_mass_out_jet, bigjet_rr, f_mass_out, maxivel


  ! ===============================
  ! 1. Инициализация описателей сетки и массивов для
  ! граничных условий
  call init_mesh  ! описатели для треугольных элементов
  call init_mesh2 ! описатели для смещенных сеток
  call init_mesh3 ! центры элементов

  call hmin_compute
!  call radiation_solverTRAB

!  if (flag_radiative == 1) call init_groups
!  call init_mesh4 ! для вытекания - точки переноса по нормали

!  call init_mesh6 ! только для задачи с аккрецией - ГУ разных типов
  call init_mesh7 ! только для задачи с аккрецией в круге - ГУ разных типов
  call init_outbound

  call strt

!  call elems_ce
  epsil = 1E-7

  errr = 1

  ! ===========================================
  ! 2. Инициализация временных параметров. Задание
  !    - начального момента времени
  !    - конечного момента времени
  !    - начального временного шага
  !    - начальное значение номера временного слоя
  call init_time

  ! ======================================
  ! 3. Задание начальных и граничных условий
  call init_ic

!  call info_read_create
  ! ========================================
  ! Цикл по временным слоям
!  call save_vtk_solution_gas
  call strt_jets ! номер ячейки - границы джета
!  call read_bufer

  timestepping: do
!	tau = 0.001

     ! ===========================
     ! пересчет временных параметров
!	if (.False.) then
    if (mod(ntime,5000)==0) then
      print *, '***************'
      print *, 'Information'
!      call save_vtk_reservation
!      call save_bufer
      print *, 'ntime =', ntime
      print *, 't_current =', t_current
      print *, 'tau =', tau
      print *, 'max SR element'
      print *, elems_center(:,nmsr)
      print *, '***************'
      print *
!	  call save_bper
    end if
!      print *, 't_current =', t_current

!      print *, 'tau =', tau

    call time_step_compute
    ntime = ntime + 1 ! номер очередного временного шага
    t_current = t_current + tau ! текущий момент времени
    nrbc = 0

    do i=1,nElems
      U_til(:,i) = U(:,i)
    end do

    call gas_hllc

    do i=1,nElems
      U_til(:,i) = U_hat(:,i)
    end do

    call magnetic
    call magnetic_kinematics
    do i=1,nElems
      U_til(:,i) = U_hat(:,i)
    end do


! !     if (flag_radiative>0.5) then
! ! !      call radiation_solverTRAB
! ! !      call radiation_kinematics
! !     end if
! !
! !     if (flag_magnetic>0.5) then
! !       call magnetic
! !       call magnetic_kinematics
! !       do i=1,nElems
! !         U_hat(:,i) = U_til(:,i)
! !       end do
! !     end if

	! ====================================================
	! сохраняем решение на текущем временном слое
    if (ntime==1 .or. t_current<tau+tau_out*floor(t_current/tau_out)) then
      print *
      print *, '***************'
      print *, 'Solution print'
      call save_vtk_solution
      call save_control_point
      call compute_flux(jet_rr, f_mass_out_jet, bigjet_rr, f_mass_out, maxivel)
      write (301,*) t_current, jet_rr, f_mass_out_jet, bigjet_rr, f_mass_out, maxivel
      print *, 'max SR element'
      print *, elems_center(:,nmsr)
      print *, 'tau = ', tau
      print *, 't_current =', t_current
      print *, '***************'
      print *
    end if


!	call compute_fluxes(f_mass_out,f_energ_out,f_mass_out_jet,f_energ_out_jet,f_mass_in,f_energ_in)

!	if (mod(ntime,1)==0) then
!    if (t_current<tau+0.01*floor(t_current/0.01)) then
!		call magn_visc_control(mag2gas_energ, dif_c, magn_decr, magn_gr)
!		write (105,*) t_current, dif_c, (4.0*pi)*mag2gas_energ, magn_decr, (4.0*pi)*mag2gas_energ+magn_decr, magn_gr
!    end if
!    if (t_current<tau+tau_out_fl*floor(t_current/tau_out_fl)) then
! 		call compute_fluxes1(jet_rr,f_mass_out,f_energ_out,f_mass_out_jet,f_energ_out_jet,f_mass_in,f_energ_in)
!		write (101,*) t_current, jet_rr, f_mass_out,f_energ_out,f_mass_out_jet,f_energ_out_jet,f_mass_in,f_energ_in
!   end if

!	if (t_current>5.1) call total_clean2

    call total_clean

!	call clean

!    if (t_current>t_final .and. flag1==1) then
!        call save_solution
!		call save_bper
!		flag1 = 0
!		print *, n_dis
!print *, 'saaave'
!    end if

	! ====================================================
!	maxxx = 0.0
!	call clean1


  do i=1,nElems
    B(:,i) = B_hat(:,i)
  end do

  do i=1,nEdges
    B_edges(:,i) = B_edges_hat(:,i)
    do j=1,2
      if (abs(B_edges(j,i))<epsil) B_edges(j,i) = 0.0
    end do
  end do
! !     end if

    do i=1,nElems
      U(:,i) = U_hat(:,i)
    end do

! !     if (flag_radiative>0.5) then
! !       intensity(i) = intensity_hat(i)
! !       Poynting(:,i) = Poynting_hat(:,i)
! !     end if

!	call correct_press_1
!	call err_compute

!	do i=1,nNodes
!		do j=1,3
!			if(abs(B_nodes(j,i))<epsil) B_nodes(j,i) = 0.0
!		end do
!	end do

!    if (mod(ntime,500)==0) then
!		write (101,*) t_current, maxxx
!		n_dis = n_dis + 1
 !   end if

     ! ==========================================
     ! Проверяем условия окончания расчета по времени
     if (is_finished_time()) &
          exit timestepping

!	call info_write

  end do timestepping

!  call stp

  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)

  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour - ddtime

  print *, 'Computation time (msec): ', ddtime

end program plasma
