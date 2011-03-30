!==============================================================
! имена!
!==============================================================

MODULE names
  IMPLICIT NONE

!! radjet2D
  integer, parameter :: precis = 8


    real (kind=precis) :: refsspeed, refMach
    real (kind=precis) :: b_scale, r_grav, init_dens
    real (kind=precis), dimension(2):: refp

    real (kind=precis) :: cleanprev = 4D0

  !!!!!basic variables================================================
  integer                  :: nElem = 0, nEdge = 0
  integer                  :: NL1 = 0, NL2 = 0, NG1 = 0, NG2 = 0
  real (kind=precis), dimension(1:5)     :: flow_elem = 0.0, flow_edge = 0.0, div_B = 0.0
  real (kind=precis), dimension(1:3)     :: dB = 0.0
  real (kind=precis), dimension(1:2)     :: xy1 = 0.0, xy2 = 0.0
  real (kind=precis), dimension(1:5,1:5) :: T = 0.0, Ti = 0.0
  real (kind=precis), dimension(1:2)     :: buf = 0.0
  real (kind=precis), dimension(1:5)     :: Ul = 0.0, Ur = 0.0 ! данные для распадника

  ! массив -- Q(i) -- локальный номер узла треугольной ячейки,
  ! стоящего после узла с локальным номером _i_
  integer, parameter, dimension(1:6):: Q = (/ 2, 3, 1, 2, 3, 1/)
  integer :: i,j,i_1,j_1
  real (kind=precis) :: A = 0.0

  integer :: nrbc
  integer :: err, errr
  real (kind=precis) :: mag2gas_energ = 0.0, magn_decr, magn_gr,dif_c
  real (kind=precis) :: energ, press, jet_rr, bigjet_rr, maxivel
  real (kind=precis) :: epsil = 10D-7
  real (kind=precis),parameter :: pi= 3.14159265358979323846
  real (kind=precis) :: lsp = 1.0
  real (kind=precis) :: gmm !gmm = 1.4 ! показатель адиабаты

  real (kind=precis) :: kkap
  real (kind=precis) :: int0
  real (kind=precis) :: p0
  real (kind=precis) :: vell
  integer :: smplflag = 0
  !!!!basic variables=================================================

  real (kind=precis) :: mag_fact = 4D0*pi

  ! problem flags
  integer :: fl_clean = 0
  integer :: flag_axial = 1
  integer :: flag_magnetic = 1
  integer :: flag_radiative = 1
  integer :: flag_nrbc = 1
  integer :: flag_debug = 1
  logical :: read_reserve = .true.

  ! jet flux variables
  real (kind=precis) :: f_mass_out,f_energ_out,f_mass_in,f_energ_in,f_mass_out_jet,f_energ_out_jet
  real (kind=precis) :: strt_jet_num = 0.0

  ! gravity parameters
  real (kind=precis) :: gr_rad = 0.1
  real (kind=precis) :: G = 0.5
!    real (kind=precis) :: G = 0.0

  ! magnetic field parameters
  real (kind=precis) :: b_sc = 3.54490770181

  ! time and timestep parameters
  real (kind=precis) :: end_time = 0.5
  real (kind=precis) :: tau_out = 0.01, tau_out_fl = 0.005 ! output timesteps
  real (kind=precis):: t_start = 0.0    ! начальные момент времени
  real (kind=precis):: t_final = 1.0    ! конечный момент времени
  real (kind=precis):: t_current = 0.0  ! текущий момент времени
  integer:: ntime  ! Номер текущего временного слоя

  ! auto timestep
  real (kind=precis) :: tau = 0.1 ! Начальное значение временного шага
  real (kind=precis) :: CFL = 0.2
  integer :: nmsr = 0 ! number of element with maximal spectral radius
  !!!!!time variables=================================================

  ! computation time
  character(LEN=10) :: ddate, dtime, dzone
  integer, dimension(1:8) :: dvalues
  integer :: dhour, dmin, dsec, dmils, ddtime

  ! jet problem parameters
  real (kind=precis) :: r_disk = 0.6
  real (kind=precis) :: r_0 = 2.5, p_0 = 0.01
  real (kind=precis) :: omega0 = 1.0
  real (kind=precis) :: rho_init = 1.0

  !!!!!mesh variables=================================================
  real (kind=precis) :: hmin
  real (kind=precis) :: r_max=0.0, z_max=0.0

  integer:: nNodes                            ! количество узлов в сетке
  real (kind=precis),allocatable:: nodes(:,:) ! массив координат узлов
  integer,           allocatable:: node_elem(:,:) !элементы, прлегающие к точке
  integer,           allocatable:: elem_nodeN(:,:) !elem_nodeN(k,j) локальный номер точки k в j-м прилегающем к ней элементе
  integer,           allocatable:: node_elem1(:,:) !элементы, прлегающие к точке
  integer,           allocatable:: node_edges_count(:) !грани, прилегающие к точке
  integer,           allocatable:: node_elems_count(:) !грани, прилегающие к точке
  integer,           allocatable:: node_elems_count1(:) !грани, прилегающие к точке
  integer,           allocatable:: node_edges(:,:) !грани, прилегающие к точке
  real (kind=precis),allocatable:: node_edges_angles(:,:) !веса граней, прилегающие к точке

  integer:: nEdges                         ! количество ячеек в сетке
  integer,           allocatable:: edges(:,:) ! список граней - номера точек грани
  integer,           allocatable:: edge_elems(:,:) ! список граней - номера соседних ячеек
  integer,           allocatable:: elem_edgeN(:,:) !elem_edgeN(k,j) локальный номер ребра k в j-м прилегающем к ней элементе
  real (kind=precis),allocatable:: edge_n(:,:) !массив фиксированых нормалей к грани, согласовано с порядком вершин в грани
  real (kind=precis),allocatable:: edge_vol(:) ! массив длин граней
  real (kind=precis),allocatable:: edge_3dsurf(:) ! массив площадей граней
  real (kind=precis),allocatable:: edge_normSign(:,:)

  integer:: nElems                         ! количество ячеек в сетке
  integer,            allocatable:: elems(:,:) ! список КЭ (ячеек)
  integer,            allocatable:: elem_edges(:,:)!ребра ячейки
  real (kind=precis), allocatable:: elems_vol(:) !, elems_vol1! массив объемов ячеек
  real (kind=precis), allocatable:: elems_3dvol(:) ! массив цилиндр. объемов ячеек
  real (kind=precis), allocatable:: elems_vol1(:)
  real (kind=precis), allocatable:: elems_center(:,:)
  real (kind=precis), allocatable:: elem_edges_center(:,:,:)
  integer,            allocatable:: elems_elems(:,:)
  ! elems_elems(k,i) -- номер ячейки, соседствующей с ячейкой _i_ через
  ! ребро с локальным номером _k_.
  ! Локальная нумерация ребер должна быть согласована с локальной нумерацией узлов. Т.е.
  ! 1-я вершина -> первое ребро -> 2-я вершина -> второе ребро -> ...
  ! т.е. ребро _i_ имеет локальные номера узлов (i,i+1) -- все, кроме последнего,
  ! последнее -- (N,1)
  ! ребро  1 -- узлы с локальными номерами (1,2)
  ! ребро  2 -- узлы с локальными номерами (2,3)
  ! ребро  3 -- узлы с локальными номерами (3,1)
  real (kind=precis), allocatable :: h(:,:)
  real (kind=precis), allocatable :: e(:,:,:), n(:,:,:)
!       n(:,nElem,nEdge)

 ! ========================================================================
  ! описатели, неодходимые для задания и вычисления граничных условий

  integer:: nElems_b ! количество граничных ячеек ! (*)
  integer,            allocatable:: bc_type(:) ! массив типов граничного условия
  integer,            allocatable:: bc_elems(:) ! (*)
  integer,            allocatable:: bc_elem_edge(:) ! локальный номер граничного ребра
  integer,            allocatable:: not_bc_node(:)
  integer,            allocatable:: is_bc_elem(:)
  integer,            allocatable:: near_ort_node(:)
  integer,            allocatable:: period(:)
  ! bc_elems(k) -- номер ячейки, с которой граничит фиктивная
  ! ячейка с _локальным_ номером k
  integer,            allocatable:: bc_edges(:)
  integer,            allocatable:: bc_edges_type(:)

  integer:: nNodes_b
  integer, allocatable:: bc_nodes(:,:) ! (*)
  ! bc_nodes(1:2,k) -- номера узлов ребра, через которое фиктивная ячейка с номером k
  ! граничит с расчетной областью
  integer,            allocatable:: bc_node(:)
  integer,            allocatable:: bc_node_type(:)

  integer,            allocatable:: bc_node_neib(:,:,:)
!  integer, dimension(1:nElems_max,1:2,1:3) :: bc_node_neib ! bc_node_neib(номер гран.точки, номер соседа из двух, параметры соседа - элемента node_elems, лок. номер ребра, глобальный номер точки)
  integer,            allocatable:: is_bc_node(:)
!  integer, dimension(1:nElems_max):: is_bc_node = 0

  integer, allocatable :: outb_nums(:)
  integer :: outb_N

  ! radiation mesh
  real (kind=precis) :: absorb0
  integer :: nAngle, nFreq, nAng1, nAng2, refl_nAng1, refl_nAng2
  integer :: nAnglgr_psi, nAnglgr_theta
  real (kind=precis) :: d_psi, d_theta,psi1,psi2,theta1,theta2,boancos
  integer, parameter :: nAnglgr_max=100
  real (kind=precis), allocatable :: psi(:)	! 1:nAnglgr_max rad_angle(1)=0, rad_angle(nAnglgr_psi+1)=2*Pi
  real (kind=precis), allocatable :: theta(:)	! 1:nAnglgr_max rad_angle(1)=0, rad_angle(nAnglgr_theta)=Pi
  real (kind=precis), allocatable :: psi_mid(:),theta_mid(:)
  real (kind=precis), allocatable :: directs(:,:,:)	! directs(3,nAnglgr_psi,nAnglgr_theta) - cylindrical coords
  real (kind=precis), allocatable :: bdir_angle(:,:,:)	! directed body angle (3,nAngle_psi, nAngle_theta) (\int \omega d\Omega)
  real (kind=precis), allocatable :: t_om(:,:,:,:)	! t_omega (3,3,nAngle_psi, nAngle_theta) (\int \omega_i \omega_j d\Omega)
  real (kind=precis), allocatable :: b_angle(:,:)	! body angle in direction (nAngle_psi, nAngle_theta) (S_bit_of_sphere)
  real (kind=precis), allocatable :: omega_hat(:,:,:)	! integral omega_hat for elem (nAngle_psi, nAngle_theta)
  integer, allocatable :: reflecters(:,:,:,:) ! number for reflected ray (nAngle_psi, nAngle_theta, nWall, nAng), nWall=1,2: 1=>Ox, 2=>Oy; nAng=1,2: 1=>nAngle_psi, 2=>nAngle_theta
  integer, allocatable :: reflecting_b(:) ! reflecting_b(nElem_b)=0,1,2: 0=>nonrefl; 1=>Ox; 2=>Oy

  !!!!!mesh variables=================================================

  !!!!!solution=======================================================
  ! gas variables
  real (kind=precis), allocatable, target:: U(:,:)    ! массив для решения
  real (kind=precis), allocatable:: U_hat(:,:) ! массив для решения
  real (kind=precis), allocatable:: U_til(:,:) ! массив для решения
  real (kind=precis), allocatable:: U_angles(:,:,:)
  real (kind=precis), allocatable:: U_hat_angles(:,:,:)
  real (kind=precis), allocatable:: V_edges(:,:)
  real (kind=precis), allocatable:: V_nodes(:,:)
  real (kind=precis), allocatable:: gas_edge_flux(:,:)

  ! magnetic field
  real (kind=precis), allocatable:: B(:,:)     ! магнитное поле
  real (kind=precis), allocatable:: B_hat(:,:) ! магнитное поле
  real (kind=precis), allocatable:: B_angles(:,:,:)
  real (kind=precis), allocatable:: B_hat_angles(:,:,:)
  real (kind=precis), allocatable:: B_edges(:,:)
  real (kind=precis), allocatable:: B_edges_hat(:,:)
  real (kind=precis), allocatable:: B_nodes(:,:)
  real (kind=precis), dimension(1:5) :: RS = 0.0
  real (kind=precis)                 :: circ, dist = 0.0
  real (kind=precis), allocatable:: gradBB(:,:)
  real (kind=precis), allocatable:: MagnEn(:)
  real (kind=precis), allocatable:: divBB(:,:)
  real (kind=precis), allocatable:: rotB(:,:) !в ячейках
  real (kind=precis), allocatable:: rotB_nodes(:) !в точках
  real (kind=precis), allocatable:: rotB_edges(:) !на ребрах

  ! radiation fields
					! центр угла по \psi находится при \phi=0
  real (kind=precis), allocatable :: intens_p(:)	! intensity_p(nElem)
  real (kind=precis), allocatable :: intensity(:)	! intensity_full(nElem)
  real (kind=precis), allocatable :: intensity_hat(:)	! intensity_full(nElem)
  real (kind=precis), allocatable :: absorp_koef(:)	! absorption koefficient(nElem) (\alpha)
  real (kind=precis), allocatable :: disper_koef(:)	! dispertion koefficient(nElem) (\beta)
  real (kind=precis), allocatable :: weaken_koef(:)	! weakening koefficient(nElem) (\kappa)
  integer, allocatable :: indic_direc(:,:,:,:,:)	! indicatris directions (nElem,nAnglgr_psi,nAnglgr_theta,1:2,1:2)
							! (nElem,nAng1,nAng2,1,:) - transient (nElem,nAng1,nAng2,2,:) - reflected;
							! (nElem,nAng1,nAng2,k,1) - n_psi (nElem,nAng1,nAng2,k,2) - n_theta;
  real (kind=precis), allocatable :: indic_koef(:,:,:,:)	! indicatris koefficients (nElem,nAnglgr_psi,nAnglgr_theta,1:2)
						!(nElem,nAng1,nAng2,1) - transient (nElem,nAng1,nAng2,2) - reflected

  integer, allocatable :: rad_raspad(:,:,:,:)
  real (kind=precis), allocatable :: raspad_koefs(:,:) ! raspad_koefs(N_zapisi,(1)svoj koef ili (2)storonnij koef)
  real (kind=precis), allocatable :: SP(:,:,:,:), RC(:,:,:)
  real (kind=precis), allocatable :: Poynting(:,:) !Pointing vector
  real (kind=precis), allocatable :: Poynting_hat(:,:) !Pointing vector
  real (kind=precis), allocatable :: radstress(:,:,:) !radiative stresses

  real (kind=precis), allocatable :: radimpulse(:,:)
  real (kind=precis), allocatable :: radenerg(:)

  ! upstream koefficients
  real (kind=precis), allocatable:: up_koef(:,:,:)
  real (kind=precis), allocatable:: upstr_nodes(:,:) !upstr_nodes(nElem_of_node,nNode)
  real (kind=precis), allocatable:: up_koef_z(:,:)
  !!!!!solution=======================================================

  !!!!!t92 variables==================================================
  integer(4) :: t92_neq, t92_nvar, t92_nonz, t92_nac, t92_nrend, t92_mxp,t92_ipd,t92_maxit,t92_iter,t92_kstep
  integer(4), dimension(4) :: t92_mode
  integer(4), dimension(3) :: t92_merr
  integer(4), allocatable :: t92_row(:),t92_h(:,:),t92_c0(:),t92_c(:), t92_r(:)
  real(8), allocatable :: t92_vt(:), t92_a0(:), t92_b(:), t92_x(:),t92_d(:),t92_diag(:), t92_rr(:), t92_v(:), t92_cna0(:), t92_barm(:), t92_a(:)
  real(8) :: t92_stab,t92_bar,t92_stpbar,t92_epsin,t92_epsout

!	integer :: t92_neq = 4
!	integer :: t92_nvar = 4
!	integer :: t92_nonz = 7
!	integer :: t92_nac = 16*t92_nonz
!	integer :: t92_nrend = t92_nac
!	integer :: t92_mxp,t92_ipd,t92_maxit,t92_iter,t92_kstep
!	integer, dimension(t92_neq+1) :: t92_row
!	integer, dimension(9,-1:t92_nvar) :: t92_h
!	real(8), dimension(t92_neq) :: t92_vt
!	integer, dimension(t92_nonz) :: t92_c0
!	integer, dimension(0:t92_nac) :: t92_c
!	integer, dimension(0:t92_nrend) :: t92_r
!	real(8), dimension(t92_nonz) :: t92_a0
!	real(8), dimension(t92_neq) :: t92_b
!	real(8), dimension(t92_nvar) :: t92_x,t92_d,t92_diag
!	real(8), dimension(t92_neq) :: t92_rr, t92_v, t92_cna0, t92_barm
!	real(8), dimension(t92_nac) :: t92_a

  !!!!!file variables=================================================
  CHARACTER(LEN=*), PARAMETER:: root = './data' ! корневой каталог данных
!  CHARACTER(LEN=*), PARAMETER:: fname_sol = 'SOLUTION' ! имя файла с решением
  CHARACTER(LEN=*), PARAMETER:: fname_sol = 'sol' ! имя файла с суперэлементом
  CHARACTER(LEN=*), PARAMETER:: fname_sol_hat = 'sol_hat' ! имя файла с решением
!  CHARACTER(LEN=*), PARAMETER:: fname_d = 'D_STAR' ! имя файла с распределением диффузии
  CHARACTER(LEN=*), PARAMETER:: delim = '-'
  CHARACTER(LEN=*), PARAMETER:: extention = '.dat'
  CHARACTER(LEN=*), PARAMETER:: extentionVTK = '.vtk'
  INTEGER, PARAMETER:: max_filename_length = 50

!  CHARACTER(LEN = max_filename_length):: current_se_name = ''
  CHARACTER(LEN = max_filename_length):: current_sol_name = ''
!  CHARACTER(LEN = max_filename_length):: current_d_name = ''
  CHARACTER(LEN = max_filename_length):: current_sol_hat_name = ''

  CHARACTER(LEN = max_filename_length):: &
      current_time_step_name = './data/last_time_step.dat'
  !!!!!file variables=================================================

CONTAINS

  ! функция возвращает имя файла с решением на временном слое t
  SUBROUTINE filename_full_sol(time, name)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name

    name = ''
    name = root // '/' // fname_sol_hat // &
         delim // 'T' // i2c(time) // extention

  END SUBROUTINE filename_full_sol

  ! функция возвращает имя файла с решением на временном слое t-1 в формате Tecplot
  SUBROUTINE filename_se_sol(time, name)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name

    name = ''
    name = root // '/' // fname_sol // &
         delim // 'T' // i2c(time) // extention

  END SUBROUTINE filename_se_sol

  ! функция возвращает имя файла с решением на временном слое t-1 в формате VTK
  SUBROUTINE filename_se_sol_vtk(time, name)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name

    name = ''
    name = root // '/' // fname_sol // &
         delim // 'T' // i2c(time) // delim // 'vtk' // extentionVTK

  END SUBROUTINE filename_se_sol_vtk

  SUBROUTINE filename_se_sol_vtk_rad(time, name)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name

    name = ''
    name = root // '/' // fname_sol // &
         delim // 'T' // i2c(time) // delim // 'vtk-rad' // extentionVTK

  END SUBROUTINE filename_se_sol_vtk_rad

  SUBROUTINE filename_se_sol_vtk_radTR(time, name, param)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name
    CHARACTER(LEN=*), INTENT(in):: param

    name = ''
    name = root // '/' // fname_sol // &
         delim // 'T' // i2c(time) // param // delim // 'vtk-rad' // extentionVTK

  END SUBROUTINE filename_se_sol_vtk_radTR

  ! функция возвращает строку, соотвествующую заданному числу
  ! длина строки фиксированна (если число содержит меньшее количество разрядов,
  ! то в начало дописываются пробелы, если больше (или число отрицательное) --
  ! то возвращается пустая строка)

  FUNCTION i2c(n) RESULT(str)
    IMPLICIT NONE

    INTEGER n
    INTEGER, PARAMETER:: r = 8 ! количество разрядов -- соотвествует 10000000

    CHARACTER(LEN=r):: str


    INTEGER, DIMENSION(0:r-1):: digit = 0
    INTEGER i
    INTEGER:: foo

    IF( (n>10**r-1).OR.(n<0)) THEN
       str = ''
    ELSE

       digit(r-1) = n/10**(r-1)
       foo = n - digit(r-1)*10**(r-1)

       DO i=r-2,0,-1
          digit(i) = foo/10**i
          foo = foo - digit(i)*10**i
       ENDDO

       str = ''

       str = ACHAR(IACHAR('0') + digit(7)) // &
             ACHAR(IACHAR('0') + digit(6)) // &
             ACHAR(IACHAR('0') + digit(5)) // &
             ACHAR(IACHAR('0') + digit(4)) // &
             ACHAR(IACHAR('0') + digit(3)) // &
             ACHAR(IACHAR('0') + digit(2)) // &
             ACHAR(IACHAR('0') + digit(1)) // &
             ACHAR(IACHAR('0') + digit(0))
    ENDIF

  END FUNCTION i2c

END MODULE names
