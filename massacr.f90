! ----------------------------------------------------------------------------------%%
!
! MASSACR 
! 
! SUMMARY: main method runs fluid dynamic simulation coupled to geochemical
!          simulation and writes selected model output to file
! 
! TO RUN: make -f theMakeFile
!		  mpirun -np 4 ./massacr
!
! ----------------------------------------------------------------------------------%%

PROGRAM main

use globals
use initialize
!use alteration
!use netcdf

implicit none
 
include 'mpif.h'
INCLUDE "IPhreeqc.f90.inc"
save

! functions within massacr.f90
interface
	
	! solves thermal energy equation
	function h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)
		use globals
		use initialize
		implicit none
		! integers
		integer :: i, j, n, ii, m=3
		! inputs 
		real(8) :: sx, sy, qx, qy, rho_in(xn,yn), phi_in(xn,yn)
		! velocity stuff
		real(8) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
		real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
		real(8) ::  velocities0(xn,2*yn)
		! matrix stuff
		real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
		! real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
		real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
		! real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
		real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
		real(8) :: kMatLong((xn-2)*(yn-2))
		real(8) :: mn(xn,yn)
		real(8) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
		real(8) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
		real(8) :: frac6_in(yn,2), temp6_in(yn,2), dt_in
	end function h_next

	! solves streamfunction vorticity equation
	function psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage,frac6_in)
		use globals
		use initialize
		implicit none
		! integers
		integer :: i, j, ii, n, m, stage
		! inputs
		real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
		real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
		! matrix stuff
		real(8) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
		real(8) :: psi_next(xn,yn)
		real(8) :: mn(xn,yn)
		! back to band
		real(8) :: aBand0(longP,4*((yn/2)-2) + 3), band_in(longP,4*((yn/2)-2) + 3)
		real(8) :: rhoLong(longP)
		real(8) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
	end function psi_next
	
	
	
	function make_band(perm_in,phi_in,permx,permy,rho_in)
		use globals
		use initialize
		implicit none
		integer :: i, j, ii, n, m
		real(8) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
		real(8) :: permx(xn,yn), permy(xn,yn), permLong(longP)
		real(8) :: permxLong(longP), permyLong(longP)
		real(8) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
		real(8) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
		real(8) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
		real(8) :: perm_long(longP)
	end function make_band
	
	
	function particles_next (trace_in, uTransport, vTransport, inval, num, num_sat)
		use globals
		use initialize
		implicit none
		integer :: i, j, ii, n, m, mm, nn, num, num_sat
		real(8) :: trace_in(5,num), particles_next(5,num)
		real(8) :: uTransport(xn,yn), vTransport(xn,yn)
		real(8) :: u_wt, v_wt, inval
		real(8) :: rando
	end function particles_next
		
	! transports solutes
	function solute_next(sol, uTransport, vTransport, seaw)
		use globals
		use initialize
		implicit none
		! integers
		integer :: i, j, ii, n, m
		! inputs
		real(8) :: sol(xn/cell,yn/cell), sol0(xn/cell,yn/cell)
		real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
		! solver stuff
		! real(8) :: uLong((xn/cell-2)*(yn/cell-2)), vLong((xn/cell-2)*(yn/cell-2))
		! real(8) :: aBand((xn/cell-2)*(yn/cell-2),5), bBand((xn/cell-2)*(yn/cell-2),5)
		! real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec((xn/cell-2)*(yn/cell-2))
		! real(8) :: sol_nextRow((xn/cell-2)*(yn/cell-2))
		real(8) :: uLong(((xn/cell)-2)*((yn/cell)-0)), vLong(((xn/cell)-0)*((yn/cell)-2))
		real(8) :: aBand(((xn/cell)-2)*((yn/cell)-0),5), bBand(((xn/cell)-0)*((yn/cell)-2),5)
		real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec(((xn/cell)-2)*((yn/cell)-0))
		real(8) :: sol_nextRow(((xn/cell)-2)*((yn/cell)-0)), sol_nextRowB(((xn/cell)-0)*((yn/cell)-2))
		real(8) :: seaw
		real(8) :: bm1(xn,yn), b0(xn,yn), bp1(xn,yn), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
		real(8) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6
	end function solute_next
	
! 	! runs geochemical alteration model for a single (coarse) cell
! 	function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
! 		use globals
! 		use initialize
! 		use alteration
! 		! declare yo shit
! 		integer :: order
! 		real(8) :: temp, timestep
! 		real(8) :: alt_next(1,167)
! 		real(8) :: alter0(1,167)
! 		real(8) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)
! 	end function alt_next

	! calculates fluid density
	function rho_next (h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(8) :: h_in(xn,yn), rho_next(xn,yn)
	end function rho_next
	
	! calculates fluid density
	function rho_one(h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(8) :: h_in, rho_one
	end function rho_one
	
	! calculates viscosity
	function visc_next(h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(8) :: h_in(xn,yn), visc_next(xn,yn)
	end function visc_next

	! calculates velocities from streamfunction values
	function velocities(psi)
		use globals
		use initialize
		implicit none
		real(8) :: velocities(xn,2*yn), psi(xn,yn)
		real(8) :: u0(xn,yn), v0(xn,yn)
	end function velocities

        ! calculates velocities from COARSE streamfunction values
	function velocitiesCoarse(psiCoarse)
		use globals
		use initialize
		implicit none
		real(8) :: velocitiesCoarse(xn/cell,2*yn/cell), psiCoarse(xn/cell,yn/cell)
		real(8) :: u0(xn/cell,yn/cell), v0(xn/cell,yn/cell)
	end function velocitiesCoarse

	! calculates partial derivative of any 1D or 2D array
	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial
	
	! calculates partial derivative of any 1D or 2D array
	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge(rows,cols)
	end function partial_edge
	
	! calculates partial derivative of any 1D or 2D array
	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge_p(rows,cols)
	end function partial_edge_p

	! writes 2D array to file
	function write_matrix ( m, n, table, filename )
		use globals
		implicit none
		integer :: m, n, j, output_status, unit0, reclen
		character ( len = * ) filename
		character ( len = 30 ) string
		real(4)  :: table(m,n) , write_matrix
	end function write_matrix

	! writes 1D array to file
	function write_vec ( n, vector, filename )
		use globals
		implicit none
		integer :: n, j, output_status, unit0
		character ( len = * ) filename 
		real(4)  :: vector(n), write_vec
	end function write_vec

end interface

!--------------DECLARE EVERYTHING 

! dependent variable arrays
real(8) :: h(xn,yn), psi(xn,yn), pside(xn,yn) ! xn rows deep & yn columns wide
real(8) :: hmat((xn*tn/(mstep*ar)),yn), psimat((xn*tn/(mstep*ar)),yn), rhomat((xn*tn/(mstep*ar)),yn)
real(8) :: velocities0(xn,2*yn)
real(8) :: phimat((xn*tn/(mstep*ar)),yn), umat((xn*tn/(mstep*ar)),yn), vmat((xn*tn/(mstep*ar)),yn), rhsmat((xn*tn/(mstep*ar)),yn)
real(8) :: permxMat((xn*tn/(mstep*ar)),yn), permyMat((xn*tn/(mstep*ar)),yn)
real(8) :: psiCoarseMat((xn*tn/(mstep*ar)),yn), uCoarseMat((xn*tn/(mstep*ar)),yn), vCoarseMat((xn*tn/(mstep*ar)),yn)
real(8) :: permmat((xn*tn/(mstep*ar)),yn)
real(8) :: u(xn,yn), v(xn,yn),  permx(xn,yn), permy(xn,yn), uTest(xn,yn), vTest(xn,yn), psiTest(xn,yn), hTest(xn,yn)


! autumn performance review
integer :: counti, countf, count_rate
real :: timeBit

! material properties
real(8) :: rho(xn,yn), visc(xn,yn)
real(8) :: rhs0(xn,yn)
integer :: unit
real(8) :: phiCoarse(xn/cell,yn/cell)
real(8) :: phi0(xn,yn), phi(xn,yn)

! netCDF & output stuff
integer :: xInt, yInt, tInt, hInt, uInt, vInt
integer :: ncid
integer :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
integer :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid
integer :: i, j, ii, m, n, jj
real(8) :: yep



! benchmark stuff
real(8) :: nusseltLocalv(xn,1), nuBar

! geochemical alteration stuff
real(8) :: alt0(1,altnum)

real(8) :: primaryShift(xn/cell,yn/cell,g_pri), secondaryShift(xn/cell,yn/cell,g_sec)


! solute transport stuff
real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
real(8) :: uCoarse(xn/cell,yn/cell), vCoarse(xn/cell,yn/cell)
real(8) :: psiCoarse(xn/cell,yn/cell), velocitiesCoarse0(xn/cell,2*yn/cell)
real(8) :: permeability0(xn,yn)

! message passing stuff
integer, parameter :: max_rows = 10000000
integer, parameter :: send_data_tag = 2001, return_data_tag = 2002
integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
integer :: num_procs, an_id, num_rows_to_receive
integer :: avg_rows_per_process, num_rows, num_rows_to_send
integer :: end_row, sender, start_row, num_rows_received
real(8) :: vector(max_rows), vector2(max_rows), partial_sum, sum
real(8) :: local_mean, global_mean
real(8) :: hLocal(xn*(yn/2)), dt_local
integer :: order

! formatted message passing arrays
real(8) :: hLong(xn*(yn/2))
real(8) :: priLong(xn*(yn/2),g_pri), priLocal(xn*(yn/2),g_pri)
real(8) :: secLong(xn*(yn/2),g_sec), secLocal(xn*(yn/2),g_sec)
real(8) :: solLong(xn*(yn/2),g_sol), solLocal(xn*(yn/2),g_sol)
real(8) :: medLong(xn*(yn/2),g_med), medLocal(xn*(yn/2),g_med)
real(8) :: priLongBit(xn*(yn/2)), priLocalBit(xn*(yn/2))
real(8) :: secLongBit(xn*(yn/2)), secLocalBit(xn*(yn/2))
real(8) :: solLongBit(xn*(yn/2)), solLocalBit(xn*(yn/2))
real(8) :: medLongBit(xn*(yn/2)), medLocalBit(xn*(yn/2))
!real(8) :: solLocal0((xn/cell)*(yn/cell),g_sol)

INTEGER(KIND=4) :: id, all=190
CHARACTER(LEN=10000) :: line
character(len=45000) :: inputz0
!real(8) :: alter(1,167)
!real(8), allocatable :: outmat(:,:)
real(8) :: outmat(4,190)
! REAL GRABS
real(8) :: temp3, timestep3, primary3(5), secondary3(114), solute3(15), medium3(7) ! important information
real(8) :: water

! MIGRATING ALTERATION STUFF



! STRINGS
character(len=25) :: s_verm_ca, s_analcime, s_phillipsite, s_clinozoisite, s_verm_na
character(len=25) :: s_diopside, s_epidote, s_minnesotaite, s_ferrite_ca, s_foshagite
character(len=25) :: s_gismondine, s_gyrolite, s_hedenbergite, s_chalcedony, s_verm_mg
character(len=25) :: s_ferrihydrite, s_lawsonite, s_merwinite, s_monticellite, s_natrolite
character(len=25) :: s_talc, s_smectite_low, s_prehnite, s_chlorite, s_rankinite 
character(len=25) :: s_scolecite, s_tobermorite_9a, s_tremolite, s_chamosite7a, s_clinochlore14a
character(len=25) :: s_clinochlore7a, s_andradite
character(len=25) :: s_saponite_ca, s_troilite, s_pyrrhotite, s_lepidocrocite, s_daphnite_7a
character(len=25) :: s_daphnite_14a, s_verm_k, s_greenalite, s_aragonite
character(len=25) :: s_siderite, s_kaolinite, s_goethite, s_dolomite, s_celadonite ! secondary
character(len=25) :: s_sio2, s_albite, s_calcite, s_mont_na, s_smectite, s_saponite ! secondary
character(len=25) :: s_stilbite, s_saponite_k, s_anhydrite, s_clinoptilolite, s_pyrite ! secondary
character(len=25) :: s_quartz, s_kspar, s_saponite_na, s_nont_na, s_nont_mg, s_nont_k ! secondary
character(len=25) :: s_nont_h, s_nont_ca, s_muscovite, s_mesolite, s_hematite, s_diaspore ! 
character(len=25) :: s_feldspar, s_pigeonite, s_augite, s_glass, s_magnetite ! primary
character(len=25) :: s_laumontite, s_mont_k, s_mont_mg, s_mont_ca
character(len=25) :: s_temp, s_timestep ! important information
character(len=25) :: s_ph, s_ca, s_mg, s_na, s_k, s_fe, s_s, s_si, s_cl, s_al, s_alk, s_co2 ! solutes
character(len=25) :: s_hco3, s_co3, s_pe
character(len=25) :: s_water, s_w, s_reactive ! medium

character(len=59000) :: L5
!character(len=6) :: L5


!path = "/home/navah/input/"
path2 = ""


!L5 = " "
L5 = "#  $Id: llnl.dat 4023 2010-02-09 21:02:42Z dlpark $" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"LLNL_AQUEOUS_MODEL_PARAMETERS" //NEW_LINE('')// &
&"-temperatures" //NEW_LINE('')// &
&"         0.0100   25.0000   60.0000  100.0000" //NEW_LINE('')// &
&"       150.0000  200.0000  250.0000  300.0000" //NEW_LINE('')// &
&"#debye huckel a (adh)" //NEW_LINE('')// &
&"-dh_a" //NEW_LINE('')// &
&"         0.4939    0.5114    0.5465    0.5995" //NEW_LINE('')// &
&"         0.6855    0.7994    0.9593    1.2180" //NEW_LINE('')// &
&"#debye huckel b (bdh)" //NEW_LINE('')// &
&"-dh_b" //NEW_LINE('')// &
&"         0.3253    0.3288    0.3346    0.3421" //NEW_LINE('')// &
&"         0.3525    0.3639    0.3766    0.3925" //NEW_LINE('')// &
&"-bdot" //NEW_LINE('')// &
&"         0.0374    0.0410    0.0438    0.0460" //NEW_LINE('')// &
&"         0.0470    0.0470    0.0340    0.0000" //NEW_LINE('')// &
&"#cco2   (coefficients for the Drummond (1981) polynomial)" //NEW_LINE('')// &
&"-co2_coefs" //NEW_LINE('')// &
&"        -1.0312              0.0012806" //NEW_LINE('')// &
&"          255.9                 0.4445" //NEW_LINE('')// &
&"      -0.001606" //NEW_LINE('')// &
&"NAMED_EXPRESSIONS" //NEW_LINE('')// &
&"#" //NEW_LINE('')// &
&"# formation of O2 from H2O " //NEW_LINE('')// &
&"# 2H2O =  O2 + 4H+ + 4e-  " //NEW_LINE('')// &
&"#" //NEW_LINE('')// &
&"	Log_K_O2" //NEW_LINE('')// &
&"	 	log_k      -85.9951" //NEW_LINE('')// &
&"		-delta_H	559.543	kJ/mol	# 	O2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-2.9 kcal/mol" //NEW_LINE('')// &
&"	        -analytic   38.0229    7.99407E-03   -2.7655e+004  -1.4506e+001  199838.45" //NEW_LINE('')// &
&"#	Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"SOLUTION_MASTER_SPECIES" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"#element species        alk     gfw_formula     element_gfw" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Al       Al+3           0.0     Al              26.9815" //NEW_LINE('')// &
 &"Alkalinity HCO3-        1.0     Ca0.5(CO3)0.5   50.05" //NEW_LINE('')// &
!&"Alkalinity	CO3-2	1.0			Ca0.5(CO3)0.5	50.05" //NEW_LINE('')// &
! &"C(-4)	CH4		0.0	CH4" //NEW_LINE('')// &
! &"C(-3)	C2H6		0.0	C2H6" //NEW_LINE('')// &
! &"C(-2)	C2H4		0.0	C2H4" //NEW_LINE('')// &
&"C        HCO3-          1.0     HCO3            12.0110" //NEW_LINE('')// &
!&"C        CO3-2          1.0     HCO3            12.0110" //NEW_LINE('')// &
! &"C(+2)	 CO		0	C" //NEW_LINE('')// &
 &"C(+4)    HCO3-          1.0     HCO3" //NEW_LINE('')// &
!&"C(+4)    CO3-2          1.0     HCO3" //NEW_LINE('')// &
&"Ca       Ca+2           0.0     Ca              40.078" //NEW_LINE('')// &
&"Cl       Cl-            0.0     Cl              35.4527" //NEW_LINE('')// &
&"Cl(-1)	 Cl-		0	Cl" //NEW_LINE('')// &
&"Cl(1)	 ClO-		0	Cl" //NEW_LINE('')// &
&"Cl(3)	 ClO2-		0	Cl" //NEW_LINE('')// &
&"Cl(5)	 ClO3-		0	Cl" //NEW_LINE('')// &
&"Cl(7)	 ClO4-		0	Cl" //NEW_LINE('')// &
&"E        e-             0.0     0.0             0.0" //NEW_LINE('')// &
&"Fe       Fe+2           0.0     Fe              55.847" //NEW_LINE('')// &
&"Fe(+2)   Fe+2           0.0     Fe" //NEW_LINE('')// &
&"Fe(+3)   Fe+3           -2.0    Fe" //NEW_LINE('')// &
&"H        H+             -1.     H               1.0079" //NEW_LINE('')// &
&"H(0)     H2             0.0     H" //NEW_LINE('')// &
&"H(+1)    H+             -1.     0.0" //NEW_LINE('')// &
&"K        K+             0.0     K               39.0983" //NEW_LINE('')// &
&"Mg       Mg+2           0.0     Mg              24.305" //NEW_LINE('')// &
&"Na       Na+            0.0     Na              22.9898" //NEW_LINE('')// &
&"O        H2O            0.0     O               15.994" //NEW_LINE('')// &
&"O(-2)    H2O            0.0     0.0" //NEW_LINE('')// &
&"O(0)     O2             0.0     O" //NEW_LINE('')// &
&"S	 SO4-2          0.0     SO4             32.066" //NEW_LINE('')// &
&"S(-2)	 HS-            1.0     S" //NEW_LINE('')// &
&"S(+2)	 S2O3-2		0.0	S" //NEW_LINE('')// &
&"S(+3)	 S2O4-2		0.0	S" //NEW_LINE('')// &
&"S(+4)	 SO3-2		0.0	S" //NEW_LINE('')// &
&"S(+5)	 S2O5-2		0.0	S" //NEW_LINE('')// &
&"S(+6)	 SO4-2          0.0     SO4" //NEW_LINE('')// &
&"S(+7)	 S2O8-2		0.0	S" //NEW_LINE('')// &
&"S(+8)	 HSO5-		0.0	S" //NEW_LINE('')// &
&"Si       SiO2         0.0     SiO2            28.0855" //NEW_LINE('')// &
!&"Si		H4SiO4	0.0	SiO2		28.0843" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"SOLUTION_SPECIES" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Al+3 =  Al+3 " //NEW_LINE('')// &
&"	-llnl_gamma	9.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Al+3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-128.681 kcal/mol" //NEW_LINE('')// &
&"-Vm  -3.3404  -17.1108  14.9917  -2.0716  2.8711 9 # supcrt" //NEW_LINE('')// &


&"Ca+2 =  Ca+2 " //NEW_LINE('')// &
&"	-llnl_gamma	6.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Ca+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-129.8 kcal/mol" //NEW_LINE('')// &
!&"	-millero -19.69 0.1058 -0.001256 1.617 -0.075 0.0008262" //NEW_LINE('')// &
&"-Vm  -0.3456  -7.252  6.149  -2.479  1.239  5  1.60  -57.1  -6.12e-3  1 # supcrt modified" //NEW_LINE('')// &


&"Cl- =  Cl- " //NEW_LINE('')// &
&"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Cl-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-39.933 kcal/mol" //NEW_LINE('')// &
!&"	-millero 16.37 0.0896 -0.001264 -1.494 0.034 -0.000621" //NEW_LINE('')// &
&"-Vm  4.465  4.801  4.325  -2.847  1.748  0  -0.331  20.16  0  1 # supcrt modified" //NEW_LINE('')// &


&"e- =  e- " //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol		e-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kJ/mol" //NEW_LINE('')// &

&"Fe+2 =  Fe+2 " //NEW_LINE('')// &
&"	-llnl_gamma	6.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Fe+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-22.05 kcal/mol" //NEW_LINE('')// &
&"-Vm  -0.3255  -9.687  1.536  -2.379  0.3033  5.5  -4.21e-2  37.96  0  1 # supcrt modified" //NEW_LINE('')// &


&"H+ =  H+ " //NEW_LINE('')// &
&"	-llnl_gamma	9.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	H+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kJ/mol" //NEW_LINE('')// &

&"HCO3- =  HCO3- " //NEW_LINE('')// &
&"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	HCO3-" //NEW_LINE('')// &
!&"-Vm  5.224  0  0  -5.85  5.126  0  0.404  110  -5.74e-3  1 # supcrt modified" //NEW_LINE('')// &
! not entirely sure about carbonate vs. bicarbonate decision here...

! &"CO3-2 =  CO3-2 " //NEW_LINE('')// &
! &"	-llnl_gamma	5.4	" //NEW_LINE('')// &
! &"	log_k 0" //NEW_LINE('')// &
! &"	-delta_H	0	kJ/mol	# 	CO3-2" //NEW_LINE('')// &
! &"-Vm  5.224  0  0  -5.85  5.126  0  0.404  110  -5.74e-3  1 # supcrt modified" //NEW_LINE('')// &

&"#	Enthalpy of formation:	-164.898 kcal/mol" //NEW_LINE('')// &
&"K+ =  K+ " //NEW_LINE('')// &
&"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	K+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-60.27 kcal/mol" //NEW_LINE('')// &
!&"	-millero 7.26 0.0892 -0.000736 2.722 -0.101 0.00151" //NEW_LINE('')// &
&"-Vm  3.322  -1.473  6.534  -2.712  9.06e-2  3.5  0  29.7  0  1 # supcrt modified" //NEW_LINE('')// &

&"Mg+2 =  Mg+2 " //NEW_LINE('')// &
&"	-llnl_gamma	8.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Mg+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-111.367 kcal/mol" //NEW_LINE('')// &
!&"	-millero -22.32 0.0868 -0.0016 2.017 -0.125 0.001457" //NEW_LINE('')// &
&"-Vm  -1.410  -8.6  11.13  -2.39  1.332  5.5  1.29  -32.9  -5.86e-3  1 # supcrt modified" //NEW_LINE('')// &

&"Na+ =  Na+ " //NEW_LINE('')// &
&"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	Na+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-57.433 kcal/mol" //NEW_LINE('')// &
!&"	-millero -3.46 0.1092 -0.000768 2.698 -0.106 0.001651" //NEW_LINE('')// &
&"-Vm  1.403  -2.285  4.419  -2.726  -5.125e-5  4.0  0.162  47.67  -3.09e-3  0.725 # sup" //NEW_LINE('')// &

&"H2O =  H2O " //NEW_LINE('')// &
&"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
&"        log_k   0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	H2O" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-68.317 kcal/mol" //NEW_LINE('')// &
&"SO4-2 =  SO4-2 " //NEW_LINE('')// &
&"	-llnl_gamma	4.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	SO4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-217.4 kcal/mol" //NEW_LINE('')// &
!&"	-millero 9.26 0.284 -0.003808 0.4348 -0.0099143 -8.4762e-05" //NEW_LINE('')// &
&"-Vm  8.0  2.51  -42.5 5.41  4.23 0 0 0 0 1 # supcrt modified" //NEW_LINE('')// &


&"SiO2 =  SiO2 " //NEW_LINE('')// &
&"	-llnl_gamma	3.0000	" //NEW_LINE('')// &
&"	log_k 0" //NEW_LINE('')// &
&"	-delta_H	0	kJ/mol	# 	SiO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-209.775 kcal/mol" //NEW_LINE('')// &
!&"-Vm  10.5  1.7  20  -2.7  0.1291 # supcrt + 2*H2O in a1" //NEW_LINE('')// &
! not sure of ion vs. sio2 here

&"2H2O =  O2 + 4H+ + 4e-  " //NEW_LINE('')// &
&"	-CO2_llnl_gamma" //NEW_LINE('')// &
&" 	log_k      -85.9951" //NEW_LINE('')// &
&"	-delta_H	559.543	kJ/mol	# 	O2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-2.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic   38.0229    7.99407E-03   -2.7655e+004  -1.4506e+001  199838.45" //NEW_LINE('')// &
&"#	Range:  0-300" //NEW_LINE('')// &
&"-Vm  5.7889  6.3536  3.2528  -3.0417  -0.3943 # supcrt" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.0000 SO4-- + 1.0000 H+  =  HS- +2.0000 O2  " //NEW_LINE('')// &
&"        -llnl_gamma           3.5    " //NEW_LINE('')// &
&"        log_k           -138.3169" //NEW_LINE('')// &
&"	-delta_H	869.226	kJ/mol	# 	HS-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-3.85 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.6251e+001 3.9525e-002 -4.5443e+004 -1.1107e+001 3.1843e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
!&"-Vm 8.2 9.2590  2.1108  -3.1618 1.1748  0 -0.3 15 0 1 # supcrt modified" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" .5000 O2 + 2.0000 HS-  = S2--  + H2O" //NEW_LINE('')// &
&"#2 HS- = S2-- +2 H+ + 2e-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           33.2673" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.21730E+02   -0.12307E-02    0.10098E+05   -0.88813E+01    0.15757E+03" //NEW_LINE('')// &
&"	-mass_balance	S(-2)2" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"#	-add_logk	Log_K_O2	0.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"2.0000 H+  + 2.0000 SO3--  = S2O3--  + O2  + H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -40.2906" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O3-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic  0.77679E+02    0.65761E-01   -0.15438E+05   -0.34651E+02   -0.24092E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

! ! new
! &"CO3-2 + H+ = HCO3-" //NEW_LINE('')// &
! &" 	-log_k	10.329" //NEW_LINE('')// &
! &"	-delta_H	-3.561	kcal" //NEW_LINE('')// &
! &"  -analytic	107.8871	0.03252849	-5151.79	-38.92561	563713.9" //NEW_LINE('')// &
! &"-Vm  8.615  0  -12.21 0  1.667  0  0  264  0  1 # supcrt modified" //NEW_LINE('')// &
!
! ! new
! &"CO3-2 + 2 H+ = CO2 + H2O" //NEW_LINE('')// &
! &" 	-log_k	16.681" //NEW_LINE('')// &
! &"	-delta_h -5.738	kcal" //NEW_LINE('')// &
! &"  -analytic	464.1965	0.09344813	-26986.16	-165.75951	2248628.9" //NEW_LINE('')// &
! &"-Vm  21.78  -49.4  -91.7  31.96 # supcrt modified" //NEW_LINE('')// &

! &" H+  + HCO3-  + H2O  = CH4 + 2.0000 O2" //NEW_LINE('')// &
! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
! &"        log_k            -144.1412" //NEW_LINE('')// &
! &"	-delta_H	863.599	kJ/mol	# 	CH4" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-21.01 kcal/mol" //NEW_LINE('')// &
! &"	-analytic    -0.41698E+02    0.36584E-01   -0.40675E+05    0.93479E+01   -0.63468E+03" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
!&"-Vm 7.7" //NEW_LINE('')// &
! unsure but unimportant


! &" CO3-2 + 10 H+ + 8 e- = CH4 + 3 H2O" //NEW_LINE('')// &
! &"        -log_k	41.071" //NEW_LINE('')// &
! &"	-delta_h -61.039 kcal" //NEW_LINE('')// &
! &"	-Vm 7.7" //NEW_LINE('')// &

! &"" //NEW_LINE('')// &
! &" 2.0000 H+  + 2.0000 HCO3-  + H2O  = C2H6  + 3.5000 O2" //NEW_LINE('')// &
! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
! &"        log_k            -228.6072" //NEW_LINE('')// &
! &"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	C2H6" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
! &"        -analytic    0.10777E+02    0.72105E-01   -0.67489E+05   -0.13915E+02   -0.10531E+04" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &

! &"" //NEW_LINE('')// &
! &" 2.000 H+  + 2.0000 HCO3-  = C2H4 + 3.0000 O2" //NEW_LINE('')// &
! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
! &"        log_k            -254.5034" //NEW_LINE('')// &
! &"	-delta_H	1446.6	kJ/mol	# 	C2H4" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	24.65 kcal/mol" //NEW_LINE('')// &
! &"        -analytic    -0.30329E+02    0.71187E-01   -0.73140E+05    0.00000E+00    0.00000E+00" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &

 &"" //NEW_LINE('')// &
! &" 1.0000 HCO3- + 1.0000 H+  =  CO +1.0000 H2O +0.5000 O2 " //NEW_LINE('')// &
! &"        -llnl_gamma           3.0    " //NEW_LINE('')// &
! &"        log_k           -41.7002" //NEW_LINE('')// &
! &"	-delta_H	277.069	kJ/mol	# 	CO" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-28.91 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.0028e+002 4.6877e-002 -1.8062e+004 -4.0263e+001 3.8031e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &

&" 1.0000 Cl- + 0.5000 O2  =  ClO-   " //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -15.1014" //NEW_LINE('')// &
&"	-delta_H	66.0361	kJ/mol	# 	ClO-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-25.6 kcal/mol" //NEW_LINE('')// &
&"        -analytic 6.1314e+001 3.4812e-003 -6.0952e+003 -2.3043e+001 -9.5128e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.0000 O2 + 1.0000 Cl-  =  ClO2-   " //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -23.108" //NEW_LINE('')// &
&"	-delta_H	112.688	kJ/mol	# 	ClO2-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-15.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic 3.3638e+000 -6.1675e-003 -4.9726e+003 -2.0467e+000 -2.5769e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.5000 O2 + 1.0000 Cl-  =  ClO3-   " //NEW_LINE('')// &
&"        -llnl_gamma           3.5    " //NEW_LINE('')// &
&"        log_k           -17.2608" //NEW_LINE('')// &
&"	-delta_H	81.3077	kJ/mol	# 	ClO3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-24.85 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.8852e+001 -4.8281e-003 -4.6779e+003 -1.0772e+001 -2.0783e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 2.0000 O2 + 1.0000 Cl-  =  ClO4-   " //NEW_LINE('')// &
&"        -llnl_gamma           3.5    " //NEW_LINE('')// &
&"        log_k           -15.7091" //NEW_LINE('')// &
&"	-delta_H	62.0194	kJ/mol	# 	ClO4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-30.91 kcal/mol" //NEW_LINE('')// &
&"        -analytic 7.0280e+001 -6.8927e-005 -5.5690e+003 -2.6446e+001 -1.6596e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&" 1.0000 H+ + 1.0000 Fe++ + 0.2500 O2  =  Fe+++ +0.5000 H2O " //NEW_LINE('')// &
&"        -llnl_gamma           9.0    " //NEW_LINE('')// &
&"        log_k           +8.4899" //NEW_LINE('')// &
&"	-delta_H	-97.209	kJ/mol	# 	Fe+3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-11.85 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.7808e+001 -1.1753e-002 4.7609e+003 5.5866e+000 7.4295e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.0000 H2O  =  H2 +0.5000 O2   " //NEW_LINE('')// &
&"	-CO2_llnl_gamma" //NEW_LINE('')// &
&"        log_k           -46.1066" //NEW_LINE('')// &
&"	-delta_H	275.588	kJ/mol	# 	H2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1 kcal/mol" //NEW_LINE('')// &
&"        -analytic 6.6835e+001 1.7172e-002 -1.8849e+004 -2.4092e+001 4.2501e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.0000 SO4-- + 1.0000 H+ + 0.5000 O2  =  HSO5-  " //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -17.2865" //NEW_LINE('')// &
&"	-delta_H	140.038	kJ/mol	# 	HSO5-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-185.38 kcal/mol" //NEW_LINE('')// &
&"        -analytic 5.9944e+001 3.0904e-002 -7.7494e+003 -2.4420e+001 -1.2094e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&" 2.0000 H+  + 2.0000 SO3--  = S2O4--  + .500 O2  + H2O" //NEW_LINE('')// &
&"        -llnl_gamma           5.0    " //NEW_LINE('')// &
&"#        log_k           -25.2075" //NEW_LINE('')// &
&"        log_k           -25.2076" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"#        -analytic  -0.15158E+05   -0.31356E+01    0.47072E+06    0.58544E+04    0.73497E+04" //NEW_LINE('')// &
&"	-analytic	-2.3172e2	2.0393e-3	-7.1011e0	8.3239e1	9.4155e-1" //NEW_LINE('')// &
&"#	changed 3/23/04, corrected to supcrt temperature dependence, GMA" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"# 2.0000 SO3--  + .500 O2  + 2.0000 H+  = S2O6--  + H2O" //NEW_LINE('')// &
&"#  H2O = .5 O2 + 2H+ + 2e- " //NEW_LINE('')// &
&"2SO3-- = S2O6-- + 2e-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           41.8289" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O6-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.14458E+03    0.61449E-01    0.71877E+04   -0.58657E+02    0.11211E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"	-add_logk  Log_K_O2	0.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 2.0000 SO3--  + 1.500 O2  + 2.0000 H+  = S2O8--  + H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           70.7489" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O8-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.18394E+03    0.60414E-01    0.13864E+05   -0.71804E+02    0.21628E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"O2 + H+ + 3.0000 HS-  = S3--  + 2.0000 H2O" //NEW_LINE('')// &
&"# 2H2O = O2 + 4H+ + 4e-" //NEW_LINE('')// &
&"#3HS- = S3-- + 3H+ + 4e-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           79.3915" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S3-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic -0.51626E+02    0.70208E-02    0.31797E+05    0.11927E+02   -0.64249E+06" //NEW_LINE('')// &
&"	-mass_balance	S(-2)3" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"#	-add_logk  Log_K_O2	1.0" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"# 3.0000 SO3--  + 4.0000 H+  = S3O6-- + .500 O2 + 2.0000 H2O" //NEW_LINE('')// &
&"# .5 O2 + 2H+ + 2e- = H2O" //NEW_LINE('')// &
&"3SO3-- + 6 H+ + 2e- = S3O6-- + 3H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -6.2316" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S3O6-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.23664E+03    0.12702E+00   -0.10110E+05   -0.99715E+02   -0.15783E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"	-add_logk	Log_K_O2	-0.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.5000 O2 + 2.0000 H+ + 4.0000 HS-  = S4--  + 3.0000 H2O" //NEW_LINE('')// &
&"#4 HS- = S4-- + 4H+ + 6e-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           125.2958" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.20875E+03    0.58133E-01    0.33278E+05   -0.85833E+02    0.51921E+03" //NEW_LINE('')// &
&"	-mass_balance	S(-2)4" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"#	-add_logk	Log_K_O2	1.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"# 4.0000 SO3-- + 6.0000 H+  = S4O6-- + 1.500 O2 + 3.0000 H2O" //NEW_LINE('')// &
&"4 SO3-- + 12 H+ + 6e- = S4O6-- + 6H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -38.3859" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S4O6-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.32239E+03    0.19555E+00   -0.23617E+05   -0.13729E+03   -0.36862E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"	-add_logk	Log_K_O2	-1.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"2.0000 O2 + 3.0000 H+  + 5.0000 HS-  = S5--  + 4.0000 H2O" //NEW_LINE('')// &
&"#5 HS- = S5-- + 5H+ + 8e-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           170.9802" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S5-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.30329E+03    0.88033E-01    0.44739E+05   -0.12471E+03    0.69803E+03" //NEW_LINE('')// &
&"	-mass_balance	S(-2)5" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"#	-add_logk	Log_K_O2	2" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"# 5.0000 SO3-- + 8.0000 H+  = S5O6-- + 2.5000 O2 + 4.0000 H2O" //NEW_LINE('')// &
&"# 2.5O2 + 10 H+ + 10e- = 5H2O" //NEW_LINE('')// &
&"5SO3-- + 18H+ + 10e- = S5O6-- + 9H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -99.4206" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S5O6-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.42074E+03    0.25833E+00   -0.43878E+05   -0.18178E+03   -0.68480E+03" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"	-add_logk	Log_K_O2	-2.5" //NEW_LINE('')// &

! &"" //NEW_LINE('')// &
! &"# 1.0000 H+  + HCO3-  + HS-  + NH3 = SCN-  + 3.0000 H2O" //NEW_LINE('')// &
! &"#        -llnl_gamma           3.5    " //NEW_LINE('')// &
! &"#        log_k            3.0070" //NEW_LINE('')// &
! &"#	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	SCN-" //NEW_LINE('')// &
! &"##	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
! &"#        -analytic  0.16539E+03    0.49623E-01   -0.44624E+04   -0.65544E+02   -0.69680E+02" //NEW_LINE('')// &
! &"##       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&" 1.0000 SO4--  =  SO3-- +0.5000 O2   " //NEW_LINE('')// &
&"        -llnl_gamma           4.5    " //NEW_LINE('')// &
&"        log_k           -46.6244" //NEW_LINE('')// &
&"	-delta_H	267.985	kJ/mol	# 	SO3-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-151.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.3771e+001 6.5102e-004 -1.3330e+004 4.7164e+000 -2.0800e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 Al+++  =  Al(OH)2+ +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -10.5945" //NEW_LINE('')// &
&"	-delta_H	98.2822	kJ/mol	# 	Al(OH)2+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-241.825 kcal/mol" //NEW_LINE('')// &
&"        -analytic 4.4036e+001 2.0168e-002 -5.5455e+003 -1.6987e+001 -8.6545e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&"" //NEW_LINE('')// &
&"2.0000 SO4-- + 1.0000 Al+++  =  Al(SO4)2-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +4.9000" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al(SO4)2-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

&" " //NEW_LINE('')// &
&"28.0000 H2O + 13.0000 Al+++  =  Al13O4(OH)24+7 +32.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           6.0    " //NEW_LINE('')// &
&"        log_k           -98.73" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al13O4(OH)24+7" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

&" " //NEW_LINE('')// &
&"2.0000 H2O + 2.0000 Al+++  =  Al2(OH)2++++ +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           5.5    " //NEW_LINE('')// &
&"        log_k           -7.6902" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al2(OH)2+4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

&" " //NEW_LINE('')// &
&"4.0000 H2O + 3.0000 Al+++  =  Al3(OH)4+5 +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           6.0    " //NEW_LINE('')// &
&"        log_k           -13.8803" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Al3(OH)4+5" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

&" " //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 Al+++  =  AlO2- +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -22.8833" //NEW_LINE('')// &
&"	-delta_H	180.899	kJ/mol	# 	AlO2-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-222.079 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.0803e+001 -3.4379e-003 -9.7391e+003 0.0000e+000 0.0000e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.0000 H2O + 1.0000 Al+++  =  AlOH++ +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.5    " //NEW_LINE('')// &
&"        log_k           -4.9571" //NEW_LINE('')// &
&"	-delta_H	49.798	kJ/mol	# 	AlOH+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-185.096 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.6224e-001 8.8816e-003 -1.8686e+003 -4.3195e-001 -2.9158e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 Al+++  =  AlSO4+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +3.0100" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	AlSO4+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &

&" " //NEW_LINE('')// &
&"1.0000 HCO3- + 1.0000 H+  =  CO2 +1.0000 H2O" //NEW_LINE('')// &
&"        -CO2_llnl_gamma" //NEW_LINE('')// &
&"        log_k           +6.3447" //NEW_LINE('')// &
&"	-delta_H	-9.7027	kJ/mol	# 	CO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-98.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.0534e+001 2.1746e-002 2.5216e+003 7.9125e-001 3.9351e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
! &"CO3-2 + 2 H+ = CO2 + H2O" //NEW_LINE('')// &
! &"        -log_k	16.681" //NEW_LINE('')// &
! &"	-delta_h -5.738	kcal" //NEW_LINE('')// &
! &"        -analytic	464.1965	0.09344813	-26986.16	-165.75951	2248628.9" //NEW_LINE('')// &
! &"        -Vm  21.78  -49.4  -91.7  31.96 # supcrt modified" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.0000 HCO3-  =  CO3-- +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.5    " //NEW_LINE('')// &
&"        log_k           -10.3288" //NEW_LINE('')// &
&"	-delta_H	14.6984	kJ/mol	# 	CO3-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-161.385 kcal/mol" //NEW_LINE('')// &
&"        -analytic -6.9958e+001 -3.3526e-002 -7.0846e+001 2.8224e+001 -1.0849e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"CO3-2 + H+ = HCO3-" //NEW_LINE('')// &
&"        -llnl_gamma           5.4    " //NEW_LINE('')// &
&"        log_k           10.3288" //NEW_LINE('')// &
&"	-delta_h -3.561	kcal" //NEW_LINE('')// &
&"        -analytic	107.8871	0.03252849	-5151.79	-38.92561	563713.9" //NEW_LINE('')// &
&"#       -Vm  8.615  0  -12.21 0  1.667  0  0  264  0  1 # supcrt modified" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.0000 HCO3- + 1.0000 Ca++  =  CaCO3 +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -7.0017" //NEW_LINE('')// &
&"	-delta_H	30.5767	kJ/mol	# 	CaCO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-287.39 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.3045e+002 5.5350e-002 -8.5056e+003 -9.1096e+001 -1.3279e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
!&"-Vm  -.2430  -8.3748  9.0417  -2.4328  -.0300 # supcrt" //NEW_LINE('')// &
! again

! &"" //NEW_LINE('')// &
! &"Ca+2 + CO3-2 = CaCO3" //NEW_LINE('')// &
! &"        -log_k	3.224" //NEW_LINE('')// &
! &"	-delta_h 3.545	kcal" //NEW_LINE('')// &
! &"        -analytic	-1228.732	-0.299440	35512.75	485.818" //NEW_LINE('')// &
! &"#       -Vm  -.2430  -8.3748  9.0417  -2.4328  -.0300 # supcrt" //NEW_LINE('')// &
!

&"" //NEW_LINE('')// &
&"1.0000 Cl- + 1.0000 Ca++  =  CaCl+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -0.6956" //NEW_LINE('')// &
&"	-delta_H	2.02087	kJ/mol	# 	CaCl+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-169.25 kcal/mol" //NEW_LINE('')// &
&"        -analytic 8.1498e+001 3.8387e-002 -1.3763e+003 -3.5968e+001 -2.1501e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 Cl- + 1.0000 Ca++  =  CaCl2" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -0.6436" //NEW_LINE('')// &
&"	-delta_H	-5.8325	kJ/mol	# 	CaCl2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-211.06 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.8178e+002 7.6910e-002 -3.1088e+003 -7.8760e+001 -4.8563e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"1.0000 HCO3- + 1.0000 Ca++  =  CaHCO3+" //NEW_LINE('')// &
! &"        -llnl_gamma           4.0    " //NEW_LINE('')// &
! &"        log_k           +1.0467" //NEW_LINE('')// &
! &"	-delta_H	1.45603	kJ/mol	# 	CaHCO3+" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-294.35 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 5.5985e+001 3.4639e-002 -3.6972e+002 -2.5864e+001 -5.7859e+000" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"-Vm  3.1911  .0104  5.7459  -2.7794  .3084 5.4 # supcrt" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"1.0000 H2O + 1.0000 Ca++  =  CaOH+ +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -12.85" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	CaOH+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &

&"#Al(OH)4-            82" //NEW_LINE('')// &
&"        Al+3 + 4H2O = Al(OH)4- + 4H+ " //NEW_LINE('')// &
&"        log_k           -22.7" //NEW_LINE('')// &
&"        delta_h 42.3 kcal" //NEW_LINE('')// &
&"        -analytical     51.578          0.0     -11168.9        -14.865         0.0" //NEW_LINE('')// &


		
&"1.0000 SO4-- + 1.0000 Ca++  =  CaSO4" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +2.1111" //NEW_LINE('')// &
&"	-delta_H	5.4392	kJ/mol	# 	CaSO4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-345.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.8618e+002 8.4084e-002 -7.6880e+003 -1.1449e+002 -1.2005e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"-Vm  2.7910  -.9666  6.1300  -2.7390  -.0010 # supcrt" //NEW_LINE('')// &

&"2.0000 H2O + 1.0000 Fe++  =  Fe(OH)2 +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -20.6" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 Fe+++  =  Fe(OH)2+ +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -5.67" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)2+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"3.0000 H2O + 1.0000 Fe+++  =  Fe(OH)3 +3.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -12" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"3.0000 H2O + 1.0000 Fe++  =  Fe(OH)3- +3.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -31" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"4.0000 H2O + 1.0000 Fe+++  =  Fe(OH)4- +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -21.6" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"4.0000 H2O + 1.0000 Fe++  =  Fe(OH)4-- +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -46" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(OH)4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 SO4-- + 1.0000 Fe+++  =  Fe(SO4)2-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +3.2137" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe(SO4)2-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 H2O + 2.0000 Fe+++  =  Fe2(OH)2++++ +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           5.5    " //NEW_LINE('')// &
&"        log_k           -2.95" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe2(OH)2+4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"4.0000 H2O + 3.0000 Fe+++  =  Fe3(OH)4+5 +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           6.0    " //NEW_LINE('')// &
&"        log_k           -6.3" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Fe3(OH)4+5" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 HCO3- + 1.0000 Fe++  =  FeCO3 +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -5.5988" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 HCO3- + 1.0000 Fe+++  =  FeCO3+ +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -0.6088" //NEW_LINE('')// &
&"	-delta_H	-50.208	kJ/mol	# 	FeCO3+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-188.748 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.7100e+002 8.0413e-002 -4.3217e+002 -7.8449e+001 -6.7948e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Fe++ + 1.0000 Cl-  =  FeCl+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -0.1605" //NEW_LINE('')// &
&"	-delta_H	3.02503	kJ/mol	# 	FeCl+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-61.26 kcal/mol" //NEW_LINE('')// &
&"        -analytic 8.2435e+001 3.7755e-002 -1.4765e+003 -3.5918e+001 -2.3064e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Fe+++ + 1.0000 Cl-  =  FeCl++" //NEW_LINE('')// &
&"        -llnl_gamma           4.5    " //NEW_LINE('')// &
&"        log_k           -0.8108" //NEW_LINE('')// &
&"	-delta_H	36.6421	kJ/mol	# 	FeCl+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-180.018 kJ/mol" //NEW_LINE('')// &
&"        -analytic 1.6186e+002 5.9436e-002 -5.1913e+003 -6.5852e+001 -8.1053e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 Cl- + 1.0000 Fe++  =  FeCl2" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -2.4541" //NEW_LINE('')// &
&"	-delta_H	6.46846	kJ/mol	# 	FeCl2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-100.37 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.9171e+002 7.8070e-002 -4.1048e+003 -8.2292e+001 -6.4108e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 Cl- + 1.0000 Fe+++  =  FeCl2+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +2.1300" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl2+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"4.0000 Cl- + 1.0000 Fe+++  =  FeCl4-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -0.79" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"4.0000 Cl- + 1.0000 Fe++  =  FeCl4--" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -1.9" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeCl4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.4108e+002 -6.0086e-003 9.7979e+003 8.4084e+001 1.5296e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 HCO3- + 1.0000 Fe++  =  FeHCO3+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +2.7200" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeHCO3+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 H2O + 1.0000 Fe++  =  FeOH+ +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -9.5" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeOH+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 H2O + 1.0000 Fe+++  =  FeOH++ +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.5    " //NEW_LINE('')// &
&"        log_k           -2.19" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeOH+2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 Fe++  =  FeSO4" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +2.2000" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	FeSO4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 Fe+++  =  FeSO4+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +1.9276" //NEW_LINE('')// &
&"	-delta_H	27.181	kJ/mol	# 	FeSO4+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-932.001 kJ/mol" //NEW_LINE('')// &
&"        -analytic 2.5178e+002 1.0080e-001 -6.0977e+003 -1.0483e+002 -9.5223e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"#1.0000 HS- + 1.0000 H+  =  H2S" //NEW_LINE('')// &
&"#        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"#        log_k           +6.99" //NEW_LINE('')// &
&"#        -analytic 1.2833e+002 5.1641e-002 -1.1681e+003 -5.3665e+001 -1.8266e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
! &"# these (above) H2S values are from " //NEW_LINE('')// &
! &"# Suleimenov & Seward, Geochim. Cosmochim. Acta, v. 61, p. 5187-5198." //NEW_LINE('')// &
! &"# values below are the original Thermo.com.v8.r6.230 data from somewhere" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 HS- + 1.0000 H+  =  H2S" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +6.9877" //NEW_LINE('')// &
&"	-delta_H	-21.5518	kJ/mol	# 	H2S" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-9.001 kcal/mol" //NEW_LINE('')// &
&"        -analytic 3.9283e+001 2.8727e-002  1.3477e+003 -1.8331e+001  2.1018e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 H+ + 1.0000 SO3--  =  H2SO3" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +9.2132" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 H+ + 1.0000 SO4--  =  H2SO4" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -1.0209" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SO4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 SiO2  =  H2SiO4-- +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -22.96" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H2SiO4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"8.0000 H2O + 4.0000 SiO2  =  H4(H2SiO4)4---- +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -35.94" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H4(H2SiO4)4-4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"8.0000 H2O + 4.0000 SiO2  =  H6(H2SiO4)4-- +2.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -13.64" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	H6(H2SiO4)4-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 Al+++  =  HAlO2 +3.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -16.4329" //NEW_LINE('')// &
&"	-delta_H	144.704	kJ/mol	# 	HAlO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-230.73 kcal/mol" //NEW_LINE('')// &
&"        -analytic 4.2012e+001 1.9980e-002 -7.7847e+003 -1.5470e+001 -1.2149e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 H+ + 1.0000 Cl-  =  HCl" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -0.67" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HCl" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 4.1893e+002 1.1103e-001 -1.1784e+004 -1.6697e+002 -1.8400e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 H+ + 1.0000 ClO-  =  HClO" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +7.5692" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HClO" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 H+ + 1.0000 ClO2-  =  HClO2" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +3.1698" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HClO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"  " //NEW_LINE('')// &
&"1.0000 H+ + 1.0000 S2O3--  =  HS2O3-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k            1.0139" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	HS2O3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SO3-- + 1.0000 H+  =  HSO3-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +7.2054" //NEW_LINE('')// &
&"	-delta_H	9.33032	kJ/mol	# 	HSO3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-149.67 kcal/mol" //NEW_LINE('')// &
&"        -analytic 5.5899e+001 3.3623e-002 -5.0120e+002 -2.3040e+001 -7.8373e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 H+  =  HSO4-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +1.9791" //NEW_LINE('')// &
&"	-delta_H	20.5016	kJ/mol	# 	HSO4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-212.5 kcal/mol" //NEW_LINE('')// &
&"        -analytic 4.9619e+001 3.0368e-002 -1.1558e+003 -2.1335e+001 -1.8051e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 SiO2 + 1.0000 H2O  =  HSiO3- +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -9.9525" //NEW_LINE('')// &
&"	-delta_H	25.991	kJ/mol	# 	HSiO3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-271.88 kcal/mol" //NEW_LINE('')// &
&"        -analytic 6.4211e+001 -2.4872e-002 -1.2707e+004 -1.4681e+001 1.0853e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 K+ + 1.0000 Cl-  =  KCl" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -1.4946" //NEW_LINE('')// &
&"	-delta_H	14.1963	kJ/mol	# 	KCl" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-96.81 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.3650e+002 3.8405e-002 -4.4014e+003 -5.4421e+001 -6.8721e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 K+ + 1.0000 H+  =  KHSO4" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +0.8136" //NEW_LINE('')// &
&"	-delta_H	29.8319	kJ/mol	# 	KHSO4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-270.54 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.2620e+002 5.7349e-002 -3.3670e+003 -5.3003e+001 -5.2576e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 K+ + 1.0000 H2O  =  KOH +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -14.46" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	KOH" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 K+  =  KSO4-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +0.8796" //NEW_LINE('')// &
&"	-delta_H	2.88696	kJ/mol	# 	KSO4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-276.98 kcal/mol" //NEW_LINE('')// &
&"        -analytic 9.9073e+001 3.7817e-002 -2.1628e+003 -4.1297e+001 -3.3779e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"4.0000 Mg++ + 4.0000 H2O  =  Mg4(OH)4++++ +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           5.5    " //NEW_LINE('')// &
&"        log_k           -39.75" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Mg4(OH)4+4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 Mg++ + 1.0000 HCO3-  =  MgCO3 +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -7.3499" //NEW_LINE('')// &
&"	-delta_H	23.8279	kJ/mol	# 	MgCO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-270.57 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.3465e+002 5.5538e-002 -8.3947e+003 -9.3104e+001 -1.3106e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Mg++ + 1.0000 Cl-  =  MgCl+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -0.1349" //NEW_LINE('')// &
&"	-delta_H	-0.58576	kJ/mol	# 	MgCl+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-151.44 kcal/mol" //NEW_LINE('')// &
&"        -analytic 4.3363e+001 3.2858e-002 1.1878e+002 -2.1688e+001 1.8403e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Mg++ + 1.0000 HCO3-  =  MgHCO3+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +1.0357" //NEW_LINE('')// &
&"	-delta_H	2.15476	kJ/mol	# 	MgHCO3+" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-275.75 kcal/mol" //NEW_LINE('')// &
&"        -analytic 3.8459e+001 3.0076e-002 9.8068e+001 -1.8869e+001 1.5187e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 Mg++  =  MgSO4" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +2.4117" //NEW_LINE('')// &
&"	-delta_H	19.6051	kJ/mol	# 	MgSO4" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1355.96 kJ/mol" //NEW_LINE('')// &
&"        -analytic 1.7994e+002 6.4715e-002 -4.7314e+003 -7.3123e+001 -8.0408e+001" //NEW_LINE('')// &
! &"#       -Range:  0-200" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 H2O + 1.0000 Na+ + 1.0000 Al+++  =  NaAlO2 +4.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -23.6266" //NEW_LINE('')// &
&"	-delta_H	190.326	kJ/mol	# 	NaAlO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-277.259 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.2288e+002 3.4921e-002 -1.2808e+004 -4.6046e+001 -1.9990e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Na+ + 1.0000 HCO3-  =  NaCO3- +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           -9.8144" //NEW_LINE('')// &
&"	-delta_H	-5.6521	kJ/mol	# 	NaCO3-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-935.885 kJ/mol" //NEW_LINE('')// &
&"        -analytic 1.6939e+002 5.3122e-004 -7.6768e+003 -6.2078e+001 -1.1984e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Na+ + 1.0000 Cl-  =  NaCl" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -0.777" //NEW_LINE('')// &
&"	-delta_H	5.21326	kJ/mol	# 	NaCl" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-96.12 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.1398e+002 3.6386e-002 -3.0847e+003 -4.6571e+001 -4.8167e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Na+ + 1.0000 HCO3-  =  NaHCO3" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +0.1541" //NEW_LINE('')// &
&"	-delta_H	-13.7741	kJ/mol	# 	NaHCO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-944.007 kJ/mol" //NEW_LINE('')// &
&"        -analytic -9.0668e+001 -2.9866e-002 2.7947e+003 3.6515e+001 4.7489e+001" //NEW_LINE('')// &
! &"#       -Range:  0-200" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 SiO2 + 1.0000 Na+ + 1.0000 H2O  =  NaHSiO3 +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -8.304" //NEW_LINE('')// &
&"	-delta_H	11.6524	kJ/mol	# 	NaHSiO3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-332.74 kcal/mol" //NEW_LINE('')// &
&"        -analytic 3.6045e+001 -9.0411e-003 -6.6605e+003 -1.0447e+001 5.8415e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 Na+ + 1.0000 H2O  =  NaOH +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           -14.7948" //NEW_LINE('')// &
&"	-delta_H	53.6514	kJ/mol	# 	NaOH" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-112.927 kcal/mol" //NEW_LINE('')// &
&"        -analytic 8.7326e+001 2.3555e-002 -5.4770e+003 -3.6678e+001 -8.5489e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"1.0000 SO4-- + 1.0000 Na+  =  NaSO4-" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           +0.8200" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	NaSO4-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 H2O  =  OH- +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           3.5    " //NEW_LINE('')// &
&"        log_k           -13.9951" //NEW_LINE('')// &
&"	-delta_H	55.8146	kJ/mol	# 	OH-" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-54.977 kcal/mol" //NEW_LINE('')// &
&"        -analytic -6.7506e+001 -3.0619e-002 -1.9901e+003 2.8004e+001 -3.1033e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"1.0000 HS-  =  S-- +1.0000 H+" //NEW_LINE('')// &
&"        -llnl_gamma           5.0    " //NEW_LINE('')// &
&"        log_k           -12.9351" //NEW_LINE('')// &
&"	-delta_H	49.0364	kJ/mol	# 	S-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	32.928 kJ/mol" //NEW_LINE('')// &
&"        -analytic 9.7756e+001 3.2913e-002 -5.0784e+003 -4.1812e+001 -7.9273e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 H+  + 2.0000 SO3--  = S2O5--  + H2O" //NEW_LINE('')// &
&"        -llnl_gamma           4.0    " //NEW_LINE('')// &
&"        log_k           9.5934" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	S2O5-2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-0 kcal/mol" //NEW_LINE('')// &
&"        -analytic 0.12262E+03    0.62883E-01   -0.18005E+04   -0.50798E+02   -0.28132E+02" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"2.0000 H+ + 1.0000 SO3--  =  SO2 +1.0000 H2O" //NEW_LINE('')// &
&"        -llnl_gamma           3.0    " //NEW_LINE('')// &
&"        log_k           +9.0656" //NEW_LINE('')// &
&"	-delta_H	26.7316	kJ/mol	# 	SO2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-77.194 kcal/mol" //NEW_LINE('')// &
&"        -analytic 9.4048e+001 6.2127e-002 -1.1072e+003 -4.0310e+001 -1.7305e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&" " //NEW_LINE('')// &
&"PHASES" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"#  1122 minerals" //NEW_LINE('')// &
! &"Afwillite" //NEW_LINE('')// &
! &"        Ca3Si2O4(OH)6 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Ca++ + 6.0000 H2O" //NEW_LINE('')// &
! &"        log_k           60.0452" //NEW_LINE('')// &
! &"	-delta_H	-316.059	kJ/mol	# 	Afwillite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-1143.31 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.8353e+001 1.9014e-003 1.8478e+004 -6.6311e+000 -4.0227e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
&"Akermanite" //NEW_LINE('')// &
&"        Ca2MgSi2O7 +6.0000 H+  =  + 1.0000 Mg++ + 2.0000 Ca++ + 2.0000 SiO2 + 3.0000 H2O" //NEW_LINE('')// &
&"        log_k           45.3190" //NEW_LINE('')// &
&"	-delta_H	-288.575	kJ/mol	# 	Akermanite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-926.497 kcal/mol" //NEW_LINE('')// &
&"        -analytic -4.8295e+001 -8.5613e-003 2.0880e+004 1.3798e+001 -7.1975e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"Al" //NEW_LINE('')// &
&"       Al +3.0000 H+ +0.7500 O2  =  + 1.0000 Al+++ + 1.5000 H2O" //NEW_LINE('')// &
&"        log_k           149.9292" //NEW_LINE('')// &
&"	-delta_H	-958.059	kJ/mol	# 	Al" //NEW_LINE('')// &
&"#	Enthalpy of formation:	0 kJ/mol" //NEW_LINE('')// &
&"        -analytic -1.8752e+002 -4.6187e-002 5.7127e+004 6.6270e+001 -3.8952e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Al2(SO4)3" //NEW_LINE('')// &
&"       Al2(SO4)3  =  + 2.0000 Al+++ + 3.0000 SO4--" //NEW_LINE('')// &
&"        log_k           19.0535" //NEW_LINE('')// &
&"	-delta_H	-364.566	kJ/mol	# 	Al2(SO4)3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-3441.04 kJ/mol" //NEW_LINE('')// &
&"        -analytic -6.1001e+002 -2.4268e-001 2.9194e+004 2.4383e+002 4.5573e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Al2(SO4)3:6H2O" //NEW_LINE('')// &
&"       Al2(SO4)3:6H2O  =  + 2.0000 Al+++ + 3.0000 SO4-- + 6.0000 H2O" //NEW_LINE('')// &
&"        log_k           1.6849" //NEW_LINE('')// &
&"	-delta_H	-208.575	kJ/mol	# 	Al2(SO4)3:6H2O" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-5312.06 kJ/mol" //NEW_LINE('')// &
&"        -analytic -7.1642e+002 -2.4552e-001 2.6064e+004 2.8441e+002 4.0691e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"################################" //NEW_LINE('')// &
! &"#      ADDITIONS BY NAVAH      #" //NEW_LINE('')// &
! &"################################" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Augite" //NEW_LINE('')// &
&"        Ca.7Fe.6Mg.7Si2O6 +4.0000 H+  =  2.0 H2O + 2.0 SiO2 + .7Ca+2 + 0.6Fe+2 + 0.7Mg+2" //NEW_LINE('')// &
&"        log_k           21.00" //NEW_LINE('')// &
&"	-delta_H	-51.8523	kJ/mol	# 	Augite" //NEW_LINE('')// &
&"	-analytic 7.84710902e+00   7.21674649e-03   1.25039649e+04  -8.82692820e+00  -8.09786954e+05" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Pigeonite" //NEW_LINE('')// &
&"        Ca1.14Fe.64Mg.22Si2O6 +4.0000 H+  =  2.0 H2O + 2.0 SiO2 + 1.14Ca+2 + 0.64Fe+2 + 0.22Mg+2" //NEW_LINE('')// &
&"        log_k           21.40" //NEW_LINE('')// &
&"	-delta_H	-51.8523	kJ/mol	# 	Pigeonite" //NEW_LINE('')// &
&"	-analytic 3.92773074e+01   1.11617261e-02   1.07613145e+04  -1.98006851e+01  -7.39527557e+05" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Plagioclase" //NEW_LINE('')// &
&"        Ca.5Na.5Al1.5Si2.5O8 +6.0000 H+  =  3.0000 H2O + 2.5 SiO2 + .5 Ca+2 + 1.5 Al+3 + 0.5Na+" //NEW_LINE('')// &
&"        log_k           14.20" //NEW_LINE('')// &
&"	-delta_H	-51.8523	kJ/mol	# 	Plagioclase" //NEW_LINE('')// &
&"	-analytic -3.80343385e+01  -7.37083665e-03   1.59944487e+04   4.95599390e+00  -1.01574822e+06" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"################################" //NEW_LINE('')// &
! &"#        END NAVAHBLOCK        #" //NEW_LINE('')// &
! &"################################" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Albite" //NEW_LINE('')// &
&"        NaAlSi3O8 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 Na+ + 2.0000 H2O + 3.0000 SiO2" //NEW_LINE('')// &
&"        log_k           2.7645" //NEW_LINE('')// &
&"	-delta_H	-51.8523	kJ/mol	# 	Albite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-939.68 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.1694e+001 1.4429e-002 1.3784e+004 -7.2866e+000 -1.6136e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Analcime" //NEW_LINE('')// &
&"        Na.96Al.96Si2.04O6:H2O +3.8400 H+  =  + 0.9600 Al+++ + 0.9600 Na+ + 2.0400 SiO2 + 2.9200 H2O" //NEW_LINE('')// &
&"        log_k           6.1396" //NEW_LINE('')// &
&"	-delta_H	-75.844	kJ/mol	# 	Analcime" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-3296.86 kJ/mol" //NEW_LINE('')// &
&"        -analytic -6.8694e+000 6.6052e-003 9.8260e+003 -4.8540e+000 -8.8780e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Andradite" //NEW_LINE('')// &
&"        Ca3Fe2(SiO4)3 +12.0000 H+  =  + 2.0000 Fe+++ + 3.0000 Ca++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
&"        log_k           33.3352" //NEW_LINE('')// &
&"	-delta_H	-301.173	kJ/mol	# 	Andradite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1380.35 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.3884e+001 -2.3886e-002 1.5314e+004 -8.1606e+000 -4.2193e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Saponite-Ca" //NEW_LINE('')// &
&"        Ca.165Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.1650 Ca++ + 0.3300 Al+++ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           26.2900" //NEW_LINE('')// &
&"	-delta_H	-207.971	kJ/mol	# Saponite-Ca" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1436.51 kcal/mol" //NEW_LINE('')// &
&"        -analytic -4.6904e+001 6.2555e-003 2.2572e+004 5.3198e+000 -1.5725e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Anhydrite" //NEW_LINE('')// &
&"        CaSO4  =  + 1.0000 Ca++ + 1.0000 SO4--" //NEW_LINE('')// &
&"        log_k           -4.3064" //NEW_LINE('')// &
&"	-delta_H	-18.577	kJ/mol	# 	Anhydrite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-342.76 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.0986e+002 -7.8823e-002 5.0969e+003 8.5642e+001 7.9594e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"-Vm 46.1 # 136.14 / 2.95" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"Phillipsite" //NEW_LINE('')// &
! &"        Na0.5K0.5AlSi3O8:H2O + 4H+ = 0.5Na+ +0.5K+ + 3H2O + Al+3 + 3SiO2" //NEW_LINE('')// &
&"        Na0.5K0.5AlSi3O8:H2O + 7H2O = 0.5Na+ +0.5K+ + Al(OH)4- + 6H2O + 3SiO2" //NEW_LINE('')// &
&"        log_k           -19.874" //NEW_LINE('')// &
&"" //NEW_LINE('')// &

&"Aragonite" //NEW_LINE('')// &
&"        CaCO3 +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 HCO3-" //NEW_LINE('')// &
&"        log_k           1.9931" //NEW_LINE('')// &
&"	-delta_H	-25.8027	kJ/mol	# 	Aragonite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-288.531 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.4934e+002 -4.8043e-002 4.9089e+003 6.0284e+001 7.6644e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
!&"-Vm 34.04" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"Calcite" //NEW_LINE('')// &
&"        CaCO3 +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 HCO3-" //NEW_LINE('')// &
&"        log_k           1.8487" //NEW_LINE('')// &
&"	-delta_H	-25.7149	kJ/mol	# 	Calcite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-288.552 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.4978e+002 -4.8370e-002 4.8974e+003 6.0458e+001 7.6464e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
!&"-Vm 36.9 cm3/mol # MW (100.09 g/mol) / rho (2.71 g/cm3)" //NEW_LINE('')// &

! &"" //NEW_LINE('')// &
!  &"Calcite" //NEW_LINE('')// &
!  &"	CaCO3 = CO3-2 + Ca+2" //NEW_LINE('')// &
!  &"	-log_k	-8.48" //NEW_LINE('')// &
!  &"	-delta_h -2.297 kcal" //NEW_LINE('')// &
!  &"	-analytic	-171.9065	-0.077993	2839.319	71.595" //NEW_LINE('')// &
!  &"	-Vm 36.9 cm3/mol # MW (100.09 g/mol) / rho (2.71 g/cm3)" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
!  &"Aragonite" //NEW_LINE('')// &
!  &"	CaCO3 = CO3-2 + Ca+2" //NEW_LINE('')// &
!  &"	-log_k	-8.336" //NEW_LINE('')// &
!  &"	-delta_h -2.589 kcal" //NEW_LINE('')// &
!  &"	-analytic	-171.9773	-0.077993	2903.293	71.595" //NEW_LINE('')// &
!  &"	-Vm 34.04" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
 &"Celadonite" //NEW_LINE('')// &
&"        KMgAlSi4O10(OH)2 +6.0000 H+  =  + 1.0000 Al+++ + 1.0000 K+ + 1.0000 Mg++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           7.4575" //NEW_LINE('')// &
&"	-delta_H	-74.3957	kJ/mol	# 	Celadonite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1394.9 kcal/mol" //NEW_LINE('')// &
&"        -analytic -3.3097e+001 1.7989e-002 1.8919e+004 -2.1219e+000 -2.0588e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Clinoptilolite-Ca" //NEW_LINE('')// &
&"        Ca1.7335Al3.45Fe.017Si14.533O36:10.922H2O +13.8680 H+  =  + 0.0170 Fe+++ + 1.7335 Ca++ + 3.4500 Al+++ + 14.5330 SiO2 + 17.8560 H2O" //NEW_LINE('')// &
&"        log_k           -7.0095" //NEW_LINE('')// &
&"	-delta_H	-74.6745	kJ/mol	# 	Clinoptilolite-Ca" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-4919.84 kcal/mol" //NEW_LINE('')// &
&"        -analytic -4.4820e+001 5.3696e-002 5.4878e+004 -3.1459e+001 -7.5491e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Clinozoisite" //NEW_LINE('')// &
&"        Ca2Al3Si3O12(OH) +13.0000 H+  =  + 2.0000 Ca++ + 3.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
&"        log_k           43.2569" //NEW_LINE('')// &
&"	-delta_H	-457.755	kJ/mol	# 	Clinozoisite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1643.78 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.8690e+001 -3.7056e-002 2.2770e+004 3.7880e+000 -2.5834e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Cronstedtite-7A" //NEW_LINE('')// &
&"        Fe2Fe2SiO5(OH)4 +10.0000 H+  =  + 1.0000 SiO2 + 2.0000 Fe++ + 2.0000 Fe+++ + 7.0000 H2O" //NEW_LINE('')// &
&"        log_k           16.2603" //NEW_LINE('')// &
&"	-delta_H	-244.266	kJ/mol	# 	Cronstedtite-7A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-697.413 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.3783e+002 -7.1026e-002 1.7752e+004 8.7147e+001 2.7707e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Daphnite-14A" //NEW_LINE('')// &
&"        Fe5AlAlSi3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Fe++ + 12.0000 H2O" //NEW_LINE('')// &
&"        log_k           52.2821" //NEW_LINE('')// &
&"	-delta_H	-517.561	kJ/mol	# 	Daphnite-14A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1693.04 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.5261e+002 -6.1392e-002 2.8283e+004 5.1788e+001 4.4137e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Daphnite-7A" //NEW_LINE('')// &
&"        Fe5AlAlSi3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Fe++ + 12.0000 H2O" //NEW_LINE('')// &
&"        log_k           55.6554" //NEW_LINE('')// &
&"	-delta_H	-532.326	kJ/mol	# 	Daphnite-7A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1689.51 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.6430e+002 -6.3160e-002 2.9499e+004 5.6442e+001 4.6035e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Dawsonite" //NEW_LINE('')// &
&"        NaAlCO3(OH)2 +3.0000 H+  =  + 1.0000 Al+++ + 1.0000 HCO3- + 1.0000 Na+ + 2.0000 H2O" //NEW_LINE('')// &
&"        log_k           4.3464" //NEW_LINE('')// &
&"	-delta_H	-76.3549	kJ/mol	# 	Dawsonite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1963.96 kJ/mol" //NEW_LINE('')// &
&"        -analytic -1.1393e+002 -2.3487e-002 7.1758e+003 4.0900e+001 1.2189e+002" //NEW_LINE('')// &
! &"#       -Range:  0-200" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Diaspore" //NEW_LINE('')// &
&"        AlHO2 +3.0000 H+  =  + 1.0000 Al+++ + 2.0000 H2O" //NEW_LINE('')// &
&"        log_k           7.1603" //NEW_LINE('')// &
&"	-delta_H	-110.42	kJ/mol	# 	Diaspore" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-238.924 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.2618e+002 -3.1671e-002 8.8737e+003 4.5669e+001 1.3850e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Dicalcium_silicate" //NEW_LINE('')// &
! &"        Ca2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 2.0000 H2O" //NEW_LINE('')// &
! &"        log_k           37.1725" //NEW_LINE('')// &
! &"	-delta_H	-217.642	kJ/mol	# 	Dicalcium_silicate" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-2317.9 kJ/mol" //NEW_LINE('')// &
! &"        -analytic -5.9723e+001 -1.3682e-002 1.5461e+004 2.1547e+001 -3.7732e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Diopside" //NEW_LINE('')// &
&"        CaMgSi2O6 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 2.0000 H2O + 2.0000 SiO2" //NEW_LINE('')// &
&"        log_k           20.9643" //NEW_LINE('')// &
&"	-delta_H	-133.775	kJ/mol	# 	Diopside" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-765.378 kcal/mol" //NEW_LINE('')// &
&"        -analytic 7.1240e+001 1.5514e-002 8.1437e+003 -3.0672e+001 -5.6880e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Dolomite" //NEW_LINE('')// &
&"        CaMg(CO3)2 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 2.0000 HCO3-" //NEW_LINE('')// &
&"        log_k           2.5135" //NEW_LINE('')// &
&"	-delta_H	-59.9651	kJ/mol	# 	Dolomite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-556.631 kcal/mol" //NEW_LINE('')// &
&"        -analytic -3.1782e+002 -9.8179e-002 1.0845e+004 1.2657e+002 1.6932e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"-Vm 64.5" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"Epidote" //NEW_LINE('')// &
&"        Ca2FeAl2Si3O12OH +13.0000 H+  =  + 1.0000 Fe+++ + 2.0000 Al+++ + 2.0000 Ca++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
&"        log_k           32.9296" //NEW_LINE('')// &
&"	-delta_H	-386.451	kJ/mol	# 	Epidote" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1543.99 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.6187e+001 -3.6436e-002 1.9351e+004 3.3671e+000 -3.0319e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Ettringite" //NEW_LINE('')// &
! &"        Ca6Al2(SO4)3(OH)12:26H2O +12.0000 H+  =  + 2.0000 Al+++ + 3.0000 SO4-- + 6.0000 Ca++ + 38.0000 H2O" //NEW_LINE('')// &
! &"        log_k           62.5362" //NEW_LINE('')// &
! &"	-delta_H	-382.451	kJ/mol	# 	Ettringite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-4193 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -1.0576e+003 -1.1585e-001 5.9580e+004 3.8585e+002 1.0121e+003" //NEW_LINE('')// &
! ! &"#       -Range:  0-200" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Fayalite" //NEW_LINE('')// &
! &"        Fe2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Fe++ + 2.0000 H2O" //NEW_LINE('')// &
! &"        log_k           19.1113" //NEW_LINE('')// &
! &"	-delta_H	-152.256	kJ/mol	# 	Fayalite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-354.119 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.3853e+001 -3.5501e-003 7.1496e+003 -6.8710e+000 -6.3310e+004" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Ferrite-Ca" //NEW_LINE('')// &
! &"        CaFe2O4 +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Fe+++ + 4.0000 H2O" //NEW_LINE('')// &
! &"        log_k           21.5217" //NEW_LINE('')// &
! &"	-delta_H	-264.738	kJ/mol	# 	Ferrite-Ca" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-363.494 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -2.8472e+002 -7.5870e-002 2.0688e+004 1.0485e+002 3.2289e+002" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Foshagite" //NEW_LINE('')// &
&"        Ca4Si3O9(OH)2:0.5H2O +8.0000 H+  =  + 3.0000 SiO2 + 4.0000 Ca++ + 5.5000 H2O" //NEW_LINE('')// &
&"        log_k           65.9210" //NEW_LINE('')// &
&"	-delta_H	-359.839	kJ/mol	# 	Foshagite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1438.27 kcal/mol" //NEW_LINE('')// &
&"        -analytic 2.9983e+001 5.5272e-003 2.3427e+004 -1.3879e+001 -8.9461e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Gismondine" //NEW_LINE('')// &
&"        Ca2Al4Si4O16:9H2O +16.0000 H+  =  + 2.0000 Ca++ + 4.0000 Al+++ + 4.0000 SiO2 + 17.0000 H2O" //NEW_LINE('')// &
&"        log_k           41.7170" //NEW_LINE('')// &
&"	-delta_H	0	      	# Not possible to calculate enthalpy of reaction	Gismondine" //NEW_LINE('')// &
&"#	Enthalpy of formation:	0 kcal/mol" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Goethite" //NEW_LINE('')// &
&"        FeOOH +3.0000 H+  =  + 1.0000 Fe+++ + 2.0000 H2O" //NEW_LINE('')// &
&"        log_k           0.5345" //NEW_LINE('')// &
&"	-delta_H	-61.9291	kJ/mol	# 	Goethite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-559.328 kJ/mol" //NEW_LINE('')// &
&"        -analytic -6.0331e+001 -1.0847e-002 4.7759e+003 1.9429e+001 8.1122e+001" //NEW_LINE('')// &
! &"#       -Range:  0-200" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Greenalite" //NEW_LINE('')// &
&"        Fe3Si2O5(OH)4 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Fe++ + 5.0000 H2O" //NEW_LINE('')// &
&"        log_k           22.6701" //NEW_LINE('')// &
&"	-delta_H	-165.297	kJ/mol	# 	Greenalite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-787.778 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.4187e+001 -3.8377e-003 1.1710e+004 1.6442e+000 -4.8290e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Gyrolite" //NEW_LINE('')// &
&"        Ca2Si3O7(OH)2:1.5H2O +4.0000 H+  =  + 2.0000 Ca++ + 3.0000 SiO2 + 4.5000 H2O" //NEW_LINE('')// &
&"        log_k           22.9099" //NEW_LINE('')// &
&"	-delta_H	-82.862	kJ/mol	# 	Gyrolite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1176.55 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.4416e+001 1.4646e-002 1.6181e+004 2.3723e+000 -1.5369e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Hedenbergite" //NEW_LINE('')// &
&"        CaFe(SiO3)2 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Fe++ + 2.0000 H2O + 2.0000 SiO2" //NEW_LINE('')// &
&"        log_k           19.6060" //NEW_LINE('')// &
&"	-delta_H	-124.507	kJ/mol	# 	Hedenbergite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-678.276 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.9473e+001 1.5288e-003 1.2910e+004 2.1729e+000 -9.0058e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Hematite" //NEW_LINE('')// &
&"        Fe2O3 +6.0000 H+  =  + 2.0000 Fe+++ + 3.0000 H2O" //NEW_LINE('')// &
&"        log_k           0.1086" //NEW_LINE('')// &
&"	-delta_H	-129.415	kJ/mol	# 	Hematite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-197.72 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.2015e+002 -6.0290e-002 1.1812e+004 8.0253e+001 1.8438e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Hillebrandite" //NEW_LINE('')// &
! &"        Ca2SiO3(OH)2:0.17H2O +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 3.1700 H2O" //NEW_LINE('')// &
! &"        log_k           36.8190" //NEW_LINE('')// &
! &"	-delta_H	-203.074	kJ/mol	# 	Hillebrandite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-637.404 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -1.9360e+001 -7.5176e-003 1.1947e+004 8.0558e+000 -1.4504e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"K-Feldspar" //NEW_LINE('')// &
&"        KAlSi3O8 +4.0000 H+  =  + 1.0000 Al+++ + 1.0000 K+ + 2.0000 H2O + 3.0000 SiO2" //NEW_LINE('')// &
&"        log_k           -0.2753" //NEW_LINE('')// &
&"	-delta_H	-23.9408	kJ/mol	# 	K-Feldspar" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-949.188 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.0684e+000 1.3111e-002 1.1671e+004 -9.9129e+000 -1.5855e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Kaolinite" //NEW_LINE('')// &
&"        Al2Si2O5(OH)4 +6.0000 H+  =  + 2.0000 Al+++ + 2.0000 SiO2 + 5.0000 H2O" //NEW_LINE('')// &
&"        log_k           6.8101" //NEW_LINE('')// &
&"	-delta_H	-151.779	kJ/mol	# 	Kaolinite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-982.221 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.6835e+001 -7.8939e-003 7.7636e+003 -1.2190e+001 -3.2354e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Larnite" //NEW_LINE('')// &
! &"        Ca2SiO4 +4.0000 H+  =  + 1.0000 SiO2 + 2.0000 Ca++ + 2.0000 H2O" //NEW_LINE('')// &
! &"        log_k           38.4665" //NEW_LINE('')// &
! &"	-delta_H	-227.061	kJ/mol	# 	Larnite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-551.74 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 2.6900e+001 -2.1833e-003 1.0900e+004 -9.5257e+000 -7.2537e+004" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Laumontite" //NEW_LINE('')// &
&"        CaAl2Si4O12:4H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 4.0000 SiO2 + 8.0000 H2O" //NEW_LINE('')// &
&"        log_k           13.6667" //NEW_LINE('')// &
&"	-delta_H	-184.657	kJ/mol	# 	Laumontite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1728.66 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.1904e+000 8.1763e-003 1.9005e+004 -1.4561e+001 -1.5851e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Lawsonite" //NEW_LINE('')// &
! &"        CaAl2Si2O7(OH)2:H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 2.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
! &"        log_k           22.2132" //NEW_LINE('')// &
! &"	-delta_H	-244.806	kJ/mol	# 	Lawsonite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-1158.1 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.3995e+001 -1.7668e-002 1.0119e+004 -8.3100e+000 1.5789e+002" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &

&"Magnesite" //NEW_LINE('')// &
&"        MgCO3 +1.0000 H+  =  + 1.0000 HCO3- + 1.0000 Mg++" //NEW_LINE('')// &
&"        log_k           2.2936" //NEW_LINE('')// &
&"	-delta_H	-44.4968	kJ/mol	# 	Magnesite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-265.63 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.6665e+002 -4.9469e-002 6.4344e+003 6.5506e+001 1.0045e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Magnetite" //NEW_LINE('')// &
&"        Fe3O4 +8.0000 H+  =  + 1.0000 Fe++ + 2.0000 Fe+++ + 4.0000 H2O" //NEW_LINE('')// &
&"        log_k           10.4724" //NEW_LINE('')// &
&"	-delta_H	-216.597	kJ/mol	# 	Magnetite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-267.25 kcal/mol" //NEW_LINE('')// &
&"        -analytic -3.0510e+002 -7.9919e-002 1.8709e+004 1.1178e+002 2.9203e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Merwinite" //NEW_LINE('')// &
! &"        MgCa3(SiO4)2 +8.0000 H+  =  + 1.0000 Mg++ + 2.0000 SiO2 + 3.0000 Ca++ + 4.0000 H2O" //NEW_LINE('')// &
! &"        log_k           68.5140" //NEW_LINE('')// &
! &"	-delta_H	-430.069	kJ/mol	# 	Merwinite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-1090.8 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -2.2524e+002 -4.2525e-002 3.5619e+004 7.9984e+001 -9.8259e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Mesolite" //NEW_LINE('')// &
&"        Na.676Ca.657Al1.99Si3.01O10:2.647H2O +7.9600 H+  =  + 0.6570 Ca++ + 0.6760 Na+ + 1.9900 Al+++ + 3.0100 SiO2 + 6.6270 H2O" //NEW_LINE('')// &
&"        log_k           13.6191" //NEW_LINE('')// &
&"	-delta_H	-179.744	kJ/mol	# 	Mesolite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-5947.05 kJ/mol" //NEW_LINE('')// &
&"        -analytic 7.1993e+000 5.9356e-003 1.4717e+004 -1.3627e+001 -9.8863e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Minnesotaite" //NEW_LINE('')// &
! &"        Fe3Si4O10(OH)2 +6.0000 H+  =  + 3.0000 Fe++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
! &"        log_k           13.9805" //NEW_LINE('')// &
! &"	-delta_H	-105.211	kJ/mol	# 	Minnesotaite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-1153.37 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -1.8812e+001 1.7261e-002 1.9804e+004 -6.4410e+000 -2.0433e+006" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"Monticellite" //NEW_LINE('')// &
! &"        CaMgSiO4 +4.0000 H+  =  + 1.0000 Ca++ + 1.0000 Mg++ + 1.0000 SiO2 + 2.0000 H2O" //NEW_LINE('')// &
! &"        log_k           29.5852" //NEW_LINE('')// &
! &"	-delta_H	-195.711	kJ/mol	# 	Monticellite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-540.8 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.5730e+001 -3.5567e-003 9.0789e+003 -6.3007e+000 1.4166e+002" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Montmor-Na" //NEW_LINE('')// &
&"        Na.33Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.3300 Mg++ + 0.3300 Na+ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           2.4844" //NEW_LINE('')// &
&"	-delta_H	-93.2165	kJ/mol	# 	Montmor-Na" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1360.69 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.9601e+000 1.1342e-002 1.6051e+004 -1.4718e+001 -1.8160e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &



&"Montmor-k" //NEW_LINE('')// &
&"        K.33Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.3300 K+ + 0.3300 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           2.1423" //NEW_LINE('')// &
&"	-delta_H	-88.184	kJ/mol	# Calculated enthalpy of reaction	Montmor-K" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1362.83 kcal/mol" //NEW_LINE('')// &
&"        -analytic 8.4757e+000 1.1219e-002 1.5654e+004 -1.6833e+001 -1.8386e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &

&"Montmor-Mg" //NEW_LINE('')// &
&"        Mg.495Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.4950 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           2.3879" //NEW_LINE('')// &
&"	-delta_H	-102.608	kJ/mol	# Calculated enthalpy of reaction	Montmor-Mg" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1357.87 kcal/mol" //NEW_LINE('')// &
&"        -analytic -6.8505e+000 9.0710e-003 1.6817e+004 -1.1887e+001 -1.8323e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &

&"Montmor-Ca" //NEW_LINE('')// &
&"        Ca.165Mg.33Al1.67Si4O10(OH)2 +6.0000 H+  =  + 0.1650 Ca++ + 0.3300 Mg++ + 1.6700 Al+++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           2.4952" //NEW_LINE('')// &
&"	-delta_H	-100.154	kJ/mol	# Calculated enthalpy of reaction	Montmor-Ca" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1361.5 kcal/mol" //NEW_LINE('')// &
&"        -analytic 6.0725e+000 1.0644e-002 1.6024e+004 -1.6334e+001 -1.7982e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &



&"Muscovite" //NEW_LINE('')// &
&"        KAl3Si3O10(OH)2 +10.0000 H+  =  + 1.0000 K+ + 3.0000 Al+++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
&"        log_k           13.5858" //NEW_LINE('')// &
&"	-delta_H	-243.224	kJ/mol	# 	Muscovite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1427.41 kcal/mol" //NEW_LINE('')// &
&"        -analytic 3.3085e+001 -1.2425e-002 1.2477e+004 -2.0865e+001 -5.4692e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Natrolite" //NEW_LINE('')// &
&"        Na2Al2Si3O10:2H2O +8.0000 H+  =  + 2.0000 Al+++ + 2.0000 Na+ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
&"        log_k           18.5204" //NEW_LINE('')// &
&"	-delta_H	-186.971	kJ/mol	# 	Natrolite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-5718.56 kJ/mol" //NEW_LINE('')// &
&"        -analytic -2.7712e+001 -2.7963e-003 1.6075e+004 1.5332e+000 -9.5765e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Nontronite-Ca" //NEW_LINE('')// &
&"        Ca.165Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.1650 Ca++ + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           -11.5822" //NEW_LINE('')// &
&"	-delta_H	-38.138	kJ/mol	# 	Nontronite-Ca" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1166.7 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.6291e+001 4.3557e-003 1.0221e+004 -1.8690e+001 -1.5427e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Nontronite-H" //NEW_LINE('')// &
&"        H.33Fe2Al.33Si3.67H2O12 +6.9900 H+  =  + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           -12.5401" //NEW_LINE('')// &
&"	-delta_H	-30.452	kJ/mol	# 	Nontronite-H" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1147.12 kcal/mol" //NEW_LINE('')// &
&"        -analytic 9.7794e+001 1.4055e-002 4.7440e+003 -4.7272e+001 -1.2103e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Nontronite-K" //NEW_LINE('')// &
&"        K.33Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 K+ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           -11.8648" //NEW_LINE('')// &
&"	-delta_H	-26.5822	kJ/mol	# 	Nontronite-K" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1167.93 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.3630e+001 4.7708e-003 1.0073e+004 -1.7407e+001 -1.5803e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Nontronite-Mg" //NEW_LINE('')// &
&"        Mg.165Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.1650 Mg++ + 0.3300 Al+++ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           -11.6200" //NEW_LINE('')// &
&"	-delta_H	-41.1779	kJ/mol	# 	Nontronite-Mg" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1162.93 kcal/mol" //NEW_LINE('')// &
&"        -analytic 5.5961e+001 1.0139e-002 8.0777e+003 -3.3164e+001 -1.4031e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Nontronite-Na" //NEW_LINE('')// &
&"        Na.33Fe2Al.33Si3.67H2O12 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 Na+ + 2.0000 Fe+++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           -11.5263" //NEW_LINE('')// &
&"	-delta_H	-31.5687	kJ/mol	# 	Nontronite-Na" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1165.8 kcal/mol" //NEW_LINE('')// &
&"        -analytic 6.7915e+001 1.2851e-002 7.1218e+003 -3.7112e+001 -1.3758e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Okenite" //NEW_LINE('')// &
! &"        CaSi2O4(OH)2:H2O +2.0000 H+  =  + 1.0000 Ca++ + 2.0000 SiO2 + 3.0000 H2O" //NEW_LINE('')// &
! &"        log_k           10.3816" //NEW_LINE('')// &
! &"	-delta_H	-19.4974	kJ/mol	# 	Okenite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-749.641 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -7.7353e+001 1.5091e-002 1.3023e+004 2.1337e+001 -1.1831e+006" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Smectite-low-Fe-Mg" //NEW_LINE('')// &
&"        Ca.02Na.15K.2Fe.29Fe.16Mg.9Al1.25Si3.75H2O12 +7.0000 H+  =  + 0.0200 Ca++ + 0.1500 Na+ + 0.1600 Fe+++ + 0.2000" // &
&" K+ + 0.2900 Fe++ + 0.9000 Mg++ + 1.2500 Al+++ + 3.7500 SiO2 + 4.5000 H2O" //NEW_LINE('')// &
&"        log_k           11.0405" //NEW_LINE('')// &
&"	-delta_H	-144.774	kJ/mol	# Smectite-low-Fe-Mg" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1352.12 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.7003e+001 6.9848e-003 1.8359e+004 -6.8896e+000 -1.6637e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Prehnite" //NEW_LINE('')// &
&"        Ca2Al2Si3O10(OH)2 +10.0000 H+  =  + 2.0000 Al+++ + 2.0000 Ca++ + 3.0000 SiO2 + 6.0000 H2O" //NEW_LINE('')// &
&"        log_k           32.9305" //NEW_LINE('')// &
&"	-delta_H	-311.875	kJ/mol	# 	Prehnite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1481.65 kcal/mol" //NEW_LINE('')// &
&"        -analytic -3.5763e+001 -2.1396e-002 2.0167e+004 6.3554e+000 -7.4967e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Pseudowollastonite" //NEW_LINE('')// &
! &"        CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2" //NEW_LINE('')// &
! &"        log_k           13.9997" //NEW_LINE('')// &
! &"	-delta_H	-79.4625	kJ/mol	# 	Pseudowollastonite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-388.9 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 2.6691e+001 6.3323e-003 5.5723e+003 -1.1822e+001 -3.6038e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Pyrite" //NEW_LINE('')// &
&"        FeS2 +1.0000 H2O  =  + 0.2500 H+ + 0.2500 SO4-- + 1.0000 Fe++ + 1.7500 HS-" //NEW_LINE('')// &
&"        log_k           -24.6534" //NEW_LINE('')// &
&"	-delta_H	109.535	kJ/mol	# 	Pyrite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-41 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.4195e+002 -8.7948e-002 -6.2911e+002 9.9248e+001 -9.7454e+000" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Pyrrhotite" //NEW_LINE('')// &
&"        FeS +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HS-" //NEW_LINE('')// &
&"        log_k           -3.7193" //NEW_LINE('')// &
&"	-delta_H	-7.9496	kJ/mol	# 	Pyrrhotite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-24 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.5785e+002 -5.2258e-002 3.9711e+003 6.3195e+001 6.2012e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Quartz" //NEW_LINE('')// &
&"        SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
&"        log_k           -3.9993" //NEW_LINE('')// &
&"	-delta_H	32.949	kJ/mol	# 	Quartz" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-217.65 kcal/mol" //NEW_LINE('')// &
&"        -analytic 7.7698e-002 1.0612e-002 3.4651e+003 -4.3551e+000 -7.2138e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"-Vm 22.67" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
! &"Rankinite" //NEW_LINE('')// &
! &"        Ca3Si2O7 +6.0000 H+  =  + 2.0000 SiO2 + 3.0000 Ca++ + 3.0000 H2O" //NEW_LINE('')// &
! &"        log_k           51.9078" //NEW_LINE('')// &
! &"	-delta_H	-302.089	kJ/mol	# 	Rankinite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-941.7 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -9.6393e+001 -1.6592e-002 2.4832e+004 3.2541e+001 -9.4630e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
&"Saponite-Mg" //NEW_LINE('')// &
&"        #Mg3Ca.165Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + .1650 Ca ++ + 3 Mg++ + " //NEW_LINE('')// &
&"#3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"	Mg3.165Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 3.1650 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           26.2523" //NEW_LINE('')// &
&"	-delta_H	-210.822	kJ/mol	# 	Saponite-Mg" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1432.79 kcal/mol" //NEW_LINE('')// &
&"        -analytic 9.8888e+000 1.4320e-002 1.9418e+004 -1.5259e+001 -1.3716e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Saponite-Na" //NEW_LINE('')// &
&"        Na.33Mg3Al.33Si3.67O10(OH)2 +7.3200 H+  =  + 0.3300 Al+++ + 0.3300 Na+ + 3.0000 Mg++ + 3.6700 SiO2 + 4.6600 H2O" //NEW_LINE('')// &
&"        log_k           26.3459" //NEW_LINE('')// &
&"	-delta_H	-201.401	kJ/mol	# 	Saponite-Na" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1435.61 kcal/mol" //NEW_LINE('')// &
&"        -analytic -6.7611e+001 4.7327e-003 2.3586e+004 1.2868e+001 -1.6493e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Scolecite" //NEW_LINE('')// &
&"        CaAl2Si3O10:3H2O +8.0000 H+  =  + 1.0000 Ca++ + 2.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
&"        log_k           15.8767" //NEW_LINE('')// &
&"	-delta_H	-204.93	kJ/mol	# 	Scolecite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-6048.92 kJ/mol" //NEW_LINE('')// &
&"        -analytic 5.0656e+001 -3.1485e-003 1.0574e+004 -2.5663e+001 -5.2769e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"SiO2(am)" //NEW_LINE('')// &
&"       SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
&"        log_k           -2.7136" //NEW_LINE('')// &
&"	-delta_H	20.0539	kJ/mol	# 	SiO2(am)" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-214.568 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.2109e+000 7.0767e-003 2.3634e+003 -3.4449e+000 -4.8591e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Siderite" //NEW_LINE('')// &
&"        FeCO3 +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HCO3-" //NEW_LINE('')// &
&"        log_k           -0.1920" //NEW_LINE('')// &
&"	-delta_H	-32.5306	kJ/mol	# 	Siderite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-179.173 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.5990e+002 -4.9361e-002 5.4947e+003 6.3032e+001 8.5787e+001" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"-Vm 29.2" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"Smectite-high-Fe-Mg" //NEW_LINE('')// &
&"#        Ca.025Na.1K.2Fe++.5Fe+++.2Mg1.15Al1.25Si3.5H2O12 +8.0000 H+  =  + 0.0250 Ca++ + 0.1000 Na+ + 0.2000 Fe+++ + 0.2000 K+ " //&
&"+ 0.5000 Fe++ + 1.1500 Mg++ + 1.2500 Al+++ + 3.5000 SiO2 + 5.0000 H2O" //NEW_LINE('')// &
&"        Ca.025Na.1K.2Fe.5Fe.2Mg1.15Al1.25Si3.5H2O12 +8.0000 H+  =  + 0.0250 Ca++ + 0.1000 Na+ + 0.2000 Fe+++ + 0.2000 K+ + 0.5000 Fe++ " //&
&"+ 1.1500 Mg++ + 1.2500 Al+++ + 3.5000 SiO2 + 5.0000 H2O " //NEW_LINE('')// &
&"        log_k           17.4200" //NEW_LINE('')// &
&"	-delta_H	-199.841	kJ/mol	# 	Smectite-high-Fe-Mg" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1351.39 kcal/mol" //NEW_LINE('')// &
&"        -analytic -9.6102e+000 1.2551e-003 1.8157e+004 -7.9862e+000 -1.3005e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"Stilbite" //NEW_LINE('')// &
&"        Ca1.019Na.136K.006Al2.18Si6.82O18:7.33H2O +8.7200 H+  =  + 0.0060 K+ + 0.1360 Na+ " //&
&"+ 1.0190 Ca++ + 2.1800 Al+++ + 6.8200 SiO2 + 11.6900 H2O" //NEW_LINE('')// &
&"        log_k           1.0545" //NEW_LINE('')// &
&"	-delta_H	-83.0019	kJ/mol	# 	Stilbite" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-11005.7 kJ/mol" //NEW_LINE('')// &
&"        -analytic -2.4483e+001 3.0987e-002 2.8013e+004 -1.5802e+001 -3.4491e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
! &"Tobermorite-9A" //NEW_LINE('')// &
! &"        Ca5Si6H6O20 +10.0000 H+  =  + 5.0000 Ca++ + 6.0000 SiO2 + 8.0000 H2O" //NEW_LINE('')// &
! &"        log_k           69.0798" //NEW_LINE('')// &
! &"	-delta_H	-329.557	kJ/mol	# 	Tobermorite-9A" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-2375.42 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -6.3384e+001 1.1722e-002 3.8954e+004 1.2268e+001 -2.8681e+006" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"Tremolite" //NEW_LINE('')// &
! &"        Ca2Mg5Si8O22(OH)2 +14.0000 H+  =  + 2.0000 Ca++ + 5.0000 Mg++ + 8.0000 H2O + 8.0000 SiO2" //NEW_LINE('')// &
! &"        log_k           61.2367" //NEW_LINE('')// &
! &"	-delta_H	-406.404	kJ/mol	# 	Tremolite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-2944.04 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 8.5291e+001 4.6337e-002 3.9465e+004 -5.4414e+001 -3.1913e+006" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &

! &"Fe" //NEW_LINE('')// &
! &"        Fe +2.0000 H+ +0.5000 O2  =  + 1.0000 Fe++ + 1.0000 H2O" //NEW_LINE('')// &
! &"        log_k           59.0325" //NEW_LINE('')// &
! &"	-delta_H	-372.029	kJ/mol		Fe" //NEW_LINE('')// &
! &"# Enthalpy of formation:	0 kcal/mol" //NEW_LINE('')// &
! &"        -analytic -6.2882e+001 -2.0379e-002 2.0690e+004 2.3673e+001 3.2287e+002" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &


&"Ferrihydrite" //NEW_LINE('')// &
&"        Fe(OH)3 + 3H+ = Fe+3 + 3H2O" //NEW_LINE('')// &
&"        log_k	3.191" //NEW_LINE('')// &
&"	delta_h	-73.374	kJ" //NEW_LINE('')// &

&"Fe3(OH)8" //NEW_LINE('')// &
&"        Fe3(OH)8 + 8H+ = 2Fe+3 + Fe+2 + 8H2O" //NEW_LINE('')// &
&"        log_k   20.222" //NEW_LINE('')// &
&"	delta_h -0      kcal" //NEW_LINE('')// &

&"Talc" //NEW_LINE('')// &
&"        Mg3Si4O10(OH)2 +6.0000 H+  =  + 3.0000 Mg++ + 4.0000 H2O + 4.0000 SiO2" //NEW_LINE('')// &
&"        log_k           21.1383" //NEW_LINE('')// &
&"-delta_H	-148.737	kJ/mol	# Calculated enthalpy of reaction	Talc" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1410.92 kcal/mol" //NEW_LINE('')// &
&"        -analytic 1.1164e+001 2.4724e-002 1.9810e+004 -1.7568e+001 -1.8241e+006" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&"Chlorite(14A)" //NEW_LINE('')// &
&"        Mg5Al2Si3O10(OH)8 + 16H+ = 5Mg+2 + 2Al+3 + 3.0 SiO2 + 12H2O" //NEW_LINE('')// &
&"        log_k           68.38" //NEW_LINE('')// &
&"delta_h -151.494 kcal" //NEW_LINE('')// &

&"" //NEW_LINE('')// &
&"Fe(OH)2" //NEW_LINE('')// &
&"        Fe(OH)2 +2.0000 H+  =  + 1.0000 Fe++ + 2.0000 H2O" //NEW_LINE('')// &
&"        log_k           13.9045" //NEW_LINE('')// &
&"	-delta_H	-95.4089	kJ/mol		Fe(OH)2" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-568.525 kJ/mol" //NEW_LINE('')// &
&"        -analytic -8.6666e+001 -1.8440e-002 7.5723e+003 3.2597e+001 1.1818e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &

! CHLORITE MINERALS

&"Chamosite-7A" //NEW_LINE('')// &
&"        Fe2Al2SiO5(OH)4 +10.0000 H+  =  + 1.0000 SiO2 + 2.0000 Al+++ + 2.0000 Fe++ + 7.0000 H2O" //NEW_LINE('')// &
&"        log_k           32.8416" //NEW_LINE('')// &
&"-delta_H	-364.213	kJ/mol	# Calculated enthalpy of reaction	Chamosite-7A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-902.407 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.5581e+002 -7.0890e-002 2.4619e+004 9.1789e+001 3.8424e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&"Clinochlore-14A" //NEW_LINE('')// &
&"        Mg5Al2Si3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Mg++ + 12.0000 H2O" //NEW_LINE('')// &
&"        log_k           67.2391" //NEW_LINE('')// &
&"-delta_H	-612.379	kJ/mol	# Calculated enthalpy of reaction	Clinochlore-14A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-2116.96 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.0441e+002 -6.2268e-002 3.5388e+004 6.9239e+001 5.5225e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&"Clinochlore-7A" //NEW_LINE('')// &
&"        Mg5Al2Si3O10(OH)8 +16.0000 H+  =  + 2.0000 Al+++ + 3.0000 SiO2 + 5.0000 Mg++ + 12.0000 H2O" //NEW_LINE('')// &
&"        log_k           70.6124" //NEW_LINE('')// &
&"-delta_H	-628.14	kJ/mol	# Calculated enthalpy of reaction	Clinochlore-7A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-2113.2 kcal/mol" //NEW_LINE('')// &
&"        -analytic -2.1644e+002 -6.4187e-002 3.6548e+004 7.4123e+001 5.7037e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&"Ripidolite-14A" //NEW_LINE('')// &
&"        Mg3Fe2Al2Si3O10(OH)8 +16.00 H+  =  + 2.00 Al+++ + 2.00 Fe++ + 3.00 Mg++ + 3.00 SiO2 + 12.00 H2O" //NEW_LINE('')// &
&"        log_k           60.9638" //NEW_LINE('')// &
&"-delta_H	-572.472	kJ/mol	# Calculated enthalpy of reaction	Ripidolite-14A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1947.87 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.8376e+002 -6.1934e-002 3.2458e+004 6.2290e+001 5.0653e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"Ripidolite-7A" //NEW_LINE('')// &
&"        Mg3Fe2Al2Si3O10(OH)8 +16.00 H+  =  + 2.00 Al+++ + 2.00 Fe++ + 3.00 Mg++ + 3.00 SiO2 + 12.00 H2O" //NEW_LINE('')// &
&"        log_k           64.3371" //NEW_LINE('')// &
&"-delta_H	-586.325	kJ/mol	# Calculated enthalpy of reaction	Ripidolite-7A" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-1944.56 kcal/mol" //NEW_LINE('')// &
&"        -analytic -1.9557e+002 -6.3779e-002 3.3634e+004 6.7057e+001 5.2489e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &

&"Fe(OH)3" //NEW_LINE('')// &
&"        Fe(OH)3 +3.0000 H+  =  + 1.0000 Fe+++ + 3.0000 H2O" //NEW_LINE('')// &
&"        log_k           5.6556" //NEW_LINE('')// &
&"	-delta_H	-84.0824	kJ/mol		Fe(OH)3" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-823.013 kJ/mol" //NEW_LINE('')// &
&"        -analytic -1.3316e+002 -3.1284e-002 7.9753e+003 4.9052e+001 1.2449e+002" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &
&"" //NEW_LINE('')// &

! &"Troilite" //NEW_LINE('')// &
! &"        FeS +1.0000 H+  =  + 1.0000 Fe++ + 1.0000 HS-" //NEW_LINE('')// &
! &"        log_k           -3.8184" //NEW_LINE('')// &
! &"	-delta_H	-7.3296	kJ/mol	# 	Troilite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-101.036 kJ/mol" //NEW_LINE('')// &
! &"        -analytic -1.6146e+002 -5.3170e-002 4.0461e+003 6.4620e+001 6.3183e+001" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"Wollastonite" //NEW_LINE('')// &
! &"        CaSiO3 +2.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 SiO2" //NEW_LINE('')// &
! &"        log_k           13.7605" //NEW_LINE('')// &
! &"	-delta_H	-76.5756	kJ/mol	# 	Wollastonite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-389.59 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 3.0931e+001 6.7466e-003 5.1749e+003 -1.3209e+001 -3.4579e+005" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"Xonotlite" //NEW_LINE('')// &
! &"        Ca6Si6O17(OH)2 +12.0000 H+  =  + 6.0000 Ca++ + 6.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
! &"        log_k           91.8267" //NEW_LINE('')// &
! &"	-delta_H	-495.457	kJ/mol	# 	Xonotlite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-2397.25 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 1.6080e+003 3.7309e-001 -2.2548e+004 -6.2716e+002 -3.8346e+002" //NEW_LINE('')// &
! &"#       -Range:  0-200" //NEW_LINE('')// &
! &"" //NEW_LINE('')// &
! &"Zoisite" //NEW_LINE('')// &
! &"        Ca2Al3(SiO4)3OH +13.0000 H+  =  + 2.0000 Ca++ + 3.0000 Al+++ + 3.0000 SiO2 + 7.0000 H2O" //NEW_LINE('')// &
! &"        log_k           43.3017" //NEW_LINE('')// &
! &"	-delta_H	-458.131	kJ/mol	# 	Zoisite" //NEW_LINE('')// &
! &"#	Enthalpy of formation:	-1643.69 kcal/mol" //NEW_LINE('')// &
! &"        -analytic 2.5321e+000 -3.5886e-002 1.9902e+004 -6.2443e+000 3.1055e+002" //NEW_LINE('')// &
! &"#       -Range:  0-300" //NEW_LINE('')// &

! MINS FROM ANOTHER DATABASE
&"Vermiculite-Na" //NEW_LINE('')// &
! &"Na0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850Na+  +0.850Al+3  -9.400H+  +3.150H4(SiO4)  -0.600H2O" //NEW_LINE('')// &
&"Na0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850Na+  +0.850Al+3  -9.400H+  +6.30H2O + 3.15SiO2  -0.600H2O" //NEW_LINE('')// &
&"  log_k  40.17  #" //NEW_LINE('')// &
&"  delta_h  -354.987   kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -6139.206  kJ/mol  07VIE" //NEW_LINE('')// &
  
  
&"Lepidocrocite" //NEW_LINE('')// &
&"FeOOH = +1.000Fe+3  -3.000H+  +2.000H2O  " //NEW_LINE('')// &
&"  log_k  0.75   #98DIA in 98CHI" //NEW_LINE('')// &
&"  delta_h  -64.26  kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -556.4  kJ/mol  " //NEW_LINE('')// &
  
&"Vermiculite-K" //NEW_LINE('')// &
&"K0.85Mg3Si3.15Al0.85O10(OH)2 = +3.000Mg+2  +0.850K+  +0.850Al+3  -9.400H+  +6.30H2O + 3.15SiO2  -0.600H2O " //NEW_LINE('')// &
&"  log_k  36.86  #" //NEW_LINE('')// &
&"  delta_h  -331.639   kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -6172.584  kJ/mol  07VIE" //NEW_LINE('')// &

&"Saponite-K" //NEW_LINE('')// &
&"K0.33Mg3Al0.33Si3.67O10(OH)2 = +3.000Mg+2  +0.330K+  +0.330Al+3  -7.320H+  +7.34H2O + 3.67SiO2  -2.680H2O" //NEW_LINE('')// &
&"  log_k  28.1   #" //NEW_LINE('')// &
&"  delta_h  -252.497   kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -6005.94   kJ/mol  07VIE" //NEW_LINE('')// &
  
&"Vermiculite-Ca" //NEW_LINE('')// &
&"Ca0.43Mg3Si3.14Al0.86O10(OH)2 = +0.430Ca+2  +3.000Mg+2  +0.860Al+3  -9.440H+  +6.28H2O + 3.14SiO2  -0.560H2O" //NEW_LINE('')// &
&"  log_k  40.68  #" //NEW_LINE('')// &
&"  delta_h  -378.219   kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -6147.38   kJ/mol  07VIE" //NEW_LINE('')// &


&"Vermiculite-Mg" //NEW_LINE('')// &
&"Mg0.43Mg3Si3.14Al0.86O10(OH)2 = +3.430Mg+2  +0.860Al+3  -9.440H+  +6.28H2O + 3.14SiO2  -0.560H2O" //NEW_LINE('')// &
&"  log_k  38.8   #" //NEW_LINE('')// &
&"  delta_h  -377.469   kJ/mol  #" //NEW_LINE('')// &
&"  # Enthalpy of formation:    -6115.45   kJ/mol  07VIE" //NEW_LINE('')// &

&"Chalcedony" //NEW_LINE('')// &
&"        SiO2  =  + 1.0000 SiO2" //NEW_LINE('')// &
&"        log_k           -3.7281" //NEW_LINE('')// &
&"	-delta_H	31.4093	kJ/mol	# Calculated enthalpy of reaction	Chalcedony" //NEW_LINE('')// &
&"#	Enthalpy of formation:	-217.282 kcal/mol" //NEW_LINE('')// &
&"        -analytic -9.0068e+000 9.3241e-003 4.0535e+003 -1.0830e+000 -7.5077e+005" //NEW_LINE('')// &
&"#       -Range:  0-300" //NEW_LINE('')// &


&""! //NEW_LINE('')// &











write(*,*) "testing..."

!--------------INITIALIZE ALL PROCESSORS

! process #0 is the root process
root_process = 0

! initialize a process
call MPI_INIT ( ierr )

! find out the process ID and how many processes were started so far
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

! print out current processor id
write(*,*) "my_id:", my_id
write(*,*) " "

! what to do if you are the master processor
if (my_id .eq. root_process) then 

!--------------DO STUFF WITH THE MASTER PROCESSOR

! initialize domain geometry
call init()




permeability0 = permeability
phiCoarse = .1
phi = .1
phi = 0.1
phi0 = phi


! fill coordinate arrays
do i = 1,(xn/cell)
	do ii = 1,(yn/cell)
		medium(i,ii,6) = x(i*cell) 
		medium(i,ii,7) = y(ii*cell)
	end do
end do

! boundary & initial condtions for flow equations
psi=0.0
frac6 = 0.0
! frac6(:,1) = 1.0e-3
! frac6(:,2) = -1.0e-3
temp6 = 280.0
psi(1,1:yn) = bcyPsi(1,1:yn)
psi(xn,1:yn) = bcyPsi(2,1:yn)
psi(1:xn,1) = bcxPsi(1:xn,1)
psi(1:xn,yn) = bcxPsi(1:xn,2)
psi = 0.0!dx*scope*1000.0 - dx*scope*0.0
rho = rho_fluid
visc = viscosity

h = param_tsw

!no slant IC, 11/15/15
do i=1,xn
	do ii=1,yn
		
		 ! (.68-(0.00133*i*1.0))
		 
		h(i,1) = param_tsw!480.0 + ( first-(0.0016*i*factor)) * dy/1.8
		
		h(i,ii) = param_tsw!480.0*(h(i,1)/h(1,1)) + (param_tsw-(480.0*(h(i,1)/h(1,1)) ))*((-y_min)+y(ii))/((-y_min))
		!h(i,ii) = h(i,ii) - (400.0 + (param_tsw-400.0)*((-y_min)-max(param_o_rhs,param_o))/((max(param_o_rhs,param_o))))
	end do
end do

! do i=1,xn
! 	do ii=1,yn
! 		if (y(ii) .ge. -max(param_o_rhs,param_o)) then
! 			h(i,ii) = 275.0
! 		end if
! 	end do
! end do


! put initial values into array for file
!hmat(1:xn,1:yn) = h
!psimat(1:xn,1:yn) = psi
!umat(1:xn,1:yn) = u
!vmat(1:xn,1:yn) = v

uTransport = 0.0
vTransport = 0.0
uCoarse = 0.0
vCoarse = 0.0


!-------RESTART STEP!

if (restart .eq. 1) then
	

		! FLUID DYNAMIC TIME SLICES
!
! 		OPEN(UNIT=10, FILE='/data/navah/ic_saturday/h_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')
! 		OPEN(UNIT=11, FILE='/data/navah/ic_saturday/psi_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')

		OPEN(UNIT=10, FILE=trim(path) // 'h.txt')
		OPEN(UNIT=11, FILE=trim(path) // 'psi.txt')

		!OPEN(UNIT=12, FILE=trim(path) // 'permeability.txt')
		!OPEN(UNIT=13, FILE=trim(path) // 'phi.txt')
		DO i = 1,yn
			READ(10,*) (hTrans(i,ii),ii=1,xn)
			READ(11,*) (psiTrans(i,ii),ii=1,xn)
			!READ(12,*) (permeabilityTrans(i,ii),ii=1,xn)
			!READ(13,*) (phiTrans(i,ii),ii=1,xn)
		END DO
		h = transpose(hTrans)
		psi = transpose(psiTrans)
		!permeability = transpose(permeabilityTrans)
		!phi = transpose(phiTrans)

!
! 		! FLUID DYNAMIC MATRICES
! 		OPEN(UNIT=14, FILE=trim(path) // 'transfer/hMat.txt')
! 		OPEN(UNIT=15, FILE=trim(path) // 'transfer/psiMat.txt')
! 		OPEN(UNIT=16, FILE=trim(path) // 'transfer/permMat.txt')
! 		OPEN(UNIT=80, FILE=trim(path) // 'transfer/uMat.txt')
! 		OPEN(UNIT=81, FILE=trim(path) // 'transfer/vMat.txt')
! 		!OPEN(UNIT=17, FILE=trim(path) // 'phiMat.txt')
! 		DO i = 1,yn*tn/(mstep*ar)
! 			!write(*,*) "check"
! 			READ(14,*) (checkMatF(i,ii),ii=1,xn)
! 		END DO
! 		hmat = transpose(checkMatF)
!
! 		DO i = 1,yn*tn/(mstep*ar)
! 			READ(15,*) (checkMatF(i,ii),ii=1,xn)
! 		END DO
! 		psimat = transpose(checkMatF)
!
! 		DO i = 1,yn*tn/(mstep*ar)
! 			READ(16,*) (checkMatF(i,ii),ii=1,xn)
! 			!READ(17,*) (checkMatF(i,ii),ii=1,xn)
! 			!phiMat = transpose(checkMatF)
! 		END DO
! 		permmat = transpose(checkMatF)
!
! 		DO i = 1,yn*tn/(mstep*ar)
! 			READ(80,*) (checkMatF(i,ii),ii=1,xn)
! 		END DO
! 		umat = transpose(checkMatF)
!
! 		DO i = 1,yn*tn/(mstep*ar)
! 			READ(81,*) (checkMatF(i,ii),ii=1,xn)
! 		END DO
! 		vmat = transpose(checkMatF)
!
! 		write(*,*) 'fluid dynamic matrices'
!
!
! 		! PRIMARY TIME SLICE
! 		OPEN(UNIT=18, FILE=trim(path) // 'checkpoint/pri_glass.txt')
! 		DO i = 1,yn/cell
! 		 	READ(18,*) (primaryTrans(i,ii,5),ii=1,xn/cell)
! 		END DO
! 		primary(:,:,n) = transpose(primaryTrans(:,:,5))
! 		write(*,*) 'primary time slices'
!
! 		! PRIMARY MATRIX/CES
! 		OPEN(UNIT=19, FILE=trim(path) // 'pri_glass.txt')
! 		DO i = 1,yn*tn/(cell*mstep*ar)
! 		 	READ(19,*) (checkMat(j,ii),ii=1,xn/cell)
! 		END DO
! 		primaryMat(:,:,5) = transpose(checkMat)
! 		write(*,*) 'primary matrices'
!
!
! 		! SOLUTION TIME SLICES
! 		OPEN(UNIT=21, FILE=trim(path) // 'checkpoint/sol_ph.txt')
! 		OPEN(UNIT=22, FILE=trim(path) // 'checkpoint/sol_alk.txt')
! 		OPEN(UNIT=23, FILE=trim(path) // 'checkpoint/sol_w.txt')
! 		OPEN(UNIT=24, FILE=trim(path) // 'checkpoint/sol_c.txt')
! 		OPEN(UNIT=25, FILE=trim(path) // 'checkpoint/sol_ca.txt')
! 		OPEN(UNIT=26, FILE=trim(path) // 'checkpoint/sol_mg.txt')
! 		OPEN(UNIT=27, FILE=trim(path) // 'checkpoint/sol_na.txt')
! 		OPEN(UNIT=28, FILE=trim(path) // 'checkpoint/sol_k.txt')
! 		OPEN(UNIT=29, FILE=trim(path) // 'checkpoint/sol_fe.txt')
! 		OPEN(UNIT=30, FILE=trim(path) // 'checkpoint/sol_s.txt')
! 		OPEN(UNIT=31, FILE=trim(path) // 'checkpoint/sol_si.txt')
! 		OPEN(UNIT=32, FILE=trim(path) // 'checkpoint/sol_cl.txt')
! 		OPEN(UNIT=33, FILE=trim(path) // 'checkpoint/sol_al.txt')
! 		OPEN(UNIT=34, FILE=trim(path) // 'checkpoint/sol_hco3.txt')
! 		OPEN(UNIT=35, FILE=trim(path) // 'checkpoint/sol_co3.txt')
! 		do n=1,15
! 			DO i = 1,yn/cell
! 			 	READ(20+n,*) (soluteTrans(i,ii,n),ii=1,xn/cell)
! 			END DO
! 			solute(:,:,n) = transpose(soluteTrans(:,:,n))
! 		end do
! 		write(*,*) 'solution time slices'
!
! 		! SOLUTION MATRICES
! 		OPEN(UNIT=36, FILE=trim(path) // 'sol_ph.txt')
! 		OPEN(UNIT=37, FILE=trim(path) // 'sol_alk.txt')
! 		OPEN(UNIT=38, FILE=trim(path) // 'sol_w.txt')
! 		OPEN(UNIT=39, FILE=trim(path) // 'sol_c.txt')
! 		OPEN(UNIT=40, FILE=trim(path) // 'sol_ca.txt')
! 		OPEN(UNIT=41, FILE=trim(path) // 'sol_mg.txt')
! 		OPEN(UNIT=42, FILE=trim(path) // 'sol_na.txt')
! 		OPEN(UNIT=43, FILE=trim(path) // 'sol_k.txt')
! 		OPEN(UNIT=44, FILE=trim(path) // 'sol_fe.txt')
! 		OPEN(UNIT=45, FILE=trim(path) // 'sol_s.txt')
! 		OPEN(UNIT=46, FILE=trim(path) // 'sol_si.txt')
! 		OPEN(UNIT=47, FILE=trim(path) // 'sol_cl.txt')
! 		OPEN(UNIT=48, FILE=trim(path) // 'sol_al.txt')
! 		OPEN(UNIT=49, FILE=trim(path) // 'sol_hco3.txt')
! 		OPEN(UNIT=50, FILE=trim(path) // 'sol_co3.txt')
! 		do n = 1,15
! 			DO i = 1,yn*tn/(cell*mstep*ar)
! 			 	READ(35+n,*) (checkMat(j,ii),ii=1,xn/cell)
! 			END DO
! 			soluteMat(:,:,n) = transpose(checkMat)
! 		end do
! 		write(*,*) 'solution matrices'
!
!
! 		! SECONDARY TIME SLICES
! 		do i = 1,g_sec/2
! 			i_unit = 50 + i
! 		        if (i < 10) then
! 					write(s_i,'(i1)') i
! 		        else
! 					write(s_i,'(i2)') i
! 		        end if
! 				OPEN(UNIT=i_unit, FILE=trim(path)//'checkpoint/sec'//trim(s_i)//'.txt')
! 				DO j = 1,yn/cell
! 				 	READ(i_unit,*) (secondaryTrans(j,ii,i),ii=1,xn/cell)
! 				END DO
! 				secondary(:,:,i) = transpose(secondaryTrans(:,:,i))
! 		end do
! 		write(*,*) 'secondary time slices'
!
! 		! SECONDARY MATRICES
! 		do i = 1,g_sec/2
! 			i_unit = i_unit+1
! 		        if (i < 10) then
! 					write(s_i,'(i1)') i
! 		        else
! 					write(s_i,'(i2)') i
! 		        end if
! 				OPEN(UNIT=i_unit, FILE=trim(path)//'sec'//trim(s_i)//'.txt')
! 				DO j = 1,yn*tn/(cell*mstep*ar)
! 				 	READ(i_unit,*) (checkMat(j,ii),ii=1,xn/cell)
! 				END DO
! 				secondaryMat(:,:,i) = transpose(checkMat)
! 		end do
! 		write(*,*) 'secondary matrices'
!
! 		! SATURATION MATRICES
! 		do i = 1,g_sec/2
! 			i_unit = i_unit+1
! 		        if (i < 10) then
! 					write(s_i,'(i1)') i
! 		        else
! 					write(s_i,'(i2)') i
! 		        end if
! 				OPEN(UNIT=i_unit, FILE=trim(path)//'sat'//trim(s_i)//'.txt')
! 				DO j = 1,yn*tn/(cell*mstep*ar)
! 				 	READ(i_unit,*) (checkMat(j,ii),ii=1,xn/cell)
! 				END DO
! 				saturationMat(:,:,i) = transpose(checkMat)
! 		end do
! 		write(*,*) 'saturation matrices'


! end restart loop
end if




permx = partial((phi/(permeability)),xn,yn,dx,dy,1)
permy = partial((phi/(permeability)),xn,yn,dx,dy,2)

outerBand = make_band(permeability,phi,permx,permy,rho)
outerBand = band(outerBand,2*((yn/2)-2) + 1,longP)

! permx = 0.0
! permy = 0.0





	
!-- DYNAMICS LOOP
	! this is the main loop that does all the solving for tn timesteps
	do j = crashstep, tn
	!do j = 2, 50


! 			! JDF
! 			! move altered cells
! 			if ( mod(j,nint(tn*cell/(xn*2.3))) .eq. 0) then
! 				write(*,*) "shift"
! 				do i =2,(xn/cell)
! 	 				primaryShift(i,:,:) = primary(i-1,:,:)
! 	 				secondaryShift(i,:,:) = secondary(i-1,:,:)
! 				end do
! 				primaryShift(1,:,5) = 9.67700
! 				secondaryShift(1,:,:) = 0.0
! 				primary = primaryShift
! 				secondary = secondaryShift
! 			end if

	! 		! solve thermal energy equation

if (j .eq. crashstep) then
	write(*,*) "STARTING STEP:" , j
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
	OPEN(UNIT=8, status = 'replace', FILE=trim(path_final) // 'dynamicStep.txt')
	write(*,*) "opened"
	write(8,*) 0
	close ( 8 )
	!velocitiesCoarse0 = 0.0
end if

	
if (mod(j,mstep) .eq. 0) then
	write(*,*) "STARTING STEP:" , j
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
	OPEN(UNIT=8, status = 'replace', FILE=trim(path_final) // 'dynamicStep.txt')
	write(*,*) "opened"
	write(8,*) j/mstep
	close ( 8 )
	!velocitiesCoarse0 = 0.0
end if



		
if (mod(j,mstep/10) .eq. 0) then
	write(*,*) "WRITING TO DYNAMIC SUB STEP"
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
		OPEN(UNIT=88, status = 'replace', FILE=trim(path_final) // 'dynamicSubStep.txt')
		write(88,*) mod(j,mstep)
		close ( 88 )
end if

! write(*,*) " "
! write(*,*) " "
! write(*,*) " "
! write(*,*) " "
! write(*,*) "j step:" , j

if (restart .ne. 1) then
	
	

		dt_bit = dt

			
		h = h_next(h, psi,rho,phi,u,v,frac6,temp6,dt_bit)
		h = h_bc(h)
		

		! short outcrop outflow condition
		do jj=2,yn-1
		do i=1,xn

			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 25.0)  ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if

			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 12.5) ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if
			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 17.5) ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if

		end do
		end do

		! inner boundaries outflow condtion
		do jj=2,yn-1
		do i=1,xn

				if ((mask(i,jj) .eq. 5.0) .and. (mask(i,jj) .eq. 5.0) .and. (u(i,jj) .gt. 0.0)) then
					h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
				end if
				if ((mask(i,jj) .eq. 12.5) .and. (u(i,jj) .gt. 0.0)) then
					h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
				end if

				! since these are also a 50.0 bc, they maybe aren't wilcock's but isothermal??
				if ((mask(i,jj) .eq. 10.0) .and. (mask(i,jj) .eq. 10.0).and. (u(i,jj) .lt. 0.0)) then
					h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
				end if
				if ((mask(i,jj) .eq. 17.5) .and. (u(i,jj) .lt. 0.0)) then
					h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
				end if

		end do
		end do
		
		! insulating boundaries during spinup phase
		if (iter .eq. 0) then
			do ii=2,yn-1
					if ((mask(f_index1,ii) .eq. 6.0) .or. (mask(f_index1,ii) .eq. 6.5) .or. (mask(f_index1,ii) .eq. 6.1) .or. (mask(f_index1,ii) .eq. 6.05)) then
						temp6(ii,2) = (4.0/3.0)*h(f_index1,ii) - (1.0/3.0)*h(f_index1+1,ii)
					end if
					if ((mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5) .or. (mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.05)) then
						temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
					end if
					
					if ((mask(f_index1,ii) .eq. 6.5)) then
						temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
					end if

					if ((mask(f_index1-1,ii) .eq. 3.5)) then
						temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
					end if
			end do
		end if
		
		! fracture temperature set by vertical heat advection after spinup phase
		if ((iter .eq. 1)) then
			
! 		temp6_mid = temp6
!
				do ii=2,yn-1
					if ((mask(f_index1-1,ii) .eq. 3.05)) then
						temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
					end if
				end do
				
				do ii=2,yn-1
					if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
						temp6(ii,1) = temp6(ii-1,1) + (dy*lambdaMat(f_index1-1,ii)/(dx*4179.0*frac6(ii,1))) * (h(f_index1-1,ii) - h(f_index1-2,ii))
					end if
				end do

! 				temp6 = temp6_mid
!
! 				do ii=2,yn-1
! 					if ((mask(f_index1-1,ii) .eq. 3.1)) then
! 						temp6_a(ii) = 0.0
! 						temp6_b(ii) = 1.0 + frac6(ii,1)
! 					end if
! 				end do
!

! 				if (mod(j,50) .eq. 0) then
!
! 				do jj=1,10000000
!
! 					temp6 = temp6_mid
! 					do ii=2,yn-1
!
! 						if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! 						end if
! ! 						if (mask(f_index1-1,ii) .eq. 3.0) then
! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! ! 						end if
! ! 						if (mask(f_index1-1,ii) .eq. 3.5) then
! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! ! 						end if
! 					end do
!
! 				end do
! ! 				write(*,*) "mod 50"
! ! 				write(*,*) temp6(:,1)
! ! 				write(*,*) " "
!
! ! ! TESTING SOMETHING
! ! do ii=2,yn-1
! ! 	if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
! ! 		temp6_mid(ii,1) = 328.0
! ! 	end if
! ! end do
!
! 				end if



				! top of fracture is no heat flow out? wilcox condition?
				do ii=yn/2,yn-1
						if ((mask(f_index1,ii) .eq. 6.5)) then
							temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
						end if

						if ((mask(f_index1-1,ii) .eq. 3.5)) then
							temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
						end if
				end do

		end if


		do ii=1,yn-1
			temp6(ii,2) = temp6(ii,1)
		end do
		
		
		
		rho = rho_next(h)
		visc = visc_next(h)
			

 

		! solve streamfunction-vorticity equation

		rhs0 = -1.0*partial((rho/rho_fluid)-1.0,xn,yn,dx,dy,1)

		frac6_last = frac6

		! find properties at top and base of fracture
		do jj=yn/2,yn-1
			if (maskP(f_index1-1,jj) .eq. 3.1) then
				h_base = temp6(jj-1,1) ! h(f_index1-1,jj)!
				y_base = y(jj)
				jj_base = jj
			end if
			if (maskP(f_index1-1,jj) .eq. 3.5) then
				h_top = param_tsw ! temp6(jj,1)!
				y_top = y(jj)
				jj_top = jj
			end if
		end do
		
		h_adjacent = sum(temp6(jj_base:jj_top,1))/(jj_top-jj_base+1)
		
! 		write(*,*) "tube temps"
! 		write(*,*) temp6(jj_base:jj_top,1)
!
! 		write(*,*) "h_adjacent"
! 		write(*,*) h_adjacent
!
! 		write(*,*) "jj_top"
! 		write(*,*) jj_top
!
! 		write(*,*) "jj_base"
! 		write(*,*) jj_base
!
! 		write(*,*) "rho_adjacent"
! 		write(*,*) rho_one(h_adjacent)
!

		if ((j .gt. spinup-1)) then
			iter = 1
		end if
		
! 		if ((j .eq. spinup-1)) then
! 			h_adjacent = 273.0
! 		end if
		
		
		if (iter .eq. 1) then
			!write(*,*) "loop of fracj 1"
			do jj=yn/2,yn-1
				if ((maskP(f_index1,jj) .eq. 6.0) .or. (maskP(f_index1,jj) .eq. 6.5) .or. (maskP(f_index1,jj) .eq. 6.1)) then
					!frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_fluid*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*y_top - rho_one(h_base)*y_base - rho_fluid*(y_top-y_base))
					frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_fluid*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*(y_top-y_top) - rho_one(h_base)*(y_base-y_top) - rho_one(param_tsw)*(y_top-y_base))
					!frac6(jj,1) = param_f_por*psi(f_index1-1,jj) + (dx*(12.0*(1.0e-16)*1.0)/(param_f_dx*param_f_dx*param_f_dx))
					
! 					write(*,*) "full frac6 jj 1"
! 					write(*,*) frac6(jj,1)
!
! 					write(*,*) "extra bit"
! 					write(*,*) (dx*(12.0*(1.0e-16)*1.0)/(param_f_dx*param_f_dx*param_f_dx))
! 					write(*,*) " "
				end if
			end do
			!write(*,*) " "

			psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
			psi = psi_bc(psi)
		
			do jj=yn/2,yn-1
				do i=1,xn
					if ((maskP(i,jj) .eq. 50.0) .and. (i .lt. f_index1)) then
						psi(i,jj+1) = maxval(frac6(:,1))
					end if
! 					if ((maskP(i,jj) .eq. 3.5)) then
! 						psi(i,jj+1) = maxval(frac6(:,1))
! 					end if
				end do
			end do

						
		end if
				

		! run with frac6 = 0 during spinup phase
		if (iter .eq. 0) then

			psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
			psi = psi_bc(psi)

		end if
				
		
			

			
			 
end if ! end if restart .ne. 1

			! get velocities from streamfunction
			velocities0 = velocities(psi)
			u = phi*velocities0(1:xn,1:yn)/(rho_fluid)!*phi
			v = phi*velocities0(1:xn,yn+1:2*yn)/(rho_fluid)!*phi
			
! 			do ii=1,yn
! 				do i=2,xn-2
! 					if ((maskP(i,ii) .eq. 6.0)) then
! 						v(i,ii) = -1.0*(frac6(ii,2) - frac6(ii,1))/(rho(i,ii)*dx)
! 					end if
! 				end do
! 			end do
			
! 			do jj=2,yn-1
! 				do i=1,xn
! 					if ((maskP(i,jj) .eq. 2.5)) then
! 						u(i,jj) = 0.0
! 						v(i,jj) = 0.0
! 					end if
! 					if ((maskP(i,jj) .eq. 7.5)) then
! 						u(i,jj) = 0.0
! 						v(i,jj) = 0.0
! 					end if
! 				end do
! 			end do
			

			
!if ((j .ge. thresh) .and. (mod(j,mstep) .eq. 0)) then

!do i = 1,cstep
! 				iso(:,:,1) = solute_next(iso(:,:,1),u,v,1.0D+00)
! 				iso(:,:,1) = iso(:,:,1)*exp(-(3.949e-12)*dt)
!
! 				iso(:,:,2) = solute_next(iso(:,:,2),u,v,1.0D+00)
!end do

if ((param_trace .eq. 1)) then
				
				isoTrace(:,:,1) = particles_next(isoTrace(:,:,1),u,v,1.0D+00,ison,particle_sat)
				
				do jj = 1,cstep
				isoTrace(3,:,1) = isoTrace(3,:,1)*exp(-(3.949e-12)*dt/cstep)
				end do
				
! 				inertTrace(:,:) = particles_next(inertTrace,u,v,10.0D+00,inertn,inert_sat)
end if
				

			
		

	! interpolate fine grid onto coarse grid
	do jj = 1,yn/cell
	do i = 1,xn/cell
		reactiveCoarse(i,jj) = reactive(i*cell,jj*cell)
		medium(i,jj,4) = reactiveCoarse(i,jj)
		!psiCoarse(i,jj) = psi(i*cell,jj*cell)
	end do
	end do
	
! 	write(*,*) "about to call velocitiesCoarse"
    velocitiesCoarse0 = velocities(psi)!velocitiesCoarse(psiCoarse)
    uCoarse = velocitiesCoarse0(1:xn/cell,1:yn/cell)
    vCoarse = velocitiesCoarse0(1:xn/cell,yn/cell+1:2*yn/cell)



				

! 	write(*,*) "about to start mstep loop"
	
	! things only done every mth timestep go here
	if (mod(j,mstep) .eq. 0) then
		
! 		! OUTER BAND THING
! 		outerBand = make_band(permeability,phi)
! 		outerBand = band(outerBand,2*(yn-2) + 1,(xn-2)*(yn-2))
! 		permx = partial_edge((phi/(permeability)),xn,yn,dx,dy,1)
! 		permy = partial_edge((phi/(permeability)),xn,yn,dx,dy,2)
		
		! make coarse grid average velocities
! 		uTransport = (uTransport/(mstep*wscale))
! 		vTransport = (vTransport/(mstep*wscale))
		!uTransport(1,:) = 0.0
		!uTransport(xn/cell,:) = 0.0

! 		write(*,*) "SOLUTE COURANT NUMBERS"
! 		!write(*,*) (dt*mstep/(cstep*dx*cell))*maxval(abs(uTransport))
! 		!write(*,*) (dt*mstep/(cstep*dy*cell))*maxval(abs(vTransport))
! 		write(*,*) dt*mstep/(cstep*dy*dy*cell*cell)

!NOT RUNNING TRANSPORT RIGHT NOW, JUST CHEMISTRY



!-- ISO
u_1d = param_f_dx*param_f_dx*param_f_dx*param_f_por*grav*22.0/(viscosity*12.0*param_h)
! cstep_int = dx*cstep/(maxval(u/phi(1,1))*mstep)
! cstep_num = 6.28e10/(cstep_int*mstep)

cstep_int = 6.28e11/tn
cstep_num = 10.0*cstep_int*mstep*maxval(abs(v/phi(1,1)))/dy

if (cstep_num .le. 1) then
	cstep_num = 1
end if

write(*,*) "max(u/phi(1,1))"
write(*,*) maxval(abs(u/phi(1,1)))
write(*,*) "max(v/phi(1,1))"
write(*,*) maxval(abs(v/phi(1,1)))
write(*,*) "u_1d"
write(*,*) u_1d
write(*,*) "cstep_int"
write(*,*) cstep_int
write(*,*) "cstep_num"
write(*,*) cstep_num

!
		!do i = 1,cstep_num
		do i = 1,cstep

! 			iso(:,:,1) = solute_next(iso(:,:,1),u/phi(1,1),v/phi(1,1),1.0D+00)
! 			iso(:,:,1) = iso(:,:,1)*exp(-(3.85e-12)*mstep*cstep_int/(cstep_num))
!
! 			iso(:,:,2) = solute_next(iso(:,:,2),u/phi(1,1),v/phi(1,1),1.0D+00)

! 			n=1 ! pH
! 			solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
  			n=2 ! alk
  			solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			!n=3 ! water?
 	 		!solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=4 ! c
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
			n=5 ! ca
			solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
			n=6 ! mg
	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=7 ! na
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=8 ! k
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=9 ! fe
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=10 ! s
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=11 ! si
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=12 ! cl
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
 			n=13 ! al
 	 		solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
! 			do jj = 1,yn/cell
! 			do ii = 1,xn/cell
! 				solute(ii,jj,13) = max(solute(ii,jj,13),1.0e-8)
! 			end do
! 			end do
		end do


		
		write(*,*) "about to stretch everything for transport to nodes"
			

		! stretch everything out
		!hLong = reshape(h(1:xn-1:cell,1:yn-1:cell), (/(xn/cell)*(yn/cell)/)) ! for cell > 1
		hLong = reshape(transpose(h(1:xn,(yn/2)+1:yn)), (/xn*(yn/2)/)) ! for cell = 1
		do i = 1,g_pri
			bit_thing = reshape(primary(:,(yn/2)+1:yn,i),(/xn, (yn/2)/))
			priLong(:,i) = reshape(transpose(bit_thing), (/xn*(yn/2)/))
		end do
		
		do i = 1,g_sec/2
			bit_thing = reshape(secondary(:,(yn/2)+1:yn,i),(/xn, (yn/2)/))
			secLong(:,i) = reshape(transpose(bit_thing), (/xn*(yn/2)/))
		end do
		
		do i = 1,g_sol
			bit_thing = reshape(solute(:,(yn/2)+1:yn,i),(/xn, (yn/2)/))
			solLong(:,i) = reshape(transpose(bit_thing), (/xn*(yn/2)/))
		end do
		
		do i = 1,g_med
			bit_thing = reshape(medium(:,(yn/2)+1:yn,i),(/xn, (yn/2)/))
			medLong(:,i) = reshape(transpose(bit_thing), (/xn*(yn/2)/))
		end do
		
	
		!--------------MESSAGE DISTRIBUTING FROM MASTER TO SLAVES
		do an_id = 1, num_procs - 1
		
			! put number of rows in vector here for hLong
			num_rows = xn*(yn/2)
			avg_rows_per_process = num_rows / (num_procs-1)
	        start_row = ( (an_id-1) * avg_rows_per_process) + 1
	        end_row = start_row + avg_rows_per_process - 1
	        if (an_id .eq. (num_procs - 1)) end_row = num_rows
	        num_rows_to_send = (end_row - start_row + 1)
		
			! send size of temperature array chunk to processor an_id
	        call MPI_SEND( num_rows_to_send, 1, MPI_INTEGER, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
			! send timestep size to processor an_id
	        call MPI_SEND( dt, 1, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
	
		
			! send temperature array chunk to processor an_id
	        call MPI_SEND( hLong(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
			! send primary array chunk to processor an_id
			do ii = 1,g_pri
				priLongBit = priLong(:,ii)
	        	call MPI_SEND( priLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
		
			! send secondary array chunk to processor an_id
			do ii = 1,g_sec/2
				secLongBit = secLong(:,ii)
	        	call MPI_SEND( secLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
		
			! send solute array chunk to processor an_id
			do ii = 1,g_sol
				solLongBit = solLong(:,ii)
	        	call MPI_SEND( solLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
			
			! send medium array chunk to processor an_id
			do ii = 1,g_med
				medLongBit = medLong(:,ii)
	        	call MPI_SEND( medLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
			
			write(*,*) "DONE SENDING TO PROCESSOR", an_id

		end do
		
		write(*,*) "DONE SENDING TO ALL PROCESSORS!!!"

		!--------------MESSAGE RECEIVING FROM SLAVE PROCESSORS
		do an_id = 1, num_procs - 1
		
			! get the size of each chunk again
			num_rows = xn*(yn/2)
			avg_rows_per_process = num_rows / (num_procs-1)
	        start_row = ( (an_id-1) * avg_rows_per_process) + 1
	        end_row = start_row + avg_rows_per_process - 1
	        if (an_id .eq. (num_procs - 1)) end_row = num_rows
	        num_rows_to_send = (end_row - start_row + 1)
		
			! receive primary chunk
			do ii = 1,g_pri
				! receive it
				call MPI_RECV( priLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				priLong(start_row:end_row,ii) = priLocal(1:num_rows_to_send,ii)
			end do
		
			! receive secondary chunk
			do ii = 1,g_sec/2
				! receive it
				call MPI_RECV( secLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				secLong(start_row:end_row,ii) = secLocal(1:num_rows_to_send,ii)
			end do
		
			! receive solute chunk
			do ii = 1,g_sol
				! receive it
				call MPI_RECV( solLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				solLong(start_row:end_row,ii) = solLocal(1:num_rows_to_send,ii)
			end do
			
			! receive medium chunk
			do ii = 1,g_med
				! receive it
				call MPI_RECV( medLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				medLong(start_row:end_row,ii) = medLocal(1:num_rows_to_send,ii)
			end do

			write(*,*) "DONE RECEIVING FROM PROCESSOR", an_id
		
		end do
		
		write(*,*) "DONE RECEIVING FROM ALL PROCESSORS!!!"
		write(*,*) "now onto stretching:"
		
		!--------------MASTER PROCESSOR SAVES OUTPUT TO BE WRITTEN TO FILE
	
		! put stretched vectors back into 2d arrays
		
		do i = 1,g_pri
			bit_thing_t = reshape(priLong(:,i),(/yn/2, xn/))
			bit_thing = transpose(bit_thing_t)
			bit_thing_t1 = reshape(bit_thing,(/xn,yn/2/))
			primary(:,(yn/2)+1:,i) = bit_thing_t1
		end do
		
		do i = 1,g_sec
			bit_thing_t = reshape(secLong(:,i),(/yn/2, xn/))
			bit_thing = transpose(bit_thing_t)
			bit_thing_t1 = reshape(bit_thing,(/xn,yn/2/))
			secondary(:,(yn/2)+1:,i) = bit_thing_t1
		end do
		
		do i = 1,g_sol
			bit_thing_t = reshape(solLong(:,i),(/yn/2, xn/))
			bit_thing = transpose(bit_thing_t)
			bit_thing_t1 = reshape(bit_thing,(/xn,yn/2/))
			solute(:,(yn/2)+1:,i) = bit_thing_t1
		end do
		
		do i = 1,g_med
			bit_thing_t = reshape(medLong(:,i),(/yn/2, xn/))
			bit_thing = transpose(bit_thing_t)
			bit_thing_t1 = reshape(bit_thing,(/xn,yn/2/))
			medium(:,(yn/2)+1:,i) = bit_thing_t1
		end do
! 		primary(:,yn/2:yn,g_pri) = reshape(priLong,(/(yn/2), xn, g_pri/))
! 		secondary(:,yn/2:yn,g_sec) = reshape(secLong,(/(yn/2), xn, g_sec/))
! 		solute(:,yn/2:yn,g_sol) = reshape(solLong,(/(yn/2), xn, g_sol/))
! 		medium(:,yn/2:yn,g_med) = reshape(medLong,(/(yn/2), xn, g_med/))
		!saturation(:,yn/2:yn,g_sec) = reshape(secLong(:,g_sec/2:),(/(yn/2), xn, g_sec/2/))
		
		write(*,*) "DONE STRETCHING!!!"
		!yep = write_matrix ( xn, yn, real(primary(:,:,5), kind = 4), trim(path) // 'glass1.txt' )

		! add timestep's output to output arrays
		if (mod(j,mstep*ar) .eq. 0) then
			 rhsmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rhs0
			 rhomat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rho
			 hmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = h
			 psimat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = psi
			 umat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = u
			 vmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = v
			 permmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permeability
			 permxMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permx
			 permyMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permy

			 
			 psiCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = psiCoarse
			 uCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = uTransport
			 vCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = vTransport
	 		! reset coarse grid velocities for next timestep
	 		uTransport = 0.0
	 		vTransport = 0.0
			 
			 primaryMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = primary
			 secondaryMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:)= secondary
			 soluteMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = solute
			 mediumMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = medium
			 isoMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = iso
			 saturationMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = saturation
			 
			 isoTraceMat(:,1+ison*(j/(mstep*ar)-1):ison*(j/(mstep*ar)),:) = isoTrace
			 inertTraceMat(:,1+inertn*(j/(mstep*ar)-1):inertn*(j/(mstep*ar))) = inertTrace

		 
! 		 ! get new porosity
 		 phiCoarse = 0.1 !medium(:,:,1) ! 1.0
		 
! 		 ! transfer new porosity to fine grid
 		 do i=1,xn/cell
 			 do ii=1,yn/cell
 				 !phi((i-1)*cell+1:(i)*cell,(ii-1)*cell+1:(ii)*cell) = phiCoarse(i,ii)
 			 end do
 		 end do
!
! 		 ! update permeability based on porosity change
 		 do i = 1,xn
 			 do ii = 1,yn
 				 !permeability(i,ii) = permeability0(i,ii)
 				 !permeability(i,ii) = permeability0(i,ii)*((phi(i,ii)/phi0(i,ii))**3)
 			 end do
 		 end do
		 
		 write(*,*) "DONE UNTIL NEXT mstep!!!"
		 
		 

	

	
	
	
	
!
! !--------------WRITE CHECKPOINT TO FILE
!
! !yep = write_vec ( 1, real(crashstep,kind=4), trim(path) // 'checkpoint/crashstep.txt' )
!
! yep = write_matrix ( xn, yn, real(permeability,kind=4), trim(path) // 'checkpoint/permeability.txt' )



if (restart .ne. 1) then



yep = write_matrix ( xn, yn, real(rhs0,kind=4), trim(path) // 'rhs.txt' )
yep = write_matrix ( xn, yn, real(psi,kind=4), trim(path) // 'psi.txt' )
yep = write_matrix ( xn, yn, real(h,kind=4), trim(path) // 'h.txt' )

yep = write_matrix ( xn, yn, real(u,kind=4), trim(path) // 'u.txt' )
yep = write_matrix ( xn, yn, real(v,kind=4), trim(path) // 'v.txt' )








!yep = write_matrix(xn,yn,real(psi,kind=4),'/data/navah/ic_saturday_400/psi_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')
!yep = write_matrix(xn,yn,real(h,kind=4),'/data/navah/ic_saturday_400/h_'// trim(param_o_string) //"_"// trim(param_o_rhs_string) //'.txt')


! 		OPEN(UNIT=10, FILE='/data/navah/ic_tsw/h_'// trim(param_tsw_string) //'.txt')
! 		OPEN(UNIT=11, FILE='/data/navah/ic_tsw/psi_ic.txt')


!
!
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,1),kind=4), trim(path) // 'checkpoint/sol_ph.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,2),kind=4), trim(path) // 'checkpoint/sol_w.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,3),kind=4), trim(path) // 'checkpoint/sol_alk.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,4),kind=4), trim(path) // 'checkpoint/sol_c.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,5),kind=4), trim(path) // 'checkpoint/sol_ca.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,6),kind=4), trim(path) // 'checkpoint/sol_mg.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,7),kind=4), trim(path) // 'checkpoint/sol_na.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,8),kind=4), trim(path) // 'checkpoint/sol_k.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,9),kind=4), trim(path) // 'checkpoint/sol_fe.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,10),kind=4), trim(path) // 'checkpoint/sol_s.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,11),kind=4), trim(path) // 'checkpoint/sol_si.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,12),kind=4), trim(path) // 'checkpoint/sol_cl.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,13),kind=4), trim(path) // 'checkpoint/sol_al.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,14),kind=4), trim(path) // 'checkpoint/sol_hco3.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(solute(:,:,15),kind=4), trim(path) // 'checkpoint/sol_co3.txt' )
!
!
! yep = write_matrix ( xn/cell, yn/cell, real(primary(:,:,1),kind=4), trim(path) // 'checkpoint/pri_feldspar.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(primary(:,:,1),kind=4), trim(path) // 'checkpoint/pri_augite.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(primary(:,:,1),kind=4), trim(path) // 'checkpoint/pri_pigeonite.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(primary(:,:,1),kind=4), trim(path) // 'checkpoint/pri_magnetite.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(primary(:,:,1),kind=4), trim(path) // 'checkpoint/pri_glass.txt' )
!
! yep = write_matrix ( xn/cell, yn/cell, real(medium(:,:,1),kind=4), trim(path) // 'checkpoint/med_phi.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(medium(:,:,2),kind=4), trim(path) // 'checkpoint/med_s_sp.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(medium(:,:,3),kind=4), trim(path) // 'checkpoint/med_v_water.txt' )
! yep = write_matrix ( xn/cell, yn/cell, real(medium(:,:,4),kind=4), trim(path) // 'checkpoint/med_kgw.txt' )
!
!
! do i = 1,g_sec/2
!         if (i < 10) then
! 			write(s_i,'(i1)') i
!         else
! 			write(s_i,'(i2)') i
!         end if
! 		yep = write_matrix(xn/cell, yn/cell, real(secondary(:,:,i),kind=4), trim(path)//'checkpoint/sec'//trim(s_i)//'.txt')
! end do
!
!
!
!
!
! write(*,*) "written big step to file"






!--------------WRITE EVERYTHING TO FILE


! yep = write_matrix ( xn, yn, real(mask, kind = 4), trim(path) // 'transfer/mask.txt' )
! yep = write_vec ( xn, real(x,kind=4), trim(path) // 'transfer/x.txt' )
! yep = write_vec ( yn, real(y,kind=4), trim(path) // 'transfer/y.txt' )
! !yep = write_vec ( tn, real(t, kind=4), 't.txt' )
! yep = write_matrix ( xn, yn*tn/(mstep*ar), real(hmat, kind = 4), trim(path) // 'transfer/hMat.txt' )
! yep = write_matrix ( xn, yn*tn/(mstep*ar), real(psimat,kind=4), trim(path) // 'transfer/psiMat.txt' )
! yep = write_matrix ( xn, yn*tn/(mstep*ar), real(umat, kind = 4), trim(path) // 'transfer/uMat.txt' )
! yep = write_matrix ( xn, yn*tn/(mstep*ar), real(vmat,kind=4), trim(path) // 'transfer/vMat.txt' )
! yep = write_matrix ( xn, yn*tn/(mstep*ar), real(permmat,kind=4), trim(path) // 'transfer/permMat.txt' )
! yep = write_matrix ( xn, yn,real(lambdaMat,kind=4), trim(path) // 'transfer/lambdaMat.txt' )

yep = write_matrix ( yn, 2, real(frac6, kind = 4), trim(path) // 'frac6.txt' )
yep = write_matrix ( yn, 2, real(temp6, kind = 4), trim(path) // 'temp6.txt' )
yep = write_matrix ( xn, yn/2, real(mask(:,(yn/2)+1:), kind = 4), trim(path) // 'mask.txt' )
yep = write_matrix ( xn, yn/2, real(maskP(:,(yn/2)+1:), kind = 4), trim(path) // 'maskP.txt' )
yep = write_vec ( xn, real(x,kind=4), trim(path) // 'x.txt' )
yep = write_vec ( yn/2, real(y(yn/2:),kind=4), trim(path) // 'y.txt' )
!yep = write_matrix (xn*tn/(mstep*ar), yn, real(rhsmat, kind = 4), trim(path) // 'rhsMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(hmat(:,(yn/2)+1:), kind = 4), trim(path) // 'hMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(psimat(:,(yn/2)+1:),kind=4), trim(path) // 'psiMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(rhomat(:,(yn/2)+1:),kind=4), trim(path) // 'rhoMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(umat(:,(yn/2)+1:), kind = 4), trim(path) // 'uMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(vmat(:,(yn/2)+1:),kind=4), trim(path) // 'vMat.txt' )
!yep = write_matrix (xn*tn/(mstep*ar), yn, real(permmat,kind=4), trim(path) // 'permMat.txt' )
! yep = write_matrix (xn*tn/(mstep*ar), yn, real(permxMat,kind=4), trim(path) // 'permxMat.txt' )
! yep = write_matrix (xn*tn/(mstep*ar), yn, real(permyMat,kind=4), trim(path) // 'permyMat.txt' )
yep = write_matrix ( xn, yn/2,real(lambdaMat(:,(yn/2)+1:),kind=4), trim(path) // 'lambdaMat.txt' )

!yep = write_matrix ( xn, yn, real(rho,kind=4), trim(path) // 'rho.txt' )
!yep = write_matrix ( xn, yn, real(visc,kind=4), trim(path) // 'visc.txt' )
yep = write_matrix ( xn, yn/2,real(permeability(:,(yn/2)+1:),kind=4), trim(path) // 'permeability.txt' )
yep = write_matrix ( xn, yn/2,real(phi(:,(yn/2)+1:),kind=4), trim(path) // 'phi.txt' )
!yep = write_matrix ( xn/cell, yn/cell, real(phiCoarse,kind=4), trim(path) // 'phiCoarse.txt' )
!yep = write_matrix ( xn, yn,real(reactive,kind=4), trim(path) // 'reactive.txt' )
!yep = write_matrix ( xn/cell, yn/cell,real(reactiveCoarse,kind=4), trim(path) // 'reactiveCoarse.txt' )


! yep = write_matrix ( 3, ison*tn/(cell*mstep*ar), real(isoTraceMat(1:3,:,1),kind=4), trim(path) // 'isopart_1.txt' )
! yep = write_matrix ( 3, inertn*tn/(cell*mstep*ar), real(inertTraceMat(1:3,:),kind=4), trim(path) // 'inert_trace.txt' )
! coarse variables
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/cell, real(uCoarseMat,kind=4), trim(path) // 'uCoarseMat.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/cell, real(vCoarseMat,kind=4), trim(path) // 'vCoarseMat.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/cell, real(psiCoarseMat,kind=4), trim(path) // 'psiCoarseMat.txt' )

!if (maxval(medium(:,:,5)) .eq. 1.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVING TIME !!!!!!!!!!!!!!!!!!!!!!!!

end if ! end write only if restart ne 1



if (maxval(medium(:,:,5)) .eq. 1.0) then



! solute concentrations
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,1),kind=4), trim(path) // 'z_sol_ph.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,3),kind=4), trim(path) // 'z_sol_w.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,2),kind=4), trim(path) // 'z_sol_alk.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,4),kind=4), trim(path) // 'z_sol_c.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,5),kind=4), trim(path) // 'z_sol_ca.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,6),kind=4), trim(path) // 'z_sol_mg.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,7),kind=4), trim(path) // 'z_sol_na.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,8),kind=4), trim(path) // 'z_sol_k.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,9),kind=4), trim(path) // 'z_sol_fe.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,10),kind=4), trim(path) // 'z_sol_s.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,11),kind=4), trim(path) // 'z_sol_si.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,12),kind=4), trim(path) // 'z_sol_cl.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,13),kind=4), trim(path) // 'z_sol_al.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,14),kind=4), trim(path) // 'sol_hco3.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(soluteMat(:,(yn/2)+1:,15),kind=4), trim(path) // 'sol_co3.txt' )

! primary minerals
! yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(primaryMat(:,(yn/2)+1:,1),kind=4), trim(path) // 'pri_feldspar.txt' )
! yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(primaryMat(:,(yn/2)+1:,2),kind=4), trim(path) // 'pri_augite.txt' )
! yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(primaryMat(:,(yn/2)+1:,3),kind=4), trim(path) // 'pri_pigeonite.txt' )
! yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(primaryMat(:,(yn/2)+1:,4),kind=4), trim(path) // 'pri_magnetite.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(primaryMat(:,(yn/2)+1:,5),kind=4), trim(path) // 'z_pri_glass.txt' )

! medium properties
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(mediumMat(:,(yn/2)+1:,1),kind=4), trim(path) // 'z_med_phi.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(mediumMat(:,(yn/2)+1:,2),kind=4), trim(path) // 'med_s_sp.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(mediumMat(:,(yn/2)+1:,3),kind=4), trim(path) // 'z_med_v_water.txt' )
!yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(mediumMat(:,(yn/2)+1:,4),kind=4), trim(path) // 'med_reactive.txt' )







write(*,*) "precipitated"
do i = 1,g_sec/2
	if (maxval(secondaryMat(:,:,i)) .gt. 0.0) then
		write(*,*) i
	end if
        if (i < 10) then
			write(s_i,'(i1)') i
        else
			write(s_i,'(i2)') i
        end if
		yep = write_matrix(xn*tn/(cell*mstep*ar),yn/2,real(secondaryMat(:,(yn/2)+1:,i),kind=4),trim(path)//'z_sec'//trim(s_i)//'.txt')
end do


write(*,*) ""
write(*,*) "supersaturated"
do i = 1,g_sec/2
	if (maxval(saturationMat(:,:,i)) .gt. 0.0) then
		write(*,*) i
	end if
        if (i < 10) then
			write(s_i,'(i1)') i
        else
			write(s_i,'(i2)') i
        end if
		!yep = write_matrix(xn*tn/(cell*mstep*ar),yn/cell,real(saturationMat(:,:,i),kind=4),trim(path)//'sat'//trim(s_i)//'.txt')
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVING TIME !!!!!!!!!!!!!!!!!!!!!!!!
!end if

!
! write(*,*) "just wrote to file dynamics, step:" , j
! !yep = write_matrix( 2, 1, real((/j*1.0, tn/), kind = 4), trim(path) // 'dynamicStep.txt' )
! OPEN(UNIT=8, status = 'replace', FILE=trim(path) // 'dynamicStep.txt')
! write(*,*) "opened"
! write(8,*) j
! close ( 8 )

end if ! end write if maxval cells on == 1

	 	end if
		! end mstep*ar loop
		
		



if (restart .eq. 1) then
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(isoMat(:,(yn/2)+1:,1),kind=4), trim(iso_path) // 'iso_c14.txt' )
yep = write_matrix ( xn*tn/(cell*mstep*ar), yn/2, real(isoMat(:,(yn/2)+1:,2),kind=4), trim(iso_path) // 'iso_control.txt' )
end if

		
		
end if 
! end mstep timestep loop, finally


end do 
! end all timestep loop




write(*,*) " "
write(*,*) "ALL DONE!"
write(*,*) tn
write(*,*) "steps"
write(*,*) tn/mstep
write(*,*) "msteps"




! what to do if you are a slave processor
else
	

	
	

	
	!--------------SLAVE PROCESSOR RECEIVES MESSAGE
	
	! message receiving has to happen every mth step
	do jj = 1, tn/mstep
		! here is a slave process, each process must receive a chunk of the h array and 
		! take the local mean, print it, send it back.
		
		! receive size of temperature array chunk
		call MPI_RECV ( num_rows_to_receive, 1 , MPI_INTEGER, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive timestep size
		call MPI_RECV ( dt_local, 1 , MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive temperature array chunk, save in local hLocal
		call MPI_RECV ( hLocal, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		num_rows_received = num_rows_to_receive

		! receive primary array chunk, save in local priLocal
		do ii = 1,g_pri
			call MPI_RECV ( priLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			priLocal(:,ii) = priLocalBit
		end do
	
		! receive secondary array chunk, save in local secLocal
		do ii = 1,g_sec/2
			call MPI_RECV ( secLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			secLocal(:,ii) = secLocalBit
		end do
	
		! receive solute chunk, save in local solLocal
		do ii = 1,g_sol
			call MPI_RECV ( solLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			solLocal(:,ii) = solLocalBit
		end do
		
		! receive medium chunk, save in local solLocal
		do ii = 1,g_med
			call MPI_RECV ( medLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			medLocal(:,ii) = medLocalBit
		end do
		
		

	
		!--------------SLAVE PROCESSOR RUNS GEOCHEMICAL MODEL

		! slave processor loops through each coarse cell
		do m=1,num_rows_to_receive

			! TRANSPORT ONLY
			!medLocal(m,5) = 0.0
			! TRANSPORT ONLY
		if (medLocal(m,5) .eq. 1.0) then
                       !write(*,*) medLocal(m,6:7)
					   
			! run the phreeqc alteration model
			!alt0 = alt_next(hLocal(m),dt_local*mstep,priLocal(m,:), &
			!		    secLocal(m,:),solLocal(m,:),medLocal(m,:))


! BEGIN MIGRATION INSANITY

primary3 = priLocal(m,:)
secondary3 = secLocal(m,:)
solute3 = solLocal(m,:)
medium3 = medLocal(m,:)
timestep3 = dt_local*mstep!/10.0
temp3 = hLocal(m)-273.0
if (temp3 .ge. 300.0) then
	temp3 = 299.0
end if
!write(*,*) "tem..."
!write(*,*) temp3

! SOLUTES TO STRINGS
write(s_ph,'(F25.10)') solute3(1)
write(s_alk,'(F25.10)') solute3(2)
write(s_water,'(F25.10)') solute3(3)
write(s_co2,'(F25.10)') solute3(4)
write(s_ca,'(F25.10)') solute3(5)
write(s_mg,'(F25.10)') solute3(6)
write(s_na,'(F25.10)') solute3(7)
write(s_k,'(F25.10)') solute3(8)
write(s_fe,'(F25.10)') solute3(9)
write(s_s,'(F25.10)') solute3(10)
write(s_si,'(F25.10)') solute3(11)
write(s_cl,'(F25.10)') solute3(12)
write(s_al,'(F25.10)') solute3(13)
write(s_hco3,'(F25.10)') solute3(14)
write(s_co3,'(F25.10)') solute3(15)

! MEDIUM TO STRINGS
write(s_w,'(F25.10)') medium3(3)!solute3(3)


!// trim(s_reactive) //"*M*(1.55e-8)*exp(-25.5/(.008314*("// trim(s_temp) //"+273.0) ))" // &
solute3(15) = (10.0**(-1.0*solute3(1)))**3

! if (solute3(13) .lt. 1e-15) then
! 	solute3(13) = 1e-8
! end if
!primary3(2) = medium3(4)*primary3(5)*(1.55e-2)*exp(-25.5/(273.0+temp3))*((solute3(15)))
!primary3(2) = medium3(4)*primary3(5)*(1.55e-6)*exp(-25.5/(273.0+temp3))*((solute3(15)**3)/solute3(13))**.33

! PRIMARIES TO STRINGS
write(s_feldspar,'(F25.10)') primary3(1)
write(s_augite,'(F25.10)') primary3(2)
write(s_pigeonite,'(F25.10)') primary3(3)
write(s_magnetite,'(F25.10)') primary3(4)
write(s_glass,'(F25.10)') primary3(5)





! SECONDARIES TO STRINGS
!secondary = 0.0
write(s_stilbite,'(F25.10)') secondary3(1)
write(s_aragonite,'(F25.10)') secondary3(2)
write(s_kaolinite,'(F25.10)') secondary3(3)
write(s_albite,'(F25.10)') secondary3(4)
write(s_saponite,'(F25.10)') secondary3(5)
write(s_celadonite,'(F25.10)') secondary3(6)
write(s_clinoptilolite,'(F25.10)') secondary3(7)
write(s_pyrite,'(F25.10)') secondary3(8)
write(s_mont_na,'(F25.10)') secondary3(9)
write(s_goethite,'(F25.10)') secondary3(10)
write(s_dolomite,'(F25.10)') secondary3(11)
write(s_smectite,'(F25.10)') secondary3(12)
write(s_saponite_k,'(F25.10)') secondary3(13)
write(s_anhydrite,'(F25.10)') secondary3(14)
write(s_siderite,'(F25.10)') secondary3(15)
write(s_calcite,'(F25.10)') secondary3(16)
write(s_quartz,'(F25.10)') secondary3(17)
write(s_kspar,'(F25.10)') secondary3(18)
write(s_saponite_na,'(F25.10)') secondary3(19) !!!!
write(s_nont_na,'(F25.10)') secondary3(20)
write(s_nont_mg,'(F25.10)') secondary3(21)
write(s_nont_k,'(F25.10)') secondary3(22)
write(s_nont_h,'(F25.10)') secondary3(23)
write(s_nont_ca,'(F25.10)') secondary3(24)
write(s_muscovite,'(F25.10)') secondary3(25)
write(s_mesolite,'(F25.10)') secondary3(26)
write(s_hematite,'(F25.10)') secondary3(27)
write(s_mont_ca,'(F25.10)') secondary3(28)
! MINERALS ADDED 10/18/2104
write(s_verm_ca,'(F25.10)') secondary3(29)
write(s_analcime,'(F25.10)') secondary3(30) 
write(s_phillipsite,'(F25.10)') secondary3(31)
write(s_diopside,'(F25.10)') secondary3(32)
write(s_epidote,'(F25.10)') secondary3(33)
write(s_gismondine,'(F25.10)') secondary3(34)
write(s_hedenbergite,'(F25.10)') secondary3(35)
write(s_chalcedony,'(F25.10)') secondary3(36)
write(s_verm_mg,'(F25.10)') secondary3(37)
write(s_ferrihydrite,'(F25.10)') secondary3(38)
write(s_natrolite,'(F25.10)') secondary3(39)
write(s_talc,'(F25.10)') secondary3(40) !!!!!!!!!
write(s_smectite_low,'(F25.10)') secondary3(41)
write(s_prehnite,'(F25.10)') secondary3(42)
write(s_chlorite,'(F25.10)') secondary3(43) !!!!!!!!
write(s_scolecite,'(F25.10)') secondary3(44)
write(s_chamosite7a,'(F25.10)') secondary3(45)
write(s_clinochlore14a,'(F25.10)') secondary3(46)
write(s_clinochlore7a,'(F25.10)') secondary3(47)
!! NEXT ROUND
write(s_saponite_ca,'(F25.10)') secondary3(48)
write(s_verm_na,'(F25.10)') secondary3(49)
write(s_pyrrhotite,'(F25.10)') secondary3(50)
write(s_magnetite,'(F25.10)') secondary3(51)
write(s_lepidocrocite,'(F25.10)') secondary3(52)
write(s_daphnite_7a,'(F25.10)') secondary3(53) !!!
write(s_daphnite_14a,'(F25.10)') secondary3(54) !!!
write(s_verm_k,'(F25.10)') secondary3(55)
write(s_mont_k,'(F25.10)') secondary3(56)
write(s_mont_mg,'(F25.10)') secondary3(57)


! OTHER INFORMATION TO STRINGS
write(s_temp,'(F25.10)') temp3
write(s_timestep,'(F25.10)') timestep3

write(s_reactive,'(F25.10)') medium3(4)

! ----------------------------------%%
! INITIAL AQUEOUS PHASE CONSITUENTS
! ----------------------------------%%

kinetics = " precipitate_only"
!kinetics = " "
write(s_pressure,'(F25.10)') 250.0 - (medium3(7)/5.0)

write(si_hematite,'(F25.10)') 1.0!-(solute3(1)*2.5) + 30.0

inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
!&"    pH " // trim(s_pH) //NEW_LINE('')// &
!&"    pe " // trim(s_pe) //NEW_LINE('')// &
&"    units   mol/kgw" //NEW_LINE('')// &
&"    temp" // trim(s_temp) //NEW_LINE('')// &
!&"    pressure 250.0" //NEW_LINE('')// &
!&"    pressure " // trim(s_pressure) //NEW_LINE('')// &
&"    Ca " // trim(s_ca) //NEW_LINE('')// &
&"    Mg " // trim(s_mg) //NEW_LINE('')// &
&"    Na " // trim(s_na) //NEW_LINE('')// &
&"    K " // trim(s_k) //NEW_LINE('')// &
&"    Fe " // trim(s_fe) //NEW_LINE('')// &
&"    S "// trim(s_s)  //NEW_LINE('')// &
&"    Si " // trim(s_si) //NEW_LINE('')// &
&"    Cl " // trim(s_cl) //NEW_LINE('')// &
&"    Al " // trim(s_al) //NEW_LINE('')// &
&"    C " // trim(s_co2) //NEW_LINE('')// &
&"    Alkalinity " // trim(s_alk) //NEW_LINE('')// &
! // "as as Ca0.5(CO3)0.5" 
&"    -water "// trim(s_w) // " # kg" //NEW_LINE('')// &

! ----------------------------------%%
! HYDROTHERMAL MINERAL CHOICES
! ----------------------------------%%


&"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
&"    Kaolinite 0.0 " // trim(s_kaolinite) // kinetics //NEW_LINE('')// &
&"    Goethite 0.0 " // trim(s_goethite) // kinetics //NEW_LINE('')// &
&"    Celadonite 0.0 " // trim(s_celadonite) // kinetics //NEW_LINE('')// &
&"    Albite 0.0 " // trim(s_albite) // kinetics //NEW_LINE('')// &
&"    Calcite 0.0 " // trim(s_calcite) // kinetics //NEW_LINE('')// & ! .135
&"    Montmor-Na 0.0 " // trim(s_mont_na) // kinetics //NEW_LINE('')// &
&"    Montmor-K 0.0 " // trim(s_mont_k) // kinetics //NEW_LINE('')// &
&"    Montmor-Mg 0.0 " // trim(s_mont_mg) // kinetics //NEW_LINE('')// &
&"    Montmor-Ca 0.0 " // trim(s_mont_ca) // kinetics //NEW_LINE('')// &
&"    Saponite-Mg 0.0 " // trim(s_saponite) // kinetics //NEW_LINE('')// &
&"    Stilbite 0.0 " // trim(s_stilbite) // kinetics //NEW_LINE('')// &
&"    Clinoptilolite-Ca 0.0 " // trim(s_clinoptilolite) // kinetics //NEW_LINE('')// & ! zeolite
&"    Pyrite 0.0 " // trim(s_pyrite) // kinetics //NEW_LINE('')// &
&"    Quartz 0.0 " // trim(s_quartz) // kinetics //NEW_LINE('')// &
&"    K-Feldspar 0.0 " // trim(s_kspar) // kinetics //NEW_LINE('')// &

! NEW MINS
 &"    Saponite-Na 0.0 " // trim(s_saponite_na) // kinetics //NEW_LINE('')// &
 &"    Nontronite-Na 0.0 " // trim(s_nont_na) // kinetics //NEW_LINE('')// &
 &"    Nontronite-Mg 0.0 " // trim(s_nont_mg) // kinetics //NEW_LINE('')// &
 &"    Nontronite-K 0.0 " // trim(s_nont_k) // kinetics //NEW_LINE('')// &
 &"    Nontronite-H 0.0 " // trim(s_nont_h) // kinetics //NEW_LINE('')// &
 &"    Nontronite-Ca 0.0 " // trim(s_nont_ca) // kinetics //NEW_LINE('')// &
&"    Muscovite 0.0 " // trim(s_muscovite) // kinetics //NEW_LINE('')// &
 &"    Mesolite 0.0 " // trim(s_mesolite) // kinetics //NEW_LINE('')// & !!!!!!!!
 &"    Anhydrite 0.0 " // trim(s_anhydrite) // kinetics //NEW_LINE('')// & ! formerly magnesite
 &"    Smectite-high-Fe-Mg 0.0 " // trim(s_smectite) // kinetics //NEW_LINE('')// &
 &"    Saponite-K 0.0 " // trim(s_saponite_k) // kinetics //NEW_LINE('')// &
  &"    Vermiculite-Na 0.0 " // trim(s_verm_na) // kinetics //NEW_LINE('')// &
  &"    Hematite 0.0 " // trim(s_hematite) // kinetics //NEW_LINE('')// &
! &"    Hematite " // trim(si_hematite) // trim(s_hematite) // kinetics //NEW_LINE('')// &

!! MINERALS ADDED 10/18/2014
  &"    Vermiculite-Ca 0.0 " // trim(s_verm_ca) // kinetics //NEW_LINE('')// &
 &"    Analcime 0.0 " // trim(s_analcime) // kinetics //NEW_LINE('')// &
 &"    Phillipsite 0.0 " // trim(s_phillipsite) // kinetics //NEW_LINE('')// &
   &"    Diopside 0.0 " // trim(s_diopside) // kinetics //NEW_LINE('')// & ! pyroxene
   &"    Epidote  0.0 " // trim(s_epidote) // kinetics //NEW_LINE('')// &
  &"    Gismondine 0.0 " // trim(s_gismondine) // kinetics //NEW_LINE('')// & ! zeolite
  &"    Hedenbergite 0.0 " // trim(s_hedenbergite) // kinetics //NEW_LINE('')// & ! pyroxene
  &"    Chalcedony 0.0 " // trim(s_chalcedony) // kinetics //NEW_LINE('')// &
  &"    Vermiculite-Mg 0.0 " // trim(s_verm_mg) // kinetics //NEW_LINE('')// &
 &"    Ferrihydrite 0.0 " // trim(s_ferrihydrite) // kinetics //NEW_LINE('')// &
  &"    Natrolite 0.0 " // trim(s_natrolite) // kinetics //NEW_LINE('')// & ! zeolite
&"    Talc 0.0 " // trim(s_talc) // kinetics //NEW_LINE('')// &
 &"    Smectite-low-Fe-Mg 0.0 " // trim(s_smectite_low) // kinetics //NEW_LINE('')// &
 &"    Prehnite 0.0 " // trim(s_prehnite) // kinetics //NEW_LINE('')// &
  &"    Chlorite(14A) 0.0 " // trim(s_chlorite) // kinetics //NEW_LINE('')// &
  &"    Scolecite 0.0 " // trim(s_scolecite) // kinetics //NEW_LINE('')// &
  &"    Chamosite-7A 0.0 " // trim(s_chamosite7a) // kinetics //NEW_LINE('')// &
  &"    Clinochlore-14A 0.0 " // trim(s_clinochlore14a) // kinetics //NEW_LINE('')// &
  &"    Clinochlore-7A 0.0 " // trim(s_clinochlore7a) // kinetics //NEW_LINE('')// &
 !! NEXT ROUND

 &"   Saponite-Ca 0.0 " // trim(s_saponite_ca) // kinetics //NEW_LINE('')// &
 &"   Pyrrhotite 0.0 " // trim(s_pyrrhotite) // kinetics //NEW_LINE('')// &
 &"   Magnetite 0.0 " // trim(s_magnetite) // kinetics //NEW_LINE('')// &
  &"   Daphnite-7a 0.0 " // trim(s_daphnite_7a) // kinetics //NEW_LINE('')// &
  &"   Daphnite-14a 0.0 " // trim(s_daphnite_14a) // kinetics //NEW_LINE('')// &
 &"   Vermiculite-K 0.0 " // trim(s_verm_k) // kinetics //NEW_LINE('')// &
 &"   Aragonite 0.0 " // trim(s_aragonite) // kinetics //NEW_LINE('')// &
! &" -force_equality"  //NEW_LINE('')// &
 &"   Lepidocrocite 0.0 " // trim(s_lepidocrocite) // kinetics //NEW_LINE('')// &



! ---------------------------------%%
! KINETIC DISSOLUTION RATE LAWS
! ----------------------------------%%
!	  &"      10 rate0 =" // s_reactive //"*M*(7.7e-11)*exp(-25.0/(.008314*TK))" //NEW_LINE('')// &
	
	
&"RATES" //NEW_LINE('')// & 

&"BGlass" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!CALC_VALUE('R(s_sp)')
! most recent
!    &"    10 rate0=M*46.5*(1.53e-5)*" // trim(s_reactive) //"*0.000001*(1e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
!    &"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &
  
!   	  &"      10 rate0 =" // s_reactive //"*M*(7.7e-11)*exp(-25.0/(.008314*TK))" //NEW_LINE('')// &

! APRIL 2016 BEST TEMP ONLY LAW FOR NOW
! &"      10 rate0 =0.01*" // trim(s_glass) //"*(7.7e-11)*exp(-25.0/(.008314*TK))" //NEW_LINE('')// &

! SAVE 04/27/16
&"    10 rate0=M*46.5*(1.53e-5)*0.1*0.001*(1e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
&"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &
! &"    10 rate0=M*46.5*(1.53e-5)*0.01*1.0*(1.0e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
! &"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &
 
 
 	  !&"      10 rate0 =TK*7.7e-17" //NEW_LINE('')// &
  
! &"    10 rate0=" // s_reactive //"*M*.000000001*200.0*23000.0*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
! &"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &

! &"    10 rate0=" // trim(s_reactive) //"*M*(1.55e-8)*exp(-25.5/(.008314*("// trim(s_temp) //"+273.0) ))" // &
! &"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &


! &"    10 rate0=" // trim(s_augite) // "*((ACT('H+')^3)/(ACT('Al+3')))^.333" //NEW_LINE('')// &
!&"    10 rate0=" // trim(s_augite) //NEW_LINE('')// &

! this!
! &"      10 rate0 =" // trim(s_reactive) //"*M*(7.7e-12)*exp(-25.0/(.008314*("// trim(s_temp) //"+273.0) ))" //NEW_LINE('')// &
!&"      10 rate0 =" // trim(s_reactive) //"*M*(7.7e-12)*exp(-25.0/(.008314*("// trim(s_temp) //"+273.0) ))" //NEW_LINE('')// &
  
!  &"      10 rate0 = (0.25)*M*(10.0^(-10))*exp(-25.0/(.008314*TK))" //NEW_LINE('')// &
! &"      10 rate0 = (.001) * (10^(-2940.0/TK)) *exp(-25.0/(.008314*TK))" //NEW_LINE('')// &
&"    20 save rate0 * time" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

! &"Plagioclase" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! !(1-SR('Plagioclase'))*
! &"    10 rate = (1-SR('Plagioclase'))*M*270.0*(1.53e-5)*1.0*(((1.58e-9)"//&
! &"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
! &"(4.78e-15)*exp(-59.0/(.008314*TK))*(ACT('H+'))^-0.57))"//NEW_LINE('')//&
! &"    20 save rate * time"//NEW_LINE('')//&
! &"-end" //NEW_LINE('')// &
!
! &"Augite" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! &"    10 rate0 = (1-SR('Augite'))*M*230.0*(1.53e-5)*1.0*(((1.58e-7)" // &
! &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //NEW_LINE('')// &
! &"    20 save rate0 * time" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
! &"Pigeonite" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! &"    10 rate0 = (1-SR('Pigeonite'))*M*236.0*(1.53e-5)*1.0*(((1.58e-7)" // &
! &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//NEW_LINE('')// &
! &"    20 save rate0 * time" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
! &"Magnetite" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! &"    10 rate0 = (1-SR('Magnetite'))*M*231.0*(1.53e-5)*1.0*(((2.57e-9)" // &
! &"*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279)+(1.66e-11)*exp(-18.6/(.008314*TK))))" //NEW_LINE('')// &
! &"    20 save rate0 * time" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &

! ----------------------------------%%
! PRIMARY (KINETIC) CONSTITUENTS
! ----------------------------------%%

!&"Use solution 1" //NEW_LINE('')// &
!&"Use equilibrium_phases 1" //NEW_LINE('')// &
&"KINETICS 1" //NEW_LINE('')// &
! &"Plagioclase" //NEW_LINE('')// &
! &"-m0 " // trim(s_feldspar) //NEW_LINE('')// &
! &"Augite" //NEW_LINE('')// &
! &"-m0 " // trim(s_augite) //NEW_LINE('')// &
! &"Pigeonite" //NEW_LINE('')// &
! &"-m0 " // trim(s_pigeonite) //NEW_LINE('')// &
! &"Magnetite" //NEW_LINE('')// &
! &"-m0 " // trim(s_magnetite) //NEW_LINE('')// &
&"BGlass" //NEW_LINE('')// &

! ! ! hole 896a
! &"-f CaO .220 SiO2 .856 Al2O3 .15 " //&
! & "FeO .15 MgO .1826 K2O .0042 " //&
! & "Na2O .0333" //NEW_LINE('')// &
! &"-m0 " // trim(s_glass) //NEW_LINE('')// &
! &"    -step " // trim(s_timestep) // " in 1" //NEW_LINE('')// &


! ! seyfried JDF
&"-f CaO .1997 SiO2 .847 Al2O3 .138 " //&
! & "Fe2O3 .149 FeO .0075 MgO .1744 K2O .002 " //&
& "Fe2O3 .149 MgO .1744 K2O .002 " //&
& "Na2O .043" //NEW_LINE('')// &
&"-m0 " // trim(s_glass) //NEW_LINE('')// &
&"    -step " // trim(s_timestep) // " in 1" //NEW_LINE('')// &

&"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &

! ----------------------------------%%
! CALCULATE POROSITY AND STUFF
! ----------------------------------%%

&"CALCULATE_VALUES" //NEW_LINE('')// &
!

&"R(sum)" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
! &"10 sum = (" //&
! &"(EQUI('stilbite')*480.19/2.15) + (EQUI('aragonite')*100.19/2.93)" // &
! &"+ (EQUI('kaolinite')*258.16/2.63) + (EQUI('albite')*263.02/2.62)" // &
! &"+ (EQUI('Saponite-Mg')*480.19/2.3) + (EQUI('celadonite')*429.02/3.0)" // &
! &"+ (EQUI('Clinoptilolite-Ca')*2742.13/2.15) + (EQUI('pyrite')*119.98/5.02)" // &
! &"+ (EQUI('montmor-na')*549.07/2.01) + (EQUI('goethite')*88.85/4.27)" // &
! &"+ (EQUI('dolomite')*180.4/2.84) + (EQUI('Smectite-high-Fe-Mg')*540.46/2.01)" // &
! &"+ (EQUI('saponite-k')*480.19/2.3) + (EQUI('anhydrite')*136.14/2.97)" // &
! &"+ (EQUI('siderite')*115.86/3.96) + (EQUI('calcite')*100.19/2.71)" // &
! &"+ (EQUI('quartz')*60.08/2.65) + (EQUI('k-feldspar')*278.33/2.56)" // &
! &"+ (EQUI('saponite-na')*480.19/2.3) + (EQUI('nontronite-na')*495.9/2.3)" // &
! &"+ (EQUI('nontronite-mg')*495.9/2.3) + (EQUI('nontronite-k')*495.9/2.3)" // &
! &"+ (EQUI('nontronite-h')*495.9/2.3) + (EQUI('nontronite-ca')*495.9/2.3)" // &
! &"+ (EQUI('muscovite')*398.71/2.81) + (EQUI('mesolite')*1164.9/2.29)" // &
! &"+ (EQUI('hematite')*159.69/5.3) + (EQUI('montmor-ca')*549.07/2.01)" // &
! &"+ (EQUI('vermiculite-ca')*504.19/2.5) + (EQUI('analcime')*220.15/2.27)" // &
! &"+ (EQUI('phillipsite')*704.93/2.2) + (EQUI('diopside')*216.55/3.3)" // &
! &"+ (EQUI('epidote')*519.3/3.41) + (EQUI('gismondine')*718.55/2.26)" // &
! &"+ (EQUI('hedenbergite')*248.08/3.56) + (EQUI('chalcedony')*60.08/2.65)" // &
! &"+ (EQUI('vermiculite-mg')*504.19/2.5) + (EQUI('ferrihydrite')*169.7/3.8)" // &
! &"+ (EQUI('natrolite')*380.22/2.23) + (EQUI('talc')*379.27/2.75)" // &
! &"+ (EQUI('smectite-low-fe-mg')*540.46/2.01) + (EQUI('prehnite')*395.38/2.87)" // &
! &"+ (EQUI('chlorite')*67.4/2.468) + (EQUI('scolecite')*392.34/2.27)" // &
! &"+ (EQUI('chamosite-7a')*664.18/3.0) + (EQUI('clinochlore-14a')*595.22/3.0)" // &
! &"+ (EQUI('clinochlore-7a')*595.22/3.0) + (EQUI('saponite-ca')*480.19/2.3)" // &
! &"+ (EQUI('vermiculite-na')*504.19/2.5) + (EQUI('pyrrhotite')*85.12/4.62)" // &
! &"+ (EQUI('magnetite')*231.53/5.15) + (EQUI('lepidocrocite')*88.85/4.08)" // &
! &"+ (EQUI('daphnite-7a')*664.18/3.2) + (EQUI('daphnite-14a')*664.18/3.2)" // &
! &"+ (EQUI('vermiculite-k')*504.19/2.5) + (EQUI('montmor-k')*549.07/2.01)" // &
! &"+ (EQUI('montmor-mg')*549.07/2.01) + (KIN('Bglass')*96.8/2.92) )" // &
&"10 sum = 5.0" //&
&"" //NEW_LINE('')// &
&"100 SAVE sum" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

!
&"R(phi)" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!&"10 phi = 1.0-(CALC_VALUE('R(sum)')/(CALC_VALUE('R(sum)')+(TOT('water')*1000.0)))" //&
&"10 phi = .351" //&
&"" //NEW_LINE('')// &
&"100 SAVE phi" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &
!
&"R(water_volume)" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!&"10 water_volume = SOLN_VOL / RHO" //NEW_LINE('')// &
&"10 water_volume = 0.3" //NEW_LINE('')// &
&"100 SAVE water_volume" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &
!
!
&"R(rho_s)" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
! &"10 rho_s = EQUI('Stilbite')*2.15 + EQUI('SiO2(am)')*2.62" //&
! &"+ EQUI('Kaolinite')*2.6 + EQUI('Albite')*2.62" // &
! &"+ EQUI('Saponite-Mg')*2.4 + EQUI('Celadonite')*3.0" // &
! &"+ EQUI('Clinoptilolite-Ca')*2.62 + EQUI('Pyrite')*4.84" // &
! &"+ EQUI('Montmor-Na')*5.3 + EQUI('Goethite')*3.8" // &
! &"+ EQUI('Dolomite')*2.84 + EQUI('Smectite-high-Fe-Mg')*2.7" // &
! &"+ EQUI('Dawsonite')*2.42 + EQUI('Anhydrite')*3.0" // &
! &"+ EQUI('Siderite')*3.96 + EQUI('Calcite')*2.71" // &
! &"+ EQUI('Quartz')*2.62 + EQUI('k-Feldspar')*2.56" // &
! &"+ KIN('Plagioclase')*2.68 + KIN('Augite')*3.4" // &
! &"+ KIN('Pigeonite')*3.38 + KIN('Magnetite')*5.15" // &
! &"+ KIN('BGlass')*2.92" //NEW_LINE('')// &
! &"20 rho_s = rho_s/ (EQUI('Stilbite') + EQUI('SiO2(am)')" //&
! &"+ EQUI('Kaolinite') + EQUI('Albite')" // &
! &"+ EQUI('Saponite-Mg') + EQUI('Celadonite')" // &
! &"+ EQUI('Clinoptilolite-Ca') + EQUI('Pyrite')" // &
! &"+ EQUI('Montmor-Na') + EQUI('Goethite')" // &
! &"+ EQUI('Dolomite') + EQUI('Smectite-high-Fe-Mg')" // &
! &"+ EQUI('Dawsonite') + EQUI('Anhydrite')" // &
! &"+ EQUI('Siderite') + EQUI('Calcite')" // &
! &"+ EQUI('Quartz') + EQUI('K-Feldspar')" // &
! &"+ KIN('Plagioclase') + KIN('Augite')" // &
! &"+ KIN('Pigeonite') + KIN('Magnetite')" // &
! &"+ KIN('BGlass'))" //NEW_LINE('')// &
! &"30 rho_s = rho_s * 1000000.0" //NEW_LINE('')// &

! ! THIS IS THE GOOD ONE
! &"10 rho_s = 1000000.0* ( (EQUI('stilbite')*480.19/2.15) + (EQUI('aragonite')*100.19/2.93)" // &
! &"+ (EQUI('kaolinite')*258.16/2.63) + (EQUI('albite')*263.02/2.62)" // &
! &"+ (EQUI('Saponite-Mg')*480.19/2.3) + (EQUI('celadonite')*429.02/3.0)" // &
! &"+ (EQUI('Clinoptilolite-Ca')*2742.13/2.15) + (EQUI('pyrite')*119.98/5.02)" // &
! &"+ (EQUI('montmor-na')*549.07/2.01) + (EQUI('goethite')*88.85/4.27)" // &
! &"+ (EQUI('dolomite')*180.4/2.84) + (EQUI('Smectite-high-Fe-Mg')*540.46/2.01)" // &
! &"+ (EQUI('saponite-k')*480.19/2.3) + (EQUI('anhydrite')*136.14/2.97)" // &
! &"+ (EQUI('siderite')*115.86/3.96) + (EQUI('calcite')*100.19/2.71)" // &
! &"+ (EQUI('quartz')*60.08/2.65) + (EQUI('k-feldspar')*278.33/2.56)" // &
! &"+ (EQUI('saponite-na')*480.19/2.3) + (EQUI('nontronite-na')*495.9/2.3)" // &
! &"+ (EQUI('nontronite-mg')*495.9/2.3) + (EQUI('nontronite-k')*495.9/2.3)" // &
! &"+ (EQUI('nontronite-h')*495.9/2.3) + (EQUI('nontronite-ca')*495.9/2.3)" // &
! &"+ (EQUI('muscovite')*398.71/2.81) + (EQUI('mesolite')*1164.9/2.29)" // &
! &"+ (EQUI('hematite')*159.69/5.3) + (EQUI('montmor-ca')*549.07/2.01)" // &
! &"+ (EQUI('vermiculite-ca')*504.19/2.5) + (EQUI('analcime')*220.15/2.27)" // &
! &"+ (EQUI('phillipsite')*704.93/2.2) + (EQUI('diopside')*216.55/3.3)" // &
! &"+ (EQUI('epidote')*519.3/3.41) + (EQUI('gismondine')*718.55/2.26)" // &
! &"+ (EQUI('hedenbergite')*248.08/3.56) + (EQUI('chalcedony')*60.08/2.65)" // &
! &"+ (EQUI('vermiculite-mg')*504.19/2.5) + (EQUI('ferrihydrite')*169.7/3.8)" // &
! &"+ (EQUI('natrolite')*380.22/2.23) + (EQUI('talc')*379.27/2.75)" // &
! &"+ (EQUI('smectite-low-fe-mg')*540.46/2.01) + (EQUI('prehnite')*395.38/2.87)" // &
! &"+ (EQUI('chlorite')*67.4/2.468) + (EQUI('scolecite')*392.34/2.27)" // &
! &"+ (EQUI('chamosite-7a')*664.18/3.0) + (EQUI('clinochlore-14a')*595.22/3.0)" // &
! &"+ (EQUI('clinochlore-7a')*595.22/3.0) + (EQUI('saponite-ca')*480.19/2.3)" // &
! &"+ (EQUI('vermiculite-na')*504.19/2.5) + (EQUI('pyrrhotite')*85.12/4.62)" // &
! &"+ (EQUI('magnetite')*231.53/5.15) + (EQUI('lepidocrocite')*88.85/4.08)" // &
! &"+ (EQUI('daphnite-7a')*664.18/3.2) + (EQUI('daphnite-14a')*664.18/3.2)" // &
! &"+ (EQUI('vermiculite-k')*504.19/2.5) + (EQUI('montmor-k')*549.07/2.01)" // &
! &"+ (EQUI('montmor-mg')*549.07/2.01) + (KIN('Bglass')*96.8/2.92) )/" // &
! &"( (EQUI('stilbite')) + (EQUI('aragonite'))" // &
! &"+ (EQUI('kaolinite')) + (EQUI('albite'))" // &
! &"+ (EQUI('Saponite-Mg')) + (EQUI('celadonite'))" // &
! &"+ (EQUI('Clinoptilolite-Ca')) + (EQUI('pyrite'))" // &
! &"+ (EQUI('montmor-na')) + (EQUI('goethite'))" // &
! &"+ (EQUI('dolomite')) + (EQUI('Smectite-high-Fe-Mg'))" // &
! &"+ (EQUI('saponite-k')) + (EQUI('anhydrite'))" // &
! &"+ (EQUI('siderite')) + (EQUI('calcite'))" // &
! &"+ (EQUI('quartz')) + (EQUI('k-feldspar'))" // &
! &"+ (EQUI('saponite-na')) + (EQUI('nontronite-na'))" // &
! &"+ (EQUI('nontronite-mg')) + (EQUI('nontronite-k'))" // &
! &"+ (EQUI('nontronite-h')) + (EQUI('nontronite-ca'))" // &
! &"+ (EQUI('muscovite')) + (EQUI('mesolite'))" // &
! &"+ (EQUI('hematite')) + (EQUI('montmor-ca'))" // &
! &"+ (EQUI('vermiculite-ca')) + (EQUI('analcime'))" // &
! &"+ (EQUI('phillipsite')) + (EQUI('diopside'))" // &
! &"+ (EQUI('epidote')) + (EQUI('gismondine'))" // &
! &"+ (EQUI('hedenbergite')) + (EQUI('chalcedony'))" // &
! &"+ (EQUI('vermiculite-mg')) + (EQUI('ferrihydrite'))" // &
! &"+ (EQUI('natrolite')) + (EQUI('talc'))" // &
! &"+ (EQUI('smectite-low-fe-mg')) + (EQUI('prehnite'))" // &
! &"+ (EQUI('chlorite')) + (EQUI('scolecite'))" // &
! &"+ (EQUI('chamosite-7a')) + (EQUI('clinochlore-14a'))" // &
! &"+ (EQUI('clinochlore-7a')) + (EQUI('saponite-ca'))" // &
! &"+ (EQUI('vermiculite-na')) + (EQUI('pyrrhotite'))" // &
! &"+ (EQUI('magnetite')) + (EQUI('lepidocrocite'))" // &
! &"+ (EQUI('daphnite-7a')) + (EQUI('daphnite-14a'))" // &
! &"+ (EQUI('vermiculite-k')) + (EQUI('montmor-k'))" // &
! &"+ (EQUI('montmor-mg')) + (KIN('BGlass')) )" // &


! don't uncomment this!
!&"30 rho_s = rho_s * 1000000.0" //NEW_LINE('')// &
&"10 rho_s = 2.93e6" //NEW_LINE('')// &
&"" //NEW_LINE('')// &
&"100 SAVE rho_s" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &
!
!
&"R(s_sp)" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!&"10 s_sp = (CALC_VALUE('R(phi)')/(1.0-CALC_VALUE('R(phi)')))*400.0/CALC_VALUE('R(rho_s)')" //&
&"10 s_sp = 1.53e-5" //&
&"" //NEW_LINE('')// &
&"100 SAVE s_sp" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &


! ----------------------------------%%
! DEFINE THE KIND OF OUTPUT
! ----------------------------------%%

! &"DUMP" //NEW_LINE('')// &
! &"    -solution 1" //NEW_LINE('')// &
! &"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -k bglass" //NEW_LINE('')// &
  &"    -ph" //NEW_LINE('')// &
  &"    -pe false" //NEW_LINE('')// &
  &"    -totals C" //NEW_LINE('')// &
  &"    -totals Ca Mg Na K Fe S Si Cl Al " //NEW_LINE('')// &
  &"    -molalities HCO3-" //NEW_LINE('')// &
  &"    -water true" //NEW_LINE('')// &
  &"    -alkalinity" //NEW_LINE('')// &
&"    -p stilbite aragonite kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
&"    -p pyrite montmor-na goethite dolomite Smectite-high-Fe-Mg saponite-k anhydrite" //NEW_LINE('')// &
&"    -p siderite calcite quartz k-feldspar saponite-na nontronite-na nontronite-mg" //NEW_LINE('')// &
&"    -p nontronite-k nontronite-h nontronite-ca muscovite mesolite hematite montmor-ca" //NEW_LINE('')// &
&"    -p vermiculate-ca analcime phillipsite diopside epidote gismondine hedenbergite" //NEW_LINE('')// &
&"    -p chalcedony vermiculite-mg ferrihydrite natrolite talc Smectite-low-Fe-Mg prehnite" //NEW_LINE('')// &
&"    -p chlorite scolecite chamosite-7a Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// &
&"    -p vermiculite-na pyrrhotite magnetite Lepidocrocite Daphnite-7A Daphnite-14A" //NEW_LINE('')// &
&"    -p vermiculite-k montmor-k montmor-mg" //NEW_LINE('')// &

&"    -s stilbite aragonite kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
&"    -s pyrite montmor-na goethite dolomite Smectite-high-Fe-Mg saponite-k anhydrite" //NEW_LINE('')// &
&"    -s siderite calcite quartz k-feldspar saponite-na nontronite-na nontronite-mg" //NEW_LINE('')// &
&"    -s nontronite-k nontronite-h nontronite-ca muscovite mesolite hematite montmor-ca" //NEW_LINE('')// &
&"    -s vermiculate-ca analcime phillipsite diopside epidote gismondine hedenbergite" //NEW_LINE('')// &
&"    -s chalcedony vermiculite-mg ferrihydrite natrolite talc Smectite-low-Fe-Mg prehnite" //NEW_LINE('')// &
&"    -s chlorite scolecite chamosite-7a Clinochlore-14A Clinochlore-7A saponite-ca" //NEW_LINE('')// &
&"    -s vermiculite-na pyrrhotite magnetite Lepidocrocite Daphnite-7A Daphnite-14A" //NEW_LINE('')// &
&"    -s vermiculite-k montmor-k montmor-mg" //NEW_LINE('')// &
  &"    -calculate_values R(phi) R(s_sp)" //NEW_LINE('')// &
  &"    -time" //NEW_LINE('')// &
&"END"


! INITIALIZE STUFF
id = CreateIPhreeqc()

! write(*,*) "ID:" , id


! IF (SetErrorFileOn(id,.TRUE.).NE.IPQ_OK) THEN
! 	CALL OutputErrorString(id)
! 	!STOP
! END IF


IF (id.LT.0) THEN
	write(*,*) "weird stop?"
	STOP
END IF


IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
	CALL OutputErrorString(id)
	write(*,*) "primary"
	write(*,*) primary3
	write(*,*) "secondary"
	write(*,*) secondary3
	write(*,*) "solute"
	write(*,*) solute3
	write(*,*) "medium"
	write(*,*) medium3
	write(*,*) "temp"
	write(*,*) temp3
	!STOP
END IF

! IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
! 	CALL OutputErrorString(id)
! 	!STOP
! END IF

! write(*,*) "we here"
!IF (LoadDatabase(id, 'l5.dat').NE.0) THEN
IF (LoadDatabaseString(id, L5).NE.0) THEN
	CALL OutputErrorString(id)
	write(*,*) "primary"
	write(*,*) primary3
	write(*,*) "secondary"
	write(*,*) secondary3
	write(*,*) "solute"
	write(*,*) solute3
	write(*,*) "medium"
	write(*,*) medium3
	write(*,*) "temp"
	write(*,*) temp3
	!STOP
END IF

! RUN INPUT
IF (RunString(id, inputz0).NE.0) THEN
	write(*,*) "issue is:" , RunString(id, inputz0)
	CALL OutputErrorString(id)
	write(*,*) "primary"
	write(*,*) primary3
	write(*,*) "secondary"
	write(*,*) secondary3
	write(*,*) "solute"
	write(*,*) solute3
	write(*,*) "medium"
	write(*,*) medium3
	write(*,*) "temp"
	write(*,*) temp3
	IF (RunString(id, inputz0).NE.0) THEN
		write(*,*) "another chance 2"
		CALL OutputErrorString(id)
	END IF
	!STOP
END IF



! ! !if errror print dump. look for other phases.
! if ((temp .GT. 35.0) .and. (temp .LT. 40.0)) then
! 	DO i=1,GetOutputStringLineCount(id)
! 		call GetOutputStringLine(id, i, line)
! 		write(*,*) trim(line)
! 	END DO
! end if
  
! NO MORE DUMPING!!!!!!!!!
! ! PRINT DUMP/OUTPUT
! DO i=1,GetOutputStringLineCount(id)
! 	call GetOutputStringLine(id, i, line)
! 	write(*,*) trim(line)
! END DO


! NOW KINDA USELESS PRINT STATEMENTS FOR WRITING TO FILES
!WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
!WRITE(*,*) "NUMBER OF LINES:"
!WRITE(*,*) GetSelectedOutputStringLineCount(id)
  
  
! OPEN FILE (don't need no file) (USEFUL)
!OPEN(UNIT=12, FILE="testMat.txt", ACTION="write", STATUS="replace") 
  
! WRITE AWAY
DO i=1,GetSelectedOutputStringLineCount(id)
	call GetSelectedOutputStringLine(id, i, line)
	! HEADER BITS YOU MAY WANT
	if (i .eq. 1) then
 	   write(12,*) trim(line)
! 	   !if ((medium3(6) .gt. 24000.0) .and. (medium3(7) .gt. -100.0)) then
! 	   write(*,*) trim(line) ! PRINT LABELS FOR EVERY FIELD (USEFUL)
! 	   !end if
	end if
	
	! MEAT
	if (i .gt. 1) then
		read(line,*) outmat(i,:)
		!!!!write(12,*) outmat(i,:) ! this writes to file, which i don't need (USEFUL)
! 		if ((medium3(6) .gt. 23000.0) .and. (medium3(7) .gt. -200.0)) then
! 		write(*,*) i
! 		write(*,*) trim(line) ! PRINT EVERY GOD DAMN LINE
! 		write(*,*) ""
! ! 		! write(*,*) solute3
! ! ! 		write(*,*) ""
! ! 		write(*,*) ""
! 		end if
	end if
END DO
  
! OUTPUT TO THE MAIN MASSACR METHOD
alt0(1,:) = outmat(3,:)
alt0(1,130:186) = outmat(2,130:186)

if (GetSelectedOutputStringLineCount(id) .ne. 3) then
	alt0(1,:) = 0.0
	write(*,*) "not = 3"
end if

! IF (RunString(id, inputz0).NE.0) THEN
! 	write(*,*) "SECOND RUNSTRING ERROR"
! 	alter(1,:) = 0.0
! END IF

! IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
! 	alter(1,:) = 0.0
! END IF
!
! IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
! 	alter(1,:) = 0.0
! END IF
!
! IF (LoadDatabase(id, "llnl.dat").NE.0) THEN
! 	alter(1,:) = 0.0
! END IF



! DESTROY INSTANCE
IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
	CALL OutputErrorString(id)
	write(*,*) "cannot be destroyed?"
	STOP
END IF

				solLocal(m,:) = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
				alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
				alt0(1,13), alt0(1,14), alt0(1,15), 0.0D+00/)
				
				!write(*,*) solLocal(m,13)

				secLocal(m,:) = (/ alt0(1,16), alt0(1,18), alt0(1,20), alt0(1,22), alt0(1,24), alt0(1,26), alt0(1,28), &
				alt0(1,30), alt0(1,32), alt0(1,34), alt0(1,36), alt0(1,38), alt0(1,40), alt0(1,42), &
				alt0(1,44), alt0(1,46), alt0(1,48), alt0(1,50), alt0(1,52), alt0(1,54), alt0(1,56), &
				alt0(1,58), alt0(1,60), alt0(1,62), alt0(1,64), alt0(1,66), alt0(1,68), alt0(1,70), &
				alt0(1,72), alt0(1,74), alt0(1,76), alt0(1,78), alt0(1,80), alt0(1,82), alt0(1,84), &
				alt0(1,86), alt0(1,88), alt0(1,90), alt0(1,92), alt0(1,94), alt0(1,96), alt0(1,98), &
				alt0(1,100), alt0(1,102), alt0(1,104), alt0(1,106), alt0(1,108), alt0(1,110), alt0(1,112), &
				alt0(1,114), alt0(1,116), alt0(1,118), alt0(1,120), alt0(1,122), alt0(1,124), &
				alt0(1,126), alt0(1,128), &
				
				alt0(1,130), alt0(1,131), alt0(1,132), alt0(1,133), alt0(1,134), alt0(1,135), alt0(1,136), &
				alt0(1,137), alt0(1,138), alt0(1,139), alt0(1,140), &
				alt0(1,141), alt0(1,142), alt0(1,143), alt0(1,144), alt0(1,145), alt0(1,146), &
				alt0(1,147), alt0(1,148), alt0(1,149), alt0(1,150), alt0(1,151), alt0(1,152), &
				alt0(1,153), alt0(1,154), alt0(1,155), alt0(1,156), alt0(1,157), alt0(1,158), &
				alt0(1,159), alt0(1,160), alt0(1,161), alt0(1,162), alt0(1,163), alt0(1,164), &
				alt0(1,165), alt0(1,166), alt0(1,167), alt0(1,168), alt0(1,169), alt0(1,170), &
				alt0(1,171), alt0(1,172), alt0(1,173), alt0(1,174), alt0(1,175), alt0(1,176), &
				alt0(1,177), alt0(1,178), alt0(1,179), alt0(1,180), alt0(1,181), alt0(1,182), &
				alt0(1,183), alt0(1,184), alt0(1,185), alt0(1,186)/)

				priLocal(m,:) = (/ 0.0*alt0(1,187), 0.0*alt0(1,187), 0.0*alt0(1,187), 0.0*alt0(1,187), alt0(1,187)/)
				
				!write(*,*) priLocal(m,5)
				
				medLocal(m,1:4) = (/ alt0(1,189), alt0(1,190), alt0(1,4), alt0(1,190)/)

			if (alt0(1,2) .lt. 1.0) then
				medLocal(m,5) = 0.0
				solLocal(m,:) = (/ solute3(1), solute3(2), solute3(3), solute3(4), solute3(5), &
				solute3(6), solute3(7), solute3(8), solute3(9), solute3(10), solute3(11), &
				solute3(12), solute3(13), solute3(14), 0.0D+00/)
			end if

			
					
			end if ! end if-cell-is-on loop
			! TRANSPORT ONLY
			!medLocal(m,5) = 1.0
			! TRANSPORT ONLY
		end do

		!--------------SLAVE PROCESSOR SENDS ALTERED MESSAGE TO MASTER PROCESSOR
	
		! send primary array chunk back to root process
		do ii = 1,g_pri
			call MPI_SEND( priLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send secondary array chunk back to root process
		do ii = 1,g_sec/2
			call MPI_SEND( secLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send solute array chunk back to root process
		do ii = 1,g_sol
			call MPI_SEND( solLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send medium array chunk back to root process
		do ii = 1,g_med
			call MPI_SEND( medLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		write(*,*) "SLAVE PROCESSOR IS DONE WITH WORK"
	
	! done with looping through coarse timesteps
	end do




! end loop through processors	
end if 




! close up shop
call MPI_FINALIZE ( ierr )


END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: computes the 2D temperature profile for the current timestep
!
! INPUTS: h(xn,yn) : temperature profile of previous timestep
!         psi(xn,yn) : 2D streamfunction array
!         rho_in(xn,yn) : 2D density array
!         flux(xn,2) : top and bottom heat flux boundary conditions
!
! RETURNS: h_next(xn,yn) : temperature profile of current timestep
!
! ----------------------------------------------------------------------------------%%

function h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)
	
use globals
use initialize
implicit none

! interface
!
! 	function velocities(psi)
! 		use globals
! 		use initialize
! 		implicit none
! 		real(8) :: velocities(xn,2*yn), psi(xn,yn)
! 		real(8) :: u0(xn,yn), v0(xn,yn)
! 	end function velocities
!
! end interface

! declare errthing

! integers
integer :: i, j, n, ii, m=3
! inputs 
real(8) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2), phi_in(xn,yn)
! velocity stuff
real(8) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
real(8) ::  velocities0(xn,2*yn)
! matrix stuff
real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
! real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
! real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(8) :: h0(xn,yn), hMid(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
real(8) :: kMatLong((xn-2)*(yn-2))
real(8) :: mn(xn,yn)
real(8) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
real(8) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
real(8) :: frac6_in(yn,2), temp6_in(yn,2), dt_in


! ! calculate velocities from streamfunction values
! velocities0 = velocities(psi)
! uf = velocities0(1:xn,1:yn)
! vf = velocities0(1:xn,yn+1:2*yn)
u = -1.0*u_in
v = -1.0*v_in
! uf = -1.0*uf!*rho_in
! vf = -1.0*vf!*rho_in


! do ii=2,yn-1
! 	do i=1,xn
! 		if ((maskP(i,ii) .eq. 2.5)) then
! 			u(i,ii) = 0.0
! 			v(i,ii) = 0.0
! 		end if
! 		if ((maskP(i,ii) .eq. 7.5)) then
! 			u(i,ii) = 0.0
! 			v(i,ii) = 0.0
! 		end if
! 	end do
! end do

! do ii=2,yn
! 	do i=1,xn
! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i,ii-1) .ne. 0.0)) then
! 			u(i,ii) = (uf(i,ii) + uf(i,ii-1))/(2.0)
! 		end if
! 	end do
! end do
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i-1,ii) .ne. 0.0)) then
! 			v(i,ii) = (vf(i,ii) + vf(i-1,ii))/(2.0)
! 		end if
! 	end do
! end do
!
!
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 9.0)) then
! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i+1,ii) - psi(i,ii))/dx
! 		end if
! 		if ((maskP(i,ii) .eq. 3.0)) then
! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i,ii) - psi(i-1,ii))/dx
! 		end if
! 	end do
! end do
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 6.0)) then
! 			v(i,ii) = (frac6_in(ii,2) - frac6_in(ii,1))/(2.0*param_f_dx)
! 		end if
! 	end do
! end do

! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 6.0)) then
! 			v(i,ii) = -1.0*(frac6_in(ii,2) - frac6_in(ii,1))/(rho_in(i,ii)*2.0*param_f_dx)
! 		end if
! 	end do
! end do

! ! fracture vertical velocity
! v(xn-1,:) = v(xn-1,:)*param_f_dx/dx






uLong = -1.0*reshape(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
vLong = -1.0*reshape(transpose(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

  
mn = h



h0 = h
hMid = h

qx = dt/(dx)
qy = dt/(dy)
sx = (2.0*dt_in*lambdaMat(1,1))/(dx*dx*rho_fluid*cp)
sy = (2.0*dt_in*lambdaMat(1,1))/(dy*dy*rho_fluid*cp)

! qxMat = 1.0*dt*4179.0*phi_in/(dx*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
! qyMat = 1.0*dt*4179.0*phi_in/(dy*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
! ! sxMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! ! syMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! sxMat = 2.0*dt*(lambdaMat)/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! syMat = 2.0*dt*(lambdaMat)/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)


qxMat = dt_in*4179.0/(dx*cp)
qyMat = dt_in*4179.0/(dy*cp)
sxMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*rho_in*cp)
syMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*rho_in*cp)


! ! print stability conditions at each timestep
! write(*,*) " "
! write(*,*) "velocity check"
! write(*,"(F10.5)") maxval(abs(u))*maxval(abs(qxMat))
! write(*,"(F10.5)") maxval(abs(v))*maxval(abs(qyMat))
! write(*,*) "conduction check"
! write(*,"(F10.5)") maxval(abs(syMat))
! write(*,*) " "
!
! write(*,*) "velocity check"
! write(*,*) maxval(abs(u))
! write(*,*) maxval(abs(v))

! vertical boundary conditions

! stretch
stretch = h0
stretch(2,:) = stretch(1,:)
stretch(xn-1,:) = stretch(xn,:)
do ii=2,yn-2
	do i=2,xn-1
		if ((mask(i,ii) .eq. 5.0) .or. (mask(i,ii) .eq. 12.5) ) then
			stretch(i,ii) = stretch(i+1,ii)
		end if
		
		if ((mask(i,ii) .eq. 3.0) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 3.1) .or. (mask(i,ii) .eq. 3.05)) then
			stretch(i,ii) = temp6_in(ii,1)
		end if
		
		if ((mask(i,ii) .eq. 10.0) .or. (mask(i,ii) .eq. 17.5)) then
			stretch(i,ii) = stretch(i-1,ii)
		end if
		
		if ((mask(i,ii) .eq. 6.0) .or. (mask(i,ii) .eq. 6.5) .or. (mask(i,ii) .eq. 6.1) .or. (mask(i,ii) .eq. 6.05)) then
			stretch(i,ii) = temp6_in(ii,2)
		end if
	end do
end do

stretchLong = reshape(stretch(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

do ii=2,yn-1
do i=2,xn-1
	! right outcrop (left boundary)
	if ((mask(i,ii) .eq. 17.5)) then
		h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 10.0)) then
		h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 6.0) .or. (mask(i,ii) .eq. 6.5) .or. (mask(i,ii) .eq. 6.1) .or. (mask(i,ii) .eq. 6.05)) then
		h(i,ii) = h(i,ii) + temp6_in(ii,2)*sxMat(i,ii)/2.0
	end if

	! left outcrop (right boundary)
	if ((mask(i,ii) .eq. 12.5)) then
		h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 5.0)) then
		h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 3.0) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 3.1) .or. (mask(i,ii) .eq. 3.05)) then
		h(i,ii) = h(i,ii) + temp6_in(ii,1)*sxMat(i,ii)/2.0
	end if



end do

if (mask(1,ii) .eq. 1.0) then
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! left
end if
if (mask(1,ii) .eq. 25.0) then
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! top left corner
end if

if (mask(xn,ii) .eq. 1.0) then
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! right
end if
if (mask(xn,ii) .eq. 50.0) then
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! top right corner
end if
end do



!h(2,2) = h(2,2) + h0(1,2)*sxMat(2,2)/2.0  ! bottom left corner
!h(xn-1,2) = h(xn-1,2) + h0(xn,2)*sxMat(xn-1,2)/2.0  ! bottom right corner
 
uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
h0Long = reshape(h0(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
sxLong = reshape(sxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
syLong = reshape(syMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
qxLong = reshape(qxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
qyLong = reshape(qyMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

! make the band
aBand = 0.0
aBand(:,1) = -sxLong/2.0 - uLong*qxLong/2.0
aBand(:,2) = 1.0+sxLong
aBand(:,3) = -sxLong/2.0 + uLong*qxLong/2.0

! bottom left corner, 2 and 3
aBand(1,1) =  0.0
aBand(1,2) = 1.0 + sxLong(1)/1.0 - uLong(1)*qxLong(1)!/2.0
aBand(1,3) = -sxLong(1)/2.0 + uLong(1)*qxLong(1)!/2.0

! top right corner, 1 and 2
aBand((xn-2)*(yn-2),1) = -sxLong((xn-2)*(yn-2))/2.0 - uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
aBand((xn-2)*(yn-2),2) = 1.0 + sxLong((xn-2)*(yn-2))/1.0 + uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
aBand((xn-2)*(yn-2),3) =  0.0



do i = 2,(xn-2)*(yn-2)-1
	
	! flow left anywhere, 2 and 3
	if (uLong(i) .lt. 0.0) then
	
		aBand(i,1) = -sxLong(i)/2.0 
		aBand(i,2) = 1.0+sxLong(i) - uLong(i)*qxLong(i)
		aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)
	
	end if
	

	! flow right anywhere, 1 and 2 
	if (uLong(i) .gt. 0.0) then
	
		aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
		aBand(i,2) = 1.0+sxLong(i) + uLong(i)*qxLong(i)
		aBand(i,3) = -sxLong(i)/2.0 
	
	end if



	! left edges, default 2 and 3
	if (((mod(i-1,xn-2).eq.0)) .or. (maskLong(i).eq.10.0) .or. (maskLong(i).eq.17.5) .or. (maskLong(i).eq.6.0) .or. (maskLong(i).eq.6.5) .or. (maskLong(i).eq.6.1) .or. (maskLong(i).eq.6.05)) then
			aBand(i,1) =  0.0
			aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
			aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)
	end if


	
	! left edge flowing to the right, 1 & 2
	if (uLong(i) .gt. 0.0) then
	
		! left edge of right outcrop
		if ((maskLong(i).eq.10.0) .or. (maskLong(i).eq.17.5) .or. (maskLong(i).eq.6.0) .or. (maskLong(i).eq.6.5) .or. (maskLong(i).eq.6.1) .or. (maskLong(i).eq.6.05)) then
				aBand(i,1) =  0.0
				aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
				aBand(i,3) = -sxLong(i)/2.0 
				uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
		end if
		
		! left edge but not uppper left corner
		if ((mod(i-1,xn-2).eq.0) .and. (maskLong(i) .ne. 25.0)) then
				aBand(i,1) =  0.0
				aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
				aBand(i,3) = -sxLong(i)/2.0 
				uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
		end if
	

	end if
	
	
	! right edge, 1 and 2 by default
	if ((mod(i,xn-2) .eq. 0) .or. (maskLong(i).eq.5.0) .or. (maskLong(i).eq.12.5) .or. (maskLong(i).eq.3.0) .or. (maskLong(i).eq.3.5) .or. (maskLong(i).eq.3.1) .or. (maskLong(i).eq.3.05)) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
			aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
			aBand(i,3) =  0.0
	end if

	
	! right edge flowing to the left, 2 and 3
	if (uLong(i) .lt. 0.0) then
	
		! right edge of left outcrop
		if ((maskLong(i).eq.5.0) .or. (maskLong(i).eq.12.5) .or. (maskLong(i).eq.3.0) .or. (maskLong(i).eq.3.5) .or. (maskLong(i).eq.3.1) .or. (maskLong(i).eq.3.05)) then
				aBand(i,1) = -sxLong(i)/2.0 
				aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
				aBand(i,3) =  0.0
				uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
		end if
		
		! right edge but not upper right corner
		if ((mod(i,xn-2) .eq. 0) .and. (maskLong(i) .ne. 25.0)) then
				aBand(i,1) = -sxLong(i)/2.0 
				aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
				aBand(i,3) =  0.0
				uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
		end if


	end if 
	
	
	! upper left corner, default 2 and 3
	if ((mod(i-1,xn-2).eq.0) .and. (maskLong(i) .eq. 25.0)) then
			aBand(i,1) =  0.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)!/2.0
	end if

	! upper right corner, default 1 and 2
	if ((mod(i,xn-2).eq.0) .and. (maskLong(i) .eq. 25.0)) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
			aBand(i,3) =  0.0
	end if
	
	! bottom right corner, 1 and 2
	if (i.eq.xn-2) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
			aBand(i,3) =  0.0
	end if


	



!

end do

do i=1,(xn-2)*(yn-2)
	! mask
	if ((maskLong(i) .eq. 0.0) .or. (maskLong(i) .eq. 600.0)) then
		aBand(i,2) = 1.0
		aBand(i,1) = 0.0
		aBand(i,3) = 0.0
	end if

end do
  
! make sure solver doesn't go out of bounds
do i=1,((yn-2)-1)
	ii = i*(xn-2)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do

!!!!!!!!!!!! THIS !!!!!!!!!!!
h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))




stretch = h0
stretch(:,2) = stretch(:,1)
do ii=2,yn-1
	do i=2,xn-1
		if ((mask(i,ii) .eq. 25.0) .or. (mask(i,ii) .eq. 50.0)  .or. (mask(i,ii) .eq. 17.5)  .or. (mask(i,ii) .eq. 12.5)  .or. (mask(i,ii) .eq. 3.5)  .or. (mask(i,ii) .eq. 6.5)) then
			stretch(i,ii) = stretch(i,ii+1)
		end if

	end do
end do

stretchLongT = reshape(transpose(stretch(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))




!h0LongT = reshape(transpose(h0(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
!h0T = transpose(h0)
sxLong = reshape(transpose(sxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
syLong = reshape(transpose(syMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
qxLong = reshape(transpose(qxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
qyLong = reshape(transpose(qyMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! horizontal boundary conditions
h(3:xn-2,2) = h(3:xn-2,2) + h0(3:xn-2,1)*syMat(3:xn-2,2)/2.0 ! bottom


! top of sediment
do ii=2,yn-1
do i=2,xn-1
	! top of sediment and any short outcrops
	if ((mask(i,ii) .eq. 50.0) .or. (mask(i,ii) .eq. 2.0) .or. (mask(i,ii) .eq. 7.0)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	end if
	
	if ((mask(i,ii) .eq. 25.0) .and. (i .ge. 3) .and. (i .le. xn-2)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top 
	end if
	if ((mask(i,ii) .eq. 25.0) .and. (i .eq. 2)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top left corner
	end if
	if ((mask(i,ii) .eq. 25.0) .and. (i .eq. xn-1)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top right corner
	end if
	
	if ((mask(i,ii) .eq. 1.0) .and. (ii .eq. 2) .and. (i .eq. 2)) then
		h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom left corner
	end if
	if ((mask(i,ii) .eq. 1.0) .and. (ii .eq. 2) .and. (i .eq. xn-1)) then
		h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom right corner
	end if
	

	if ((mask(i,ii) .eq. 12.5) .or. (mask(i,ii) .eq. 17.5) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 6.5)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	end if


end do
end do

!genTrans = transpose(h)
!h_nextRow = reshape(genTrans(2:yn-1,2:xn-1), (/(xn-2)*(yn-2)/))
h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! make the band
bBand = 0.0
bBand(:,1) = -syLong/2.0 - vLong(:)*qyLong/2.0
bBand(:,2) = 1.0+syLong
bBand(:,3) = -syLong/2.0 + vLong*qyLong/2.0

! bottom left corner
bBand(1,1) =  0.0
bBand(1,2) = 1.0 + syLong(1)/1.0 - vLong(1)*qyLong(1)!/2.0
bBand(1,3) = -syLong(1)/2.0 + vLong(1)*qyLong(1)!/2.0

! top right corner
bBand((xn-2)*(yn-2),1) = -syLong((xn-2)*(yn-2))/2.0 - vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
bBand((xn-2)*(yn-2),2) = 1.0 + syLong((xn-2)*(yn-2))/1.0 + vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
bBand((xn-2)*(yn-2),3) =  0.0
do i = 2,(xn-2)*(yn-2)-1
	
		
	! flow going down anywhere, 2 and 3
	if (vLong(i) .lt. 0.0) then
		
		bBand(i,1) = -syLong(i)/2.0 
		bBand(i,2) = 1.0+syLong(i) - vLong(i)*qyLong(i)
		bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)

	end if
	
	! flow coming up anywhere, 1 and 2
	if (vLong(i) .gt. 0.0) then
		
		bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
		bBand(i,2) = 1.0+syLong(i) + vLong(i)*qyLong(i)
		bBand(i,3) = -syLong(i)/2.0 
		
	end if


	!!!!! TOP EDGES
	
	
	! bottom rrow, default 2 and 3 !! 
	if (mod(i-1,yn-2) .eq. 0) then
			bBand(i,1) =  0.0
			bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
			bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)
	end if
	
	! bottom rrow, if flow is coming up (1 and 2)
	if (vLong(i) .gt. 0.0) then
		
		! bottom row but not bottom right corner
		if ((mod(i-1,yn-2) .eq. 0) .and. (i .ne. (xn-2)*(yn-2)-(yn-2)+1)) then
				bBand(i,1) =  0.0
				bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
				bBand(i,3) = -syLong(i)/2.0 
				h_nextRow(i) = h_nextRow(i) + vLong(i)*qyLong(i)*stretchLongT(i)
		end if

	end if
	

	! last/top edge, default 1 and 2
	if ((maskLongT(i).eq.25.0) .or. (maskLongT(i).eq.50.0) .or. (maskLongT(i).eq.2.0) .or. (maskLongT(i).eq.7.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
			bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
			bBand(i,3) =  0.0
			
	end if
	
	if ((maskLongT(i).eq.12.5).or.(maskLongT(i).eq.17.5).or.(maskLongT(i).eq.3.5).or.(maskLongT(i).eq.6.5)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0
	end if
	
	
	! top edges if flow coming down (2, 3)
	if (vLong(i) .lt. 0.0) then
		
		if ((maskLongT(i).eq.50.0) .or. (maskLongT(i).eq.2.0) .or. (maskLongT(i).eq.7.0)) then
				bBand(i,1) = -syLong(i)/2.0 
				bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)
		end if
		
		! top but not top left or top right corner
		if ((maskLongT(i).eq.25.0) .and. (i .gt. yn-2) .and. (i .le. (xn-2)*(yn-2)-(yn-2))) then
				bBand(i,1) = -syLong(i)/2.0 
				bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)
			
		end if
	
		if ((maskLongT(i).eq.12.5).or.(maskLongT(i).eq.17.5).or.(maskLongT(i).eq.3.5).or.(maskLongT(i).eq.6.5)) then
				bBand(i,1) = -syLong(i)/2.0 
				bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)!/2.0
		end if
		
	end if 
	
	
	
	
	! upper left corner
	if ((i .le. (yn-2)) .and. (maskLongT(i) .eq. 25.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0
	end if
	
	! upper right corner
	if ((i .gt. (yn-2)*(xn-2)-(yn-2)) .and. (maskLongT(i) .eq. 25.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0
			
	end if
	

	
	! bottom right corner
	if (i .eq. (yn-2)*(xn-2)-(yn-2)+1) then
			bBand(i,1) =  0.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)!/2.0
	end if




end do

do i=1,(xn-2)*(yn-2)
	! mask
	if ((maskLongT(i) .eq. 0.0) .or. (maskLongT(i) .eq. 600.0)) then
		bBand(i,2) = 1.0
		bBand(i,1) = 0.0
		bBand(i,3) = 0.0
	end if

end do
 
! make sure solver doesn't go out of bounds
do i=1,((xn-2)-1)
	ii = i*(yn-2)
	bBand(ii,3) = 0.0
	bBand(ii+1,1) = 0.0
end do

h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))
h_next(2:xn-1,2:yn-1) = transpose(reshape(h_nextRow, (/yn-2, xn-2/))) 
hMid = h_next
!
!
!

! ! something awful right here
! ! possible not-doing-upwind correction?
! do ii=2,yn-1
! do i = 2,xn-1
! 	if ( (hMid(i,ii).gt.h0(i-1,ii)).and.(hMid(i,ii).gt.h0(i+1,ii)).and. &
! 	(hMid(i,ii).gt.h0(i,ii-1)).and.(hMid(i,ii).gt.h0(i,ii+1))) then
! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
! 	end if
!
! 	if ( (hMid(i,ii).lt.h0(i-1,ii)).and.(hMid(i,ii).lt.h0(i+1,ii)).and. &
! 	(hMid(i,ii).lt.h0(i,ii-1)).and.(hMid(i,ii).lt.h0(i,ii+1))) then
! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
! 	end if
!
!
! end do
! end do




! check out how this equation is converging to steady-state
!write(*,*) "deltaT"
!write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-h_next(2:xn-1,2:yn-1))/h_next(2:xn-1,2:yn-1)))

return

end function h_next


! ----------------------------------------------------------------------------------%%
!
! PSI_NEXT
!
! SUMMARY: computes the 2D streamfunction array of the current timestep
!
! INPUTS: h(xn,yn) : temperature profile
!         rhs0(xn,yn) : right hand side of streamfunction-vorticity equation
!         psi(xn,yn) : 2D streamfunction array of previous timestep
!         top_in(xn,1) : permeable upper boundary
!         rho_in(xn,yn) : 2D density array
!
! RETURNS: psi_next(xn,yn): 2D streamfunction array for current timestep
!
! ----------------------------------------------------------------------------------%%

function psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage, frac6_in)

use globals
use initialize
implicit none

! interface
!
! 	function partial(array,rows,cols,d1,d2,dim)
! 		use globals
! 		use initialize
! 		implicit none
! 		integer :: rows, cols, dim, i, j, ii, jj
! 		real(8) :: array(rows,cols), d1, d2, d
! 		real(8) :: partial(rows,cols)
! 	end function partial
!
! end interface

! declare errthing

! integers
integer :: i, j, ii, n, m, stage
! inputs
real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
! matrix stuff
real(8) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
real(8) :: psi_next(xn,yn)
real(8) :: mn(xn,yn)
! back to band
real(8) :: aBand0(longP,2*((yn/2)-2) + 1), band_in(longP,2*((yn/2)-2) + 1)
real(8) :: rhoLong(longP)
real(8) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
real(8) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
real(8) :: psi_f

psi_next = 0.0

mn = psi

m = 2*((yn/2)-2) + 1

rhs1 = rhs0

rho_in = rho_fluid

!phi_in = 1.0

permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
do ii=2,yn-1
do i=2,xn-1
	permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
	permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)
	
	if (maskP(i,ii) .eq. 5.0) then
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
	end if

	permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
	permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)

	if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
		permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	end if
	if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1) .or. (maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	end if
end do
end do

do ii=2,yn-1
do i=2,xn-1
	if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
		permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
	end if
end do
end do


! do ii=yn/2,yn-1
! do i=2,xn-1
! 	if ((y(ii) .ge. sed(i)-dy) .and. (any(maskP(i,:) .eq. 50.0))) then
! 		permy_top(i,ii) = phi_in(i,ii) / ((1e-16 + 1e-16) / 2.0)
! 	end if
! end do
! end do

do ii=yn/2,yn
	if (maskP(1,ii) .eq. 1.0) then
	 ! left
	 rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
	end if
	
	if (maskP(xn,ii) .eq. 1.0) then
	! right
	 rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
	end if

	! bottom left corner
	if (maskP(1,ii) .eq. 100.0) then
	 ! left
		rhs1(2,ii) = rhs1(2,ii) + psi(2,ii-1)*permy_bottom(2,ii)/(dy*dy)
		rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
	end if
	
	! bottom right corner
	if (maskP(xn,ii) .eq. 100.0) then
	 ! left
		rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn-1,ii-1)*permy_bottom(xn-1,ii)/(dy*dy)
		rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
	end if

	
end do




 ! mask boundary conditions
 ! vertical boundary conditions
do ii=yn/2,yn-1
do i=2,xn-1

	if ((maskP(i,ii) .eq. 17.5)) then
	    rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 10.0)) then
	    !rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
		rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i+1,ii))*permx_left(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
	    rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,2)*permx_left(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 12.5)) then
	    rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 5.0)) then
		!rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
		rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i-1,ii))*permx_right(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
		rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
	end if
	
	if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1)) then
		rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
		!rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*phi_in(i,ii)*2.0*viscosity/(grav*rho_fluid*(rho_fluid*perm_in(i,ii) + (rho_fluid*param_f_dx*param_f_dx/12.0) ))
		!rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*(phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/12.0)*rho_fluid) / 2.0))
	end if
	

end do
end do




!
! ! at the bottom of everything
do ii=yn/2,yn-1
do i=3,xn-2
	if ((maskP(i,ii) .eq. 100.0) .or. (maskP(i,ii) .eq. 3.01) .or. (maskP(i,ii) .eq. 6.01)) then
		rhs1(i,ii) = rhs1(i,ii) + psi(i,ii-1)*permy_bottom(i,ii)/(dy*dy)
	end if
end do
end do

! do ii=yn/2,yn
! do i=2,xn-1
!
! if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 25.0) .or. (maskP(i,ii) .eq. 12.5) .or. (maskP(i,ii) .eq. 17.5)) then
! 	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
! end if
!
! end do
! end do

do ii=yn/2,yn
do i=2,xn-1

if ((maskP(i,ii) .eq. 25.0) .or. (maskP(i,ii) .eq. 12.5) .or. (maskP(i,ii) .eq. 17.5) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
end if

end do
end do

do ii=yn/2,yn
do i=2,xn-1

if ((maskP(i,ii) .eq. 50.0)) then
	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
end if

end do
end do

! ! CORNER 2
! do ii=yn/2,yn
! do i=2,xn-1
!
! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i-1,ii) .ne. 50.0)) then
! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
! end if
!
! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i+1,ii) .ne. 50.0)) then
! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
! end if
!
! end do
! end do


do ii=1,yn
do i=1,xn
	if ((maskP(i,ii) .eq. 0.0) .or. (maskP(i,ii) .eq. 600.0)) then
		rhs1(i,ii) = 0.0
	end if
end do
end do


uVec = reshape(transpose(rhs1( 2:xn-1 , (yn/2)+2:yn-1 )),(/longP/))

psi_next = 0.0



 

! THIS IS WHERE THE BAND IS MADE
aband0 = band_in

! use the banded solver here
psi_nextRow = solve(aBand0,uVec,2*((yn/2)-2) + 1,longP)
psi_next(2:xn-1,(yn/2)+2:yn-1) = transpose(reshape(psi_nextRow, (/(yn/2)-2, xn-2/)))




!write(*,*) "deltaPSI"
!write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))


return

end function psi_next




! ----------------------------------------------------------------------------------%%
!
! MAKE BAND
!
! SUMMARY: only make the big matrix every once in a while
!
!
! ----------------------------------------------------------------------------------%%





function make_band(perm_in,phi_in,permx,permy,rho_in)
	
	
	use globals
	use initialize
	implicit none
	
	
	interface

	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial
	
	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge(rows,cols)
	end function partial_edge
	
	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge_p(rows,cols)
	end function partial_edge_p
	
	end interface
	
	
	integer :: i, j, ii, n, m
	real(8) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
	real(8) :: permx(xn,yn), permy(xn,yn), permLong(longP)
	real(8) :: permxLong(longP), permyLong(longP)
	real(8) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
	real(8) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
	real(8) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
	real(8) :: perm_long(longP)
	
!	phi_in = 1.0
	
	rho_in = rho_fluid

	permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	do ii=2,yn-1
	do i=2,xn-1
		permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)
		
		if (maskP(i,ii) .eq. 5.0) then
			permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
		end if

		permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
		permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)
	
		if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
			permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
		end if
		if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1) .or. (maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
			permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
		end if
	end do
	end do
	
	do ii=2,yn-1
	do i=2,xn-1
		if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
			permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
		end if
	end do
	end do


	permx_left_long = reshape(transpose(permx_left(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permx_right_long = reshape(transpose(permx_right(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permy_bottom_long = reshape(transpose(permy_bottom(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permy_top_long = reshape(transpose(permy_top(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	
	
! 	do ii=2,yn-1
! 	do i=2,xn-1
! 		if (maskP(i,ii) .eq. 50.0) then
! 			permy_top(:,ii) = perm_in(:,ii)! phi_in(i,ii) / ((perm_in(i,ii) + perm_in(i,ii+1)) / 2.0)
! 		end if
! 	end do
! 	end do
	
	! make the band
	innerBand = 0.0
	make_band = 0.0
	m = 2*((yn/2)-2) + 1
	

	! trial (7 days) inner band
	innerBand = 0.0
		do i = 1,longP
			! diagonal
			innerBand(i,(m+1)/2) = (permx_right_long(i) + permx_left_long(i))/(dx*dx) + (permy_top_long(i) + permy_bottom_long(i))/(dy*dy)
			
			
			! off-diagonals
			! all but left edge
			if ((i .gt. ((yn/2)-2)) .and. (maskPLongT(i) .ne. 10.0) .and. (maskPLongT(i) .ne. 17.5) .and. (maskPLongT(i) .ne. 6.0) .and. (maskPLongT(i) .ne. 6.5) .and. (maskPLongT(i) .ne. 6.1) .and. (maskPLongT(i) .ne. 6.05) .and. (maskPLongT(i) .ne. 6.01)) then
				innerBand(i,1) = -permx_left_long(i)/(dx*dx)
			end if
			! all but right edge
			if ((i .le. (longP)-((yn/2)-2)) .and. (maskPLongT(i) .ne. 5.0) .and. (maskPLongT(i) .ne. 12.5) .and. (maskPLongT(i) .ne. 3.0) .and. (maskPLongT(i) .ne. 3.5) .and. (maskPLongT(i) .ne. 3.1) .and. (maskPLongT(i) .ne. 3.05) .and. (maskPLongT(i) .ne. 3.01)) then
				innerBand(i,m) = -permx_right_long(i)/(dx*dx)
			end if
		
		
			! all but bottom
			if ((maskPLongT(i) .ne. 100.0) .or. (maskPLongT(i) .ne. 3.01) .or. (maskPLongT(i) .ne. 6.01)) then
				innerBand(i,(m+1)/2-1) = -permy_bottom_long(i)/(dy*dy)
			end if
		
			! all but top
			if ((maskPLongT(i) .ne. 50.0) .and. (maskPLongT(i) .ne. 25.0) .and. (maskPLongT(i) .ne. 12.5) .and. (maskPLongT(i) .ne. 17.5) .and. (maskPLongT(i) .ne. 3.5) .and. (maskPLongT(i) .ne. 6.5)) then
				innerBand(i,(m+1)/2+1) = -permy_top_long(i)/(dy*dy)
			end if
			
! 			if (maskPLongT(i) .eq. 2.0) then
! 				innerBand(i,m) = 0.0 ! skip right
! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
! 			end if
!
! 			if (maskPLongT(i) .eq. 7.0) then
! 				innerBand(i,1) = 0.0 ! skip left
! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
! 			end if
		
		end do

	do i = 1,longP
	 
		! mask
		if ((maskPLongT(i) .eq. 0.0) .or. (maskPLongT(i) .eq. 600.0)) then
			innerBand(i,:) = 0.0
			innerBand(i,(m+1)/2) = 1.0!(2.0)/(permLong(i)*dx*dx) + (2.0)/(permLong(i)*dy*dy)
		end if

	end do

	!write(*,*) innerBand(1,:)
	
	make_band = innerBand
	
return
end function make_band








! ----------------------------------------------------------------------------------%%
!
! PARTICLES NEXT
!
! SUMMARY: transports particles by euler or 4th order runge kutta
!
! INPUTS: 
!
!
! RETURNS: 
!
! ----------------------------------------------------------------------------------%%

function particles_next (trace_in, uTransport, vTransport, inval, num, num_sat)
	
use globals
use initialize
implicit none

! declare errthing

! integers
integer :: i, j, ii, n, m, mm, nn, num, num_sat
! inputs
real(8) :: trace_in(5,num), particles_next(5,num)
real(8) :: uTransport(xn,yn), vTransport(xn,yn)
real(8) :: u_wt, v_wt, inval
real(8) :: rando
real(8) :: r_i, r_ii


do mm=1,cstep

! delete particles out of bound

! do n = 1,num
! 	if ((trace_in(3,n) .ne. 0.0)) then
! 		i = floor(trace_in(1,n)/dx) + 1
! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
! 		if ((i .lt. 1) .or. (i .ge. xn) .or.  (ii .lt. 1) .or. (ii .ge. yn) )then
! 			trace_in(:,n) = 0.0
! 		end if
! 	end if
! end do
!
! do n = 1,num
! 	if ((trace_in(3,n) .ne. 0.0)) then
! 		i = floor(trace_in(1,n)/dx) + 1
! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
! 		if (maskP(i,ii) .eq. 0.0) then
! 			trace_in(:,n) = 0.0
! 		end if
!
! 		if ((maskP(i,ii+1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii+1)) .lt. dy/2.0)) then
! 			trace_in(:,n) = 0.0
! 		end if
!
! 		if ((maskP(i,ii-1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii-1)) .lt. dy/2.0)) then
! 			trace_in(:,n) = 0.0
! 		end if
!
! 	end if
! end do


! put in new particles in inflow cells

do nn = 1, xn
	
		if (vTransport(twentyfives(1,nn),twentyfives(2,nn)) .lt. 0.0) then
			
		! generate a new particle so the cell is at saturation
		do j=1,num_sat
			m = 0
			n = 0
			! find an empty slot to put the particle in
			do while (m .eq. 0)
				n = n +1
				if ((trace_in(3,n) .eq. 0.0)) then
					call random_number(rando) 
					trace_in(1,n) = dx*rando + (x(twentyfives(1,nn)) - dx/2)
					trace_in(2,n) = dy*rando + (y(twentyfives(2,nn)) - dy/2)
					trace_in(3,n) = inval
					
					trace_in(4,n) = twentyfives(1,nn)
					trace_in(5,n) = twentyfives(2,nn)
					m = 1
					
				end if
			end do
		end do
		
		end if

end do


do nn = yn/2, yn
	
		if (uTransport(fives(1,nn),fives(2,nn)) .lt. 0.0) then
			
		! generate a new particle so the cell is at saturation
		do j=1,num_sat
			m = 0
			n = 0
			! find an empty slot to put the particle in
			do while (m .eq. 0)
				n = n +1
				if ((trace_in(3,n) .eq. 0.0)) then
					call random_number(rando) 
					trace_in(1,n) = dx*rando + (x(fives(1,nn)) - dx/2)
					trace_in(2,n) = dy*rando + (y(fives(2,nn)) - dy/2)
					trace_in(3,n) = inval
					
					trace_in(4,n) = fives(1,nn)
					trace_in(5,n) = fives(2,nn)
					m = 1
				end if
			end do
		end do
		
		end if
		
		
		if (uTransport(tens(1,nn),tens(2,nn)) .gt. 0.0) then
			
		! generate a new particle so the cell is at saturation
		do j=1,num_sat
			m = 0
			n = 0
			! find an empty slot to put the particle in
			do while (m .eq. 0)
				n = n +1
				if ((trace_in(3,n) .eq. 0.0)) then
					call random_number(rando) 
					trace_in(1,n) = dx*rando + (x(tens(1,nn)) - dx/2)
					trace_in(2,n) = dy*rando + (y(tens(2,nn)) - dy/2)
					trace_in(3,n) = inval
					
					trace_in(4,n) = tens(1,nn)
					trace_in(5,n) = tens(2,nn)
					m = 1
				end if
			end do
		end do
		
		end if

end do


! do nn = 1, xn
!
! 		if (vTransport(fifties(1,nn),fifties(2,nn)) .lt. 0.0) then
!
! 		! generate a new particle so the cell is at saturation
! 		do j=1,num_sat
! 			m = 0
! 			n = 0
! 			! find an empty slot to put the particle in
! 			do while (m .eq. 0)
! 				n = n +1
! 				if ((trace_in(3,n) .eq. 0.0)) then
! 					call random_number(rando)
! 					trace_in(1,n) = dx*rando + (x(fifties(1,nn)) - dx/2)
! 					trace_in(2,n) = dy*rando + (y(fifties(2,nn)) - dy/2)
! 					trace_in(3,n) = inval
!
! 					trace_in(4,n) = fifties(1,nn)
! 					trace_in(5,n) = fifties(2,nn)
! 					m = 1
! 				end if
! 			end do
! 		end do
!
! 		end if
!
! end do



! advect all of the particles

do n = 1,num

	if ((trace_in(3,n) .ne. 0.0)) then
! 		i = floor(trace_in(1,n)/dx) + 1
! 		i = min(i,xn-1)
! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
! 		ii = min(ii,yn-1)
		i = trace_in(4,n)
		ii = trace_in(5,n)


					u_wt = ( uTransport(i,ii) * (y(ii+1) - trace_in(2,n)) * (x(i+1) - trace_in(1,n)) + &
						   uTransport(i+1,ii) * (y(ii+1) - trace_in(2,n)) * (trace_in(1,n) - x(i)) + &
						   uTransport(i+1,ii+1) * (trace_in(2,n) - y(ii)) * (trace_in(1,n) - x(i)) + &
						   uTransport(i,ii+1) * (trace_in(2,n) - y(ii)) * (x(i+1) - trace_in(1,n)) ) / (dx*dy) 
			
   					v_wt = ( vTransport(i,ii) * (y(ii+1) - trace_in(2,n)) * (x(i+1) - trace_in(1,n)) + &
   						   vTransport(i+1,ii) * (y(ii+1) - trace_in(2,n)) * (trace_in(1,n) - x(i)) + &
   						   vTransport(i+1,ii+1) * (trace_in(2,n) - y(ii)) * (trace_in(1,n) - x(i)) + &
   						   vTransport(i,ii+1) * (trace_in(2,n) - y(ii)) * (x(i+1) - trace_in(1,n)) ) / (dx*dy) 

! 					u_wt = uTransport(i,ii)
! 					v_wt = vTransport(i,ii)

		
		! advect the thing
		trace_in(1,n) = trace_in(1,n) + u_wt*dt*mstep/(cstep)!/1000.0
		trace_in(2,n) = trace_in(2,n) + v_wt*dt*mstep/(cstep)!/1000.0
		
		if ((trace_in(1,n).gt.x_max) .or. (trace_in(1,n).lt.x_min) .or. (trace_in(2,n).gt.y_max) .or. (trace_in(2,n).lt.y_min)) then
			trace_in(:,n) = 0.0
		end if
		
! 		i = 1.0*floor(trace_in(1,n)/dx) + 1
! 		i = 1.0*min(i,xn-1)
! 		trace_in(4,n) = i
!
! 		ii = yn - floor(-1.0*trace_in(2,n)/dy) - 1
! 		ii = 1.0*min(ii,yn-1)
! 		trace_in(5,n) = ii

		trace_in(4,n) = min(floor(trace_in(1,n)/dx) + 1,xn-1)
		trace_in(5,n) = min(yn - floor(-1.0*trace_in(2,n)/dy) - 1,yn-1)
		i = trace_in(4,n)
		ii = trace_in(5,n)
		
		if (maskP(i,ii) .eq. 0.0) then
			trace_in(:,n) = 0.0
		end if
		
		if ((maskP(i,ii+1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii+1)) .lt. dy/2.0)) then
			trace_in(:,n) = 0.0
		end if
		
		if ((maskP(i,ii-1) .eq. 0.0) .and. (abs(trace_in(2,n) - y(ii-1)) .lt. dy/2.0)) then
			trace_in(:,n) = 0.0
		end if
!
! 		if (permeability(i,ii) .le. 1e-13) then
! 			trace_in(:,n) = 0.0
! 		end if
!
	end if
end do


end do	

particles_next = trace_in

return


end function particles_next


! ----------------------------------------------------------------------------------%%
!
! SOLUTE_NEXT
!
! SUMMARY: transports solutes on coarse mesh
!
! INPUTS: sol(xn/cell,yn/cell) : 2D array of initial solute concentrations
!         uTransport(xn/cell,yn/cell) : lateral velocities (coarse mesh)
!         vTransport(xn/cell,yn/cell) : vertical velocities (coarse mesh)
!
! RETURNS: solute_next(xn/cell,yn/cell): 2D array of solute concentrations
!
! ----------------------------------------------------------------------------------%%

function solute_next (sol, uTransport, vTransport, seaw)
	
use globals
use initialize
implicit none

! declare errthing

! integers
integer :: i, j, ii, n, m
! inputs
real(8) :: sol(xn/cell,yn/cell), sol0(xn/cell,yn/cell)
real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
! solver stuff
! real(8) :: uLong((xn/cell-2)*(yn/cell-2)), vLong((xn/cell-2)*(yn/cell-2))
! real(8) :: aBand((xn/cell-2)*(yn/cell-2),5), bBand((xn/cell-2)*(yn/cell-2),5)
! real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec((xn/cell-2)*(yn/cell-2))
! real(8) :: sol_nextRow((xn/cell-2)*(yn/cell-2))
real(8) :: uLong(((xn/cell)-2)*((yn/cell)-0)), vLong(((xn/cell)-0)*((yn/cell)-2))
real(8) :: aBand(((xn/cell)-2)*((yn/cell)-0),5), bBand(((xn/cell)-0)*((yn/cell)-2),5)
real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec(((xn/cell)-2)*((yn/cell)-0))
real(8) :: sol_nextRow(((xn/cell)-2)*((yn/cell)-0)), sol_nextRowB(((xn/cell)-0)*((yn/cell)-2))
real(8) :: seaw
real(8) :: bm1(xn,yn), b0(xn,yn), bp1(xn,yn), correction, sigma1, sigma2, sigma1a, sigma1b, sigma2a, sigma2b
real(8) :: sigma3, sigma4, sigma3a, sigma3b, sigma4a, sigma4b, sigma5, sigma6



do i = 1,xn
	do j = 1,yn
		if ((maskP(i,j) .eq. 0.0)) then
			sol(i,j) = seaw
		end if
	end do
end do

sol(1,:) = (4.0/3.0)*sol(2,:) - (1.0/3.0)*sol(3,:)
sol(xn,:) = (4.0/3.0)*sol(xn-1,:) - (1.0/3.0)*sol(xn-2,:)

sol0 = sol
solute_next = sol

! ! cstep_num = 6.28e10/(mstep)
! ! cstep_int = dx*cstep_num/(maxval(uTransport)*mstep)
!
! cstep_int = 6.28e11/tn
! cstep_num = 10.0*cstep_int*mstep*maxval(abs(vTransport))/dy
!
!
!
! ! cstep_int = 6.28e11/tn
! ! cstep_num = cstep_int*mstep*maxval(u/phi(1,1))/dx
!
! if (cstep_num .le. 1) then
! 	cstep_num = 1
! end if

!
! write(*,*) "cstep_int"
! write(*,*) cstep_int
! write(*,*) "cstep_num"
! write(*,*) cstep_num

! qx = cstep_int*mstep/(cstep_num*dx)
! qy = cstep_int*mstep/(cstep_num*dy)

qx = dt*mstep/(dx*ar)
dy = dy*mstep/(dy*ar)

! uLong = reshape(uTransport(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
! !! transpose coarse needed!
! vLong = reshape(transpose(vTransport(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))

uTransport(1,:) = 0.0
uTransport(xn,:) = 0.0

! write(*,*) qx*maxval(abs(uTransport))
! write(*,*) qy*maxval(abs(vTransport))



do i = 2,xn-1
	do j = yn/2,yn
		if ((maskP(i,j) .eq. 5.0) .or. (maskP(i,j) .eq. 120.5)) then
			
			if (uTransport(i,j) .lt. 0.0) then
				sol(i+1,j) = seaw
				uTransport(i+1,j) = uTransport(i,j)
			else
				sol(i+1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i-1,j)
			end if
			
		end if

		if ((maskP(i,j) .eq. 10.0) .or. (maskP(i,j) .eq. 170.5)) then
			if (uTransport(i,j) .gt. 0.0) then
				sol(i-1,j) = seaw
				uTransport(i-1,j) = uTransport(i,j)
			else
				sol(i-1,j) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i+1,j)
			end if
		end if
	end do
end do


! bm1 = (qx*uTransport/2.0)*( (qx*uTransport) + 1.0)
! b0 = 1.0 - ( (qx*uTransport) * (qx*uTransport) )
! bp1 = (qx*uTransport/2.0)*( (qx*uTransport) - 1.0)
!
! do i = 2,xn-1
! 	do j = yn/2,yn
! 		if (maskP(i,j) .ge. 1.0) then
! 			solute_next(i,j) = bm1(i-1,j)*sol(i-1,j) + b0(i,j)*sol(i,j) + bp1(i+1,j)*sol(i+1,j)
! 		end if
! 	end do
! end do



do i = 2,xn-1
	do j = yn/2,yn-1
		!if ( (solute_next(i,j).gt.sol0(i,j)) .and. (solute_next(i,j).gt.sol0(i+1,j)) .and. (solute_next(i,j).gt.sol0(i-1,j)) ) then
			
			
		
			if (uTransport(i,j) .gt. 0.0) then
				! upwind including LHS value
				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i,j) - sol0(i-1,j) )
				
! 				! correction loop: sort of a mess
! 				if (i .gt. 2) then
! 					if (maskP(i-2,j) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i+1,j) - sol0(i,j))/dx
! 						sigma1b = (sol0(i,j) - sol0(i-1,j))/dx
!
! 						if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = sign(1.0D+00,sigma1a)*minval((/abs(sigma1a), abs(sigma1b)/))
! 						end if
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i,j) - sol0(i-1,j))/dx
! 						sigma2b = (sol0(i-1,j) - sol0(i-2,j))/dx
!
! 						if (sigma2a*sigma2b .gt. 0.0) then
! 							sigma2 = sign(1.0D+00,sigma2a)*minval((/abs(sigma2a), abs(sigma2b)/))
! 						end if
!
! 						correction = (uTransport(i,j)*(dt*mstep/cstep)*0.5/dx) * (sigma1 - sigma2) * (dx - uTransport(i,j)*(dt*mstep/cstep))
! 						solute_next(i,j) = solute_next(i,j) - correction
!
! 					end if
! 				end if

				! correction loop: sort of a mess
				if (i .gt. 2) then
					if (maskP(i-2,j) .ne. 0.0) then
						sigma1 = 0.0
						sigma1a = (sol0(i+1,j) - sol0(i,j))/dx
						sigma1b = 2.0*(sol0(i,j) - sol0(i-1,j))/dx
						
						!if (sigma1a*sigma1b .gt. 0.0) then
							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
							if (sigma1 .eq. abs(sigma1a)) then
								sigma1 = sign(1.0D+00,sigma1a)*sigma1
							end if
							if (sigma1 .eq. abs(sigma1b)) then
								sigma1 = sign(1.0D+00,sigma1b)*sigma1
							end if
							!end if
						
						sigma3 = 0.0
						sigma3a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
						sigma3b = (sol0(i,j) - sol0(i-1,j))/dx
						
						!if (sigma3a*sigma3b .gt. 0.0) then
						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
						if (sigma3 .eq. abs(sigma3a)) then
							sigma3 = sign(1.0D+00,sigma3a)*sigma3
						end if
						if (sigma3 .eq. abs(sigma3b)) then
							sigma3 = sign(1.0D+00,sigma3b)*sigma3
						end if
							!end if
						
						
						! choosing sigma5
						sigma5 = 0.0
						!if (sigma1*sigma3 .gt. 0.0) then
							sigma5 = sign(1.0D+00,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
							!end if



						

						sigma2 = 0.0
						sigma2a = (sol0(i,j) - sol0(i-1,j))/dx
						sigma2b = 2.0*(sol0(i-1,j) - sol0(i-2,j))/dx

						!if (sigma2a*sigma2b .gt. 0.0) then
						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
						if (sigma2 .eq. abs(sigma2a)) then
							sigma2 = sign(1.0D+00,sigma2a)*sigma2
						end if
						if (sigma2 .eq. abs(sigma2b)) then
							sigma2 = sign(1.0D+00,sigma2b)*sigma2
						end if
							!end if
						
						sigma4 = 0.0
						sigma4a = 2.0*(sol0(i,j) - sol0(i-1,j))/dx
						sigma4b = (sol0(i-1,j) - sol0(i-2,j))/dx

						!if (sigma4a*sigma4b .gt. 0.0) then
						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
						if (sigma4 .eq. abs(sigma4a)) then
							sigma4 = sign(1.0D+00,sigma4a)*sigma4
						end if
						if (sigma4 .eq. abs(sigma4b)) then
							sigma4 = sign(1.0D+00,sigma4b)*sigma4
						end if
							!end if
						
						
						! choosing sigma6
						sigma6 = 0.0
						!if (sigma2*sigma4 .gt. 0.0) then
							sigma6 = sign(1.0D+00,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
							!end if
						
! 						write(*,*) "sigma5"
! 						write(*,*) sigma5
! 						write(*,*) "sigma6"
! 						write(*,*) sigma6
						correction = (uTransport(i,j)*(mstep/ar)*0.5/dx) * (sigma5 - sigma6) * (dx - uTransport(i,j)*(mstep/ar))
						solute_next(i,j) = solute_next(i,j) - correction

					end if
				end if
				! end correction loop
				
			end if
			
			if (uTransport(i,j) .lt. 0.0) then
				! upwind including RHS value
				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i+1,j) - sol0(i,j) )
				
! 				! correction loop: sort of a mess
! 				if (i .lt. xn-1) then
! 					if (maskP(i+2,j) .ne. 0.0) then
! 						sigma1 = 0.0
! 						sigma1a = (sol0(i+2,j) - sol0(i+1,j))/dx
! 						sigma1b = (sol0(i+1,j) - sol0(i,j))/dx
!
! 						if (sigma1a*sigma1b .gt. 0.0) then
! 							sigma1 = sign(1.0D+00,sigma1a)*minval((/abs(sigma1a), abs(sigma1b)/))
! 						end if
!
! 						sigma2 = 0.0
! 						sigma2a = (sol0(i+1,j) - sol0(i,j))/dx
! 						sigma2b = (sol0(i,j) - sol0(i-1,j))/dx
!
! 						if (sigma2a*sigma2b .gt. 0.0) then
! 							sigma2 = sign(1.0D+00,sigma2a)*minval((/abs(sigma2a), abs(sigma2b)/))
! 						end if
!
! 						correction = (uTransport(i,j)*(dt*mstep/cstep)*0.5/dx) * (sigma1 - sigma2) * (dx - uTransport(i,j)*(dt*mstep/cstep))
! 						solute_next(i,j) = solute_next(i,j) - correction
! 					end if
! 				end if

				! correction loop: sort of a mess
				if (i .lt. xn-1) then
					if (maskP(i+2,j) .ne. 0.0) then
						sigma1 = 0.0
						sigma1a = (sol0(i+2,j) - sol0(i+1,j))/dx
						sigma1b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
						
						!if (sigma1a*sigma1b .gt. 0.0) then
							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
							if (sigma1 .eq. abs(sigma1a)) then
								sigma1 = sign(1.0D+00,sigma1a)*sigma1
							end if
							if (sigma1 .eq. abs(sigma1b)) then
								sigma1 = sign(1.0D+00,sigma1b)*sigma1
							end if
							!end if
						
						sigma3 = 0.0
						sigma3a = 2.0*(sol0(i+2,j) - sol0(i+1,j))/dx
						sigma3b = (sol0(i+1,j) - sol0(i,j))/dx
						
						!if (sigma3a*sigma3b .gt. 0.0) then
						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
						if (sigma3 .eq. abs(sigma3a)) then
							sigma3 = sign(1.0D+00,sigma3a)*sigma3
						end if
						if (sigma3 .eq. abs(sigma3b)) then
							sigma3 = sign(1.0D+00,sigma3b)*sigma3
						end if
							!end if
						
						
						! choosing sigma5
						sigma5 = 0.0
						if (sigma1*sigma3 .gt. 0.0) then
							sigma5 = sign(1.0D+00,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
						end if



						

						sigma2 = 0.0
						sigma2a = (sol0(i+1,j) - sol0(i,j))/dx
						sigma2b = 2.0*(sol0(i+1,j) - sol0(i,j))/dx

						!if (sigma2a*sigma2b .gt. 0.0) then
						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
						if (sigma2 .eq. abs(sigma2a)) then
							sigma2 = sign(1.0D+00,sigma2a)*sigma2
						end if
						if (sigma2 .eq. abs(sigma2b)) then
							sigma2 = sign(1.0D+00,sigma2b)*sigma2
						end if
							!end if
						
						sigma4 = 0.0
						sigma4a = 2.0*(sol0(i+1,j) - sol0(i,j))/dx
						sigma4b = (sol0(i+1,j) - sol0(i,j))/dx

						!if (sigma4a*sigma4b .gt. 0.0) then
						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
						if (sigma4 .eq. abs(sigma4a)) then
							sigma4 = sign(1.0D+00,sigma4a)*sigma4
						end if
						if (sigma4 .eq. abs(sigma4b)) then
							sigma4 = sign(1.0D+00,sigma4b)*sigma4
						end if
							!end if
						
						
						! choosing sigma6
						sigma6 = 0.0
						if (sigma2*sigma4 .gt. 0.0) then
							sigma6 = sign(1.0D+00,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
						end if
						
! 						write(*,*) "sigma5"
! 						write(*,*) sigma5
! 						write(*,*) "sigma6"
! 						write(*,*) sigma6
						correction = (uTransport(i,j)*(mstep/ar)*0.5/dx) * (sigma6 - sigma5) * (dx - uTransport(i,j)*(mstep/ar))
						solute_next(i,j) = solute_next(i,j) - correction
					end if
				end if
				! end correction loop
				
			end if
				
			!end if

! 		!if ( (solute_next(i,j).lt.sol0(i,j)) .and. (solute_next(i,j).lt.sol0(i+1,j)) .and. (solute_next(i,j).lt.sol0(i-1,j)) ) then
! 			if (uTransport(i,j) .gt. 0.0) then
! 				! upwind including LHS value
! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i,j) - sol0(i-1,j) )
! 			end if
!
! 			if (uTransport(i,j) .lt. 0.0) then
! 				! upwind including RHS value
! 				solute_next(i,j) = sol0(i,j) - qx*uTransport(i,j)*( sol0(i+1,j) - sol0(i,j) )
! 			end if
! 			!end if
		
	end do
end do


sol = solute_next



!
! bm1 = (qy*vTransport/2.0)*( (qy*vTransport) + 1.0)
! b0 = 1.0 - ( (qy*vTransport) * (qy*vTransport) )
! bp1 = (qy*vTransport/2.0)*( (qy*vTransport) - 1.0)
!

do i = 1,xn
	do j = yn/2,yn
		if ((maskP(i,j) .eq. 25.0) .or. (maskP(i,j) .eq. 51.0) .or. (maskP(i,j) .eq. 120.5) .or. (maskP(i,j) .eq. 170.5) ) then

			if (vTransport(i,j+1) .lt. 0.0) then
				sol(i,j+1) = seaw
			else
				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
			end if

		end if

		if ((maskP(i,j) .eq. 50.0)) then

! 			if (vTransport(i,j+1) .lt. 0.0) then
				sol(i,j+1) = seaw
! 			else
! 				sol(i,j+1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j-1)
! 			end if

		end if
		
		
		
		if (maskP(i,j) .eq. 100.0) then
			
				sol(i,j-1) = (4.0/3.0)*sol(i,j) - (1.0/3.0)*sol(i,j+1)

		end if
		
		
	end do
end do

! do i = 1,xn
! 	do j = yn/2,yn-1
! 		if (maskP(i,j) .ge. 1.0) then
! 			solute_next(i,j) = bm1(i,j-1)*sol(i,j-1) + b0(i,j)*sol(i,j) + bp1(i,j+1)*sol(i,j+1)
! 		end if
! 	end do
! end do

sol0 = sol


do i = 1,xn
	do j = yn/2,yn-1
		!if ( (solute_next(i,j).gt.sol0(i,j)) .and. (solute_next(i,j).gt.sol0(i,j+1)) .and. (solute_next(i,j).gt.sol0(i,j-1)) ) then
			

				
			if (vTransport(i,j) .gt. 0.0) then
				! upwind including bottom value
				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j) - sol0(i,j-1) )
!
			! correction loop: sort of a mess
				if (maskP(i,j-2) .ne. 0.0) then
						sigma1 = 0.0
						sigma1a = (sol0(i,j+1) - sol0(i,j))/dy
						sigma1b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
						
						!if (sigma1a*sigma1b .gt. 0.0) then
							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
							if (sigma1 .eq. abs(sigma1a)) then
								sigma1 = sign(1.0D+00,sigma1a)*sigma1
							end if
							if (sigma1 .eq. abs(sigma1b)) then
								sigma1 = sign(1.0D+00,sigma1b)*sigma1
							end if
							!end if
						
						sigma3 = 0.0
						sigma3a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
						sigma3b = (sol0(i,j) - sol0(i,j-1))/dy
						
						!if (sigma3a*sigma3b .gt. 0.0) then
						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
						if (sigma3 .eq. abs(sigma3a)) then
							sigma3 = sign(1.0D+00,sigma3a)*sigma3
						end if
						if (sigma3 .eq. abs(sigma3b)) then
							sigma3 = sign(1.0D+00,sigma3b)*sigma3
						end if
							!end if
						
						
						! choosing sigma5
						sigma5 = 0.0
						if (sigma1*sigma3 .gt. 0.0) then
							sigma5 = sign(1.0D+00,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
						end if



						

						sigma2 = 0.0
						sigma2a = (sol0(i,j) - sol0(i,j-1))/dy
						sigma2b = 2.0*(sol0(i,j-1) - sol0(i,j-2))/dy

						!if (sigma2a*sigma2b .gt. 0.0) then
						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
						if (sigma2 .eq. abs(sigma2a)) then
							sigma2 = sign(1.0D+00,sigma2a)*sigma2
						end if
						if (sigma2 .eq. abs(sigma2b)) then
							sigma2 = sign(1.0D+00,sigma2b)*sigma2
						end if
							!end if
						
						sigma4 = 0.0
						sigma4a = 2.0*(sol0(i,j) - sol0(i,j-1))/dy
						sigma4b = (sol0(i,j-1) - sol0(i,j-2))/dy

						!if (sigma4a*sigma4b .gt. 0.0) then
						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
						if (sigma4 .eq. abs(sigma4a)) then
							sigma4 = sign(1.0D+00,sigma4a)*sigma4
						end if
						if (sigma4 .eq. abs(sigma4b)) then
							sigma4 = sign(1.0D+00,sigma4b)*sigma4
						end if
							!end if
						
						
						! choosing sigma6
						sigma6 = 0.0
						if (sigma2*sigma4 .gt. 0.0) then
							sigma6 = sign(1.0D+00,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
						end if
						
! 						write(*,*) "sigma5"
! 						write(*,*) sigma5
! 						write(*,*) "sigma6"
! 						write(*,*) sigma6
						correction = (vTransport(i,j)*(mstep/ar)*0.5/dy) * (sigma5 - sigma6) * (dy - vTransport(i,j)*(mstep/ar))
						solute_next(i,j) = solute_next(i,j) - correction
				end if
				! end correction loop
				
				
			end if
			
			if (vTransport(i,j) .lt. 0.0) then
				! upwind including top value
				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j+1) - sol0(i,j) )
				
				! correction loop: sort of a mess
				if (j .lt. yn-1) then
					if (maskP(i,j-2) .ne. 0.0) then
						sigma1 = 0.0
						sigma1a = (sol0(i,j+2) - sol0(i,j+1))/dy
						sigma1b = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
						
						!if (sigma1a*sigma1b .gt. 0.0) then
							sigma1 = minval((/abs(sigma1a), abs(sigma1b)/))
							if (sigma1 .eq. abs(sigma1a)) then
								sigma1 = sign(1.0D+00,sigma1a)*sigma1
							end if
							if (sigma1 .eq. abs(sigma1b)) then
								sigma1 = sign(1.0D+00,sigma1b)*sigma1
							end if
							!end if
						
						sigma3 = 0.0
						sigma3a = 2.0*(sol0(i,j+2) - sol0(i,j+1))/dy
						sigma3b = (sol0(i,j+1) - sol0(i,j))/dy
						
						!if (sigma3a*sigma3b .gt. 0.0) then
						sigma3 = minval((/abs(sigma3a), abs(sigma3b)/))
						if (sigma3 .eq. abs(sigma3a)) then
							sigma3 = sign(1.0D+00,sigma3a)*sigma3
						end if
						if (sigma3 .eq. abs(sigma3b)) then
							sigma3 = sign(1.0D+00,sigma3b)*sigma3
						end if
							!end if
						
						
						! choosing sigma5
						sigma5 = 0.0
						if (sigma1*sigma3 .gt. 0.0) then
							sigma5 = sign(1.0D+00,sigma1)*maxval((/abs(sigma1), abs(sigma3)/))
						end if



						

						sigma2 = 0.0
						sigma2a = (sol0(i,j+1) - sol0(i,j))/dy
						sigma2b = 2.0*(sol0(i,j) - sol0(i,j-1))/dy

						!if (sigma2a*sigma2b .gt. 0.0) then
						sigma2 = minval((/abs(sigma2a), abs(sigma2b)/))
						if (sigma2 .eq. abs(sigma2a)) then
							sigma2 = sign(1.0D+00,sigma2a)*sigma2
						end if
						if (sigma2 .eq. abs(sigma2b)) then
							sigma2 = sign(1.0D+00,sigma2b)*sigma2
						end if
							!end if
						
						sigma4 = 0.0
						sigma4a = 2.0*(sol0(i,j+1) - sol0(i,j))/dy
						sigma4b = (sol0(i,j) - sol0(i,j-1))/dy

						!if (sigma4a*sigma4b .gt. 0.0) then
						sigma4 = minval((/abs(sigma4a), abs(sigma4b)/))
						if (sigma4 .eq. abs(sigma4a)) then
							sigma4 = sign(1.0D+00,sigma4a)*sigma4
						end if
						if (sigma4 .eq. abs(sigma4b)) then
							sigma4 = sign(1.0D+00,sigma4b)*sigma4
						end if
							!end if
						
						
						! choosing sigma6
						sigma6 = 0.0
						if (sigma2*sigma4 .gt. 0.0) then
							sigma6 = sign(1.0D+00,sigma2)*maxval((/abs(sigma2), abs(sigma4)/))
						end if
						
! 						write(*,*) "sigma5"
! 						write(*,*) sigma5
! 						write(*,*) "sigma6"
! 						write(*,*) sigma6
						correction = (vTransport(i,j)*(mstep/ar)*0.5/dy) * (sigma6 - sigma5) * (dy - vTransport(i,j)*(mstep/ar))
						solute_next(i,j) = solute_next(i,j) - correction

					end if
				end if
				! end correction loop
					
					
			end if
				
			!end if

! 		!if ( (solute_next(i,j).lt.sol0(i,j)) .and. (solute_next(i,j).lt.sol0(i,j+1)) .and. (solute_next(i,j).lt.sol0(i,j-1)) ) then
! 			if (vTransport(i,j) .gt. 0.0) then
! 				! upwind including bottom value
! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j) - sol0(i,j-1) )
! 			end if
!
! 			if (vTransport(i,j) .lt. 0.0) then
! 				! upwind including top value
! 				solute_next(i,j) = sol0(i,j) - qy*vTransport(i,j)*( sol0(i,j+1) - sol0(i,j) )
! 			end if
!
!
!
! 			!end if
	end do
end do




! do i = 2,xn-1
! 	do j = 2,yn-1
! 		if ((solute_next(i,j).gt.sol0(i,j)).and.(solute_next(i,j).gt.sol0(i+1,j)).and.(solute_next(i,j).gt.sol0(i-1,j)).and.(solute_next(i,j).gt.sol0(i,j-1)).and.(solute_next(i,j) .gt. sol0(i,j+1))) then
! 			solute_next(i,j) = sol0(i,j)
! 		end if
!
! 		if ((solute_next(i,j).lt.sol0(i,j)).and.(solute_next(i,j).lt.sol0(i+1,j)).and.(solute_next(i,j).lt.sol0(i-1,j)).and.(solute_next(i,j).lt.sol0(i,j-1)).and.(solute_next(i,j) .lt. sol0(i,j+1))) then
! 			solute_next(i,j) = sol0(i,j)
! 		end if
! 	end do
! end do

!
! do i = 1,xn
! 	do j = 1,yn
! 		if ((maskP(i,j) .eq. 0.0)) then
! 			solute_next(i,j) = seaw
! 		end if
! 	end do
! end do

! ! do left, right first i guess?
! solute_next(2:xn-1,:) = solute_next(2:xn-1,:)*bm1(2:xn-1,:) + solute_next(2:xn-1,:)*b0(2:xn-1,:) + solute_next(2:xn-1,:)*bp1(2:xn-1,:)

! put left, right, 10, and 5 bcs back in, i think

! ! ADVECTION STEP
! !
! !
! do ii=1,yn
! 	if ((uTransport(2,ii) .gt. 0.0) .and. (maskP(1,ii) .eq. 1.0)) then
! 		sol(2,i) = sol(2,ii) + uTransport(2,ii)*sol(1,ii)*qx
! 	end if
! 	if ((uTransport(xn-1,ii) .lt. 0.0) .and. (maskP(xn,ii) .eq. 1.0)) then
! 		sol(xn-1,ii) = sol(xn-1,ii) - uTransport(xn-1,ii)*sol(xn,ii)*qx
! 	end if
! end do
!
! do i = 1, xn
! 	do ii = 1,yn
! 		if ((maskP(i,ii) .eq. 5.0) .or. (maskP(i,ii) .eq. 12.5)) then
! 			if (uTransport(i,ii) .lt. 0.0) then
! 				sol(i,ii) = sol(i,ii) - seaw*uTransport(i,ii)*qx
! 			end if
! 		end if
!
! 		if ((maskP(i,ii) .eq. 10.0) .or. (maskP(i,ii) .eq. 17.5)) then
! 			if (uTransport(i,ii) .ge. 0.0) then
! 				sol(i,ii) = sol(i,ii) + seaw*uTransport(i,ii)*qx
! 			end if
! 		end if
! 	end do
! end do
!
!
! vec = reshape(sol(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
!
! ! MAKE THE BAND
! aBand = 0.0
!
!
! aBand(:,1) =  - uLong*qx/2.0
! aBand(:,2) = 1.0
! aBand(:,3) = uLong*qx/2.0
!
! aBand(1,1) = 0.0
! aBand(1,2) = 1.0 - uLong((xn-2)*(yn-0))*qx
! aBand(1,3) = uLong((xn-2)*(yn-0))*qx
!
! aBand((xn-2)*(yn-2),2) = 1.0 + uLong((xn-2)*(yn-0))*qx
! aBand((xn-2)*(yn-2),1) = - uLong((xn-2)*(yn-0))*qx
! aBand((xn-2)*(yn-2),3) =  0.0
!
!
! do i = 2,(xn-2)*(yn-0)-1
!
! 	! flow left
! 	if ((uLong(i) .lt. 0.0)) then
! 		aBand(i,1) =  0.0
! 		aBand(i,2) = 1.0 - uLong(i)*qx
! 		aBand(i,3) = uLong(i)*qx
! 	end if
!
! 	! flow right
! 	if ((uLong(i) .ge. 0.0)) then
! 		aBand(i,1) = - uLong(i)*qx
! 		aBand(i,2) = 1.0 + uLong(i)*qx
! 		aBand(i,3) =  0.0
! 	end if
!
! 	!left edge
! 	if ((maskLongU(i) .eq. 10.0)) then
! 		aBand(i,1) =  0.0
! 	end if
! 	! right edge
! 	if ((maskLongU(i) .eq. 5.0)) then
! 		aBand(i,3) =  0.0
! 	end if
!
!
!
! 	!left edge
! 	if ((mod(i-1,xn-2) .eq. 0)) then
! 		aBand(i,1) =  0.0
! 	end if
! 	! right edge
! 	if ((mod(i,xn-2) .eq. 0)) then
! 		aBand(i,3) =  0.0
! 	end if
!
! if ((i .gt. (xn-2)) .and. (i .lt. yn*(xn-2)-(xn-2))) then
! 	if ((maskLongU(i) .eq. 0.0) .and. (maskLongU(i-(xn-2)) .eq. 0.0) .and. (maskLongU(i+(xn-2)) .eq. 0.0)) then
! 		aBand(i,1) = 0.0
! 		aBand(i,2) = 1.0
! 		aBand(i,3) = 0.0
! 	end if
! end if
!
!
!
! end do
!
!
!
! !!!!!!!!!!!! THIS !!!!!!!!!!!
! sol_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),vec,(xn-2)*(yn-0))
! sol(2:xn-1,1:yn) = reshape(sol_nextRow, (/xn-2, yn-0/))
!
! solute_next(2:xn-1,1:yn) = reshape(sol_nextRow, (/xn-2, yn-0/))
! sol(2:xn-1,1:yn) = solute_next(2:xn-1,1:yn)
!
!
! ! top bottom boundary conditions
! do i = 1, xn
! 	do ii = 1,yn
! 		if ((mask(i,ii) .eq. 25.0) .or. (mask(i,ii) .eq. 50.0) .or. (mask(i,ii) .eq. 12.5) .or. (mask(i,ii) .eq. 17.5)) then
! 			if (vTransport(i,ii) .lt. 0.0) then
! 				sol(i,ii) = sol(i,ii) - vTransport(i,ii)*seaw*qy
! 		 	end if
! 		end if
! 	end do
! end do
!
!
! sol_nextRowB = reshape(transpose(sol(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))
!
!
!
! ! MAKE THE BAND
! bBand = 0.0
!
! bBand(:,1) =  - vLong*qy/2.0
! bBand(:,2) = 1.0
! bBand(:,3) = vLong*qy/2.0
!
!
! bBand(1,1) = 0.0
! bBand(1,2) = 1.0 - vLong((xn-0)*(yn-2))*qy
! bBand(1,3) = vLong((xn-0)*(yn-2))*qy
!
! bBand((xn-0)*(yn-2),2) = 1.0 + vLong((xn-0)*(yn-2))*qy
! bBand((xn-0)*(yn-2),1) = - vLong((xn-0)*(yn-2))*qy
! bBand((xn-0)*(yn-2),3) =  0.0
!
!
! do i = 2,(xn-0)*(yn-2)-1
!
! 	! flow down
! 	if ((vLong(i) .lt. 0.0)) then
! 		bBand(i,2) = 1.0 - vLong(i)*qy
! 		bBand(i,1) =  0.0
! 		bBand(i,3) = vLong(i)*qy
! 	end if
!
! 	! flow up
! 	if ((vLong(i) .ge. 0.0)) then
! 		bBand(i,2) = 1.0 + vLong(i)*qy
! 		bBand(i,1) = - vLong(i)*qy
! 		bBand(i,3) =  0.0
! 	end if
!
! 	! top edge, flow down
! 	if ((maskLongTV(i) .eq. 25.0) .and. (vLong(i) .lt. 0.0)) then
! 		bBand(i,3) = 0.0
! 	end if
!
! 	! top edge, flow down
! 	if ((maskLongTV(i) .eq. 50.0) .and. (vLong(i) .lt. 0.0)) then
! 		bBand(i,3) = 0.0
! 	end if
!
! 	! top edge, flow down
! 	if ((maskLongTV(i) .eq. 12.5) .and. (vLong(i) .lt. 0.0)) then
! 		bBand(i,3) = 0.0
! 	end if
!
! 	! top edge, flow down
! 	if ((maskLongTV(i) .eq. 17.5) .and. (vLong(i) .lt. 0.0)) then
! 		bBand(i,3) = 0.0
! 	end if
!
! 	! bottom edge, has to flow up...
! 	if ( (maskLongTV(i) .eq. 100.0) .and. (vLong(i) .gt. 0.0)) then
! 		bBand(i,1) = 0.0
! 	end if
!
! 	if (maskLongTV(i) .eq. 0.0) then
! 		bBand(i,2) = 1.0
! 		bBand(i,1) = 0.0
! 		bBand(i,3) = 0.0
! 	end if
!
!
! end do
!
!
! sol_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),sol_nextRowB,(xn-0)*(yn-2))
! solute_next(1:xn,2:yn-1) = transpose(reshape(sol_nextRow, (/(yn-2), (xn-0)/)))
!
! ! do i = 1, xn
! ! 	do ii = 1,yn
! ! 		if (solute_next(i,ii) .gt. 1.0) then
! ! 			solute_next(i,ii) = 1.0
! ! 		end if
! !
! ! 	end do
! ! end do
!
! do i = 1, xn
! do ii = 1,yn
! 	if ((maskP(i,ii) .eq. 25.0) .and. (uTransport(i,ii) .le. 0.0)) then
!
! 	end if
! end do
! end do
! !solute_next = 0.0

return

end function solute_next





! ! ----------------------------------------------------------------------------------%%
! !
! ! ALT_NEXT
! !
! ! SUMMARY: solves for equilibrium at a single grid cell using PHREEQC
! !
! ! INPUTS: temp : temperature of grid cell
! !         timestep : time elapsed
! !         primaryList(5) : amounts of primary minerals
! !         secondaryList(16) : amounts of secondary minerals
! !         soluteList(11) : concentrations of solutes
! !
! ! RETURNS: alt_next(1,altnum): returns everything from PHREEQC in a big pile
! !          and it gets parsed in the main method's geochem loop
! !
! ! ----------------------------------------------------------------------------------%%
!
! function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
! use globals
! use initialize
! use alteration
! implicit none
!
! interface
!
! end interface
!
! ! declare errthing
!
! integer :: order
! real(8) :: temp, timestep
! real(8) :: alt_next(1,167)
! real(8) :: alter0(1,167)
! real(8) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)
!
! ! use the alteration module
! alter0 = alter(temp-272.9, timestep, primaryList, secondaryList, soluteList, mediumList)
!
! ! rename it for a reason that i now forget
! alt_next = alter0
!
! end function alt_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

function rho_next(h_in)
	
use globals
use initialize
implicit none
  
! declare errthing
integer :: i,j
real(8) :: h_in(xn,yn), rho_next(xn,yn)


do i=1,xn
	do j = 1,yn
		rho_next(i,j) = rho_fluid*(1.0 - alpha*(h_in(i,j)-273.0))
	end do 
end do

return

end function rho_next



! ----------------------------------------------------------------------------------%%
!
! RHO_ONE
!
!
! ----------------------------------------------------------------------------------%%

function rho_one(h_in)
	
use globals
use initialize
implicit none
  
! declare errthing
integer :: i,j
real(8) :: h_in, rho_one


		rho_one = rho_fluid*(1.0 - alpha*((h_in-273.0)))


return

end function rho_one


! ----------------------------------------------------------------------------------%%
!
! VISC_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

function visc_next(h_in)
	
use globals
use initialize
implicit none
  
! declare errthing
integer :: i,j
real(8) :: h_in(xn,yn), visc_next(xn,yn)


do i=1,xn
	do j = 1,yn
		!visc_next(i,j) = (2.44e-5)*10.0**(247.8/(h_in(i,j)-140.0))
		visc_next(i,j) = .001!.001!.0002 !.00350*exp(-(h_in(i,j)-273.16)/35.0) + .0002
		! from nist.gov
	end do 
end do

return

end function visc_next




! ----------------------------------------------------------------------------------%%
!
! VELOCITIES
!
! SUMMARY : computes the darcy velocity (specific discharge) from the streamfunction
!           using finite difference partial derivatives
!
! INPUTS : psi(xn,yn) : 2D streamfunction array of current timestep
!
! RETURNS : velocities(xn,2*yn) : both u and v velocities in one matrix
!
! ----------------------------------------------------------------------------------%%


function velocities(psi)
	
use globals
use initialize
implicit none

interface

	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge(rows,cols)
	end function partial_edge
	
	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial_edge_p(rows,cols)
	end function partial_edge_p

end interface

! declare errthing
integer :: i,ii
real(8) :: velocities(xn,2*yn), psi(xn,yn)
real(8) :: u0(xn,yn), v0(xn,yn)

u0 = partial_edge_p(psi,xn,yn,dx,dy,2)
v0 = -partial_edge_p(psi,xn,yn,dx,dy,1)

velocities(1:xn,1:yn) = u0
velocities(1:xn,yn+1:2*yn) = v0

return
end function velocities


! calculates velocities from COARSE streamfunction values
function velocitiesCoarse(psiCoarse)
	use globals
	use initialize
	implicit none
	integer :: i,ii
	real(8) :: velocitiesCoarse(xn/cell,2*yn/cell), psiCoarse(xn/cell,yn/cell)
	real(8) :: u0(xn/cell,yn/cell), v0(xn/cell,yn/cell)

interface

	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial

end interface

u0 = partial(psiCoarse,xn/cell,yn/cell,dx*cell,dy*cell,2)
v0 = -partial(psiCoarse,xn/cell,yn/cell,dx*cell,dy*cell,1)

! do i =1,xn
! 	do ii = 1,yn
! 		if (mask(i,ii) .eq. 50.0) then
! 			u0(i,ii+1) = 0.0
! 		end if
! 	end do
! end do

velocitiesCoarse(1:xn/cell,1:yn/cell) = u0
velocitiesCoarse(1:xn/cell,yn/cell+1:2*yn/cell) = v0

return
end function velocitiesCoarse



! ----------------------------------------------------------------------------------%%
!
! PARTIAL
!
! SUMMARY : versatile function for solving for second-order accurate partial 
!           derivatives of 1D or 2D arrays with respect to specified dimension
!           
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%

function partial(array,rows,cols,d1,d2,dim)
	
use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols

! figure out which direction derivative goes (dx or dy)

partial = 0.0

if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1
	
	! compute edges beforehand
	partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
	partial(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)
	
	do i = 2,rows-1
		do j = 1,cols
			partial(i,j) = (array(i+1,j) - array(i-1,j))/(2.0*d)
		end do
	end do
	
	do i = 2,rows-1
		do j = 1,cols
			if ((maskP(i,j) .eq. 3.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 3.1) .or. (maskP(i,j) .eq. 3.05) .or. (maskP(i,j) .eq. 3.01) .or. (mask(i,j) .eq. 3.05)) then
				partial(i,j) = (array(i,j) - array(i-1,j))/d
				!partial(i,j) = (3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j))/(2.0*d)
			end if
			if ((maskP(i,j) .eq. 6.0) .or. (maskP(i,j) .eq. 6.5) .or. (maskP(i,j) .eq. 6.1) .or. (maskP(i,j) .eq. 6.05) .or. (maskP(i,j) .eq. 6.01) .or. (mask(i,j) .eq. 6.05)) then
				partial(i,j) = (array(i+1,j) - array(i,j))/d
				!partial(i,j) = (-3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j))/(2.0*d)
			end if
		end do
	end do
	
end if


if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
	
	! compute edges beforehand 
	partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
end if

! ! use central difference method ignoring edges (already done)
! do i = 2-jj,rows-1+jj
!     do j = 2-ii,cols-1+ii
!     	partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
! 	end do
! end do


! if (dim .eq. 1) then
! 	do i=2,rows-1
! 		do j=2,cols-1
!
! 			if (maskP(i,j) .eq. 5.0) then
! 				partial(i,j) = (array(i,j) - array(i-1,j))/d
! 				partial(i+1,j) = (array(i+1,j) - array(i,j))/d
! 			end if
!
! 			if (maskP(i,j) .eq. 10.0) then
! 				partial(i,j) = (array(i+1,j) - array(i,j))/d
! 				partial(i-1,j) = (array(i,j) - array(i-1,j))/d
! 			end if
!
!
! 		end do
! 	end do
!
!
! end if

! ! if jj=1, dim=2, cols=yn, rows=xn

! ! compute edges beforehand 
! partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
! partial(:,yn) = ( 3.0*array(:,yn) - 4.0*array(:,yn-1) + array(:,yn-2) ) / (2.0*d)


! do i = 1,xn
!     do j = 2,yn-1
!     	partial(i,j) = (array(i,j+1)-array(i,j-1))/(2.0*d)
! 	end do
! end do






! ! if ii=1, dim=1, cols=yn, rows=xn

! ! compute edges beforehand 
! partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
! partial(xn,:) = ( 3.0*array(xn,:) - 4.0*array(xn-1,:) + array(xn-2,:) ) / (2.0*d)


! do i = 2,xn-1
!     do j = 1,yn
!     	partial(i,j) = (array(i+1,j)-array(i-1,j))/(2.0*d)
! 	end do
! end do

return
end function partial







! ----------------------------------------------------------------------------------%%
!
! PARTIAL_EDGE
!
! SUMMARY : versatile function for solving for second-order accurate partial 
!           derivatives of 1D or 2D arrays with respect to specified dimension
!           
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%

function partial_edge(array,rows,cols,d1,d2,dim)
	
use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial_edge(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols


! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1
end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
end if
partial_edge = 0.0
! use central difference method ignoring edges (already done)
do i = 2-jj,rows-1+jj
    do j = 2-ii,cols-1+ii
		if ((mask(i,j) .ne. 0.0)) then
    		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
		end if
	end do
end do

! do i = 2,rows-1
!     do j = 2,cols-1
! 		if ((mask(i,j) .ne. 0.0) .or. (mask(i+jj,j) .ne. 0.0) .or. (mask(i-jj,j-ii) .ne. 0.0)) then
!     		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
! 		end if
! 	end do
! end do

do i = 2,rows-1
    do j = 2,cols-1
		if ((mask(i,j) .ne. 0.0) .or. (mask(i-1,j) .ne. 0.0) .or. (mask(i+1,j) .ne. 0.0) .or. (mask(i,j-1) .ne. 0.0) .or. (mask(i,j+1) .ne. 0.0) ) then
    		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
		end if
	end do
end do


! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	! compute edges beforehand
	partial_edge(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
	partial_edge(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)
	
	! inner edges
	do j=2,cols-1
		do i=2,rows-1
			if ((mask(i-1,j) .eq. 12.5)) then
				partial_edge(i,j) =( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
			end if
			if ((mask(i-1,j) .eq. 5.0) .and. (mask(i-1,j) .eq. 5.0)) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d) 
			end if
			
			if ((mask(i+1,j) .eq. 17.5)) then
				partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)
			end if
			if ((mask(i+1,j) .eq. 10.0) .and. (mask(i+1,j) .eq. 10.0)) then
				partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d) 
			end if
			
			
			if (mask(i,j-1) .eq. 12.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
			end if 
			
			if (mask(i,j-1) .eq. 17.5) then
				partial_edge(i,j) =( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d) 
			end if 


		end do
	end do
	
	
	do j=1,cols
		do i=2,rows-1
			
			if (mask(i-1,j-1) .eq. 12.5) then
				partial_edge(i,j) = partial_edge(i-1,j)
			end if

			if (mask(i+1,j-1) .eq. 17.5) then
				partial_edge(i,j) = partial_edge(i+1,j)
			end if
		
		end do
	end do
	
	
end if




if (dim .eq. 2) then
	! compute edges beforehand 
	partial_edge(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial_edge(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
	
	
	! inner edges
	do j=2,cols
		do i=1,rows
			if ((mask(i,j-1).eq.25.0) .or. (mask(i,j-1).eq.12.5)  .or. (mask(i,j-1).eq.17.5)) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if
		end do
	end do
	
	do j=2,cols-1
		do i=2,rows-1
			if ((mask(i,j-1).eq.50.0) .and. (mask(i,j-1).eq.50.0) .and. (mask(i,j-1).eq.50.0)) then
				partial_edge(i,j) =  ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)!
			end if
		end do
	end do
	


	! JANUARY DELETION...
	do j=2,cols
		do i=2,rows-1

			if (mask(i-1,j) .eq. 12.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if

			if (mask(i+1,j) .eq. 17.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if

		end do
	end do
	
	do j=2,cols
		do i=2,rows-1
 
			if ((mask(i-1,j-1) .eq. 12.5)) then
				partial_edge(i,j) = partial_edge(i,j-1)
			end if

			if ((mask(i+1,j-1) .eq. 17.5)) then
				partial_edge(i,j) = partial_edge(i,j-1)
			end if
			

		end do
	end do
	

	
end if





return
end function partial_edge










! ----------------------------------------------------------------------------------%%
!
! partial_edge_p
!
! SUMMARY : versatile function for solving for second-order accurate partial 
!           derivatives of 1D or 2D arrays with respect to specified dimension
!           
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%


function partial_edge_p(array,rows,cols,d1,d2,dim)
	
use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial_edge_p(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols

partial_edge_p = 0.0

! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1
	
	do i=2,rows-1
		do j=1,cols
			if (maskP(i,j) .ne. 0.0) then
				partial_edge_p(i,j) = (array(i+1,j) - array(i-1,j)) / (2.0*d)
			end if
		end do
	end do
	
	do i=2,rows-1
		do j=1,cols
			if ((maskP(i,j) .eq. 5.0) .or. (maskP(i,j) .eq. 2.5)) then
				partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / (1.0*d)
			end if
		end do
	end do
	
	do i=1,rows-1
		do j=1,cols
! 			if (maskP(i,j) .eq. 6.0) then
! 				partial_edge_p(i,j) = 0.0
! 			end if
			if ((maskP(i,j) .eq. 6.0) .or. (maskP(i,j) .eq. 6.5) .or. (maskP(i,j) .eq. 6.1) .or. (maskP(i,j) .eq. 6.05) .or. (maskP(i,j) .eq. 6.01)) then
				partial_edge_p(i,j) = (array(i+1,j) - array(i,j)) / d
			end if
			if ((maskP(i,j) .eq. 3.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 3.1) .or. (maskP(i,j) .eq. 3.05) .or. (maskP(i,j) .eq. 3.01)) then
				partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / d
			end if
		end do
	end do
!
	! fracture vertical velocity
	!partial_edge_p(xn-1,:) = (array(xn,:) - array(xn-1,:)) / param_f_dx

	
end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
	
	do i=1,rows
		do j=2,cols-1
			if ((maskP(i,j) .ne. 0.0)) then
				partial_edge_p(i,j) = (array(i,j+1) - array(i,j-1)) / (2.0*d)
			end if
		end do
	end do
	
	do i=1,rows
		do j=2,cols-1
			if ((maskP(i,j) .eq. 50.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 6.5)) then
				partial_edge_p(i,j) = (array(i,j) - array(i,j-1)) / (1.0*d)
			end if
		end do
	end do
	


end if




return
end function partial_edge_p




! ----------------------------------------------------------------------------------%%
!
! WRITE_VEC
!
! SUMMARY: Write linspace-style vector to file
!
! INPUTS: n : dimension
!         vector : vector with data
!         filename : file name
!
! RETURNS: write_vec
!
! ----------------------------------------------------------------------------------%%

function write_vec ( n, vector, filename )
use globals
  implicit none
  integer :: n, j, output_status, unit0
  character ( len = * ) filename 
  real(4)  :: vector(n), write_vec



  unit0 = get_unit ()
  open ( unit = unit0, file = filename, status = 'replace', iostat = output_status )
  if ( output_status /= 0 ) then
    write ( *, '(a,i8)' ) 'COULD NOT OPEN OUTPUT FILE "' // &
      trim ( filename ) // '" USING UNIT ', unit0
    unit0 = -1
    stop
  end if
  

  if ( 0 < n ) then
    do j = 1, n
      write ( unit0, '(2x,g24.16)' ) vector(j)
    end do

  end if


  close ( unit = unit0 )
  write_vec = 1.0
  return
end function write_vec




! ----------------------------------------------------------------------------------%%
!
! WRITE_MATRIX
!
! SUMMARY: Write 2d array to file
!
! INPUTS: m,n : 2d dimensinons
!         table : 2d array with data
!         filename : file name
!
! RETURNS: write_matrix
!
! ----------------------------------------------------------------------------------%%

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0, reclen
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix

  
  
  INQUIRE(iolength=reclen)table
  unit0 = get_unit ()
  open ( unit = unit0, file = filename, &
    status = 'replace', iostat = output_status, buffered='YES', buffercount=500)

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) 'Could not open the output file "' // &
      trim ( filename ) // '" on unit ', unit0
    unit0 = -1
    stop
  end if
  
  

!
 	write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
! 	!write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!     do j = 1, n
!       write ( unit0, string) table(1:m,j)
!     end do

    do j = 1, n
      write ( unit0, 400) table(1:m,j)
    end do
400 FORMAT(<m>g24.16)


  close ( unit = unit0 )
  write_matrix = 2.0
  return
end function write_matrix