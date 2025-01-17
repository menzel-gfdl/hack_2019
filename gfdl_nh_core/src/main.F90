program main
use iso_fortran_env, only: output_unit, error_unit
use netcdf
use nh_core_mod, only: Riem_Solver_c, allocate_data
implicit none

!Namelist parameters.
integer :: num_iter = 1000
character(len=128) :: riem_solver_c_file = "Riem_Solver_c_data.nc"

namelist /nh_core_nml/ num_iter

integer :: nx=0, ny=0, nz=0 
integer :: ng = 3
integer :: ms
real :: dt, akap, cp, ptop, scale_m, p_fac, a_imp

integer :: npx, npy
integer :: isc, iec, jsc, jec
integer :: isd, ied, jsd, jed


integer :: i, j, k, iter


real, allocatable, dimension(:,:) :: hs, ws
real, allocatable, dimension(:,:,:) :: cappa, w3, pt, q_con, delp, gz, pef, gz_out, pef_out

real, allocatable, dimension(:,:) :: hs_old, ws_old
real, allocatable, dimension(:,:,:) :: cappa_old, w3_old, pt_old, q_con_old, delp_old, gz_old



logical :: fexist, opened
integer :: unit, io_status
integer :: info, myrank
integer :: start_time, end_time, cr, cm
real    :: rate, total_time
integer :: ncid, did_axis, vid
real :: gz_diff !< Difference between actual value and benchmark.
real :: pef_diff !< Difference between actual value and benchmark.
real, parameter :: tolerance = 1.e-8

call system_clock(count_rate=cr)
call system_clock(count_max=cm)
rate = real(cr)
inquire(file="input.nml", exist=fexist)
if (fexist) then
  do unit = 103, 512
    inquire(unit, opened=opened)
    if (.not. opened) then
      exit
    endif
  enddo
  open(unit, file="input.nml", iostat=io_status)
  if (io_status .gt. 0) then
    call error("failed at open input.nml")
  endif
  read(unit, nh_core_nml, iostat=io_status)
  if (io_status .gt. 0) then
    call error("failed at read nh_core_nml.")
  endif
  close(unit)
endif


  !--- check if input file exist. If exist, read grid size from the file.
inquire(file=trim(riem_solver_c_file), exist=fexist)
if (.not. fexist) then
  call error("file "//trim(riem_solver_c_file)// " does not exist.")
endif
write(output_unit, *) "reading riem_solver_c data from file "//trim(riem_solver_c_file)
call handle_error(nf90_open(trim(riem_solver_c_file), nf90_nowrite, ncid))
call handle_error(nf90_inq_varid(ncid, 'npx', vid))
call handle_error(nf90_get_var(ncid, vid, npx))
call handle_error(nf90_inq_varid(ncid, 'npy', vid))
call handle_error(nf90_get_var(ncid, vid, npy))
call handle_error(nf90_inq_varid(ncid, 'km', vid))
call handle_error(nf90_get_var(ncid, vid, nz))
call handle_error(nf90_inq_varid(ncid, 'ng', vid))
call handle_error(nf90_get_var(ncid, vid, ng))

nx = npx-1; ny = npy-1
write(output_unit, *) "nx, ny, nz", nx, ny, nz, ng


isc = 1
iec = nx
jsc = 1
jec = ny
isd = isc - ng
ied = iec + ng
jsd = jsc - ng
jed = jec + ng

allocate(cappa(isd:ied,jsd:jed,nz))
allocate(hs(isd:ied,jsd:jed))
allocate(w3(isd:ied,jsd:jed,nz))
allocate(pt(isd:ied,jsd:jed,nz))
allocate(q_con(isd:ied,jsd:jed,nz))
allocate(delp(isd:ied,jsd:jed,nz))
allocate(gz(isd:ied,jsd:jed,nz+1))
allocate(ws(isd:ied,jsd:jed))
allocate(gz_out(isd:ied,jsd:jed,nz+1))
allocate(pef(isd:ied,jsd:jed,nz+1))
allocate(pef_out(isd:ied,jsd:jed,nz+1))

allocate(cappa_old(isd:ied,jsd:jed,nz))
allocate(hs_old(isd:ied,jsd:jed))
allocate(w3_old(isd:ied,jsd:jed,nz))
allocate(pt_old(isd:ied,jsd:jed,nz))
allocate(q_con_old(isd:ied,jsd:jed,nz))
allocate(delp_old(isd:ied,jsd:jed,nz))
allocate(gz_old(isd:ied,jsd:jed,nz+1))
allocate(ws_old(isd:ied,jsd:jed))

call allocate_data(isc,iec,jsc,jec,nz)


pef = 1.e8
call handle_error(nf90_inq_varid(ncid, 'ms', vid))
call handle_error(nf90_get_var(ncid, vid, ms))
call handle_error(nf90_inq_varid(ncid, 'dt', vid))
call handle_error(nf90_get_var(ncid, vid, dt))
call handle_error(nf90_inq_varid(ncid, 'akap', vid))
call handle_error(nf90_get_var(ncid, vid, akap))
!call handle_error(nf90_inq_varid(ncid, 'cappa', vid))
!call handle_error(nf90_get_var(ncid, vid, cappa))
call handle_error(nf90_inq_varid(ncid, 'cp', vid))
call handle_error(nf90_get_var(ncid, vid, cp))
call handle_error(nf90_inq_varid(ncid, 'ptop', vid))
call handle_error(nf90_get_var(ncid, vid, ptop))
call handle_error(nf90_inq_varid(ncid, 'hs', vid))
call handle_error(nf90_get_var(ncid, vid, hs))
call handle_error(nf90_inq_varid(ncid, 'w3', vid))
call handle_error(nf90_get_var(ncid, vid, w3))
call handle_error(nf90_inq_varid(ncid, 'pt', vid))
call handle_error(nf90_get_var(ncid, vid, pt))
!call handle_error(nf90_inq_varid(ncid, 'q_con', vid))
!call handle_error(nf90_get_var(ncid, vid, q_con))
call handle_error(nf90_inq_varid(ncid, 'delp', vid))
call handle_error(nf90_get_var(ncid, vid, delp))
call handle_error(nf90_inq_varid(ncid, 'gz', vid))
call handle_error(nf90_get_var(ncid, vid, gz))
call handle_error(nf90_inq_varid(ncid, 'ws', vid))
call handle_error(nf90_get_var(ncid, vid, ws))
call handle_error(nf90_inq_varid(ncid, 'p_fac', vid))
call handle_error(nf90_get_var(ncid, vid, p_fac))
call handle_error(nf90_inq_varid(ncid, 'a_imp', vid))
call handle_error(nf90_get_var(ncid, vid, a_imp))
call handle_error(nf90_inq_varid(ncid, 'scale_m', vid))
call handle_error(nf90_get_var(ncid, vid, scale_m))
call handle_error(nf90_inq_varid(ncid, 'gz_out', vid))
call handle_error(nf90_get_var(ncid, vid, gz_out))
call handle_error(nf90_inq_varid(ncid, 'pef', vid))
call handle_error(nf90_get_var(ncid, vid, pef_out))


!cappa_old = cappa
hs_old = hs
w3_old = w3
pt_old = pt
!q_con_old = q_con
delp_old = delp
gz_old = gz
ws_old = ws

call system_clock(start_time)
!$ACC DATA copyin(ms,dt,isc,iec,jsc,jec,npx,npy,nz,ng,akap,cappa,cp, &
!$ACC             ptop,hs,w3,pt,q_con,delp, gz, ws, p_fac, a_imp, &
!$ACC             scale_m,gz_old) copyout(gz) copy(pef)
do iter = 1, num_iter
!$acc parallel loop collapse(3) present(gz,gz_old)
  do k = 1, nz+1
    do j = jsd, jed
      do i= isd, ied
        gz(i,j,k) = gz_old(i,j,k)
      enddo
    enddo
  enddo
  call Riem_Solver_c(ms, dt, isc, iec, jsc, jec, nz, ng, akap, &
                     cappa, cp, ptop, hs, w3, pt, q_con, &
                     delp, gz, pef, ws, p_fac, a_imp, scale_m)
enddo
!$ACC END DATA
call system_clock(end_time)

!Check the results.
write(output_unit, '(a,es18.10)') "Sum of gz_out  = ", sum(gz_out)
write(output_unit, '(a,es18.10)') "Sum of gz      = ", sum(gz)  
write(output_unit, '(a,es18.10)') "Sum of pef_out = ", sum(pef_out)
write(output_unit, '(a,es18.10)') "Sum of pef     = ", sum(pef)  

gz_diff = maxval(abs(gz-gz_out))
write(output_unit, '(a,es18.10)') "max diff of gz = ", gz_diff
if (gz_diff .gt. tolerance) then
  call error("gz diff greater than tolerance.")
endif
pef_diff = maxval(abs((pef(isc-1:iec+1,jsc-1:jec+1,:)- &
                  pef_out(isc-1:iec+1,jsc-1:jec+1,:))/pef(isc-1:iec+1,jsc-1:jec+1,:)))
write(output_unit, '(a,es18.10)') "max relative diff of pef= ", pef_diff
if (pef_diff .gt. tolerance) then
  call error("pef diff greater than tolerance.")
endif
write(output_unit, *) "Elapsed time (s) = ", (end_time-start_time)/rate


contains


subroutine error(errmsg)
  character(len=*), intent(in) :: errmsg

  write(error_unit, "(a)") "Error: "//trim(errmsg)
  call abort()
end subroutine error


!> @brief Simple error handle for netCDF.
subroutine handle_error(error_flag)
  integer, intent(in) :: error_flag

  if (error_flag .ne. nf90_noerr) then
    call error(nf90_strerror(error_flag))
  endif
end subroutine handle_error


end program main
