!
! Argonne Leadership Computing Facility
! GFDL application kernel for Held Suarez benchmark
! Vitali Morozov (morozov@anl.gov)
! Aug 11, 2010
!
program main
use iso_fortran_env, only: error_unit, output_unit, real32, real64
use mpi
#ifdef _OPENACC
use openacc
#endif
use fyppm_mod, only: fp, fyppm

implicit none

real(kind=fp), parameter :: expected_gsum = 9.6326547367e+10_fp !< Benchmark result.
real(kind=fp), parameter :: tolerance = 1.e-10_fp !< Tolerance for test to pass.
integer, parameter :: nx = 96 !<Grid x size.
integer, parameter :: ny = 96 !< Grid y size.
integer, parameter :: km = 91 !< Grid z size.
integer :: ierror !< MPI return code.
real(kind=fp) :: gsum !< Global sum, used to compare to benchmark.
integer :: myrank !< Current MPI rank.
integer :: num_ranks !< Number of MPI ranks.
integer :: num_devices !< Number of devices.
real(kind=real64) :: time1, time2 !< Time stamps for main loop.

integer, dimension(2) :: layout
integer   :: i, j, k, iter
integer :: ifirst, ilast, isd, ied
integer :: jfirst, jlast, jsd, jed, js, je
integer :: jord, npy, ng
real(kind=fp), allocatable, dimension(:,:,:) :: c
real(kind=fp), allocatable, dimension(:,:,:) :: q
real(kind=fp), allocatable, dimension(:,:,:) :: wk
real(kind=fp), allocatable, dimension(:,:,:) :: area
real(kind=fp), allocatable, dimension(:,:,:) :: al, bl, br, dq
real(kind=fp), allocatable, dimension(:,:) :: dya
real(kind=fp), allocatable, dimension(:,:,:) :: flux, dm
real(kind=fp) :: ppm_limiter = 2.0_fp
integer :: nxc, nyc, ipos, jpos
integer :: istart, iend, jstart, jend

!Initialize MPI and if using OpenACC associate one device per rank.
call mpi_init(ierror)
call mpi_error_check(ierror)
call mpi_comm_rank(mpi_comm_world, myrank, ierror)
call mpi_error_check(ierror)
call mpi_comm_size(mpi_comm_world, num_ranks, ierror)
call mpi_error_check(ierror)
#ifdef _OPENACC
num_devices = acc_get_num_devices(acc_device_nvidia)
#else
num_devices = 0
#endif
if (num_devices .gt. 0) then
  if (num_ranks .ne. num_devices) then
    call error("This code requires that the number of MPI equals the number of GPUs.")
  endif
#ifdef _OPENACC
  call acc_set_device_num(myrank, acc_device_nvidia)
#endif
  if (myrank .eq. 0) then
    write(output_unit, "(a)") "MPI rank to device mapping:"
    do i = 0, num_ranks-1
      write(output_unit, "(a,i2.2,a,i2.2)") "MPI rank ", i, " using device ", i
    enddo
  endif
else
  if (myrank .eq. 0) then
    write(output_unit, "(a,i2.2,a)") "No devices found, using ", num_ranks, " MPI ranks."
  endif
endif
call mpi_barrier(mpi_comm_world, ierror)
call mpi_error_check(ierror)

!Divide grid up amongst ranks.
call define_layout(nx, ny, num_ranks, layout)
nxc = nx/layout(1)
nyc = ny/layout(2)
ipos = mod(myrank, layout(1))
jpos = myrank/layout(1)
isd = ipos*nxc + 1
ied = isd + nxc - 1
jsd = jpos*nyc + 1
jed = jsd + nyc - 1
ifirst = isd
ilast  = ied
jfirst = jsd
jlast  = jed
js = jsd
je = jed
ng = 3
npy = ny + 1
if (myrank .eq. 0) then
  write(output_unit, "(a,i,a)") "Number of grid points per rank in x-direction: ", nxc
  write(output_unit, "(a,i,a)") "Number of grid points per rank in y-direction: ", nyc
  write(output_unit, "(a,i,a)") "Number of grid points per rank in z-direction: ", km
endif

!Allocate and initialize buffers.
allocate(c(isd:ied, js:je+1, km))
allocate(dya(isd:ied, js:je))
allocate(q(ifirst:ilast, jfirst-ng:jlast+ng, km))
allocate(area(ifirst:ilast, jfirst:jlast+1, km))
allocate(wk(ifirst:ilast, jfirst:jlast+1, km))
allocate(flux(isd:ied, js:je+1, km))
allocate(dm(isd:ied, js-2:je+2, km))
allocate(al(ifirst:ilast, jfirst-1:jlast+2, km))
allocate(bl(ifirst:ilast, jfirst-1:jlast+1, km))
allocate(br(ifirst:ilast, jfirst-1:jlast+1, km))
allocate(dq(ifirst:ilast, jfirst-3:jlast+2, km))
do k = 1, km
  do j = jfirst-ng, jlast+ng
    do i = ifirst, ilast
      q(i,j,k) = i + j + k
    enddo
  enddo
enddo
do k = 1, km
  do j = js, je+1
    do i = isd, ied
      c(i,j,k) = 1.0_fp + i*1.e-2_fp + j*2.e-4_fp + k*3._fp*1.e-6_fp
    enddo
  enddo
enddo
do k = 1, km
  do j = jfirst, jlast+1
    do i = ifirst, ilast
      wk(i,j,k) = 1.0_fp + k*1.e-2_fp + j*2.e-4_fp + i*3._fp*1.e-6_fp
      area(i,j,k) = 12.0_fp + k*1.e-2_fp + j*2.e-4_fp + i*3._fp*1.e-6_fp
    enddo
  enddo
enddo
do j = js, je
  do i = isd, ied
    dya(i,j) = 2.0_fp + i*1.e-2_fp + j*2.e-4_fp
  enddo
enddo

!Get time before main loop.
call cpu_time(time1)

!!$ACC enter data copyin(wk,area,c,q,dya,dm,flux)
!$ACC DATA copy(wk) copyin(area,c,q,dya,dm,flux,ppm_limiter,jord,jfirst,jlast,npy) copyout(al,bl,br,dq)
do iter = 1, 10000
  if (mod(iter, 10) .lt. 6) then
    jord = 9
  else
    jord = 12
  endif
!$omp parallel do
!$acc kernels
  do k = 1, km
    flux(:,:,k) = 1.0_fp
    dm(:,:,k) = 1.0_fp
  enddo
!$acc end kernels
  call fyppm(c, q, flux, jord, ifirst, ilast, jfirst, jlast, npy, dm, ppm_limiter, &
             dya, isd, ied, jsd, jed, js, je, ng, km, al, bl, br, dq)
!$omp parallel do
!$acc kernels loop collapse(3)
  do k = 1, km
    do j = jfirst, jlast+1
      do i = ifirst, ilast
        wk(i,j,k) = wk(i,j,k) + flux(i,j,k)/area(i,j,k)
      enddo
    enddo
  enddo
!$acc end kernels
enddo
!$ACC END DATA
!!!$ACC exit data copyout(wk)

!Get time after main loop.
call cpu_time(time2)

!Check against the benchmark.
istart = ifirst
iend = ilast
jstart = jfirst
jend = jlast
if (jlast .eq. ny) then
  jend = jend + 1
endif
gsum = global_sum(wk(istart:iend,jstart:jend,:))
if (myrank .eq. 0) then
  write(output_unit, "(a,es18.10)") "Using tolerance = ", tolerance
  write(output_unit, "(a,es18.10)") "Expected result = ", expected_gsum
  write(output_unit, "(a,es18.10)") "Actual result   = ", gsum
  if (abs(expected_gsum - gsum)/expected_gsum .gt. tolerance) then
    call error("Global sum does not match expected output to with tolerance.")
  else
    write(output_unit, "(a)") "Test passed successfully."
  endif
  write(output_unit, "(a,es18.10)") "Elapsed walltime (s) = ", time2 - time1
endif

!Clean up.
deallocate(c)
deallocate(dya)
deallocate(q)
deallocate(area)
deallocate(wk)
deallocate(flux)
deallocate(dm)
deallocate(al)
deallocate(bl)
deallocate(br)
deallocate(dq)
call mpi_finalize(ierror)


contains


!> @brief Abort if an error occurs.
subroutine error(errmsg, code)
  character(len=*), intent(in) :: errmsg !< Error message.
  integer, intent(in), optional :: code !< Return code. Defaults to one.

  integer :: err

  write(error_unit, "(a)") "Error: "//trim(errmsg)
  if (present(code)) then
    call mpi_abort(mpi_comm_world, code, err)
  else
    call mpi_abort(mpi_comm_world, 1, err)
  endif
end subroutine error


!> @brief Trap MPI errors.
subroutine mpi_error_check(ierror)
  integer, intent(in) :: ierror !< Return code from MPI call.

  character(len=mpi_max_error_string) :: msg
  integer :: resultlen
  integer :: e

  if (ierror .ne. mpi_success) then
    call mpi_error_string(ierror, msg, resultlen, e)
    call error(msg, ierror)
  endif
end subroutine mpi_error_check


!> @brief Calculate domain layout.
subroutine define_layout(isz, jsz, ndivs, layout)
  integer, intent(in) :: isz
  integer, intent(in) :: jsz
  integer, intent(in) :: ndivs 
  integer, dimension(2), intent(out) :: layout

  layout(1) = 1
  layout(2) = ndivs
  if (mod(jsz, layout(2)) .ne. 0) then
    call error("the y size of the grid must be evenly divisible by the number of MPI ranks.")
  endif
end subroutine define_layout


!> @brief Compute a global sum.
function global_sum(data)
  real(kind=fp), dimension(:,:,:), intent(in) :: data !< Data to be summed.

  real(kind=fp) :: sum_data
  real(kind=fp) :: global_sum
  integer :: err
  integer :: s
  character(len=128) :: msg

  if (fp .eq. real64) then
    s = mpi_real8
  elseif (fp .eq. real32) then
    s = mpi_real4
  else
    call error("floating point kind must be either real64 or real32.")
  endif
  sum_data = sum(data)
  call mpi_allreduce(sum_data, global_sum, 1, s, mpi_sum, mpi_comm_world, err)
  call mpi_error_check(err)
end function global_sum


end program main
