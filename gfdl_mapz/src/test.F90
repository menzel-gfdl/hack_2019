module test
use iso_fortran_env, only: real64
implicit none


public :: map_level_k
public :: in_layer
public :: with_fractional_parts
public :: positive_definite_limiter
public :: standard_ppm_limiter
public :: full_monotonicity_limiter
integer, parameter :: fp = real64


contains


!> @brief Determines the k bounds for pe2 within pe1.
subroutine map_level_k(is, ie, js, je, km, kn, pe1, pe2, k1, k2)
  integer, intent(in) :: is
  integer, intent(in) :: ie
  integer, intent(in) :: js
  integer, intent(in) :: je
  integer, intent(in) :: km
  integer, intent(in) :: kn
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe1
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe2
  integer, dimension(is:, js:, :), intent(inout) :: k1
  integer, dimension(is:, js:, :), intent(inout) :: k2

  integer :: i
  integer :: j

  do j = js, je
    do i = is, ie
      call bounds_old(km, kn, pe1(i,j,:), pe2(i,j,:), k1(i,j,:), k2(i,j,:))
    enddo
  enddo
end subroutine map_level_k


subroutine bounds_old(km, kn, array, val, lower, upper)
  integer, intent(in) :: km
  integer, intent(in) :: kn
  real(kind=fp), dimension(:), intent(in) :: array
  real(kind=fp), dimension(:), intent(in) :: val
  integer, dimension(:), intent(inout) :: lower
  integer, dimension(:), intent(inout) :: upper

  integer :: k0
  integer :: k
  integer :: l
  integer :: m

  k0 = 1
  do k = 1, kn
    do l = k0, km
      if (val(k) .ge. array(l) .and. val(k) .le. array(l+1)) then
        lower(k) = l
        if (val(k+1) .le. array(l+1)) then
          upper(k) = l
          k0 = l
        else
          do m = l+1, km
            if (val(k+1) .le. array(m+1)) then
              upper(k) = m
              k0 = m
              exit
            endif
          enddo
        endif
        exit
      endif
    enddo
  enddo
end subroutine bounds_old


subroutine in_layer(is, ie, js, je, km, kn, pe1, pe2, k1, k2, dp1, q4, q2)
  integer, intent(in) :: is
  integer, intent(in) :: ie
  integer, intent(in) :: js
  integer, intent(in) :: je
  integer, intent(in) :: km
  integer, intent(in) :: kn
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe1
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe2
  integer, dimension(is:, js:, :), intent(in) :: k1
  integer, dimension(is:, js:, :), intent(in) :: k2
  real(kind=fp), dimension(is:, js:, :), intent(in) :: dp1
  real(kind=fp), dimension(:, is:, js:, :), intent(in) :: q4
  real(kind=fp), dimension(is:, js:, :), intent(inout) :: q2

  real(kind=fp) :: pl
  real(kind=fp) :: pr
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  real, parameter :: r3 = 1./3.

!$acc kernels
!$acc loop collapse(3)
  do k = 1, kn
    do j = js, je
      do i = is, ie
        l = k1(i,j,k)
        if (k2(i,j,k) .eq. l) then
          pl = (pe2(i,j,k) - pe1(i,j,l))/dp1(i,j,l)
          pr = (pe2(i,j,k+1) - pe1(i,j,l))/dp1(i,j,l)
          q2(i,j,k) = q4(2,i,j,l) + 0.5*(q4(4,i,j,l) + q4(3,i,j,l) - q4(2,i,j,l)) &
                      *(pr + pl) - q4(4,i,j,l)*r3*(pr*(pr + pl) + pl*pl)
        endif
      enddo
    enddo
  enddo
!$acc end kernels
end subroutine in_layer


subroutine with_fractional_parts(is, ie, js, je, km, kn, pe1, pe2, k1, k2, dp1, q4, q2)
  integer, intent(in) :: is
  integer, intent(in) :: ie
  integer, intent(in) :: js
  integer, intent(in) :: je
  integer, intent(in) :: km
  integer, intent(in) :: kn
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe1
  real(kind=fp), dimension(is:, js:, :), intent(in) :: pe2
  integer, dimension(is:, js:, :), intent(in) :: k1
  integer, dimension(is:, js:, :), intent(in) :: k2
  real(kind=fp), dimension(is:, js:, :), intent(in) :: dp1
  real(kind=fp), dimension(:, is:, js:, :), intent(in) :: q4
  real(kind=fp), dimension(is:, js:, :), intent(inout) :: q2

  real(kind=fp) :: pl
  real(kind=fp) :: qsum
  real(kind=fp) :: dp
  real(kind=fp) :: esl
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: m
  integer :: u
  real, parameter :: r3 = 1./3.
  real, parameter :: r23 = 2./3.

!$acc kernels
!$acc loop collapse(3)
  do k = 1, kn
    do j = js, je
      do i = is, ie
        l = k1(i,j,k)
        u = k2(i,j,k)
        if (u .gt. l) then
          !Fractional area.
          pl = (pe2(i,j,k) - pe1(i,j,l))/dp1(i,j,l)
          qsum = (pe1(i,j,l+1) - pe2(i,j,k))*(q4(2,i,j,l) + 0.5*(q4(4,i,j,l) + &
                 q4(3,i,j,l) - q4(2,i,j,l))*(1. + pl) - q4(4,i,j,l)* &
                 (r3*(1. + pl*(1. + pl))))
          do m = l+1, u-1
            !Whole layer.
            qsum = qsum + dp1(i,j,m)*q4(1,i,j,m)
          enddo
          !Fractional area.
          dp = pe2(i,j,k+1) - pe1(i,j,u)
          esl = dp/dp1(i,j,u)
          qsum = qsum + dp*(q4(2,i,j,u) + 0.5*esl*(q4(3,i,j,u) - q4(2,i,j,u) + &
                 q4(4,i,j,u)*(1. - r23*esl)))
          q2(i,j,k) = qsum/(pe2(i,j,k+1) - pe2(i,j,k))
        endif
      enddo
    enddo
  enddo
!$acc end kernels
end subroutine with_fractional_parts


!> @brief Limiter that keeps the data positive definite.
subroutine positive_definite_limiter(i1, i2, j1, j2, k1, k2, a4)
  integer, intent(in) :: i1
  integer, intent(in) :: i2
  integer, intent(in) :: j1
  integer, intent(in) :: j2
  integer, intent(in) :: k1
  integer, intent(in) :: k2
  real, dimension(:, i1:, j1:, k1:), intent(inout) :: a4

  integer :: i
  integer :: j
  integer :: k
  real, parameter :: r12 = 1./12.

!$acc kernels
!$acc loop collapse(3)
  do k = k1, k2
    do j = j1, j2
      do i = i1, i2
        if (abs(a4(3,i,j,k) - a4(2,i,j,k)) .lt. -a4(4,i,j,k)) then
          if (a4(1,i,j,k) + 0.25*(a4(3,i,j,k) - a4(2,i,j,k))**2/a4(4,i,j,k) + a4(4,i,j,k)*r12 .lt. 0) then
            if (a4(1,i,j,k) .lt. a4(3,i,j,k) .and. a4(1,i,j,k) .lt. a4(2,i,j,k)) then
              a4(3,i,j,k) = a4(1,i,j,k)
              a4(2,i,j,k) = a4(1,i,j,k)
              a4(4,i,j,k) = 0.
            elseif (a4(3,i,j,k) .gt. a4(2,i,j,k)) then
              a4(4,i,j,k) = 3.*(a4(2,i,j,k) - a4(1,i,j,k))
              a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
            else
              a4(4,i,j,k) = 3.*(a4(3,i,j,k) - a4(1,i,j,k))
              a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
            endif
          endif
        endif
      enddo
    enddo
  enddo
!$acc end kernels
end subroutine positive_definite_limiter


!> @brief Standard PPM limiter.
subroutine standard_ppm_limiter(i1, i2, j1, j2, k1, k2, dm, a4)
  integer, intent(in) :: i1
  integer, intent(in) :: i2
  integer, intent(in) :: j1
  integer, intent(in) :: j2
  integer, intent(in) :: k1
  integer, intent(in) :: k2
  real, dimension(i1:, j1:, k1:), intent(in) :: dm
  real, dimension(:, i1:, j1:, k1:), intent(inout) :: a4

  integer :: i
  integer :: j
  integer :: k
  real :: da1
  real :: da2
  real :: a6da

!$acc kernels
!$acc loop collapse(3)
  do k = k1, k2
    do j = j1, j2
      do i = i1, i2
        if (dm(i,j,k) .eq. 0.) then
          a4(2,i,j,k) = a4(1,i,j,k)
          a4(3,i,j,k) = a4(1,i,j,k)
          a4(4,i,j,k) = 0.
        else
          da1 = a4(3,i,j,k) - a4(2,i,j,k)
          da2 = da1*da1
          a6da = a4(4,i,j,k)*da1
          if (a6da .lt. -da2) then
            a4(4,i,j,k) = 3.*(a4(2,i,j,k) - a4(1,i,j,k))
            a4(3,i,j,k) = a4(2,i,j,k) - a4(4,i,j,k)
          elseif (a6da .gt. da2) then
            a4(4,i,j,k) = 3.*(a4(3,i,j,k) - a4(1,i,j,k))
            a4(2,i,j,k) = a4(3,i,j,k) - a4(4,i,j,k)
          endif
        endif
      enddo
    enddo
  enddo
!$acc end kernels
end subroutine standard_ppm_limiter


!> @brief Full monotonicity limiter.
subroutine full_monotonicity_limiter(i1, i2, j1, j2, k1, k2, dm, a4)
  integer, intent(in) :: i1
  integer, intent(in) :: i2
  integer, intent(in) :: j1
  integer, intent(in) :: j2
  integer, intent(in) :: k1
  integer, intent(in) :: k2
  real, dimension(i1:, j1:, k1:), intent(in) :: dm
  real, dimension(:, i1:, j1:, k1:), intent(inout) :: a4

  real :: qmp
  integer :: i
  integer :: j
  integer :: k

!$acc kernels
!$acc loop collapse(3)
  do k = k1, k2
    do j = j1, j2
      do i = i1, i2
        qmp = 2.*dm(i,j,k)
        a4(2,i,j,k) = a4(1,i,j,k) - sign(min(abs(qmp), abs(a4(2,i,j,k) - a4(1,i,j,k))), qmp)
        a4(3,i,j,k) = a4(1,i,j,k) + sign(min(abs(qmp), abs(a4(3,i,j,k) - a4(1,i,j,k))), qmp)
        a4(4,i,j,k) = 3.*(2.*a4(1,i,j,k) - (a4(2,i,j,k) + a4(3,i,j,k)))
      enddo
    enddo
  enddo
!$acc end kernels
end subroutine full_monotonicity_limiter


!> @brief Find indices in array that bracket value.
subroutine bounds(km, array, val, lower, upper)
  integer, intent(in) :: km
  real(kind=fp), dimension(:), intent(in) :: array
  real(kind=fp), intent(in) :: val
  integer, intent(out) :: lower
  integer, intent(out) :: upper

  integer :: mid

  lower = 1
  upper = km
  do while (.true.)
    mid = (upper + lower)/2
    if (val .gt. array(mid)) then
      lower = mid
    elseif (val .lt. array(mid)) then
      upper = mid
    else
      lower = mid - 1
      upper = mid
      exit
    endif
    if (upper - lower .eq. 1) then
      exit
    elseif (lower .eq. upper) then
      stop
    endif
  enddo
end subroutine bounds


end module test
