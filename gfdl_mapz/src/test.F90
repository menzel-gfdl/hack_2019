module test
use iso_fortran_env, only: real64
implicit none


public :: map_level_k
public :: in_layer
public :: with_fractional_parts
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

!!!$acc kernels
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
!!!$acc end kernels
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

!!!$acc kernels
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
!!!$acc end kernels
end subroutine with_fractional_parts


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
