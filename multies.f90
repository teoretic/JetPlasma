module multies

use names

implicit none

contains

pure function vect_mult2d (a,b) result (c)

  real (kind=precis),dimension(1:2), Intent(In) :: a,b
  real (kind=precis) :: c

  c = a(1)*b(2)-a(2)*b(1)

end function vect_mult2d


pure function vect_mult (a,b) result (c)

  real (kind=precis),dimension(1:3), Intent(In) :: a,b
  real (kind=precis),dimension(1:3) :: c

  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)

end function vect_mult

pure function TriangleVol(a,b) result (S)
  real (kind=precis), dimension(1:2), Intent(In):: a,b
  real (kind=precis) :: S

  S = 0.5*abs(a(1)*b(2)-b(1)*a(2))
end function TriangleVol

pure function TriangleVol3(a,b,c) result (S)
  real (kind=precis), dimension(1:3), Intent(In):: a,b,c
  real (kind=precis), dimension(1:3) :: d,e,f
  real (kind=precis) :: S

  d = b-a
  e = c-a
  f(1) = d(2)*e(3)-e(2)*d(3)
  f(2) = d(3)*e(1)-e(3)*d(1)
  f(3) = d(1)*e(2)-e(1)*d(2)
  S = 0.5*norm2_3(f)
end function TriangleVol3

pure function distance(a,b) result (c)
  real (kind=precis), dimension(1:2), Intent(In):: a,b
  real (kind=precis) :: c

  c = sqrt(sum((b-a)**2))
end function distance

pure function distance3(a,b) result (c)
  real (kind=precis), dimension(1:3), Intent(In):: a,b
  real (kind=precis) :: c

  c = sqrt(sum((b-a)**2))
end function distance3


pure function signum(numb) result(sig)

  integer, Intent(In) :: numb
  integer :: sig

  sig = 1
  if (numb<0) sig = -1

end function signum

pure function signumr(numb) result(sig)

  real(kind=precis), Intent(In) :: numb
  real(kind=precis) :: sig

  sig = 1.
  if (numb<0._precis) sig = -1.

end function signumr

pure function norm2(a) result (c)
  real (kind=precis), dimension(1:2), Intent(In):: a
  real (kind=precis) :: c

  c = sqrt(sum(a**2))
end function norm2

pure function norm2_3(a) result (c)
  real (kind=precis), dimension(1:3), Intent(In):: a
  real (kind=precis) :: c

  c = sqrt(sum(a**2))
end function norm2_3

pure function pintr (a,t1,t2,t3) result (pin)
  real (kind=precis), dimension(1:2), Intent(In):: a,t1,t2,t3
  logical :: pin

  pin = (signumr(vect_mult2d(a-t1,t2-t1))<1E-14 .and. signumr(vect_mult2d(a-t2,t3-t2))<1E-14 .and. signumr(vect_mult2d(a-t3,t1-t3))<1E-14)

!  pin = Abs(TriangleVol(t2-t1,t3-t1)-(TriangleVol(a-t1,t2-t1)+TriangleVol(a-t2,t3-t2)+TriangleVol(a-t3,t1-t3)))<1E-4

!print *, vect_mult2d(a-t1,t2-t1), vect_mult2d(a-t2,t3-t2), vect_mult2d(a-t3,t1-t3)
end function

end module multies