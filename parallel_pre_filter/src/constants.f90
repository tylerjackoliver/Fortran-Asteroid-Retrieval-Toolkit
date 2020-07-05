module constants

  use precision_kinds
  
  implicit none

  double precision :: pi = 4.0_dp*datan(1.0_dp)

  double precision :: mu = 1.32712440018e11 ! km^3/s^-2 

  double precision :: au = 149597870.7000 ! km

end module constants
