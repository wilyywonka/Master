module PrecMod

  ! Module that helps with global declaration of precision

  use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64
  implicit none

  integer, parameter :: wpf = real64
  integer, parameter :: wpi = int64

  
end module PrecMod