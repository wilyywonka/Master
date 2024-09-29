module PrecMod

  ! Module that helps with global declaration of precision

  use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64
  implicit none

  integer, parameter :: wpf = real32
  integer, parameter :: wpi = int32

  
end module PrecMod