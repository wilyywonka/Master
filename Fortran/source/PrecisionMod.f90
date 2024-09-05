module PrecMod

  ! Module that helps with global declaration of precision
  use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64
  implicit none

  ! Integer, Parameter :: I4B = Selected_int_kind(9)
  ! Integer, Parameter :: I2B = Selected_int_kind(4)
  ! Integer, Parameter :: I1B = Selected_int_kind(2)
  ! Integer, Parameter :: SP  = Kind(1.0)
  ! Integer, Parameter :: DP  = Kind(1.0D0)
  ! Integer, Parameter :: SPC = Kind((1.0,1.0))
  ! Integer, Parameter :: DPC = Kind((1.0D0,1.0D0))
  ! Integer, Parameter :: LGT = Kind(.True.)

  integer, parameter :: wpf = real32
  integer, parameter :: wpi = int32

  
end module PrecMod