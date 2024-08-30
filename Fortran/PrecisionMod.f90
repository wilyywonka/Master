module PrecMod

  ! Module that helps with global declaration of precision
  ! Most of this module is taken from Ingve Simonsen, who again modelled it after the Numerical Recipes nrtype

  implicit none
  
  Integer, Parameter :: I4B = Selected_int_kind(9)
  Integer, Parameter :: I2B = Selected_int_kind(4)
  Integer, Parameter :: I1B = Selected_int_kind(2)
  Integer, Parameter :: SP  = Kind(1.0)
  Integer, Parameter :: DP  = Kind(1.0D0)
  Integer, Parameter :: SPC = Kind((1.0,1.0))
  Integer, Parameter :: DPC = Kind((1.0D0,1.0D0))
  Integer, Parameter :: LGT = Kind(.True.)

  integer, parameter :: wpf = SP
  integer, parameter :: wpi = I2B

  
end module PrecMod