module TypeModule
  implicit none
  
  type :: paramType
    

  end type paramType

  contains

  subroutine Initialize(param_type)
    implicit none
    type(paramType), intent(inout) :: param_type

    
  end subroutine Initialize
  
end module TypeModule