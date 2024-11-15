Subroutine GetColor_typShape(this, color)
    Implicit None
    Class(typShape), intent(in) :: this
    Character(len=:), Allocatable :: color
    
    color = this%color
End Subroutine

Function IsFilled_typShape(this) result(filled)
    Implicit None
    Class(typShape), intent(in) :: this
    Logical :: filled
    filled = this%filled

End Function