!! -------------------------------------------------------------------------- !!
Program app1
!! -------------------------------------------------------------------------- !!
    use glbApp1
!! -------------------------------------------------------------------------- !!
Implicit None
!! -------------------------------------------------------------------------- !!
    Type(typShape) :: shape1
    Character(len=:), Allocatable :: shape1Color
    Logical :: shape1Filled

    shape1%color = "red"
    shape1%filled = .TRUE.

    Call shape1%GetColor( shape1Color )
    shape1Filled = shape1%IsFilled()

    write(*,*) "Shape1 Color : ", shape1Color
    write(*,*) "Shape1 Filled: ", shape1Filled

!! -------------------------------------------------------------------------- !!
End Program
!! -------------------------------------------------------------------------- !!
