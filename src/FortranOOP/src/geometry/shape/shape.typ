!!... Class: Shape
Type typShape

    !!... Color
    Character(len=:), Allocatable :: color

    !!... Does the color is filled ?
    Logical :: filled = .FALSE.

Contains

    !!... Get Color 
    Procedure :: GetColor => GetColor_typShape

    !!... Get Filled
    Procedure :: IsFilled => IsFilled_typShape
    
End Type