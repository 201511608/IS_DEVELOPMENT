Visual basic
 'Dim declare in Meamory
 Dim x as integer
 Dim y as integer = 100
 Dim x,y as integer
 
  Dim anArray() As Integer = {1, 3, 5, 7, 9}
  
  Dim fname, lname, fullname, greetings As String
  fname = "Rowan"
	  
 Const maxval As Long = 4999
 Public Const message As String = "HELLO"   ' works in class decleration
 Private Const piValue As Double = 3.1415

 Console.WriteLine("Hello World")
 Console.WriteLine("declaring on the day of: {0}", da)
 Console.Write(c & "muzafar" & c & "haha " & da & " Date")
   Console.WriteLine("Area = " & Str(area))
  Console.WriteLine("Length: {0}", length)
 
 
 Console.ReadKey()
 Console.ReadLine()  ' Enter it moves

 Dim message As String
 message = Console.ReadLine
 Dim message2 As String = Console.ReadLine


 'sub no return value
'function return value

'shared invoked without creating an object
Shared Sub Main()   End Sub
'static invoked with out creating an object 
static Sub Main() End Sub


'Global
Enum Colors
    red = 1
    orange = 2
End Enum
Console.Write("colour  " & Colors.red)

MsgBox("User name is" & Colors.red)


Class Box
Public length As Double
Public breadth As Double   
Public height As Double
End Class


Class Box
Public length As Double
Public breadth As Double   
Public height As Double
End Class

If (a <= 20) Then
   c= c+1
End If

For a = 10 To 20
    Console.WriteLine("value of a: {0}", a)
Next

 Dim anArray() As Integer = {1, 3, 5, 7, 9}
 Dim arrayItem As Integer
 For Each arrayItem In anArray
       Console.WriteLine(arrayItem)
 Next

While a < 20
   Console.WriteLine("value of a: {0}", a)
   a = a + 1
End While


Dim grade As Char
grade = "B"

Select grade
          Case "A"
              Console.WriteLine("Excellent!")
          Case "B", "C"
              Console.WriteLine("Well done")
          Case Else
              Console.WriteLine("Invalid grade")
 End Select
 
 
 Continue Do