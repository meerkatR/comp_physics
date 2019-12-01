program circle_area
 implicit None
 real(8)::  radius, circum, area
 real(8)::  PI=3.1415926
 integer:: model_n=1
 print *, 'enter a radius:'
 read *, radius
 circum=2.0*PI*radius
 area=radius*radius*PI
 print *, 'radius=', radius
 print *, 'circumference=', circum
 print *, 'area=', area
 print *, 'inter=', model_n
end program circle_area
