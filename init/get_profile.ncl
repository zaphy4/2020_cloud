; read namelist
line = asciiread("../namelist",-1,"string")
nt   = tofloat(str_get_field(line(1),2,"=, "))
nz   = toint(str_get_field(line(2),2,"=, "))
ps   = tofloat(str_get_field(line(3),2,"=, "))
ptop = tofloat(str_get_field(line(4),2,"=, "))
grid = str_get_field(line(5),2,"=, ")
stat = str_get_field(line(6),2,"=, ")
;print(line)
;print(nt+" "+nz+" "+ps+" "+ptop+" "+grid+" "+stat)


; read file
fdir = "/data2/eraint/1.5x1.5/Tpl/"
f    = addfile(fdir+"Tpl.201808.6hr.grib","r")
temp = f->T_GDS0_ISBL(0,:,{38},{127})
lev  = f->lv_ISBL1

fdir = "/data2/eraint/1.5x1.5/Qpl/"
f    = addfile(fdir+"Qpl.201808.6hr.grib","r")
qv = f->Q_GDS0_ISBL(0,:,{38},{127})


; -- make pressure level
splev = new( toint(nz+1),"float")
if (grid .eq."'const'") then

  dz = (ptop-ps)/nz
  splev(0) = ps
  do i = 0, nz-1
    splev(i+1) = dz + splev(i)
  end do

else if (grid .eq."'stret'") then


end if
end if
; print(splev) ; check


; -- 1D interpolation
; function
tt_f = int2p(lev, temp, splev, 1) ; 0 is log interpolation


; code
;tt = temp(0) + (X-X1 )*(y2-y1)/(x2-x1)
tt = new( nz+1, float)



else

do i = 0, nz
aa = closest_val(splev(i), lev)

if (splev(i).eq.lev(aa)) then 
tt(i) = temp(aa)

else if (splev(i).lt.min(lev))then      ; extrapolation
tt(i) = temp(aa) + (splev(i)-lev(aa))*(temp(aa+1)-temp(aa))/(lev(aa+1)-lev(aa))
print("aa")
else if (splev(i).gt.max(lev))then ; extrapolation, just adjust former gradient 
tt(i) = temp(aa) + (splev(i)-lev(aa))*(temp(aa)-temp(aa-1))/(lev(aa)-lev(aa-1))


else if (splev(i).gt.lev(aa)) then ; interpolation
tt(i) = temp(aa) + (splev(i)-lev(aa))*(temp(aa+1)-temp(aa))/(lev(aa+1)-lev(aa))

else if (splev(i).lt.lev(aa)) then ; interpolation
tt(i) = temp(aa-1) + (splev(i)-lev(aa-1))*(temp(aa)-temp(aa-1))/(lev(aa)-lev(aa-1))

end if
end if
end if
end if
end if

end do


end if

print("fun:"+tt_f+" made:"+tt)





