fdir = "/data2/eraint/1.5x1.5/Tpl/"
f    = addfile(fdir+"Tpl.201808.6hr.grib","r")
temp = f->T_GDS0_ISBL(0,:,{38},{127})
lev  = f->lv_ISBL1

fdir = "/data2/eraint/1.5x1.5/Qpl/"
f    = addfile(fdir+"Qpl.201808.6hr.grib","r")
qv = f->Q_GDS0_ISBL(0,:,{38},{127})


print(lev+" "+temp+" "+qv)

