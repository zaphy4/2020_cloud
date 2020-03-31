#!/bin/ksh

cat > namelist << EOF
&variable
    nt   = 100,
    nz   = 20, 
    ps   = 1000.,
    ptop = 10.,
    stat = 'ideal'          ! -- 'ideal' or 'real' case
/
EOF

gfortran file_io.f90 temp_ideal.f90 driver.f90 -o a.out
./a.out
