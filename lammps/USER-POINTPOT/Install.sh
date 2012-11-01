# install/uninstall package files in LAMMPS

if (test $1 = 1) then
    cp pair_point_pot.* ..
elif (test $1 = 0) then
    rm ../pair_point_pot.*
fi
