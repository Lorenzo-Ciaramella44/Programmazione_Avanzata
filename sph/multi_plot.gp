# Terminale for PNG output
set terminal pngcairo size 800,600
load "param.txt"
# Loop through data files
do for [i=0:a] {
    set output sprintf("frame_%03d.png", i)
    plot "particles_step_".i.".txt" u 1:4
}

# 1 x [m]
# 2 rho [kg/m^3]
# 3 velocity [m/s]
# 4 pressure [Pa]
# 5 accel [m/s^2]
# 6 energy [J/kg]
# 7 du_dt [J/(kg·s)]
# 8 h [m]
