set terminal png
set output '1_output.png'
set yrange [0:1050]
# Legge due colonne da a.txt
plot 'final.txt' using 1:2 with points title 'Data t < 1Gyr', \
     'final.txt' using 1:3 with points title 'Data t < 1.5Gyr', \
     'final.txt' using 1:4 with points title 'Data t < 2Gyr', \
     'final.txt' using 1:5 with points title 'Data t < 4Gyr', \
     'final.txt' using 1:6 with points title 'Data t < 8Gyr', \
     'final.txt' using 1:7 with points title 'Data t < 20Gyr', \
     'result.txt' using 1:3:4 with yerrorbars title 'Theoric Points', \
     'result.txt' using 2:5:6 with yerrorbars title 'Experimental Points'