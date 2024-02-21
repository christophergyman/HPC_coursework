# This gnuplot script plots the results of all the averages of sheer advected u[i][j]
# it is assumed that the data is plotted is in a file called average.dat
#to use this file 

# Set the output file format and name
set terminal pngcairo enhanced font 'arial,10' size 800, 600
set output 'average.png'

# Set plot title and axis labels
set title "Line Graph"
set xlabel "X Axis"
set ylabel "Y Axis"

# Plot data from the text file
plot 'average.dat' using 1:2 with linespoints title "Data"
