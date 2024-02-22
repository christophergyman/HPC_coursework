# HPC_coursework
This is my submission for my HPC coursework

## Final images breakdown:
1. 'initial.png' --> Initial advection2D.c without modifications
2. 'modifyingCalc.png' --> After the initial values were changed 
3. 'final.png' --> After the vertical sheer was added
4. 'average.png' --> Task4 average advected distance compared to y values

## How to run with OpenMP?
Go into the '/src' folder and run the following commands

Assuming your running linux run the following commands in the '/src' directory:

```
gcc -fopenmp -o advection2D -std=c99 advection2D.c -lm
./advection2D
gnuplot plot_final.sh
gnuplot plot_average.sh
open final.png
open plot_average.sh
```

## How to run without OpenMP?
Go into the '/src' folder and run the following commands

Assuming your running linux run the following commands in the '/src' directory:

```
gcc -o advection2D -std=c99 advection2D.c -lm
./advection2D
gnuplot plot_final.sh
gnuplot plot_average.sh
open final.png
open plot_average.sh
```
