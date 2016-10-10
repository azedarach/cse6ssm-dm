set terminal x11

set title "CSE6SSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='CSE6SSM_rgflow.dat'

plot for [i=2:248+1] filename using 1:(column(i)) title columnhead(i)
