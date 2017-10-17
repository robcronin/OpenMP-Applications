reset
set term postscript eps

set xlabel "Size"
set ylabel "Speed-up"

set title "Size vs Speed-up for Gauss Elimination"
set output "gauss-size.ps"
gauss="gauss_size.dat"
plot gauss w l

set title "Size vs Speed-up for Sieve of Eratosthenes"
set output "sieve-size.ps"
sieve="sieve_size.dat"
plot sieve w l

set xlabel "Number of Processes"
set title "Processes vs Speed-up"
set output "procs.ps"
gauss ="gauss_procs.dat"
sieve ="sieve_procs.dat"
plot gauss w l, sieve w l
