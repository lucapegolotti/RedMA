#!/bash/bin

for bfunctions in zernike chebyshev
do
    for nmax in 0 1 2 3
    do
        for hmesh in 0.80 0.70 0.60 0.50 0.40
        do
            killall Ross* ; mpirun -n 2 ./RossEthierSteinmanBenchmark $bfunctions $nmax $hmesh
        done
    done
done
