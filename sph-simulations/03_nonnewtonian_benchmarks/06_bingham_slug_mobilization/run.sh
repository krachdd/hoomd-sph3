python3 create_slug_tube.py

rm -f mobilization_results.txt

# g sweep as multiples of the analytic g_crit = 2*tau_y*l/(rho*r*L) = tau_y/6
for TAUY in 3 6; do
    GC=$(python3 -c "print($TAUY/6.0)")
    for F in 0.4 0.6 0.8 0.9 1.1 1.25 1.5 2.0; do
        G=$(python3 -c "print($F*$GC)")
        python3 run_mobilization.py slugtube_init.gsd $TAUY $G 30000
    done
done

python3 postprocess.py
