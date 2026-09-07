python3 create_packing.py 2 16

mpirun -np 6 python3 run_invasion.py packing_init.gsd 0  0.030 8000
mpirun -np 8 python3 run_invasion.py packing_init.gsd 15 0.030 25000

python3 postprocess.py packing 0 15
