mpirun -np 1 ./create_input_geometies_couette_flow.py 20
mpirun -np 1 ./create_input_geometies_couette_flow.py 30
mpirun -np 1 ./create_input_geometies_couette_flow.py 50
mpirun -np 1 ./create_input_geometies_couette_flow.py 100

mpirun -np 8 ./run_couette_flow.py 20 couette_flow_20_28_17_vs_5e-05_init.gsd
mpirun -np 8 ./run_couette_flow.py 30 couette_flow_30_38_17_vs_3.3333333333333335e-05_init.gsd
mpirun -np 8 ./run_couette_flow.py 50 couette_flow_50_58_17_vs_2e-05_init.gsd
mpirun -np 8 ./run_couette_flow.py 100 couette_flow_100_108_17_vs_1e-05_init.gsd
