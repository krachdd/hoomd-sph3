mpirun -np 1 ./create_input_geometries_couette_flow_nn.py 20
mpirun -np 1 ./create_input_geometries_couette_flow_nn.py 30
mpirun -np 1 ./create_input_geometries_couette_flow_nn.py 50
mpirun -np 1 ./create_input_geometries_couette_flow_nn.py 100

mpirun -np 8 ./run_couette_flow_nn.py 20 couette_flow_nn_20_28_17_vs_5e-05_init.gsd 20001
mpirun -np 8 ./run_couette_flow_nn.py 30 couette_flow_nn_30_38_17_vs_3.3333333333333335e-05_init.gsd 20001
mpirun -np 8 ./run_couette_flow_nn.py 50 couette_flow_nn_50_58_17_vs_2e-05_init.gsd 20001
mpirun -np 8 ./run_couette_flow_nn.py 100 couette_flow_nn_100_108_17_vs_1e-05_init.gsd 20001
