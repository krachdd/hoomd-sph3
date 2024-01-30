mpirun -np 1 ./create_input_geometries_channel_flow.py 20
mpirun -np 1 ./create_input_geometries_channel_flow.py 30
mpirun -np 1 ./create_input_geometries_channel_flow.py 50
mpirun -np 1 ./create_input_geometries_channel_flow.py 100


# mpirun -np 8 ./run_channel_flow.py 20 channel_flow_20_32_32_vs_5e-05_init.gsd 10001
# mpirun -np 8 ./run_channel_flow.py 30 channel_flow_30_42_42_vs_3.3333333333333335e-05_init.gsd 10001
mpirun -np 8 ./run_channel_flow.py 50 channel_flow_50_62_62_vs_2e-05_init.gsd 10001
# mpirun -np 8 ./run_channel_flow.py 100 channel_flow_100_112_112_vs_1e-05_init.gsd 20001


mpirun -np 8 ./run_channel_flow_TV.py 30 channel_flow_30_42_42_vs_3.3333333333333335e-05_init.gsd 10001
mpirun -np 8 ./run_channel_flow_TV.py 50 channel_flow_50_62_62_vs_2e-05_init.gsd 10001
