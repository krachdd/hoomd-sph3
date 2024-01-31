mpirun -np 1 ./create_input_geometries_tube.py 20 
mpirun -np 1 ./create_input_geometries_tube.py 30 
mpirun -np 1 ./create_input_geometries_tube.py 50 
mpirun -np 1 ./create_input_geometries_tube.py 100 


mpirun -np 8 ./run_tube_TV.py 20 poiseuille_flow_20_32_32_vs_5e-05_init.gsd 10001
# mpirun -np 8 ./run_tube_TV.py 50 poiseuille_flow_50_62_62_vs_2e-05_init.gsd 10001

# mpirun -np 8 ./run_tube.py 20 poiseuille_flow_20_32_32_vs_5e-05_init.gsd 10001
# mpirun -np 8 ./run_tube.py 30 poiseuille_flow_30_42_42_vs_3.3333333333333335e-05_init.gsd 10001
# mpirun -np 8 ./run_tube.py 50 poiseuille_flow_50_62_62_vs_2e-05_init.gsd 10001
# mpirun -np 8 ./run_tube.py 100 poiseuille_flow_100_112_112_vs_1e-05_init.gsd 20001
