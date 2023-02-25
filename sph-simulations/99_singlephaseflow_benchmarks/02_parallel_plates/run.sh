mpirun -np 1 ./create_input_geometries_parallel_plates.py 20
mpirun -np 1 ./create_input_geometries_parallel_plates.py 30
mpirun -np 1 ./create_input_geometries_parallel_plates.py 50
mpirun -np 1 ./create_input_geometries_parallel_plates.py 100

mpirun -np 8 ./run_parallel_plates.py 20 parallel_plates_20_32_17_vs_5e-05_init.gsd
mpirun -np 8 ./run_parallel_plates.py 30 parallel_plates_30_42_17_vs_3.3333333333333335e-05_init.gsd
mpirun -np 8 ./run_parallel_plates.py 50 parallel_plates_50_62_17_vs_2e-05_init.gsd
mpirun -np 8 ./run_parallel_plates.py 100 parallel_plates_100_112_17_vs_1e-05_init.gsd

