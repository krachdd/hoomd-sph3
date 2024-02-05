# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=20 --reynolds=10
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=30 --reynolds=10
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=50 --reynolds=1
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=50 --reynolds=10
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=50 --reynolds=100
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=50 --reynolds=400
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=100 --reynolds=10
# mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=10

mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=1
mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=100
mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=400
mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=1000
mpirun -np 1 ./create_input_geometries_ldc.py  --resolution=200 --reynolds=10000


#mpirun -np 8  ./run_ldc.py --resolution=20 --reynolds=10 --steps=20001 --initgsd="liddrivencavity_28_28_17_vs_0.05_re_10.0_init.gsd"
#mpirun -np 8  ./run_ldc.py --resolution=30 --reynolds=10 --steps=20001 --initgsd="liddrivencavity_38_38_17_vs_0.03333333333333333_re_10.0_init.gsd"
#mpirun -np 8  ./run_ldc.py --resolution=50 --reynolds=10 --steps=20001 --initgsd="liddrivencavity_58_58_17_vs_0.02_re_10.0_init.gsd"
#mpirun -np 8  ./run_ldc.py --resolution=100 --reynolds=10 --steps=20001 --initgsd="liddrivencavity_108_108_17_vs_0.01_re_10.0_init.gsd"
#mpirun -np 8 ./run_ldc_resolution.py --resolution=50 --reynolds=10 --steps=20001 --initgsd="liddrivencavity_58_58_17_vs_0.02_re_10.0_init.gsd"
