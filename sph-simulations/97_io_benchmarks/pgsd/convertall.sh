mpirun -np 2 create_parallel_gsd_from_raw.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > log2.log
mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_init.gsd c2.gsd
mpirun -np 3 create_parallel_gsd_from_raw.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > log3.log
mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_init.gsd c3.gsd
mpirun -np 4 create_parallel_gsd_from_raw.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > log4.log
mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_init.gsd c4.gsd
mpirun -np 8 create_parallel_gsd_from_raw.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > log8.log
mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_init.gsd c8.gsd


python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py c2.gsd 
python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py c3.gsd 
python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py c4.gsd 
python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py c8.gsd
