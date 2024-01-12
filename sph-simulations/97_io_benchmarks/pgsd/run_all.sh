mpirun -np 2 create_parallel_gsd_from_raw.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > testcreate.log

mpirun -np 2 run_parallel_gsd.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > rlog2.log
mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_run_gsd.gsd rc2.gsd
python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py rc2.gsd 

# mpirun -np 3 run_parallel_gsd.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > rlog3.log
# mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_run_gsd.gsd rc3.gsd
# python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py rc3.gsd 

# mpirun -np 4 run_parallel_gsd.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > rlog4.log
# mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_run_gsd.gsd rc4.gsd
# python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py rc4.gsd 

# mpirun -np 8 run_parallel_gsd.py input_SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05.txt > rlog8.log
# mv SPHERE3D_bcc_NX100_NY100_NZ100_PHI0.5_DIAMETER0.0007815926417967722_VSIZE1e-05_run_gsd.gsd rc8.gsd
# python3 ../../../helper_modules/gsd2vtu/spfgsd2vtu.py rc8.gsd


