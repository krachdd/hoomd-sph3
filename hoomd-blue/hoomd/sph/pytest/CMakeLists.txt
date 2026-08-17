#-----------------------------------------------------------
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de
# ----------------------------------------------------------

set(files __init__.py
          conftest.py
          test_kernels.py
          test_eos.py
          test_lattice.py
          test_singlephase.py
          test_strain_rate.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph/pytest
       )

copy_files_to_build("${files}" "sph-pytest" "*.py")
