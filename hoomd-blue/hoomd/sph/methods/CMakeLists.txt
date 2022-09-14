#-----------------------------------------------------------
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de
# ----------------------------------------------------------

set(files __init__.py
          methods.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph/methods
       )

copy_files_to_build("${files}" "sph-methods" "*.py")
