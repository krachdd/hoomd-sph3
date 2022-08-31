set(files __init__.py
          sphmodel.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph/sphmodel
       )

copy_files_to_build("${files}" "sph-sphmodel" "*.py")
