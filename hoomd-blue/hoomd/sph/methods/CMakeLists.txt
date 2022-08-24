set(files __init__.py
          methods.py
#           rattle.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph/methods
       )

copy_files_to_build("${files}" "sph-methods" "*.py")
