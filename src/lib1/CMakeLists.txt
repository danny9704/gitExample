SET(ITEM_PREFIX "LIB")
SET(ITEM_NAME   "lib1")

set(${ITEM_PREFIX}_SRC_${ITEM_NAME}
    ##.... Wrapper Module of EaFort
    src/glbLib1.F90
    src/pkgLib1.F90
)

##... Compile
AutoCompile(
    IS_COMPILE_LIB      True
    IS_DYNAMIC_LIB      "${COMPILE_DYNAMIC_LIB}"
    INSTALL_NAME        "${${ITEM_PREFIX}_NAME_${ITEM_NAME}}"
    INSTALL_SRC         "${${ITEM_PREFIX}_SRC_${ITEM_NAME}}"
    INSTALL_PATH        "${${ITEM_PREFIX}_PATH_${ITEM_NAME}}"
    INSTALL_DEPEND      "${${ITEM_PREFIX}_DEPEND_${ITEM_NAME}}"
    INSTALL_EXTLIB      "${${ITEM_PREFIX}_EXTLIB_${ITEM_NAME}}"
    INSTALL_EXTINC      "${${ITEM_PREFIX}_EXTINC_${ITEM_NAME}}"
    INSTALL_EXTFLAG     "${${ITEM_PREFIX}_EXTFLAG_${ITEM_NAME}}"
)
