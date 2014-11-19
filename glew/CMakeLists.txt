CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include )
PROJECT( glew )
ADD_LIBRARY( ${PROJECT_NAME} STATIC include/GL/glew.h src/glew.c )
# Default preprocessors
#set_default_cs_preproc( ${PROJECT_NAME} )
# Add custom prepocessor definitions
set_property( TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS GLEW_STATIC )
