# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "3dSEBox_withVTK_autogen"
  "CMakeFiles\\3dSEBox_withVTK_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\3dSEBox_withVTK_autogen.dir\\ParseCache.txt"
  )
endif()
