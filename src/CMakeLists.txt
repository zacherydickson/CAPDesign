cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
 
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
 
# build app with seqan3
add_executable (CAPDesign CAPDesign.cpp)
target_link_libraries (CAPDesign seqan3::seqan3)
