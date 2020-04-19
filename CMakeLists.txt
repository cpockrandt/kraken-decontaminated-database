cmake_minimum_required (VERSION 3.2)
project (decontam_tools CXX)

set (CXX11_FOUND "1")
set (CMAKE_PREFIX_PATH "include/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "include/seqan/include")

find_package(SeqAn REQUIRED)

include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})
add_definitions(-DSEQAN_DISABLE_VERSION_CHECK="YES")

find_package (SeqAn3 REQUIRED
		     HINTS "include/seqan3/build_system")
find_package(OpenMP REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -fopenmp")

add_executable (decontam_conterminator decontam_conterminator.cpp)
target_link_libraries(decontam_conterminator ${SEQAN_LIBRARIES})

add_executable (decontam_kmer decontam_kmer.cpp include/KMC/kmc_api/kmc_file.cpp include/KMC/kmc_api/kmer_api.cpp include/KMC/kmc_api/mmer.cpp)
target_link_libraries (decontam_kmer seqan3::seqan3)
target_link_libraries (decontam_kmer pthread)
