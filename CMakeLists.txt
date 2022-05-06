# CMakeList.txt : CMake project for VocalTractLabApi, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/$<CONFIG>)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/$<CONFIG>)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/$<CONFIG>)

add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)

set(BackendSources 
 src/VocalTractLabBackend/AnatomyParams.cpp 
 src/VocalTractLabBackend/Dsp.cpp 
 src/VocalTractLabBackend/F0EstimatorYin.cpp 
 src/VocalTractLabBackend/GeometricGlottis.cpp 
 src/VocalTractLabBackend/Geometry.cpp 
 src/VocalTractLabBackend/GesturalScore.cpp 
 src/VocalTractLabBackend/Glottis.cpp 
 src/VocalTractLabBackend/GlottisFactory.cpp 
 src/VocalTractLabBackend/IirFilter.cpp 
 src/VocalTractLabBackend/ImpulseExcitation.cpp 
 src/VocalTractLabBackend/LfPulse.cpp 
 src/VocalTractLabBackend/Matrix2x2.cpp 
 src/VocalTractLabBackend/PoleZeroPlan.cpp 
 src/VocalTractLabBackend/Sampa.cpp 
 src/VocalTractLabBackend/SegmentSequence.cpp 
 src/VocalTractLabBackend/Signal.cpp 
 src/VocalTractLabBackend/Speaker.cpp
 src/VocalTractLabBackend/Splines.cpp 
 src/VocalTractLabBackend/StaticPhone.cpp 
 src/VocalTractLabBackend/Surface.cpp 
 src/VocalTractLabBackend/Synthesizer.cpp 
 src/VocalTractLabBackend/TdsModel.cpp 
 src/VocalTractLabBackend/TimeFunction.cpp 
 src/VocalTractLabBackend/TlModel.cpp 
 src/VocalTractLabBackend/TriangularGlottis.cpp 
 src/VocalTractLabBackend/Tube.cpp 
 src/VocalTractLabBackend/TwoMassModel.cpp 
 src/VocalTractLabBackend/VocalTract.cpp 
 src/VocalTractLabBackend/VoiceQualityEstimator.cpp 
 src/VocalTractLabBackend/VowelLf.cpp 
 src/VocalTractLabBackend/XmlHelper.cpp 
 src/VocalTractLabBackend/XmlNode.cpp
)

project("VocalTractLabApi")
add_library(VocalTractLabApi SHARED)
target_include_directories (VocalTractLabApi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/VocalTractLabApi)
target_include_directories (VocalTractLabApi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(VocalTractLabApi PRIVATE ${BackendSources} 
src/VocalTractLabApi/VocalTractLabApi.cpp 
src/VocalTractLabApi/VocalTractLabApi.def)


project("VocalTractLabBackend")
add_library(VocalTractLabBackend STATIC)
target_include_directories (VocalTractLabBackend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(VocalTractLabBackend PUBLIC ${BackendSources})


include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  "VtlApiTests"
  "test/VtlApiTests.cpp"
 )

target_include_directories (VtlApiTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/VocalTractLabApi)

target_link_libraries(
  VtlApiTests
  gtest_main
  VocalTractLabApi
)

add_executable(
  "VtlBackendTests"
  "test/VtlBackendTests.cpp"
 )

target_include_directories (VtlBackendTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/VocalTractLabBackend)

target_link_libraries(
  VtlBackendTests
  gtest_main
  VocalTractLabBackend
)

include(GoogleTest)
gtest_discover_tests(VtlApiTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
gtest_discover_tests(VtlBackendTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
