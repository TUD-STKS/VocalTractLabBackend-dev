# CMakeList.txt : CMake project for VocalTractLabApi, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)

project("VocalTractLabApi")
add_library(VocalTractLabApi SHARED)
set_target_properties(VocalTractLabApi PROPERTIES PUBLIC_HEADER include/VocalTractLabApi/VocalTractLabApi.h)
target_include_directories (VocalTractLabApi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/VocalTractLabApi)
target_include_directories (VocalTractLabApi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/Backend)
target_sources(VocalTractLabApi PRIVATE src/VocalTractLabApi.def src/AnatomyParams.cpp src/Dsp.cpp src/F0EstimatorYin.cpp src/GeometricGlottis.cpp src/Geometry.cpp src/GesturalScore.cpp src/Glottis.cpp src/IirFilter.cpp src/ImpulseExcitation.cpp src/LfPulse.cpp src/Matrix2x2.cpp src/PoleZeroPlan.cpp src/Sampa.cpp src/SegmentSequence.cpp src/Signal.cpp src/Splines.cpp src/StaticPhone.cpp src/Surface.cpp src/Synthesizer.cpp src/TdsModel.cpp src/TimeFunction.cpp src/TlModel.cpp src/TriangularGlottis.cpp src/Tube.cpp src/TwoMassModel.cpp src/VocalTract.cpp src/VocalTractLabApi.cpp src/VoiceQualityEstimator.cpp src/VowelLf.cpp src/XmlHelper.cpp src/XmlNode.cpp)


project("VocalTractLabBackend")
add_library(VocalTractLabBackend STATIC)
target_include_directories (VocalTractLabBackend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/Backend)
target_sources(VocalTractLabBackend PUBLIC src/AnatomyParams.cpp src/Dsp.cpp src/F0EstimatorYin.cpp src/GeometricGlottis.cpp src/Geometry.cpp src/GesturalScore.cpp src/Glottis.cpp src/IirFilter.cpp src/ImpulseExcitation.cpp src/LfPulse.cpp src/Matrix2x2.cpp src/PoleZeroPlan.cpp src/Sampa.cpp src/SegmentSequence.cpp src/Signal.cpp src/Splines.cpp src/StaticPhone.cpp src/Surface.cpp src/Synthesizer.cpp src/TdsModel.cpp src/TimeFunction.cpp src/TlModel.cpp src/TriangularGlottis.cpp src/Tube.cpp src/TwoMassModel.cpp src/VocalTract.cpp src/VoiceQualityEstimator.cpp src/VowelLf.cpp src/XmlHelper.cpp src/XmlNode.cpp)



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

target_link_libraries(
  VtlApiTests
  gtest_main
  VocalTractLabApi
)

add_executable(
  "VtlBackendTests"
  "test/VtlBackendTests.cpp"
 )

target_link_libraries(
  VtlBackendTests
  gtest_main
  VocalTractLabBackend
)

include(GoogleTest)
gtest_discover_tests(VtlApiTests)
gtest_discover_tests(VtlBackendTests)