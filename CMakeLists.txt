# CMakeList.txt : CMake project for VocalTractLabApi, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("VocalTractLabApi")

add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
add_library (VocalTractLabApi SHARED VocalTractLabApi.def AnatomyParams.cpp Dsp.cpp F0EstimatorYin.cpp GeometricGlottis.cpp Geometry.cpp GesturalScore.cpp Glottis.cpp GlottisFactory.cpp IirFilter.cpp ImpulseExcitation.cpp LfPulse.cpp Matrix2x2.cpp PoleZeroPlan.cpp PoleZeroPlan.h Sampa.cpp SegmentSequence.cpp Signal.cpp Speaker.cpp Splines.cpp StaticPhone.cpp Surface.cpp Synthesizer.cpp TdsModel.cpp TimeFunction.cpp TlModel.cpp TriangularGlottis.cpp Tube.cpp TwoMassModel.cpp VocalTract.cpp VocalTractLabApi.cpp VoiceQualityEstimator.cpp VowelLf.cpp XmlHelper.cpp XmlNode.cpp)
