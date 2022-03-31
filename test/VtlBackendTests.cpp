#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include "VocalTractLabBackend/GeometricGlottis.h"
#include "VocalTractLabBackend/GesturalScoreWithHistory.h"

#include "VocalTractLabBackend/Sampa.h"

TEST(BackendTest, GesturalScoreWithHistory)
{
    const auto vt = make_unique<VocalTract>();
    const auto g = make_unique<GeometricGlottis>();

	GesturalScoreWithHistory gs(vt.get(), g.get());
    gs.initTestScore();

    gs.addClosingGesture(GesturalScore::VOWEL_GESTURE, "a", 0, 2, true);


    
    std::cout << "End of test." << std::endl;
}

TEST(BackendTest, Sampa) 
{
    std::cout << "I am a Test!" << std::endl;
    EXPECT_FALSE(Sampa::isVowel("C"));
    EXPECT_TRUE(Sampa::isVowel("a"));
}