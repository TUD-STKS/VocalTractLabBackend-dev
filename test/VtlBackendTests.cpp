#include <gtest/gtest.h>

#include <iostream>

#include "VocalTractLabBackend/Sampa.h"

TEST(BackendTest, Sampa) 
{
    std::cout << "I am a Test!" << std::endl;
    EXPECT_FALSE(Sampa::isVowel("C"));
    EXPECT_TRUE(Sampa::isVowel("a"));
}