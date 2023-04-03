
#include "../main.cpp"
#include "gtest/gtest.h"
#define TEST(test_suite_name, test_name) GTEST_TEST(test_suite_name, test_name)

namespace {

    TEST(Arbitrary, oneTransaction){EXPECT_EQ(task1({{1,2},{3,5}}), {1, 0, 1})



}




