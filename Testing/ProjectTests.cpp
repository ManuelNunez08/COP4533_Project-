
#include "../main.cpp"
#include "gtest/gtest.h"
#define TEST(test_suite_name, test_name) GTEST_TEST(test_suite_name, test_name)

namespace {
    void compareArrays(vector<int>& A, vector<int>& B) {
        ASSERT_EQ(A.size(), B.size()) << "Vectors x and y are of unequal length";

        for (int i = 0; i < A.size(); ++i) {
            EXPECT_EQ(A[i], B[i]) << "Vectors x and y differ at index " << i;
        }
    }

    TEST(Arbitrary, oneTransaction){
        vector<int> A = {1,0,1};
        vector<vector<int>> x = {{1,2},{3,5}};
        vector<int> B = task1(x);
        compareArrays(A, B);

    }

}



int runTests (){
    return RUN_ALL_TESTS();
}




