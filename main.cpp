#include <iostream>

int main() {
    //std::cout << "Hello, World!" << std::endl;
    int x[3][3] = {{1,2,3}, {3,5,8}, {2,3,4}};

    /* we are essentially finding the greatest percentage change in our initial investment of an arbitrary amount A.
     For the three days and three companies presented above if we start with am amount A invested in company 1,
     by day 2 we have an amount 3A. If we then sell and buy company 2, by day three we have an amount 3A* 4/5 = 12A/5.

     to calculate the transactions that maximize profit we have to consider at every day if we should buy stock
     (assuming we don't have any stock ), whether to sell our current stock (if we have any), and whether to buy stock
     (assuming we just sold our stock or don't have stock)*/


    return 0;
}
