#include <iostream>
#include <vector>

using namespace std;


vector<int> task1(vector<vector<int>>& A) {
    int m = A.size();
    int n = A[0].size();
    int max_profit = INT_MIN;  // Set to negative infinity to always return a transaction, even if the profit is negative
    vector<int> best_transaction(3);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                int profit = A[i][k] - A[i][j];
                if (profit > max_profit) {
                    max_profit = profit;
                    best_transaction[0] = i;
                    best_transaction[1] = j;
                    best_transaction[2] = k;
                }
            }
        }
    }
    return best_transaction;
}

vector<int> task2(vector<vector<int>>& A) {
    int m = A.size();
    int n = A[0].size();
    int maxProfit = INT_MIN;
    int stockIndex = 0;
    int buyDay = 0;
    int sellDay = 0;
    vector<int> best_transaction(3);

    for (int i = 0; i < m; i++) { // Iterating through each stock
        // Keeping track of the minimum of each stock
        int minPrice = A[i][0];
        int minDay = 0;

        for (int j = 1; j < n; j++) { // Iterating through each day of each stock
            int profit = A[i][j] - minPrice;
            if (profit > maxProfit) { // If the profit is greater than the max profit, this is the new transaction
                maxProfit = profit;
                stockIndex = i;
                buyDay = minDay;
                sellDay = j;
            }
            if (A[i][j] < minPrice) { // Checking for a new buy day
                minPrice = A[i][j];
                minDay = j;
            }
        }
    }
    best_transaction[0] = stockIndex;
    best_transaction[1] = buyDay;
    best_transaction[2] = sellDay;
    return best_transaction;
}

int main() {

//    int m, n;
//    cin >> m >> n;
//    vector<vector<int>> A(m, vector<int>(n));
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < n; j++) {
//            cin >> A[i][j];
//        }
//    }
    vector<vector<int>> A = {{6,5,4,1}};
    vector<int> task2Result = task2(A);

    cout << task2Result[0] << " " << task2Result[1] << " " << task2Result[2] << endl;

    int x[3][3] = {{1,2,3}, {3,5,9}, {2,3,4}};

    /* we are essentially finding the greatest percentage change in our initial investment of an arbitrary amount A.
     For the three days and three companies presented above if we start with am amount A invested in company 1,
     by day 2 we have an amount 3A. If we then sell and buy company 2, by day three we have an amount 3A* 4/5 = 12A/5.

     to calculate the transactions that maximize profit we have to consider at every day if we should buy stock
     (assuming we don't have any stock ), whether to sell our current stock (if we have any), and whether to buy stock
     (assuming we just sold our stock or don't have stock)


     Observations:
     1. We are not interested in the magnitude of a stock price, we are interested in its change over the days in question.

     2. A set of stock prices for a given company, such as {2,4,1,6}, can be represented through an array specifying changes in
        the stock price; {1,2,0.25, 6}, since day 2 has a stock price worth twice as much as day 1, day 3 has a stock price
        worth a fourth as much as day 2, etc. We call this array a derivative array.

    3. With one company, one set of stock prices (as listed above), and no limits on transactions. We can maximize profit by buying
        on all days n where day n + 1 has a value in the derivative array greater than 1 and selling on all day's where the
        derivative array has a value smaller than 1. We are indifferent on day's where the stock price equals 1.

    2. With a set of companies, a set of corresponding stock prices, and no limits on transactions, we can also create a
     derivative array. For the 2d-array shown above we have {{1,1,1}, {3,2.5,3}, {2/3,0.6,4/9}}. We can then maximize profit
     by buying stock for a company m on day n if they have the greatest n+1 derivative array value assuming a value greater than 1 exists.
     If no value greater than 1 exists we would buy nothing on day n.


     */




    return 0;
}
