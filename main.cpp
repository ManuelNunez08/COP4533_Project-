#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;


/*This brute force algorithm find the most profitable transaction possible for m companies and their respective stock
 * prices over n days. Given a single transaction only involves one company, the algorithm considers a transaction for each
 * company with a purchase of stock at day x and a sale of said stock at day y, where 0 <= x < y <= n. After finding the
 * maximum profit for a transaction pertaining to a given company's stock, the profit is compared to maximum profit
 * found for all other companies considered beforehand and updated if necessary. The algorithm has a complexity of O(m*n^2)*/
vector<int> task1(vector<vector<int>>& A) {
    int m = A.size();
    int n = A[0].size();
    /*initialize initial maximum profit to the lowest possible stock price minus the maximum possible stock price
      This is given by the assignment specs as 0 - 10^5 */
    int max_profit = INT_MIN;
    // Create an array to hold output data: Stock index, buy date, and sell date
    vector<int> best_transaction(3);
    //for every company
    for (int i = 0; i < m; i++) {
        // for every day
        for (int j = 0; j < n; j++) {
            // for every possible day of sale after the day i being considered
            for (int k = j + 1; k < n; k++) {
                int profit = A[i][k] - A[i][j];
                // update maximum profit for a given company if better combination of days is found
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


/*The algorithm below finds the most profitable transaction possible for m companies and their respective stock
 * prices over n days. The algorithm takes a greedy approach, and consequently decreases the solution's complexity,
 * by recognizing that a transaction which yields maximum profit always employs the buy date with the minimum stock price
 * before the sell date in question. This leads to a variable which tracks and updates the minimum price as it's encountered
 * while parsing the stock prices for a given company. For each company, a profit is then only considered once for every day
 * by only comparing its value to the value of the minimum stock price before it. This exploit allows for a complexity of O(m*n). */

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


/*Given a 1 dimensional array A holding the prices of stock over n days for a given company and a series of other parameters
 * this function finds the minimum price of a stock preceding the day specified by parameter i. Additionally,
 * the function updates a parameter "delta" to always hold the distance between index i and the index holding its
 * respective minimum (the stock preceding it of a minimum price). Using both delta and its returned value, the function
 * ultimately updates the table B which task 3A's Dynamic Programming Solution is based on. Table B is of the same size
 * as the 2 dimensional array passed into task 3A holding stock prices over a series of days and holds, in each entry [i][j],
 * a pair of the form (minimum stock price prior to j for company i, distance to minimum prior stock price).
 * */
int findMin( vector<int>& A, vector<pair <int, int>>& B, int i, int &delta) {
    // if an entry in B holding the sought after minimum is yet to be filled out.
    if (B[i].first == -1){
        // if we are considering the first index, no prior minimum exists so B values are invalid.
        if (i == 0) {
            B[i].first = INT_MAX;
            B[i].second = -1;
            delta = 1; //update delta for A[1] to have a distance of 1 from its minimum from A[0]
            return A[i]; //return A[i] since A[1] will always have a minimum of A[0] before it.
        }
        /* computes both potential minimums before index i:
            1. other: what the minimum is for i-1 (considering values before i -1)
            2. current: considering the value at i -1*/
        int other = findMin(A, B, i - 1, delta);
        int current  = A[i - 1];
        // new minimum is found
        if (current < other){
            delta = 1; //reset delta
            B[i].first = current;
            B[i].second = delta;
            delta++;
            return current;
        }
        // continue with current minimum
        else{
            B[i].first = other;
            B[i].second = delta;
            delta ++;
            return other;
        }
    }
    //if the minimum has already been computed O(1) operation.
    else {
        return B[i].first;
    }
}


/*findMax is also a part of Task 3A. It is given a 1 dimensional array A holding the prices of stock
 * over n days for a given company. It is also passed references to buy and sell dates which it updates as it finds new
 * transactions that maximize profit as well as a parameter i specifying the day being considered. The remainder of
 * parameters only serve to call the findMin function. findMax ultimately compares the profit of selling at i to the profit
 * at selling at day i - 1. It does so using the findMin function. The sell and buy date are updated as new max profits are found. */
int findMax( vector<int>& A,  vector<pair <int, int>>& B, int i,int &delta, int &buy, int &sell) {
    //selling at day 0 is invalid so a minimum value is given such that its never chosen in a max comparison.
    if (i == 0){
        return INT_MIN;
    }
    else {
        /* compute both potential maximum profits:
         * 1. current = selling at i with the minimum price before i.
         * 2. other = selling at an index j before i with the minimum price before j.
         */
        int current = A[i] - findMin(A, B, i , delta);
        int other = findMax(A, B, i - 1, delta, buy, sell);
        // update sell and buy if sell date yielding maximum profit is found
        if (current > other && B[i].second != 0){
            sell = i;
            buy = i - B[i].second;
            return current;
        }
        else {
            return other;
        }
    }
}

/*The algorithm below finds the most profitable transaction possible for m companies and their respective stock
 * prices over n days using a top-down dynamic programming approach. We first present the following definitions:
 *      - We define, for a company m and a day n, min(m,n) as the minimum stock price for m at a day < n.
 *      - We define, for a company m and a day n, delta(m,n) as the index distance to the minimum stock price for m at a day < n.
 *
 * The following algorithm recursively computes the maximum profit of a given transaction while building a table B with entries
 * of the form (min, delta) for each m and n. Using the table the buy day corresponding to the maximum profit found can be outputted.
 *
 *  A recursive formulation to the problem is given in the assignment document.
 *
*/
vector<int> task3A(vector<vector<int> >& A) {
    int m = A.size();
    vector<int> best_transaction(3);
    int maxProfit = INT_MIN;
    //temporary buy and sell days
    int buy  = 0;
    int sell  = 0;
    //variable to keep track of the distance of an index from its minimum.
    int delta = 0;
    for (int i = 0; i < m; i++){
        delta = 0;
        vector<pair <int, int>>Empty (A[0].size(), {-1,0});
        // find the maximum profit for a given company
        int profit = findMax(A[i], Empty, A[0].size() - 1, delta, buy, sell);
        // if profit is greater than for other companies considered update the values off best transaction.
        if (profit > maxProfit){
            best_transaction[0] = i;
            best_transaction[1] = buy;
            best_transaction[2] = sell;
            maxProfit = profit;
        }
    }
    return best_transaction;
}



/*The algorithm below finds the most profitable transaction possible for m companies and their respective stock
* prices over n days using a bottom up dynamic programming approach. Unlike other algorithms above, the algorithm
 * considers the maximum price of a stock and only computes profits for sell dates corresponding to the maximum stock prices.
 * This is done through preprocessing . The algorithm has a complexity of O(m*n)*/
vector<int> task3B(vector<vector<int> >& A) {
    //PREPROCESSING
    /*Create a secondary 2 dimensional array where each entry stores the max value to the right of it. The array can be
     * built in O(n^2) time by starting on the last index of each subarray and considering for each element either the
     * element to its right or the element that its right element found the largest to its right. */
    int m = A.size();
    int n = A[0].size();

    vector<vector<int> > B(m, vector<int>(n));
    for (int i = 0; i < m; i++) {
        for (int j = n-1; j >= 0 ; j--) {
            if(j == n - 1){
                B[i][j] = -1;
            }
            else if ( A[i][j+1] > B[i][j+1] ) {
                B[i][j] = A[i][j + 1];
            }
            else {
                B[i][j] = B[i][j + 1];
            }
        }
    }
    // IMPLEMENTATION
    int stockIndex = 0;
    int buyDay = 0;
    int sellDay = 0;
    int maxProfit = INT_MIN;
    vector<int> best_transaction(3);

    for (int i = 0; i < m; i++) { // Iterating through each stock
        // Keeping track of the minimum and maximum price of each stock
        int soldAt = -1;
        for (int j = 0; j < n; j++) {
            //identify date sold at
            if (A[i][j] == soldAt){
                sellDay = j;
            }
            // consider buying in every day except the last
            if (j < n-1) {
                // Iterating through each day of each stock
                int profit = B[i][j] - A[i][j];
                if (profit > maxProfit) {
                    maxProfit = profit;
                    stockIndex = i;
                    buyDay = j;
                    soldAt = B[i][j];
                    //sell day has to be found by comparing all entries to the
                }
            }
        }
    }
    best_transaction[0] = stockIndex;
    best_transaction[1] = buyDay;
    best_transaction[2] = sellDay;
    return best_transaction;
}

//vector<int> task4(vector<vector<int> >& A, int k) {
//    /*For each day,
//     *  if we are currently holding stock
//     *      if we sell it
//     *          buy another or don't
//     *      if we don't sell it
//     *  if we are not currently holding stock
//     *     if we buy stock (k-1)
//     *     if we don't*/
//
//    int numTransactions  = k;
//
//    for(int i  = 0; i < A[0].size(); i++){
//        for(int j = 0; j < A.size(); j++){
//            A[j][i]
//        }
//    }
//
//
//
//}
//
//int task4Recurse(vector<vector<int> >& A,int day, int owning, int k) {
//    if (day  == 0) {
//        if (owning != -1){
//            return A[owning][day];
//        }
//    }
//    else {
//        for(int i  = 0; i < A.size(); i ++){
//            if (owning != -1){
//                return A[owning][day];
//            }
//        }
//    }
//
//}


bool Equal(vector<int>& A, vector<int>& B){
    if (A.size() != B.size()){
        return false;
    }

    for (int i = 0; i < A.size(); i++){
        if (A[i] != B[i]) {return false;}
    }

    return true;
}



int main() {

    vector<vector<int>> A(5, vector<int>(5));
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < A.size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                A[j][k] = rand() % 1000;
            }
        }
        vector<int> task1Result = task1(A);
        vector<int> task2Result = task2(A);
        vector<int> task3aResult = task3A(A);
        vector<int> task3bResult = task3B(A);
        if (!Equal(task1Result, task2Result) || !Equal(task1Result, task3aResult) ){
            for (int x = 0; x < A.size(); x++) {
                for (int y = 0; y < A[0].size(); y++) {
                    cout << A[x][y] << " ";
                }
                cout << endl;
            }
            cout << endl;

            cout << "Task 1 Result: " << task1Result[0] << " " << task1Result[1] << " " << task1Result[2] << endl;

            cout << "Task 2 Result: " << task2Result[0] << " " << task2Result[1] << " " << task2Result[2] << endl;

            cout << "Task 3a Result: " << task3aResult[0] << " " << task3aResult[1] << " " << task3aResult[2] << endl;

            cout << "Task 3b Result: " << task3bResult[0] << " " << task3bResult[1] << " " << task3bResult[2] << endl;

            cout << endl;
        }
        else {
            cout << "Passed Test #" << i << endl;
        }




    }



//    vector<vector<int>> A = {{5,1,3,9,7}};
//    for (int i = 0; i < A.size(); i++) {
//        for (int j = 0; j < A[i].size(); j++) {
//            cout << A[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    vector<int> task1Result = task1(A);
//
//    cout << "Task 1 Result: " << task1Result[0] << " " << task1Result[1] << " " << task1Result[2] << endl;
//
//    vector<int> task2Result = task2(A);
//
//    cout << "Task 2 Result: " << task2Result[0] << " " << task2Result[1] << " " << task2Result[2] << endl;
//
//    vector<int> task3aResult = task3A(A);
//
//    cout << "Task 3a Result: " << task3aResult[0] << " " << task3aResult[1] << " " << task3aResult[2] << endl;
//
//    vector<int> task3bResult = task3B(A);
//
//
//    cout << "Task 3b Result: " << task3bResult[0] << " " << task3bResult[1] << " " << task3bResult[2] << endl;
//    cout << endl;





    return 0;
}
