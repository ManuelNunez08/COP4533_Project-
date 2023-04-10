#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>
using namespace std;
using namespace std::chrono;


//task 1 functions
vector<int> task1(vector<vector<int>>& A);
//task 2 functions
vector<int> task2(vector<vector<int>>& A);
//task 3A functions
int findMin( vector<int>& A, vector<pair <int, int>>& B, int i, int &delta);
int findMax( vector<int>& A,  vector<pair <int, int>>& B, int i,int &delta, int &buy, int &sell);
vector<int> task3A(vector<vector<int> >& A);
//task3B functions
vector<int> task3B(vector<vector<int> >& A);
// time complexity analysis of tasks 1-3 functions
void Plot1();
void Plot2();
//testing functions task 1 -3


//task 4 functions
vector<vector<int>> task4K2(vector<vector<int> >& A);
vector<vector<int>> task4AnyK(vector<vector<int>>& A, int k, int startCol);


//task 6 functions
vector<vector<int>> findIndices( vector<vector<vector<pair<int,int>>>>& T, vector<vector<int>>& A,vector<vector<int>>& transactions, int k, int startDay, int lastDay, int lastCompany);
int findBest(vector<vector<int> >& A, pair<int,int>& holding, vector<vector<vector<pair<int,int>>>>& recorded,  int day, int k);
vector<vector<int>> task6(vector<vector<int> >& A, int k);

//time complexity analysis for tasks 4-6


//helper/testing functions
bool Equal(vector<int>& A, vector<int>& B);
void testProblem2(int days, int companies, int k,vector<string>& methods);

int main() {

    vector<string> test2 = {"4","6"};
    testProblem2(5,3,2,test2);

    return 0;
}

bool Equal(vector<int>& A, vector<int>& B){
    if (A.size() != B.size()){
        return false;
    }

    for (int i = 0; i < A.size(); i++){
        if (A[i] != B[i]) {return false;}
    }

    return true;
}


vector<vector<int>> findIndices( vector<vector<vector<pair<int,int>>>>& T, vector<vector<int>>& A,vector<vector<int>>& transactions, int k, int startDay, int lastDay, int lastCompany){

    vector<pair<int, int>> maxesK;
    int currMax = INT_MIN;
    int max = INT_MIN;
    int sellDay = -1;
    int buyDay = -1;
    int company = -1;
    for (int i = startDay; i >= 0; i--) {
        for (int j = 0 ; j <  T.size() ; j++) {
            pair<int, int> current = T[j][i][k - 1];
            if (current.first >= max  && startDay == A[0].size() -1) {
                max = current.first;
                sellDay = i;
                buyDay  = current.second;
                company = j;
            }
            else if (startDay != A[0].size() -1 )
            {
                int lastMax = T[lastCompany][lastDay][k].first;
                int lookFor = A[lastCompany][lastDay] - A[lastCompany][startDay];
                if(current.first == lastMax - lookFor) {
                    max = current.first;
                    sellDay = i;
                    buyDay = current.second;
                    company = j;
                }
            }
        }
    }
    if (max > 0 && k > 0){
        transactions.push_back({company, buyDay, sellDay});
        if (k >= 1) {
            findIndices(T, A, transactions, k - 1, buyDay, sellDay, company);
        }
    }
    return transactions;
}
int findBest(vector<vector<int> >& A, pair<int,int>& holding, vector<vector<vector<pair<int,int>>>>& recorded,  int day, int k){
    // if day 0 is reached, or we've run out of transactions and are not currently holding stock.
    if (day == 0 || (k == 0 && holding.first == -1)){
        // if we are holding stock, it must be day 0 for us to be in this method, so sell whatever you have
        if (holding.first != -1){
            recorded[holding.first][holding.second][k].second = day;
            return (A[holding.first][holding.second] - A[holding.first][day]);
        }
            // else, either because its day 0 and we have no stock or because it's an arbitrary day, but we have no stock or
            // transactions left.
        else
        {
            return 0;
        }
    }
        // else if we are not holding anything, consider buying something or skipping the current day.
    else if(holding.first == -1){
        int maxProfit = INT_MIN;
        for(int i = 0; i < A.size(); i++){
            int profit  = INT_MIN;
            if(recorded[i][day][k-1].first != -1){
                profit = recorded[i][day][k-1].first;
            }
            else {
                pair<int, int> sold = {i, day};
                profit = findBest(A, sold, recorded, day - 1, k - 1);
                recorded[i][day][k-1].first = profit;
            }
            if(profit > maxProfit){
                maxProfit = profit;
            }
        }
        return max(maxProfit, findBest(A, holding, recorded, day-1, k));
    }
        // else if we are holding stock consider selling it, or to continue holding it.
    else{
        pair<int, int> clearHold = {-1,-1};
        int current = (A[holding.first][holding.second] - A[holding.first][day]) + findBest(A,clearHold,recorded, day, k);
        int other = findBest(A, holding,recorded, day-1, k) ;
        if (current > other){
            recorded[holding.first][holding.second][k].second = day;
        }
        return max(current, other);
    }
}

vector<vector<int>> task6(vector<vector<int>>& A, int k){
    vector<vector<vector<pair<int,int>>>>recorded(A.size() , vector<vector<pair<int,int>>>(A[0].size(), vector<pair<int,int>>(k, {-1,-1})));
    pair<int, int> holding = {-1, -1};

    int max  = findBest(A, holding, recorded,  A[0].size() - 1, k);
    for (int i = 0; i < recorded.size(); i++) {
        for (int j = 0; j < recorded[0].size(); j++) {
            cout << "(";
            for (int k = 0; k < recorded[0][0].size(); k++) {
                cout << "(" << recorded[i][j][k].first << ", " << recorded[i][j][k].second << ")";
            }
            cout << ")      ";
        }
        cout << endl;
    }
    cout << "Task 6 profit :" << max << endl;
    vector<vector<int>> transactions;
    return findIndices(recorded,A,transactions, k,A[0].size()-1,-1,-1);
}









void testProblem2(int days, int companies, int k, vector<string>& methods){
    vector<vector<int>> A(companies, vector<int>(days));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < companies; j++) {
            for (int k = 0; k < days; k++) {
                A[j][k] = rand() % 10;
            }
        }
//        // Formatting for the Matrix for readability
//        cout << "  ";
//        for(int i = 0; i < days;i++){
//            cout << " " << "\033[4m" << i << "\033[0m";
//        }
//        cout << endl;
//        for (int i = 0; i < companies; i++) {
//            cout << i << "| ";
//            for (int j = 0; j < days; j++) {
//                cout << A[i][j] << " ";
//            }
//            cout << endl;
//        }
//        cout << endl;
//        cout << endl;

        if (methods.begin(), methods.end(), "4") {
            vector<vector<int>> task4bResult = task4AnyK(A, k, 0);
            int profit = 0;
            cout << endl;
            cout << "Task 4 Result: " << endl;
            for (int l = 0; l < task4bResult.size(); l++) {
                int stock = task4bResult[l][0];
                int buyDay = task4bResult[l][1];
                int sellDay = task4bResult[l][2];
                cout << stock << " " << buyDay << " " << sellDay << endl;
                profit += A[stock][sellDay] - A[stock][buyDay];
            }
            cout << "Profit: " << profit << endl << endl;
        }

        if (methods.begin(), methods.end(), "6") {
            vector<vector<int>> task6Result = task6(A, k);
            cout << endl;
            cout << "Task 6 Result: " << endl;
            int profit = 0;
            for (int l = 0; l < task6Result.size(); l++) {
                int stock = task6Result[l][0];
                int buyDay = task6Result[l][1];
                int sellDay = task6Result[l][2];
                cout << stock << " " << buyDay << " " << sellDay << endl;
                profit += A[stock][sellDay] - A[stock][buyDay];
            }
            cout << "Profit: " << profit << endl;
        }

        cout << "-----------------------------------------------------------------------" << endl;
    }
}











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
        for (int j = 0; j < n-1; j++) {
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
    vector<int> best_transaction(3);

    for (int i = 0; i < m; i++) { // Iterating through each stock
        // Keeping track of the minimum of each stock
        int minPrice = A[i][0];
        int minDay = 0;

        for (int j = 1; j < n; j++) { // Iterating through each day of each stock
            int profit = A[i][j] - minPrice;
            if (profit > maxProfit) { // If the profit is greater than the max profit, this is the new transaction
                maxProfit = profit;
                best_transaction[0] = i;
                best_transaction[1] = minDay;
                best_transaction[2] = j;
            }
            if (A[i][j] < minPrice) { // Checking for a new buy day
                minPrice = A[i][j];
                minDay = j;
            }
        }
    }
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
    /* Create a secondary 2 dimensional array where each entry stores the max value to the right of it. The array can be
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
vector<vector<int>> task4K2(vector<vector<int> >& A) {
    int m = A.size();
    int n = A[0].size();
    int max_profit = 0;
    vector<vector<int>> best_transactions(2,vector<int>(3, -1)); // best_transactions 2d array with k transactions.

    // ATTEMPT AT K=2
    for (int i = 0; i < m; i++) { // Stock loop for FIRST transaction
        for (int j = 0; j < n; j++) { // Buyday loop for FIRST transaction
            for (int k = j + 1; k < n; k++) { // SellDay loop for FIRST transaction
                for(int z = 0; z < m; z++) { // Stock loop for SECOND transaction
                    for (int l = k; l < n; l++) { //Buy day loop for SECOND transaction
                        for (int o = l + 1; o < n; o++) { // SellDay loop for SECOND transaction
                            int profit1 = A[i][k] - A[i][j]; // Transaction 1: i j k
                            int profit2 = A[z][o] - A[z][l]; // Transaction 2: z l o
                            int total_profit = profit1 + profit2;
                            if (total_profit > max_profit) {
                                max_profit = total_profit;
                                best_transactions[0][0] = i;
                                best_transactions[0][1] = j;
                                best_transactions[0][2] = k;
                                best_transactions[1][0] = z;
                                best_transactions[1][1] = l;
                                best_transactions[1][2] = o;
                            }
                        }
                    }
                }
                int profit4 = A[i][k] - A[i][j]; // Transaction 1: i j k
                // One transaction
                if (profit4 > max_profit) {
                    max_profit = profit4;
                    best_transactions[0][0] = i;
                    best_transactions[0][1] = j;
                    best_transactions[0][2] = k;
                    best_transactions[1][0] = -1;
                    best_transactions[1][1] = -1;
                    best_transactions[1][2] = -1;
                }
            }
        }
    }
    return best_transactions;
}
vector<vector<int>> task4AnyK(vector<vector<int>>& A, int k, int startCol) {
    int m = A.size();
    int n = A[0].size();
    int max_profit = 0;
    vector<vector<int>> best_transactions; // best_transactions 2d array with k transactions.
    vector<vector<int>> tempBest_transactions; // temporary best_transactions 2d array that will hold the transactions returned by the recursive calls

    for (int i = 0; i < m; i++) { // Stock loop for transactions
        for (int j = startCol; j < n-1; j++) { // Buyday loop for transactions
            for (int l = j + 1; l < n; l++) { // SellDay loop for transactions
                if(k > 1){
                    tempBest_transactions = task4AnyK(A, k-1, l);
                }
                int profit1 = A[i][l] - A[i][j]; // Transaction 1: Stock: i BuyDay: j SellDay: l
                int profit2 = 0; // Max profit from the best transactions from the recursive call;
                for(int z = 0; z<tempBest_transactions.size(); z++){ // Looping through best transactions from recursive call to add their profits
                    int stock = tempBest_transactions[z][0];
                    int buyDay = tempBest_transactions[z][1];
                    int sellDay = tempBest_transactions[z][2];
                    profit2 += A[stock][sellDay] - A[stock][buyDay];
                }
                int total_profit = profit1 + profit2;
                if (total_profit > max_profit) {
                    max_profit = total_profit;
                    best_transactions.clear(); // Clears the best transactions 2d array because we found new better transactions
                    best_transactions.push_back({i,j,l}); // Adding the transaction from this call to the beginning of the best_transactions list;
                    for(int z = 0; z<tempBest_transactions.size(); z++){
                        best_transactions.push_back(tempBest_transactions[z]);  // Adding the transactions from the recursive call to the best list
                    }
                }
            }
        }
    }
    return best_transactions;
}

void Plot1() {

    vector<int> days = {1000, 2000, 3000, 4000, 5000};
    vector<int> companies = {1000, 2000, 3000, 4000, 5000};

    //PLOT 1: variable days (n) and fixed companies(m = 2000)
    vector<vector<int>> Plot1 (days.size(), vector<int>(4));
    for (int i = 0; i < days.size(); i++ ){
        int numDays = days[i];
        vector<vector<int>> A(2000, vector<int>(numDays));
        for (int j = 0; j < A.size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                A[j][k] = rand() % 10000;
            }
        }

        //Task 1
        auto startTask1 = high_resolution_clock::now();
        task1(A);
        auto stopTask1 = high_resolution_clock::now();
        auto durationTask1 = duration_cast<microseconds>(stopTask1 - startTask1);
        Plot1[i][0] = durationTask1.count();

        //Task 2
        auto startTask2 = high_resolution_clock::now();
        task2(A);
        auto stopTask2 = high_resolution_clock::now();
        auto durationTask2 = duration_cast<microseconds>(stopTask2 - startTask2);
        Plot1[i][1] = durationTask2.count();

        //Task 3A
        auto startTask3A = high_resolution_clock::now();
        task3A(A);
        auto stopTask3A = high_resolution_clock::now();
        auto durationTask3A = duration_cast<microseconds>(stopTask3A - startTask3A);
        Plot1[i][2] = durationTask3A.count();

        //Task 3B
        auto startTask3B = high_resolution_clock::now();
        task3B(A);
        auto stopTask3B = high_resolution_clock::now();
        auto durationTask3B = duration_cast<microseconds>(stopTask3B - startTask3B);
        Plot1[i][3] = durationTask3B.count();

    }

    for (int j = 0; j < Plot1.size(); j++) {
        for (int k = 0; k <Plot1[0].size(); k++) {
            cout << Plot1[j][k] << " ";
        }
        cout << endl;
    }

}

void Plot2() {

    vector<int> days = {1000, 2000, 3000, 4000, 5000};
    vector<int> companies = {1000, 2000, 3000, 4000, 5000};

    //PLOT 2: fixed days (n = 2000) and variable companies
    vector<vector<int>> Plot2 (companies.size(), vector<int>(4));
    for (int i = 0; i < companies.size(); i++ ){
        int numCompanies = companies[i];
        vector<vector<int>> A(numCompanies, vector<int>(2000));
        for (int j = 0; j < A.size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                A[j][k] = rand() % 10000;
            }
        }

        //Task 1
        auto startTask1 = high_resolution_clock::now();
        task1(A);
        auto stopTask1 = high_resolution_clock::now();
        auto durationTask1 = duration_cast<microseconds>(stopTask1 - startTask1);
        Plot2[i][0] = durationTask1.count();

        //Task 2
        auto startTask2 = high_resolution_clock::now();
        task2(A);
        auto stopTask2 = high_resolution_clock::now();
        auto durationTask2 = duration_cast<microseconds>(stopTask2 - startTask2);
        Plot2[i][1] = durationTask2.count();

        //Task 3A
        auto startTask3A = high_resolution_clock::now();
        task3A(A);
        auto stopTask3A = high_resolution_clock::now();
        auto durationTask3A = duration_cast<microseconds>(stopTask3A - startTask3A);
        Plot2[i][2] = durationTask3A.count();

        //Task 3B
        auto startTask3B = high_resolution_clock::now();
        task3B(A);
        auto stopTask3B = high_resolution_clock::now();
        auto durationTask3B = duration_cast<microseconds>(stopTask3B - startTask3B);
        Plot2[i][3] = durationTask3B.count();

    }

    for (int j = 0; j < Plot2.size(); j++) {
        for (int k = 0; k <Plot2[0].size(); k++) {
            cout << Plot2[j][k] << " ";
        }
        cout << endl;
    }

}


//    for (int x = 0; x < trans.size(); x++) {
//        cout << "Buy: " << trans[x][0] << " | Sell: " << trans[x][1] <<  " | Company: " << trans[x][2] <<  " | Profit:" << trans[x][3];
//        cout << endl;
//    }
//    cout << endl;
//
//
//    for (int x = 0; x < table.size(); x++) {
//        for (int y  = 0; y < table[0].size(); y++  ) {
//            cout << "( " << table[x][y].first << ", " << table[x][y].second << " )";
//        }
//        cout << endl;
//    }
//    cout << endl;

// create a vector to store computed recursive values. Memoization is incorporated into the solution through this vector.
//vector<vector<vector<int>>>recorded(days , vector<vector<int>>(days, vector<int>(k, -1)));
//vector to store transactions
//    cout << "count: " << count << endl;


//    cout << "table: " << endl;
//    for (int x = 0; x < recorded.size(); x++) {
//        for (int y  = 0; y < recorded[0].size(); y++  ) {
//            cout << recorded[x][y] << ", " ;
//        }
//        cout << endl;
//    }
//    cout << endl;
//cout << "Max Profit: " << max << endl;

// vector<pair<int, int>> maxesK;
//    int currMax = INT_MIN;
//    int max = INT_MIN;
//    int day = -1;
//    int company = -1;
//    int considerBefore = T[0].size();
//    for (int tran = k -1; tran >= 0; tran--) {
//        for (int i = 0; i < T.size(); i++) {
//            for (int j = 0; j < considerBefore; j++) {
//               if (T[i][j][tran] > max){
//                   max = T[i][j][tran];
//                   day = j;
//                   company = i;
//               }
//            }
//        }
//
//        if (max > 0) {
//            maxesK.push_back({day, company});
//        }
//        considerBefore = day;
//        max = INT_MIN;
//    }
//
//
//    for (int i = 0; i < maxesK.size(); i++){
//        cout << "(" << maxesK[i].first << ", " << maxesK[i].second <<  ")" << endl;
//    }
//
//    vector<vector<int>> transactions;
//
//    for (int i = maxesK.size() - 1; i >= 0; i--) {
//        int lookUntil;
//        if (i == maxesK.size() - 1) {
//            lookUntil = 0;
//        } else {
//            lookUntil = maxesK[i + 1].first;
//        }
//        int min = INT_MAX;
//        int buy = -1;
//        for (int j = maxesK[i].first - 1; j >= lookUntil; j--) {
//            if (A[maxesK[i].second][j] < min) {
//                min = A[maxesK[i].second][j];
//                buy = j;
//            }
//        }
//        transactions.push_back({ maxesK[i].second, buy, maxesK[i].first});
//    }
//
//    for (int i = 0; i < transactions.size(); i++){
//        cout << transactions[i][0] << " " << transactions[i][1] << " " << transactions[i][2] << endl;
//    }