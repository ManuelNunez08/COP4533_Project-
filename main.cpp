#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>
using namespace std;
using namespace std::chrono;


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
        for (int j = startCol; j < n; j++) { // Buyday loop for transactions
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

//int bestProfit(vector<vector<pair<int, int>>>& T, vector<int>& recorded, vector<vector<int>>& optimal, int i, int k ){
//    if (i == 0 || k == 0 ){
//        return 0;
//    }
//    else if (recorded[i] != -1){
//        return max(recorded[i], bestProfit(T, recorded, optimal, i - 1, k));
//    }
//    else {
//        int maxProfit = 0;
//        optimal.push_back({-1,-1,-1});
//        for (int j = i - 1; j >= 0; j--) {
//            int profit = T[j][i].second + bestProfit(T, recorded,optimal,  j, k - 1);
//            if (profit > maxProfit) {
//                maxProfit = profit;
//                optimal[optimal.size()-1][0] = j;
//                optimal[optimal.size()-1][1] = i;
//                optimal[optimal.size()-1][2] = T[j][i].first;
//            }
//        }
//        recorded[i] = maxProfit; // record the maximum profit of a transaction with a sell day of i
//        int next = bestProfit(T, recorded, optimal, i - 1, k); // calculate the maximum profit of a transaction with a sell day of i -1
//
//        if (maxProfit <= next){
//            optimal.pop_back();
//            return next;
//        }
//        else{
//            return maxProfit;
//        }
//    }
//}
int bestProfit(vector<vector<pair<int, int>>>& T, vector<int>& recorded, int i, int k ){
    if (i == 0 || k == 0 ){
        return 0;
    }
    else if (recorded[i] != -1){
        return max(recorded[i], bestProfit(T, recorded,  i - 1, k));
    }
    else {
        int maxProfit = 0;
        for (int j = i - 1; j >= 0; j--) {
            int profit = T[j][i].second + bestProfit(T, recorded,  j, k - 1);
            if (profit > maxProfit) {
                maxProfit = profit;
            }
        }
        recorded[i] = maxProfit; // record the maximum profit of a transaction with a sell day of i
        return max(maxProfit, bestProfit(T, recorded, i-1, k));
    }
}
/* /*The algorithm below finds the most profitable combination of at most k transactions for m companies and their respective
 * stock prices over n days using a dynamic programming approach that incorporates both bottom up and top down techniques.
 * The algorithm first computes an nxn table and stores for each entry [i][j] the most profitable transaction starting at day i and
 * ending at day j (where i < j). It then uses recursion to traverse the nxn array to find the maximum profit possible
 * from a max of k transactions by summing entries whose respective i and j indexes don't overlap with one another since
 * only pone piece of stock can be held at once.
 *
 * The preprocessing element of the algorithm has a complexity of O(m*n*n)
 * Using memoization, the recursive aspect of the algorithm has a complexity of O(n*n)
 *
 * The total time complexity is of the order O(m*n*n). */
int task6(vector<vector<int> >& A, int k) {
    int companies = A.size();
    int days = A[0].size();
    vector<vector<int>> trans;
    /* The 2 dimensional table vector contains pairs storing the max profit for a transaction
     * started (buy) on day i and finished (sell) on day j as well as the company through which it's generated  */
    vector<vector<pair<int, int>>> table(days, vector<pair<int, int>>(days,{-1, -1}));
    for (int i  = days - 1; i > 0; i--){
        for (int k = i -1; k >= 0; k --){
            vector<int> tran = {-1,-1,-1,INT_MIN};
            for (int j  = 0; j < companies; j++){
                int profit = A[j][i] - A[j][k];
                // update a transaction covering a given interval if a better transaction covering that interval is found.
                if ( profit > tran[3]) {
                    tran[0] = k;
                    tran[1] = i;
                    tran[2] = j;
                    tran[3] = profit;
                }
            }
            trans.push_back(tran);
            table[tran[0]][tran[1]] = {tran[2], tran[3]};
        }
    }

    for (int x = 0; x < trans.size(); x++) {
        cout << "Buy: " << trans[x][0] << " | Sell: " << trans[x][1] <<  " | Company: " << trans[x][2] <<  " | Profit:" << trans[x][3];
        cout << endl;
    }
    cout << endl;


    for (int x = 0; x < table.size(); x++) {
        for (int y  = 0; y < table[0].size(); y++  ) {
            cout << "( " << table[x][y].first << ", " << table[x][y].second << " )";
        }
        cout << endl;
    }
    cout << endl;

    // create a vector to store computed recursive values. Memoization is incorporated into the solution through this vector.
    vector<int>recorded(days - 1, -1);
    //vector to store transactions
    vector<vector<int>> optimal;
    int max =  bestProfit(table, recorded, table.size()-1, k);
    cout << "transactions: " << endl;
    for (int x = 0; x < optimal.size(); x++) {
        for (int y  = 0; y < optimal[0].size(); y++  ) {
            cout << optimal[x][y] << ", " ;
        }
        cout << endl;
    }
    cout << endl;
    cout << "Max Profit: " << max << endl;

    return max;



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
    //vector<vector<int>> A = {{1, 9, 8 ,7,6}, {10, 6, 5 ,2, 1} };
    vector<vector<int>> A(5, vector<int>(5));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < A.size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                A[j][k] = rand() % 10;
            }
        }
        // Formatting for the Matrix for readability
        cout << "  ";
        for(int i = 0; i < A.size();i++){
            cout << " " << "\033[4m" << i << "\033[0m";
        }
        cout << endl;
        for (int i = 0; i < A.size(); i++) {
            cout << i << "| ";
            for (int j = 0; j < A[0].size(); j++) {
                cout << A[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        task6(A, 2);
        cout << endl;
        vector<vector<int>> task4Result = task4K2(A);
        vector<vector<int>> task4bResult = task4AnyK(A, 2, 0);
        int profit = 0;
        cout << "Task 4K2 Result: " << endl;
        for (int l = 0; l < task4Result.size(); l++) {
            cout << "Transaction " << l + 1 << ": Stock: " << task4Result[l][0] << " | BuyDay: " << task4Result[l][1] << " | SellDay:"
                 << task4Result[l][2] << endl;
            int stock = task4bResult[l][0];
            int buyDay = task4bResult[l][1];
            int sellDay = task4bResult[l][2];
            profit += A[stock][sellDay] - A[stock][buyDay];
        }
        cout << "Profit: " << profit << endl;
        cout << endl;
        cout << "Task 4AnyK Result: " << endl;
        profit = 0;
        for (int l = 0; l < task4bResult.size(); l++) {
            cout << "Transaction: " << l + 1 << ": Stock: " << task4bResult[l][0] << " | BuyDay: " << task4bResult[l][1] << " | SellDay:"
                 << task4bResult[l][2] << endl;
            int stock = task4bResult[l][0];
            int buyDay = task4bResult[l][1];
            int sellDay = task4bResult[l][2];
            profit += A[stock][sellDay] - A[stock][buyDay];
        }
        cout << "Profit: " << profit << endl;
        cout << "-----------------------------------------------------------------------" << endl;
    }
        return 0;
}













//    vector<vector<int>> A(4, vector<int>(5));
//    for (int i = 0; i < 1000; i++) {
//        for (int j = 0; j < A.size(); j++) {
//            for (int k = 0; k < A[0].size(); k++) {
//                A[j][k] = rand() % 1000;
//            }
//        }
//        vector<int> task1Result = task1(A);
//        vector<int> task2Result = task2(A);
//        vector<int> task3aResult = task3A(A);
//        vector<int> task3bResult = task3B(A);
//        if (!Equal(task1Result, task2Result) || !Equal(task1Result, task3aResult) ){
//            for (int x = 0; x < A.size(); x++) {
//                for (int y = 0; y < A[0].size(); y++) {
//                    cout << A[x][y] << " ";
//                }
//                cout << endl;
//            }
//            cout << endl;
//
//            cout << "Task 1 Result: " << task1Result[0] << " " << task1Result[1] << " " << task1Result[2] << endl;
//
//            cout << "Task 2 Result: " << task2Result[0] << " " << task2Result[1] << " " << task2Result[2] << endl;
//
//            cout << "Task 3a Result: " << task3aResult[0] << " " << task3aResult[1] << " " << task3aResult[2] << endl;
//
//            cout << "Task 3b Result: " << task3bResult[0] << " " << task3bResult[1] << " " << task3bResult[2] << endl;
//
//            cout << endl;
//        }
//        else {
//            cout << "Passed Test #" << i << endl;
//        }
//
//
//
//
//    }



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