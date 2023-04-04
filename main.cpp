#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

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


int findMin( vector<int>& A, vector<pair <int, int>>& B, int i, int &delta, int &buy, int &sell) {

    if (B[i].first == -1){

        if (i == 0) {
            B[i].first = INT_MAX;
            B[i].second = -1;
            delta = 1;
            return A[i];
        }

        int other = findMin(A, B, i - 1, delta, buy, sell);
        int current  = A[i - 1];


        if (current < other){
            delta = 1;
            B[i].first = current;
            B[i].second = delta;
            delta++;
            return current;
        }
        else{
            B[i].first = other;
            B[i].second = delta;
            delta ++;
            return other;
        }

    }
    else {
        return B[i].first;
    }
}



int findMax( vector<int>& A,  vector<pair <int, int>>& B, int i,int &delta, int &buy, int &sell) {
    if (i == 0){
        return INT_MIN;
    }
    else {
        int current = A[i] - findMin(A, B, i , delta, buy, sell);
        int other = findMax(A, B, i - 1, delta, buy, sell);

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

vector<int> task3A(vector<vector<int> >& A) {
    /*
    * starting from the right, for every company find
    *     OPT(b, s) = max{ (A[i][j]- min(A[i][j] for 0<j<A[0].size()), (OPT(b, s-1) ))
     *
    */

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
        int profit = findMax(A[i], Empty, A[0].size() - 1, delta, buy, sell);
        if (profit > maxProfit){
            best_transaction[0] = i;
            best_transaction[1] = buy;
            best_transaction[2] = sell;
            maxProfit = profit;
        }
    }

    return best_transaction;
}




vector<int> task3B(vector<vector<int> >& A) {
    //PREPROCESSING
    /*Create a secondary 2 dimensional array where each entry stores the max value to the right of it. The array can be
     * built in O(n^2) time by starting on the last index of each subarray and considering for each element either the
     * element to its right or the element that its right element found the largest to the right. */
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

vector<int> task4(vector<vector<int> >& A, int k) {

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
