#include <iostream>
#include <vector>
#include <omp.h>
#include <limits>
#include <algorithm>

using namespace std;

// Constant for maximum possible cost
const long long INF = numeric_limits<long long>::max();

class TSPSolver {
private:
    int N; // Number of cities
    vector<vector<long long>> travelCost; // Cost matrix for travel between cities
    vector<vector<int>> precedenceConstraints; // Precedence constraints between cities
    vector<vector<long long>> dp; // DP table for memoization
    vector<vector<int>> path; // Path table to reconstruct the tour
    
    bool checkPrecedence(int currentSet, int nextCity) {
        for (int i = 0; i < N; ++i) {
            if ((currentSet & (1 << i)) && precedenceConstraints[nextCity][i]) {
                return false; // Violation of precedence constraints
            }
        }
        return true;
    }

    long long lowerBound(int currentSet, int currentCity) {
        long long lb = 0;
        for (int i = 0; i < N; ++i) {
            if (!(currentSet & (1 << i))) { // City not visited
                long long minCost = INF;
                for (int j = 0; j < N; ++j) {
                    if (i != j) {
                        minCost = min(minCost, travelCost[i][j]);
                    }
                }
                lb += minCost;
            }
        }
        return lb;
    }

    long long tsp(int currentSet, int currentCity, long long upperBound) {
        if (currentSet == (1 << N) - 1) {
            return travelCost[currentCity][0]; // Return to start city
        }

        if (dp[currentSet][currentCity] != -1) {
            return dp[currentSet][currentCity]; // Return memoized result
        }

        long long result = INF;

        long long lb = lowerBound(currentSet, currentCity);
        if (lb + travelCost[currentCity][0] >= upperBound) {
            return result; // Prune this branch
        }

        #pragma omp parallel for schedule(dynamic) reduction(min: result)
        for (int nextCity = 0; nextCity < N; ++nextCity) {
            if (!(currentSet & (1 << nextCity)) && checkPrecedence(currentSet, nextCity)) {
                int nextSet = currentSet | (1 << nextCity); // Mark city as visited
                long long nextCost = travelCost[currentCity][nextCity];

                if (nextCost == INF) continue;

                long long cost = tsp(nextSet, nextCity, upperBound - nextCost);
                if (cost != INF) {
                    cost += nextCost;
                    #pragma omp critical
                    {
                        if (cost < result) {
                            result = cost;
                            path[currentSet][currentCity] = nextCity; // Track path
                        }
                    }
                }
            }
        }

        dp[currentSet][currentCity] = result;
        return result;
    }

    long long greedyUpperBound() {
        long long cost = 0;
        int currentCity = 0;
        vector<bool> visited(N, false);
        visited[0] = true;

        for (int count = 1; count < N; ++count) {
            long long minCost = INF;
            int nextCity = -1;
            for (int i = 0; i < N; ++i) {
                if (!visited[i] && travelCost[currentCity][i] < minCost) {
                    minCost = travelCost[currentCity][i];
                    nextCity = i;
                }
            }
            if (nextCity != -1) {
                cost += minCost;
                visited[nextCity] = true;
                currentCity = nextCity;
            }
        }
        cost += travelCost[currentCity][0]; // Return to the start city
        return cost;
    }

public:
    TSPSolver(int N, const vector<vector<long long>>& travelCost,
              const vector<vector<int>>& precedenceConstraints)
        : N(N), travelCost(travelCost), precedenceConstraints(precedenceConstraints) {
        dp.assign(1 << N, vector<long long>(N, -1)); // Initialize DP table
        path.assign(1 << N, vector<int>(N, -1)); // Initialize path tracking table
    }

    long long solve() {
        long long upperBound = greedyUpperBound(); // Initial upper bound
        long long minCost = tsp(1, 0, upperBound); // Start from city 0

        if (minCost == INF) {
            cout << "No valid tour found!" << endl;
        } else {
            cout << "Minimum Cost: " << minCost << endl;
            printPath(); // Print the optimal path
        }

        return minCost;
    }

    void printPath() {
        int currentSet = 1;
        int currentCity = 0;
        cout << "Tour: " << currentCity;

        while (currentSet != (1 << N) - 1) {
            int nextCity = path[currentSet][currentCity];
            if (nextCity == -1) break; // Path tracking failed
            cout << " -> " << nextCity;
            currentSet |= (1 << nextCity);
            currentCity = nextCity;
        }

        cout << " -> 0" << endl; // Return to start city
    }
};

int main() {
    int N;
    cout << "Enter the number of cities: ";
    cin >> N;

    vector<vector<long long>> travelCost(N, vector<long long>(N));
    cout << "Enter the travel cost matrix (use " << INF << " for no direct path):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cin >> travelCost[i][j];
        }
    }

    vector<vector<int>> precedenceConstraints(N, vector<int>(N));
    cout << "Enter the precedence constraints matrix (1 if city i must be visited before city j, else 0):" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cin >> precedenceConstraints[i][j];
        }
    }

    TSPSolver solver(N, travelCost, precedenceConstraints);
    solver.solve();

    return 0;
}
