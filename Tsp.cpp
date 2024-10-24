#include <iostream>
#include <vector>
#include <omp.h>
#include <limits>
#include <algorithm>

// Constant for maximum possible cost
const long long INF = std::numeric_limits<long long>::max(); 

/**
 * @class TSPSolver
 * @brief Solves the Traveling Salesman Problem with precedence constraints using Dynamic Programming and Bounding.
 */
class TSPSolver {
private:
    int N; // Number of cities
    std::vector<std::vector<long long>> travelCost; // Cost matrix for travel between cities
    std::vector<std::vector<int>> precedenceConstraints; // Precedence constraints between cities
    std::vector<std::vector<long long>> dp; // DP table for memoization
    std::vector<std::vector<int>> path; // Path table to reconstruct the tour
    
    /**
     * @brief Check if the precedence constraints are satisfied for the next city.
     * @param currentSet A bitmask representing the set of visited cities.
     * @param nextCity The city we want to visit next.
     * @return True if constraints are satisfied, otherwise false.
     */
    bool checkPrecedence(int currentSet, int nextCity) {
        for (int i = 0; i < N; ++i) {
            if ((currentSet & (1 << i)) && precedenceConstraints[nextCity][i]) {
                return false; // Violation of precedence constraints
            }
        }
        return true;
    }

    /**
     * @brief Lower Bound Heuristic to prune suboptimal paths.
     * @param currentSet A bitmask representing the set of visited cities.
     * @param currentCity The current city being evaluated.
     * @return The calculated lower bound for the remaining cities.
     */
    long long lowerBound(int currentSet, int currentCity) {
        long long lb = 0;
        for (int i = 0; i < N; ++i) {
            if (!(currentSet & (1 << i))) { // City not visited
                long long minCost = INF;
                for (int j = 0; j < N; ++j) {
                    if (i != j) {
                        minCost = std::min(minCost, travelCost[i][j]);
                    }
                }
                lb += minCost;
            }
        }
        return lb;
    }

    /**
     * @brief Dynamic Programming with bounding to solve the TSP.
     * @param currentSet A bitmask representing the set of visited cities.
     * @param currentCity The current city.
     * @param upperBound The current upper bound for the solution cost.
     * @return The minimum cost to complete the tour starting from the current city.
     */
    long long tsp(int currentSet, int currentCity, long long upperBound) {
        if (currentSet == (1 << N) - 1) {
            return travelCost[currentCity][0]; // Return to start city
        }

        if (dp[currentSet][currentCity] != -1) {
            return dp[currentSet][currentCity]; // Return memoized result
        }

        long long result = INF;

        // Only prune if the lower bound guarantees the path is worse
        long long lb = lowerBound(currentSet, currentCity);
        if (lb + travelCost[currentCity][0] >= upperBound) {
            return result; // Prune this branch
        }

        // Explore all next cities in parallel
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

    /**
     * @brief Generate an upper bound using a greedy heuristic.
     * @return The greedy upper bound cost.
     */
    long long greedyUpperBound() {
        long long cost = 0;
        int currentCity = 0;
        std::vector<bool> visited(N, false);
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
    /**
     * @brief Constructor to initialize the TSP solver with the number of cities and cost/precedence matrices.
     * @param N Number of cities.
     * @param travelCost The cost matrix for traveling between cities.
     * @param precedenceConstraints Matrix that defines precedence constraints between cities.
     */
    TSPSolver(int N, const std::vector<std::vector<long long>>& travelCost,
              const std::vector<std::vector<int>>& precedenceConstraints)
        : N(N), travelCost(travelCost), precedenceConstraints(precedenceConstraints) {
        dp.assign(1 << N, std::vector<long long>(N, -1)); // Initialize DP table
        path.assign(1 << N, std::vector<int>(N, -1)); // Initialize path tracking table
    }

    /**
     * @brief Solve the TSP problem with precedence constraints and return the minimum cost.
     * @return The minimum cost of the tour.
     */
    long long solve() {
        long long upperBound = greedyUpperBound(); // Initial upper bound
        long long minCost = tsp(1, 0, upperBound); // Start from city 0

        if (minCost == INF) {
            std::cout << "No valid tour found!" << std::endl;
        } else {
            std::cout << "Minimum Cost: " << minCost << std::endl;
            printPath(); // Print the optimal path
        }

        return minCost;
    }

    /**
     * @brief Print the tour path.
     */
    void printPath() {
        int currentSet = 1;
        int currentCity = 0;
        std::cout << "Tour: " << currentCity;

        while (currentSet != (1 << N) - 1) {
            int nextCity = path[currentSet][currentCity];
            if (nextCity == -1) break; // Path tracking failed
            std::cout << " -> " << nextCity;
            currentSet |= (1 << nextCity);
            currentCity = nextCity;
        }

        std::cout << " -> 0" << std::endl; // Return to start city
    }
};

// Example usage
int main() {
    // Cost matrix for 4 cities
    std::vector<std::vector<long long>> travelCost = {
        {INF, 10, 15, 20},
        {10, INF, 35, 25},
        {15, 35, INF, 30},
        {20, 25, 30, INF}
    };

    // Precedence constraints matrix (1 if city i must be visited before city j)
    std::vector<std::vector<int>> precedenceConstraints = {
        {0, 1, 0, 0},  // City 0 must be visited before City 1
        {0, 0, 0, 0},
        {0, 0, 0, 1},  // City 2 must be visited before City 3
        {0, 0, 0, 0}
    };

    // Create the TSP solver
    TSPSolver solver(4, travelCost, precedenceConstraints);

    // Solve the TSP problem and print the result
    solver.solve();

    return 0;
}
