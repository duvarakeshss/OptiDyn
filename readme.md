# Traveling Salesman Problem with Precedence Constraints (TSP-PC)

This project solves the Traveling Salesman Problem with precedence constraints using dynamic programming and bounding (DPBB). It leverages OpenMP for parallelism to speed up the exploration of possible solutions and uses a greedy upper bound heuristic to assist in bounding.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Usage](#usage)
- [Explanation of Code](#explanation-of-code)
- [Example Output](#example-output)

## Introduction

The **Traveling Salesman Problem with Precedence Constraints** (TSP-PC) is a variant of the classic TSP problem where certain cities must be visited before others. The goal is to find the shortest possible tour that visits every city exactly once while adhering to these precedence constraints.

### Problem Definition

Given:

- `N` cities and a cost matrix `travelCost[N][N]` where `travelCost[i][j]` is the cost of traveling from city `i` to city `j`.
- A precedence constraints matrix `precedenceConstraints[N][N]`, where `precedenceConstraints[i][j] = 1` means city `i` must be visited before city `j`.

The objective is to find the optimal tour starting and ending at city 0, visiting all cities, and satisfying the precedence constraints, while minimizing the total travel cost.

## Features

- **Dynamic Programming with Bounding (DPBB)**: Utilizes memoization to store the minimum cost for each subset of visited cities, with bounding to prune unpromising paths.
- **Parallelization**: Uses OpenMP to explore multiple cities in parallel, speeding up the solution process.
- **Greedy Upper Bound**: A heuristic upper bound is calculated to improve the efficiency of bounding.
- **Path Reconstruction**: After solving for the minimum cost, the optimal tour path is reconstructed and printed.

## Requirements

To compile and run this project, you'll need:

- A C++ compiler that supports C++11 or later (e.g., g++)
- OpenMP for parallelism (optional, but recommended for performance)

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/duvarakeshss/TSP-BoundQuest
   ```
2. Navigate to the project directory:

   ```bash
   cd TSP-BoundQuest
   ```
3. Compile the code using a C++ compiler:

   ```bash
   g++ -fopenmp -o tsp tsp.cpp
   ```
4. Run the compiled program:

   ```bash
   ./tsp
   ```

## Usage

You can modify the number of cities, the travel cost matrix, and the precedence constraints directly in the code.

### Example Problem

In the current code setup, we have 4 cities and the following travel cost matrix:
