#include "City.h"
#include "DataGenerator.h"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

double calculateDistance(City &city1, City &city2) {
    return sqrt(pow(city1.getX() - city2.getX(), 2) + pow(city1.getY() - city2.getY(), 2));
}

double checkEcoWaste(std::vector<std::vector<double>> citiesGraph, const std::vector<int>& road, int numberOfCities, double batteryCapacity, int ecoWasteFactor) {
    double totalDistance = 0.0;
    double batteryLevel = batteryCapacity;

    for (int i = 0; i < numberOfCities - 1; ++i) {
        totalDistance += citiesGraph[road[i]][road[i + 1]];
        batteryLevel -= citiesGraph[road[i]][road[i + 1]];
        if (batteryLevel < 0) {
            totalDistance += ecoWasteFactor;
            batteryLevel = batteryCapacity;
        }
    }

    totalDistance += citiesGraph[road.back()][road.front()];
    batteryLevel -= citiesGraph[road.back()][road.front()];

    if (batteryLevel < 0)
        totalDistance += ecoWasteFactor;

    return totalDistance;
}

std::pair<double, std::vector<int>> traditionalGreedyTSP(std::vector<City>& cities, double batteryCapacity, int ecoWasteFactor) {

    int numberOfCities = cities.size();
    std::vector<bool> visited(numberOfCities, false);
    std::vector<int> road;
    int currentCity = 0;
    double currentBatteryCapacity = batteryCapacity;
    double totalDistance = 0.0;

    road.push_back(currentCity);
    visited[currentCity] = true;

    for (int i = 0; i < numberOfCities - 1; i++) {
        int nextCity;
        double minCost = DBL_MAX;

        for (int j = 0; j < numberOfCities; j++) {
            if (!visited[j] && j != currentCity) {
                double distance = calculateDistance(cities[currentCity], cities[j]);
                double newBatteryLevel = currentBatteryCapacity - distance;
                double travelCost = distance;

                if(newBatteryLevel < 0)
                    travelCost += ecoWasteFactor;

                if (travelCost < minCost) {
                    minCost = travelCost;
                    nextCity = j;
                }
            }
        }

        currentBatteryCapacity -= calculateDistance(cities[currentCity], cities[nextCity]);
        if (currentBatteryCapacity < 0)
            currentBatteryCapacity = batteryCapacity;
        
        totalDistance += minCost;

        road.push_back(nextCity);
        visited[nextCity] = true;
        currentCity = nextCity;
    }

    std::clog << "Road: ";
    for (int city : road)
        std::clog << city << " -> ";
    road.push_back(0);
    std::clog << road.back() << std::endl;

    double distance = calculateDistance(cities[road.at(road.size()-1)], cities[road.at(road.size()-2)]);
    totalDistance += distance;
    currentBatteryCapacity -= distance;

    if(currentBatteryCapacity < 0)
        totalDistance += ecoWasteFactor;

    std::clog << "Distance: " << totalDistance << std::endl;
    return std::make_pair(totalDistance, road);
}

std::pair<double, std::vector<int>> randomisedGreedyTSP(std::vector<City>& cities, double batteryCapacity, int ecoWasteFactor) {
    int numberOfCities = cities.size();
    std::vector<bool> visited(numberOfCities, false);
    std::vector<int> road;
    double currentBatteryCapacity = batteryCapacity;
    double totalDistance = 0.0;

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, numberOfCities - 1);
    int currentCity = dis(gen);
    int firstCity = currentCity;

    road.push_back(currentCity);
    visited[currentCity] = true;

    for (int i = 0; i < cities.size() - 1; i++) {
        int nextCity;
        double minCost = DBL_MAX;

        for (int j = 0; j < cities.size(); j++) {
            if (!visited[j] && j != currentCity) {
                double distance = calculateDistance(cities[currentCity], cities[j]);
                double newBatteryLevel = currentBatteryCapacity - distance;
                double travelCost = distance;

                if(newBatteryLevel < 0)
                    travelCost += ecoWasteFactor;

                if (travelCost < minCost) {
                    minCost = travelCost;
                    nextCity = j;
                }
            }
        }

        currentBatteryCapacity -= calculateDistance(cities[currentCity], cities[nextCity]);
        if (currentBatteryCapacity < 0)
            currentBatteryCapacity = batteryCapacity;
        
        totalDistance += minCost;
        road.push_back(nextCity);
        visited[nextCity] = true;
        currentCity = nextCity;
    }     

    std::clog << "Road: ";
    for (int city : road)
        std::clog << city << " -> ";
    road.push_back(firstCity);
    std::clog << road.back() << std::endl;

    double distance = calculateDistance(cities[road.at(road.size()-1)], cities[road.at(road.size()-2)]);
    totalDistance += distance;
    currentBatteryCapacity -= distance;

    if(currentBatteryCapacity < 0)
        totalDistance += ecoWasteFactor;

    std::clog << "Distance: " << totalDistance << std::endl;
    return std::make_pair(totalDistance, road);
}

std::pair<double, std::vector<int>> bruteForceTSP(std::vector<City>& cities, std::vector<std::vector<double>> citiesGraph, double batteryCapacity, int ecoWasteFactor) {
    int numberOfCities = cities.size();
    std::vector<int> answerRoad;
    std::vector<int> minRoad;
    double minCost = DBL_MAX;

    std::vector<int> order(numberOfCities);
    for (int i = 0; i < numberOfCities; i++)
        order[i] = i;

    do {
        double currCost = checkEcoWaste(citiesGraph, order, numberOfCities, batteryCapacity, ecoWasteFactor);
        if (currCost < minCost) {
            minCost = currCost;
            minRoad = order;
        }
    } while (std::next_permutation(order.begin(), order.end()));

    answerRoad = minRoad;

    std::clog << "Order of cities: ";
    for (int city : answerRoad)
        std::clog << city << " ";
    std::clog << "\nDistance: " << minCost << '\n';
    std::clog << '\n';

    return std::make_pair(minCost, answerRoad);
}

std::vector<int> generateRandomTour(int numberOfCities) {

    std::vector<int> tour;
    for (int i = 0; i < numberOfCities; i++)
        tour.push_back(i);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(tour.begin(), tour.end(), std::default_random_engine(seed));

    return tour;
}

std::pair<double, std::vector<int>> tabuSearchTSP(const std::vector<City>& cities, std::vector<std::vector<double>> citiesGraph, double batteryCapacity, int ecoWasteFactor) {

    int maxIterations = 1000; // 100000
    int numberOfCities = cities.size();
    int tabuListSize = numberOfCities; // 5

    std::vector<int> bestTour = generateRandomTour(numberOfCities);
    std::vector<int> currentTour = bestTour;
    std::vector<std::vector<int>> tabuList;

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::vector<std::vector<int>> neighborhood;

        for (int i = 1; i < numberOfCities - 1; ++i) {
            for (int j = i + 1; j < numberOfCities; ++j) {
                std::vector<int> neighbor = currentTour;
                std::swap(neighbor[i], neighbor[j]);
                neighborhood.push_back(neighbor);
            }
        }

        std::vector<int> bestNeighbor;
        double minNeighborEcoWaste = DBL_MAX;

        for (const auto& neighbor : neighborhood) {
            if (std::find(tabuList.begin(), tabuList.end(), neighbor) == tabuList.end()) {
                double neighborEcoWaste = checkEcoWaste(citiesGraph, neighbor, numberOfCities, batteryCapacity,
                                                        ecoWasteFactor);
                if (neighborEcoWaste < minNeighborEcoWaste) {
                    minNeighborEcoWaste = neighborEcoWaste;
                    bestNeighbor = neighbor;
                }
            }
        }

        if (!bestNeighbor.empty()) {
            currentTour = bestNeighbor;
            tabuList.push_back(currentTour);
            if (tabuList.size() > tabuListSize) {
                tabuList.erase(tabuList.begin());
            }

            double minEcoWaste = checkEcoWaste(citiesGraph, bestTour, numberOfCities, batteryCapacity, ecoWasteFactor);
            double currentTourEcoWaste = checkEcoWaste(citiesGraph, currentTour, numberOfCities, batteryCapacity,
                                                       ecoWasteFactor);
            if (currentTourEcoWaste < minEcoWaste) {
                bestTour = currentTour;
            }
        }
    }

    double bestTourEcoWaste = checkEcoWaste(citiesGraph, bestTour, numberOfCities, batteryCapacity, ecoWasteFactor);

    std::clog << "Optimized Tour: ";
    for (int city : bestTour) {
        std::clog << city << " ";
    }
    std::clog << std::endl;

    std::clog << "Total Cost: " << bestTourEcoWaste << std::endl;

    return std::make_pair(bestTourEcoWaste, bestTour);
}

int main() {

    std::string outputFileName = "testy4.csv";
    std::ofstream outputFile(outputFileName, std::ios::app);
    if (!outputFile.is_open()) {
        std::cerr << "Nie można otworzyć pliku wynikowego.\n";
        return -1;
    }

    outputFile << "numberOfCities" << ";" << "method" << ";" << "iteration" << ";" << "time" << ";" << "ecoWaste" << std::endl;

    for(int m = 4; m <= 13; m++) {
        for(int p = 0; p < 5; p++) {

            std::ifstream dataFile;
            std::string nameOfFile;
            std::string fileName;

            // Opening file and creating objects
            nlohmann::json dataFromFile;
            std::vector<City> cities;
            int leftConstraint = 0, rightConstraint = 10000;
            DataGenerator generator;

            generator.generateData("data" + std::to_string(m) + "_" + std::to_string(p+1), m, leftConstraint, rightConstraint);
            nameOfFile = "./data/data" + std::to_string(m) + "_" + std::to_string(p+1) + ".json";
            dataFile.open(nameOfFile);
            if (!dataFile.is_open()) {
                std::cerr << "Unable to open JSON file." << std::endl;
                return 1;
            }
            std::cout << "Udalo sie " + std::to_string(m) + "." + std::to_string(p) + "\n";

            // Importing data
            dataFile >> dataFromFile;
            dataFile.close();
            for (const auto& oneCity : dataFromFile) {
                City city(oneCity["id"], oneCity["x"], oneCity["y"]);
                cities.push_back(city);
            }

            // Variables needed for TSP
            std::vector<std::vector<double>> citiesGraph(cities.size(), std::vector<double>(cities.size()));
            int ecoWasteFactor = cities.size();
            double batteryCapacity = 0;

            // Generating battery capacity
            for (int i = 0; i < cities.size(); i++)
                for (int j = 0; j < cities.size(); j++)
                    if (i != j)
                        citiesGraph[i][j] = sqrt(pow((cities.at(i).getX() - cities.at(j).getX()), 2) + pow((cities.at(i).getY() - cities.at(j).getY()), 2));
                    else
                        citiesGraph[i][j] = 0;

            for (int i = 0; i < citiesGraph[0].size(); i++)
                batteryCapacity += citiesGraph[0][i];
            batteryCapacity = ceil(batteryCapacity * 0.3);

            // Algorithms
            auto startBF = std::chrono::high_resolution_clock::now();
            std::pair<double, std::vector<int>> BF = bruteForceTSP(cities, citiesGraph, batteryCapacity, ecoWasteFactor);
            auto endBF = std::chrono::high_resolution_clock::now();
            auto durationBF = std::chrono::duration_cast<std::chrono::nanoseconds>(endBF - startBF).count();
            outputFile << m << ";" << "BruteForce" << ";" << p << ";" << durationBF << ";" << BF.first << std::endl;

            auto startG1 = std::chrono::high_resolution_clock::now();
            std::pair<double, std::vector<int>> G1 = traditionalGreedyTSP(cities, batteryCapacity, ecoWasteFactor);
            auto endG1 = std::chrono::high_resolution_clock::now();
            auto durationG1 = std::chrono::duration_cast<std::chrono::nanoseconds>(endG1 - startG1).count();
            outputFile << m << ";" << "Greedy" << ";" << p << ";" << durationG1 << ";" << G1.first << std::endl;

            auto startG2 = std::chrono::high_resolution_clock::now();
            std::pair<double, std::vector<int>> G2 = randomisedGreedyTSP(cities, batteryCapacity, ecoWasteFactor); auto endG2 = std::chrono::high_resolution_clock::now();
            auto durationG2 = std::chrono::duration_cast<std::chrono::nanoseconds>(endG2 - startG2).count();
            outputFile << m << ";" << "RandomisedGreedy" << ";" << p << ";" << durationG2 << ";" << G2.first << std::endl;

            auto startTS = std::chrono::high_resolution_clock::now();
            std::pair<double, std::vector<int>> TS = tabuSearchTSP(cities, citiesGraph, batteryCapacity, ecoWasteFactor);
            auto endTS = std::chrono::high_resolution_clock::now();
            auto durationTS = std::chrono::duration_cast<std::chrono::nanoseconds>(endTS - startTS).count();
            outputFile << m << ";" << "TabuSearch" << ";" << p << ";" << durationTS << ";" << TS.first << std::endl;
        }
    }
    outputFile.close();
}