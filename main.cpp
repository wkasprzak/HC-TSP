#include "City.h"
#include "DataGenerator.h"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <direct.h>
#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <vector>

double calculateDistance(City &city1, City &city2) {
    return sqrt(pow(city1.getX() - city2.getX(), 2) + pow(city1.getY() - city2.getY(), 2));
}

void traditionalGreedyTSP(std::vector<City>& cities, double batteryCapacity, int ecoWasteFactor) {
    
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
}

void randomisedGreedyTSP(std::vector<City>& cities, double batteryCapacity, int ecoWasteFactor) {
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
}

void bruteForceTSP(std::vector<City>& cities, std::vector<std::vector<double>> citiesGraph, double batteryCapacity, int ecoWasteFactor) {
    
    int numberOfCities = cities.size();

    std::vector<int> answerRoad;
    std::vector<int> minRoad;

    double minCost = DBL_MAX;
    double currentBatteryCapacity = batteryCapacity;
    
    std::vector<int> order(numberOfCities);
    for (int i = 0; i < numberOfCities; i++)
        order[i] = i;

    do {
        double currCost = 0.0;
        for (int i = 0; i < numberOfCities - 1; ++i) {
            currCost += citiesGraph[order[i]][order[i + 1]];
            currentBatteryCapacity -= citiesGraph[order[i]][order[i + 1]];
            if(currentBatteryCapacity < 0) {
                currCost += ecoWasteFactor;
                currentBatteryCapacity = batteryCapacity;
            }
        }
        currCost += citiesGraph[order[numberOfCities - 1]][order[0]];
        currentBatteryCapacity -= citiesGraph[order[numberOfCities - 1]][order[0]];
        if(currentBatteryCapacity < 0) {
            currCost += ecoWasteFactor;
            currentBatteryCapacity = batteryCapacity;
        }

        if (currCost < minCost) {
            minCost = currCost;
            minRoad = order;
        }
    } while (next_permutation(order.begin() + 1, order.end()));
    answerRoad = minRoad;

    std::clog << "Order of cities: ";
    for (int city : answerRoad)
        std::clog << city << " ";
    std::clog << "\nDistance: " << minCost << std::endl;
    std::clog << std::endl;
}

std::vector<int> generateRandomTour(int numberOfCities) {
    std::vector<int> tour;
    tour.reserve(numberOfCities);

    for (int i = 0; i < numberOfCities; i++)
        tour.push_back(i);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(tour.begin(), tour.end(), std::default_random_engine(seed));

    return tour;
}

double checkEcoWaste(std::vector<std::vector<double>> citiesGraph, const std::vector<int>& tour, int numberOfCities, double batteryCapacity, int ecoWasteFactor) {
    double totalDistance = 0.0;
    double ecologicalFootprint = 0.0;
    double batteryLevel = batteryCapacity;

    for (int i = 0; i < numberOfCities - 1; ++i) {

        totalDistance += citiesGraph[tour[i]][tour[i + 1]];
        batteryLevel -= citiesGraph[tour[i]][tour[i + 1]];
        if (batteryLevel < 0) {
            ecologicalFootprint += ecoWasteFactor; 
            batteryLevel = batteryCapacity;
        }
    }

    totalDistance += citiesGraph[tour.back()][tour.front()];
    batteryLevel -= citiesGraph[tour.back()][tour.front()];

    if (batteryLevel < 0) {
        ecologicalFootprint += ecoWasteFactor;
        batteryLevel = batteryCapacity;
    }
    return totalDistance + ecologicalFootprint;
}

void tabuSearchTSP(const std::vector<City>& cities, std::vector<std::vector<double>> citiesGraph, double batteryCapacity, int ecoWasteFactor) {

    int tabuListSize = 5;
    int maxIterations = 1000;
    int numberOfCities = cities.size();

    std::vector<int> currentTour = generateRandomTour(numberOfCities);
    std::vector<int> bestTour = currentTour;
    std::vector<std::vector<int>> tabuList;

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::vector <std::vector<int>> neighborhood;

        for (int i = 0; i < numberOfCities - 1; ++i) {
            for (int j = i + 1; j < numberOfCities; ++j) {
                std::vector<int> neighbor = currentTour;
                std::swap(neighbor[i], neighbor[j]);
                neighborhood.push_back(neighbor);
            }
        }

        std::vector<int> bestNeighbor;
        double minNeighborEcoWaste = DBL_MAX;

        for (const auto& neighbor : neighborhood) {
            if (find_if(tabuList.begin(), tabuList.end(), [&neighbor](const auto& element) {
                return element == neighbor;
            }) == tabuList.end()) {
                double neighborEcoWaste = checkEcoWaste(citiesGraph, neighbor, numberOfCities, batteryCapacity,
                                                        ecoWasteFactor);
                if (neighborEcoWaste <= minNeighborEcoWaste) {
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
}

int main() {
    
    // Generating or getting data - query
    std::ifstream dataFile;
    std::string nameOfFile;
    std::string fileName;
    char option;

    // Opening file and creating objects
    nlohmann::json dataFromFile;
    std::vector<City> cities;
    for(int i = 4; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            unsigned int n;
            int leftConstraint = 0, rightConstraint = 1000;
            DataGenerator generator;
            //std::clog << "\nPlease enter the data in order: \nFile name, number of cities, left constraint, right constraint" << std::endl;
            //std::cin >> fileName >> n >> leftConstraint >> rightConstraint;
            generator.generateData("data" + std::to_string(i) + "_" + std::to_string(j+1), i, leftConstraint, rightConstraint);
            nameOfFile = "./data/data" + std::to_string(i) + "_" + std::to_string(j+1) + ".json";
            dataFile.open(nameOfFile);
            std::cout << "Udalo sie " + std::to_string(i) + "." + std::to_string(j) + "\n";

            // Importing data
            dataFile >> dataFromFile;
            dataFile.close();
            for (const auto& oneCity : dataFromFile) {
                City city(oneCity["id"], oneCity["x"], oneCity["y"]);
                cities.push_back(city);
            }

            for (const City& city : cities)
                std::clog << "City: " << city.getId() << " x = " << city.getX() << " y = " << city.getY() << std::endl;


            // Printing cities
//            std::clog << "\nWould you like to see the cities? (y/n) " << std::endl;
//            std::cin >> option;
//            if (option == 'y')

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

            for (int i = 0; i < citiesGraph.size(); i++) {
                for (int j = 0; j < citiesGraph[i].size(); j++)
                    std::clog << citiesGraph[i][j] << "\t";
                std::clog << std::endl;
            }

            for (int i = 0; i < citiesGraph[0].size(); i++)
                batteryCapacity += citiesGraph[0][i];
            batteryCapacity = ceil(batteryCapacity * 0.3);

            // USUNĄĆ !!!!
            //batteryCapacity = 25;

            // Algorithms
            std::chrono::steady_clock::time_point beginG1 = std::chrono::steady_clock::now();
            traditionalGreedyTSP(cities, batteryCapacity, ecoWasteFactor);
            std::chrono::steady_clock::time_point endG1 = std::chrono::steady_clock::now();
            std::clog << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (endG1 - beginG1).count() << "[ns]" << std::endl;

            std::chrono::steady_clock::time_point beginG2 = std::chrono::steady_clock::now();
            randomisedGreedyTSP(cities, batteryCapacity, ecoWasteFactor);
            std::chrono::steady_clock::time_point endG2 = std::chrono::steady_clock::now();
            std::clog << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (endG2 - beginG2).count() << "[ns]" << std::endl;

            std::chrono::steady_clock::time_point beginBF = std::chrono::steady_clock::now();
            bruteForceTSP(cities, citiesGraph, batteryCapacity, ecoWasteFactor);
            std::chrono::steady_clock::time_point endBF = std::chrono::steady_clock::now();
            std::clog << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (endBF - beginBF).count() << "[ns]" << std::endl;

            std::chrono::steady_clock::time_point beginTS = std::chrono::steady_clock::now();
            tabuSearchTSP(cities, citiesGraph, batteryCapacity, ecoWasteFactor);
            std::chrono::steady_clock::time_point endTS = std::chrono::steady_clock::now();
            std::clog << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (endTS - beginTS).count() << "[ns]" << std::endl;
        }
    }
//
//    std::clog << "Would you like to generate data and use it? (y/n) " << std::endl;
//    std::cin >> option;
//    switch(option) {
//        case 'y':
//        case 'Y':
//
//        break;
//        case 'n':
//        case 'N':
//            std::clog << "\nWhat file would you like to use?" << std::endl;
//            std::cin >> fileName;
//            nameOfFile = "./data/" + fileName + ".json";
//            dataFile.open(nameOfFile);
//        break;
//        default:
//            std::clog << "Choice error" << std::endl;
//            return 1;
//    }
//
//    if (!dataFile.is_open()) {
//        std::cerr << "Unable to open JSON file." << std::endl;
//        return 1;
//    }




}