#include "DataGenerator.h"
#include "nlohmann/json.hpp"

#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <random>
#include <chrono>

void DataGenerator::generateData(const std::string& dataFile, unsigned int numberOfCities, int leftConstraint, int rightConstraint) {

    nlohmann::json dataArray;

    if (leftConstraint > rightConstraint)
        std::swap(leftConstraint, rightConstraint);

    // Use current time as a seed
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(static_cast<unsigned>(seed));
    std::uniform_int_distribution<int> dist(leftConstraint, rightConstraint);

    // Generating data
    for(int i = 0; i < numberOfCities; i++) {
        nlohmann::json generatedData;
        generatedData["id"] = i;
        generatedData["x"] = dist(gen);
        generatedData["y"] = dist(gen);
        dataArray.push_back(generatedData);
    }

    // Saving data into file
    std::string filePath = "./data/" + dataFile + ".json";
    std::ofstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open JSON file." << std::endl;
        return;
    }
    file << dataArray.dump(4);
    file.close();
}