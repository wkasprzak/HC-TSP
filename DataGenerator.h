#ifndef _DATAGENERATOR_H_
#define _DATAGENERATOR_H_

#include <string>

class DataGenerator {

public:
    void generateData(const std::string& dataFile, unsigned int numberOfCities, int leftConstraint, int rightConstraint);
};

#endif