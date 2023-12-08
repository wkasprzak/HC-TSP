#include "City.h"

City::City(int id, double x, double y) {
    this->id = id;
    this->x = x;
    this->y = y;
}

int City::getId() const {
    return id;
}

double City::getX() const {
    return x;
}

double City::getY() const {
    return y;
}