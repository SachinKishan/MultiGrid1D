#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "color.h"

void parseCSV(const std::string& filename, const double value) {
    std::string searchValue = std::to_string(value);
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }

        if (!row.empty() && row[0] == searchValue) {
            std::cout << line << std::endl;
            double val = std::stod(row[0]);
            std::cout << val + 4;
        }
    }

    file.close();
}

color matchColorToWavelength(const std::string& filename, const int value) {
    std::string searchValue = std::to_string(value);
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return Black;
    }

    std::string line;
    std::vector<std::string> row;

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        row.clear();

        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }

        if (!row.empty() && row[0] == searchValue) {
            //std::cout << line << std::endl;
            //double val = std::stod(row[0]);
            //std::cout << val + 4;
            break;
        }
    }

    file.close();
    return color(std::stod(row[1]), std::stod(row[2]), std::stod(row[3]));
}


inline color sRGBtoCIEXYZ(color rgb)
{


    double var_R = (rgb.r / 255);
    double var_G = (rgb.g / 255);
    double var_B = (rgb.b / 255);
    if (var_R > 0.04045) var_R = pow((var_R + 0.055) / 1.055, 2.4);
    else                   var_R = var_R / 12.92;
    if (var_G > 0.04045) var_G = pow((var_G + 0.055) / 1.055, 2.4);
    else                   var_G = var_G / 12.92;
    if (var_B > 0.04045) var_B = pow((var_B + 0.055) / 1.055, 2.4);
    else                   var_B = var_B / 12.92;

    var_R = var_R * 100;
    var_G = var_G * 100;
    var_B = var_B * 100;

    double X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
    double Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
    double Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;
    return color(X, Y, Z);
}
double clamp(double min, double max, double in)
{
    if (in < min)in = min;
    if (in > max)in = max;
    return in;
}
color CIEXYZtosRGB(color XYZ)
{
    //normalize
    double var_X = XYZ.r;
    double var_Y = XYZ.g;
    double var_Z = XYZ.b;

    double var_R = var_X * 3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
    double var_G = var_X * -0.9689 + var_Y * 1.8758 + var_Z * 0.0415;
    double var_B = var_X * 0.0557 + var_Y * -0.2040 + var_Z * 1.0570;

    if (var_R > 0.0031308) var_R = 1.055 * pow(var_R, (1 / 2.4)) - 0.055;
    else                     var_R = 12.92 * var_R;
    if (var_G > 0.0031308) var_G = 1.055 * pow(var_G, (1 / 2.4)) - 0.055;
    else                     var_G = 12.92 * var_G;
    if (var_B > 0.0031308) var_B = 1.055 * pow(var_B, (1 / 2.4)) - 0.055;
    else                     var_B = 12.92 * var_B;


    double sR = clamp(0, 1, var_R);
    double sG = clamp(0, 1, var_G);
    double sB = clamp(0, 1, var_B);
    return color(sR, sG, sB);
}
