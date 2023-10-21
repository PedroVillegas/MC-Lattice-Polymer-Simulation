#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <fstream>
#include <string>


void WriteDataToFile(std::string filename, const std::vector<double>& data, uint32_t dataWidth);
uint32_t GetMaxSampleSizeFromFile(std::string filename);
uint32_t GetMaxWalkLengthFromFile(std::string filename);

struct Step
{
    int dx;
    int dy;

    bool operator == (const Step &other) const 
    {
        return dx == other.dx && dy == other.dy;
    }
};

struct Site
{
    int x;
    int y;

    bool operator == (const Site &other) const 
    {
        return x == other.x && y == other.y;
    }
};

extern std::array<Step, 4> STEPS;

uint32_t Index(uint32_t width, uint32_t row_idx, uint32_t col_idx);
float Randf01();
Site GetNeighbour(Site site, Step step);

struct SAW_Analysis
{
    std::vector<double> samples;
    std::vector<double> weights;
};
