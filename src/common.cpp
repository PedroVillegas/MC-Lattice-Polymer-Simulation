#include "common.h"


std::array<Step, 4> STEPS = { {{1,0}, {0,1}, {-1,0}, {0,-1}} };

uint32_t Index(uint32_t width, uint32_t row_idx, uint32_t col_idx)
{
    return width * row_idx + col_idx;
}

float Randf01()
{
    return ((double)std::rand()) / RAND_MAX;
}

Site GetNeighbour(Site site, Step step)
{
    Site neighbour { site.x + step.dx, site.y + step.dy };
    return neighbour;
}

uint32_t GetMaxSampleSizeFromFile(std::string filename)
{
    std::ifstream file;
    file.open(filename);

    std::string firstLine;
    std::getline(file, firstLine);

    // std::string::substr(offset, upto)
    // std::atof(const char* string) converts scientific notation to double
    uint32_t maxSampleSize = std::atof(firstLine.substr(0, firstLine.find(" ")).c_str());

    return maxSampleSize;
}

uint32_t GetMaxWalkLengthFromFile(std::string filename)
{
    std::ifstream file;
    file.open(filename);

    int numLines = 0;
    std::string unused;
    while (std::getline(file, unused))
        ++numLines;

    return numLines;
}

void WriteDataToFile(std::string filename, const std::vector<double>& data, uint32_t dataWidth)
{
    std::ofstream file;
    file.open(filename);
    
    for (uint32_t i = 0; i < dataWidth * dataWidth; i++)
    {
        if (i % dataWidth == 0 && i != 0)
            file << "\n";
        file << data[i] << " ";
    }

    file.close();
}