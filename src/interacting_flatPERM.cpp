#include "common.h"


SAW_Analysis FlatPERMInteractingSAW(uint32_t max_size, uint32_t max_tours)
{
    const uint32_t MS  = max_size + 1;
    const uint32_t MS2 = MS * MS;

    std::cout << "MS, MS2: " << MS << ", " << MS2 << std::endl;
    std::cout << "Max Tours: " << max_tours << std::endl;

    // We can use a 1D vector of size MS2 and access [row][col] via std::vector[width * row_idx + col_idx]
    std::vector<double>                 s(MS2, 0.f);
    std::vector<double>                 w(MS2, 0.f);
    std::vector<uint8_t>                copy(MS, 0);
    std::vector<double>                 weight(MS, 0);
    std::vector<uint8_t>                interactions(MS, 0); // There can only be at most 3 interactions at any given path vertex
    std::vector<Site>                   emptyNeighbours;
    std::vector<Site>                   path;

    uint32_t tours = 0;
    uint32_t n = 0;
    uint32_t m = 0;

    // Initialize walk at origin
    uint8_t atmosphere = 1;
    copy[0] = 1;
    weight[0] = 1;
    s[0] += 1.f;
    w[0] += weight[0];
    Site site { 0,0 };
    path.push_back(site);

    while (tours < max_tours)
    {
        // # If maximal length has been reached or the atmosphere is zero: Dont grow.
        if (n == max_size || atmosphere == 0)
        {
            copy[n] = 0;
        }
        else
        {
            // Pruning/enrichment by comparing the target weight with the averaged weight estimate.
            float w_nm = w[Index(MS,n,m)];
            float ratio = weight[n] / (w_nm / s[0]);
            float p = std::fmod(ratio, 1.f);
            float r = Randf01();
            // Generate a random number between 0 and 1
            if (r < p)
                copy[n] = std::floor(ratio) + 1;
            else
                copy[n] = std::floor(ratio);
            weight[n] = w_nm / s[0];
        }

        if (copy[n] == 0)
        {
            while (n > 0 && copy[n] == 0)
            {
                // Delete interactions before retracing a step since interactions are
                // recorded after each next step is taken.
                if (interactions[n] > 0)
                {
                    m -= interactions[n];
                    interactions[n] = 0;
                }
                // Delete last site of walk.
                path.pop_back(); 
                site = path.back();
                n -= 1;
            }
        }

        if (n == 0 && copy[0] == 0)
        {
            if (tours % 1000 == 0)
                std::cout << "Tour " << tours << " / " << max_tours << std::endl;
            // Start new tour.
            tours += 1;
            // Start new walk with step size zero.
            atmosphere = 1;
            copy[0] = 1;
            s[0] += 1;
            w[0] += weight[0];
            site = Site { 0,0 };
            path.clear();
            path.push_back(site);
        }
        else
        {
            emptyNeighbours.clear();
            // Create list of neighboring unoccupied sites, determine the atmosphere a.
            for (Step& step : STEPS)
            {
                Site new_step = GetNeighbour(site, step);
                // std::find returns path.end() if new_step is not found
                if (std::find(path.begin(), path.end()-1, new_step) == path.end()-1)
                {
                    // path does not contain new_step
                    emptyNeighbours.push_back(new_step);
                }
            }

            atmosphere = emptyNeighbours.size();
            // If the walk cannot continue, reject entire walk and exit loop.
            if (atmosphere > 0)
            {
                copy[n] -= 1;
                // Randomly select element from emptyNeighbours
                uint8_t random_index = std::floor(Randf01() * atmosphere);
                site = emptyNeighbours[random_index];
                path.push_back(site);
                n++;
                // Check for interactions between new site and all other visited sites, ignore
                // the last 2 steps before new site.

                for (uint32_t i = 0; i < path.size()-2; i++)
                {
                    Site visited = path[i];
                    // Check for adjacency
                    Step difference = { visited.x - site.x, visited.y - site.y };
                    if (std::find(STEPS.begin(), STEPS.end(), difference) != STEPS.end())
                    {
                        interactions[n] += 1;
                        m += 1;
                    }
                }
                weight[n] = weight[n-1] * atmosphere;     
                s[Index(MS,n,m)] += 1.f;
                w[Index(MS,n,m)] += weight[n];
            }
        }
    }
    s[0] -= 1.f;
    w[0] -= 1.f;

    for (uint32_t i = 0; i < MS2; i++)
    {
        w[i] /= s[0];
    }

    return SAW_Analysis { s, w };
}

int main()
{

    std::string outputSamplesCurrent    = "../data/interacting/interacting_saw_samples.ser";
    std::string outputCountsCurrent     = "../data/interacting/interacting_saw_est.ser";
    std::string outputSamplesBest       = "../data/interacting/interacting_saw_best_samples.ser";
    std::string outputCountsBest        = "../data/interacting/interacting_saw_best_est.ser";

    std::srand(time(0));

    SAW_Analysis interactingSAW;
    
    uint32_t MS = 100;
    uint32_t MT = 100;
    std::cout << "Running flatPERM Algorithm with Max Size { " << MS << " } and Max Tours { " << MT << " }.\n";
    auto start = std::chrono::high_resolution_clock::now();
    interactingSAW = FlatPERMInteractingSAW(MS, MT);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Successfully executed in " << duration.count() << "s" << std::endl;
    std::cout << "\n";

    // Only update 'best' data files when more samples AND longer walks have been evaluated
    if (MT > GetMaxSampleSizeFromFile(outputSamplesBest) && 
        MS + 1 > GetMaxWalkLengthFromFile(outputCountsBest))
    {
        std::cout << "New best data generated, updating records...\n";
        WriteDataToFile(outputSamplesBest, interactingSAW.samples, MS + 1);
        WriteDataToFile(outputCountsBest, interactingSAW.weights, MS + 1);
    }
    
    std::cout << "Recording samples to 'interacting_saw_samples.ser'\n";
    std::cout << "Recording counts to 'interacting_saw_est.ser'\n";
    // Always write data generated per run
    WriteDataToFile(outputSamplesCurrent, interactingSAW.samples, MS + 1);
    WriteDataToFile(outputCountsCurrent, interactingSAW.weights, MS + 1);

    return 0;
}