#include "common.h"


SAW_Analysis FlatPERM_CactusSAW(uint32_t max_size, uint32_t max_tours)
{
    const uint32_t MS  = max_size + 1;
    const uint32_t MS2 = MS * MS;
    const uint32_t max_steps = 3 * max_size + 2;

    std::cout << "MS, MS2: " << MS << ", " << MS2 << std::endl;
    std::cout << "Max Tours: " << max_tours << std::endl;

    // We can use a 1D vector of size max_size + 1 and access [row][col] via std::vector[width * row_idx + col_idx]
    std::vector<double>                     s(MS2, 0.f);
    std::vector<double>                     w(MS2, 0.f);
    std::vector<uint8_t>                    copy(max_steps, 0);
    std::vector<double>                     weight(max_steps, 0);
    std::vector<uint8_t>                    spike_locations(MS, 0);
    std::vector<Site>                       spikes;
    std::vector<Site>                       emptyNeighbours;
    std::vector<Site>                       path;

    uint32_t tours = 0;
    uint32_t spine = 0; // Size of the cactus walk excluding spike steps
    uint32_t n = 0; // Size of the cactus walk including spike steps
    uint32_t m = 0;

    // Initialize walk at origin
    uint8_t atmosphere = 1;
    copy[0] = 1;
    weight[0] = 1;
    s[0] += 1.f;
    w[0] += weight[0];
    Site site { 0,0 };
    path.push_back(site);
    uint32_t max_growth = 0;

    while (tours < max_tours)
    {
        if (spine == std::min(max_growth + 1, max_size) || atmosphere == 0)
        {
            copy[n] = 0;
        }
        else
        {
            // Pruning/enrichment by comparing the target weight with the averaged weight estimate.
            float w_spinem = w[Index(MS,spine,m)];
            float ratio = weight[n] / (w_spinem / s[0]);
            float p = std::fmod(ratio, 1.f);
            float r = Randf01();
            // Generate a random number between 0 and 1
            if (r < p)
                copy[n] = std::floor(ratio) + 1;
            else
                copy[n] = std::floor(ratio);
            weight[n] = w_spinem / s[0];
        }

        if (copy[n] == 0)
        {
            while (n > 0 && copy[n] == 0)
            {
                if (spike_locations[spine] > 0)
                {
                    // Site has a spike, remove spike and adjust spine size accordingly
                    m -= 1;
                    spike_locations[spine] -= 1;
                    spine += 2;
                    spikes.pop_back();
                    spikes.pop_back();
                }
                // Delete last site of walk.
                path.pop_back(); 
                site = path.back();
                spine -= 1;
                n -= 1;
            }
        }

        if (n == 0 && copy[0] == 0)
        {
            // Start new tour.
            tours += 1;
            if (tours % 500 == 0)
            {
                std::cout << "Tour " << tours << " / " << max_tours << std::endl;
                std::cout << tours << " tours ran with spine length " << std::min(max_growth + 1, max_size) << "\n";
                max_growth++;
            }
                
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
            std::vector<Site> neighbours { 
                GetNeighbour(site, STEPS[0]),
                GetNeighbour(site, STEPS[1]),
                GetNeighbour(site, STEPS[2]),
                GetNeighbour(site, STEPS[3]) 
            };

            for (auto& site : neighbours)
            {
                if (std::find(path.begin(), path.end() - 2, site) == path.end() - 2)
                {
                    if (std::find(spikes.begin(), spikes.end(), site) == spikes.end())
                        emptyNeighbours.push_back(site);
                }
            }

            atmosphere = emptyNeighbours.size();
            // If the walk cannot continue, reject entire walk and exit loop.
            if (atmosphere > 0)
            {
                copy[n] -= 1;
                // Randomly select element from emptyNeighbours
                site = emptyNeighbours[std::floor(Randf01() * atmosphere)];
                path.push_back(site);
                spine++;
                n++;

                // Prevent index being out of bound for occupied_sites
                if (spine > 1)
                {
                    if (site == path.end()[-3])
                    {
                        // Spike has formed, add to list of spikes and adjust the spine accordingly
                        spikes.push_back(path.end()[-2]);
                        spikes.push_back(path.end()[-1]);
                        m += 1;
                        spine -= 2;
                        spike_locations[spine] += 1;
                    }
                }
                weight[n] = weight[n-1] * atmosphere;     
                s[Index(MS,spine,m)] += 1.f;
                w[Index(MS,spine,m)] += weight[n];
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
    std::srand(time(0));

    std::string outputSamplesCurrent    = "../data/cactus/cactus_saw_samples.ser";
    std::string outputCountsCurrent     = "../data/cactus/cactus_saw_est.ser";
    std::string outputSamplesBest       = "../data/cactus/cactus_saw_best_samples.ser";
    std::string outputCountsBest        = "../data/cactus/cactus_saw_best_est.ser";

    SAW_Analysis cactusSAW;
    
    uint32_t MS = 50;
    uint32_t MT = 100000;
    std::cout << "Running flatPERM Algorithm with Max Size { " << MS << " } and Max Tours { " << MT << " }.\n";
    auto start = std::chrono::high_resolution_clock::now();
    cactusSAW = FlatPERM_CactusSAW(MS, MT);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Successfully executed in " << duration.count() << "s" << std::endl;
    std::cout << "\n";

    // Only update 'best' data files when more samples AND longer walks have been evaluated
    if (MT >= GetMaxSampleSizeFromFile(outputSamplesBest) && 
        MS + 1 >= GetMaxWalkLengthFromFile(outputCountsBest))
    {
        std::cout << "New best data generated, updating records...\n";
        WriteDataToFile(outputSamplesBest, cactusSAW.samples, MS + 1);
        WriteDataToFile(outputCountsBest, cactusSAW.weights, MS + 1);
    }
    
    std::cout << "Recording samples to 'cactus_saw_samples.ser'\n";
    std::cout << "Recording counts to 'cactus_saw_est.ser'\n";
    // Always write data generated per run
    WriteDataToFile(outputSamplesCurrent, cactusSAW.samples, MS + 1);
    WriteDataToFile(outputCountsCurrent, cactusSAW.weights, MS + 1);

    return 0;
}
