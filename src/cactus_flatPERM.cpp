#include "common.h"


SAW_Analysis FlatPERM_CactusSAW(const uint32_t maxSize, const uint32_t maxTours)
{
    const uint32_t MS  = maxSize + 1;
    const uint32_t MS2 = MS * MS;
    const uint32_t max_steps = 3 * maxSize + 2;

    // We can use a 1D vector of size maxSize + 1 and access [row][col] via std::vector[width * row_idx + col_idx]
    std::vector<double> s(MS2, 0.f);
    std::vector<double> w(MS2, 0.f);
    std::vector<double> weight(max_steps, 0);
    std::vector<uint8_t> copy(max_steps, 0);
    std::vector<uint8_t> spikeLocations(MS, 0);
    std::vector<Site> spikes;
    std::vector<Site> emptyNeighbours;
    std::vector<Site> path;

    // n is the length of the cactus walk including spike steps
    // spine is the length of the cactus walk excluding spike steps
    uint32_t tours = 0;
    uint32_t spine = 0;
    uint32_t n = 0;
    uint32_t m = 0;

    // Initialize walk at origin
    Site site { 0,0 };
    path.push_back(site);
    copy[0] = 1;
    weight[0] = 1;
    s[0] += 1.f;
    w[0] += weight[0];
    uint8_t atmosphere = 1;
    uint32_t lengthToGrowTo = 0;
    const uint32_t growthConstant = 5000;

    while (tours < maxTours)
    {
        if (spine == std::min(lengthToGrowTo + 1, maxSize) || atmosphere == 0)
        {
            copy[n] = 0;
        }
        else
        {
            // Pruning/enrichment by comparing the target weight with the normalized weight estimate
            double w_sm = w[Index(MS,spine,m)];
            // We only want to normalize w_sm by the number of tours that had a max walk length of min(lengthToGrowTo, maxSize)
            double tours_s = (spine == maxSize ? tours+1 - (growthConstant * maxSize) : tours+1 % growthConstant);
            double normalizedWeight = w_sm / tours_s;

            float ratio = weight[n] / normalizedWeight;
            float p = std::fmod(ratio, 1.f);
            float r = Randf01();
            // Generate a random number between 0 and 1
            if (r < p)
                copy[n] = std::floor(ratio) + 1;
            else
                copy[n] = std::floor(ratio);
            weight[n] = normalizedWeight;
        }

        if (copy[n] == 0)
        {
            while (n > 0 && copy[n] == 0)
            {
                if (spikeLocations[spine] > 0)
                {
                    // Site has a spike, remove spike and adjust spine size accordingly
                    m -= 1;
                    spikeLocations[spine] -= 1;
                    spine += 2;
                    spikes.pop_back();
                    spikes.pop_back();
                }
                // Delete last site of walk
                path.pop_back(); 
                site = path.back();
                spine -= 1;
                n -= 1;
            }
        }

        if (n == 0 && copy[0] == 0)
        {
            // Start new tour
            tours += 1;
            if (tours % growthConstant == 0)
            {
                std::cout << "Tour " << tours << " / " << maxTours << std::endl;
                std::cout << tours << " tours ran with spine length " << std::min(lengthToGrowTo + 1, maxSize) << "\n";
                lengthToGrowTo++;
            }
                
            // Start new walk with step size zero
            site = Site { 0,0 };
            path.clear();
            path.push_back(site);
            copy[0] = 1;
            s[0] += 1;
            w[0] += weight[0];
            atmosphere = 1;
        }
        else
        {
            emptyNeighbours.clear();
            // Create list of neighboring unoccupied sites, determine the atmosphere a
            std::vector<Site> neighbours { 
                GetNeighbour(site, STEPS[0]),
                GetNeighbour(site, STEPS[1]),
                GetNeighbour(site, STEPS[2]),
                GetNeighbour(site, STEPS[3]) 
            };

            for (auto& site : neighbours)
            {   
                if (path.size() < 3)
                {
                    emptyNeighbours.push_back(site);
                    continue;
                }

                if (std::find(path.begin(), path.end() - 2, site) == path.end() - 2)
                {
                    if (std::find(spikes.begin(), spikes.end(), site) == spikes.end())
                        emptyNeighbours.push_back(site);
                }
            }

            atmosphere = emptyNeighbours.size();
            // If the walk cannot continue, reject entire walk and exit loop
            if (atmosphere > 0)
            {
                copy[n] -= 1;
                // Each site can only be revisited once - one spike per site
                // Randomly select element from emptyNeighbours
                uint8_t randomIndex = rand() % atmosphere;
                site = emptyNeighbours[randomIndex];
                path.push_back(site);
                spine += 1;
                n += 1;

                // Prevent index being out of bound for occupied_sites
                if (path.size() >= 3)
                {
                    if (site == path[path.size() - 3])
                    {
                        // Spike has formed, add to list of spikes and adjust the spine accordingly
                        spikes.push_back(path[path.size() - 3]);
                        spikes.push_back(path[path.size() - 2]);
                        m += 1;
                        spine -= 2;
                        spikeLocations[spine] += 1;
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

    spine = 0;
    for (uint32_t i = 0; i < MS2; i++)
    {
        if (i % MS == 0)
            spine++;
        
        if (spine < 2)
            w[i] /= s[0];
        else
            w[i] /= s[0] - (growthConstant * spine);
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
    
    uint32_t MS;
    uint32_t MT;

    std::cout << "Max Length: ";
    std::cin >> MS;
    std::cout << "Max Tours: ";
    std::cin >> MT;

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
