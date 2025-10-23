#include "Constants.h"

int main() {
    std::cout << "\n path=" << path << "\n PERCENT_SAVE=" << PERCENT_SAVE << "\n TOTAL_NUMBER_SAMPLE=" << TOTAL_NUMBER_SAMPLE << "\n OBSERVABLES=";
    for (const auto& obs : OBSERVABLES) {
        std::cout << " " << obs;
    }

    for (int i = 0; i < Ls.size(); ++i) {
        for (double mu : mus) {
            for (double K : Ks) {
                double OUTPUT_DISTANCE = OUTPUT_DISTANCE_SWEEPS * Ls[i] * Bs[i];
                double lambda_r = lambda_r_prefactor / (Ls[i] * Bs[i]);
                mcSampler Sampler(Ls[i], Bs[i], K, g, mu, lambda_r, lambda_w, TOTAL_NUMBER_SAMPLE, OUTPUT_DISTANCE, PERCENT_SAVE, OBSERVABLES, path, data_block_avg);
                Sampler.init_simu();
                Sampler.N_new_sample();
            }
        }
        std::cout << "\n";
    }
    return 0;
}
