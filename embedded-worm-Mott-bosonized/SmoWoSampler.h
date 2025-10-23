#pragma once
#include "PhysSystemSmoWo.h"
#include "Imports.h"

namespace fs = std::filesystem;


class mcSampler {
private:
    int L;
    int B;
    double K;
    double g;
    double mu;
    double lambda_r;
    double lambda_w;
    std::vector<double> Theta_save;

    PhysSystem Phys;
    double running_output_distance;
    double OUTPUT_DISTANCE;
    double TOTAL_NUMBER_SAMPLE;
    int PERCENT_SAVE;
    std::vector<std::string> OBSERVABLES;
    std::string path;
    int data_block_avg;
    int number_events;
    int number_sample;
    std::string dir_root;
    std::string t_init;

    std::unordered_map<std::string, std::vector<std::vector<double>> > data_vector;
    std::unordered_map<std::string, std::vector<double> > data_scalar;
    std::unordered_map<std::string, std::vector<int64_t> > data_int;
    std::unordered_map<std::string, std::string> dic_files;


public:
    inline std::tm localtime_xp(std::time_t timer) {
        std::tm bt{};
#if defined(__unix__)
        localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
        localtime_s(&bt, &timer);
#else
        static std::mutex mtx;
        std::lock_guard<std::mutex> lock(mtx);
        bt = *std::localtime(&timer);
#endif
        return bt;
    }

    inline std::string time_stamp(const std::string& fmt = "%Y_%m_%d_%H_%M_%S")
    {
        auto bt = localtime_xp(std::time(0));
        char buf[64];
        return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
    }


    //constructor
    mcSampler(int L, int B, double K, double g, double mu, double lambda_r, double lambda_w, double TOTAL_NUMBER_SAMPLE, double OUTPUT_DISTANCE, int PERCENT_SAVE, std::vector<std::string> OBSERVABLES, std::string path, int data_block_avg) :
        L{ L },
        B{ B },
        K{ K },
        g{ g },
        mu{ mu },
        lambda_r{ lambda_r },
        lambda_w{ lambda_w },
        data_block_avg{ data_block_avg },
        OUTPUT_DISTANCE{ OUTPUT_DISTANCE },
        TOTAL_NUMBER_SAMPLE{ TOTAL_NUMBER_SAMPLE },
        PERCENT_SAVE{ PERCENT_SAVE },
        OBSERVABLES{ OBSERVABLES },
        path{ path },
        Phys(L, B, K, g, mu, lambda_r, lambda_w), running_output_distance(OUTPUT_DISTANCE), number_events(0), number_sample(0)
    {
        t_init = time_stamp();

        char filepath[200];
        sprintf_s(filepath, sizeof(filepath), "/K_%g/g_%g/mu_%g", K, g, mu);

        dir_root = path + filepath;
        fs::create_directories(dir_root);   //can update the run_id if necessary
    }

    void new_sample() {
        while (true) {
            number_events += 1;

            auto [new_distance, event_type, event_neighbor] = Phys.get_next_event();

            if (new_distance > running_output_distance) {
                Phys.do_lift(running_output_distance);
                Theta_save[Phys.worm_size()] += running_output_distance;
                running_output_distance = OUTPUT_DISTANCE;
                break;
            }
            else {
                running_output_distance -= new_distance;

                Phys.do_lift(new_distance);
                Theta_save[Phys.worm_size()] += new_distance;
                Phys.update_lift_var(event_type, event_neighbor);
            }
        }
    }

    void N_new_sample() {
        auto start = std::chrono::high_resolution_clock::now();
        int percent = 1;

        while (number_sample < TOTAL_NUMBER_SAMPLE) {
            new_sample();
            if (Phys.head_label == Phys.tail_label) { //physical field configuration
                Phys.update_varphi();
                for (const auto& obs : OBSERVABLES) {
                    if (obs == "field") {
                        data_vector["field"].push_back(Phys.field());
                    }
                    if (obs == "C_2kf") {
                        data_scalar["C_2kf"].push_back(Phys.C_2kf());
                    }
                    if (obs == "N_space") {
                        data_scalar["N_space"].push_back(Phys.N_space());
                    }
                    if (obs == "N_time") {
                        data_scalar["N_time"].push_back(Phys.N_time());
                    }
                    if (obs == "kappa") {
                        data_scalar["kappa"].push_back(Phys.kappa());
                    }
                    if (obs == "rho_s") {
                        data_scalar["rho_s"].push_back(Phys.rho_s());
                    }
                    if (obs == "Cx_phi") {
                        data_vector["Cx_phi"].push_back(Phys.Cx_phi());
                    }
                    if (obs == "Ct_phi") {
                        data_vector["Ct_phi"].push_back(Phys.Ct_phi());
                    }
                    if (obs == "algotime") {
                        data_int["algotime"].push_back(Phys.algotime());
                    }
                }
                number_sample++;
            }

            if (number_sample / TOTAL_NUMBER_SAMPLE > percent * PERCENT_SAVE / 100.0) {
                store();
                percent++;
            }
        }
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        std::cout << " 100% (" << std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() << "s)";
        store();
    }

    void init_simu() {
        std::cout << "\n L=" << L << " B=" << B << " K = " << K << " g = " << g << " mu = " << mu << " OUTPUT_DISTANCE=" << OUTPUT_DISTANCE;
        Phys.init_config();
        for (const auto& obs : OBSERVABLES) {
            data_vector[obs].clear();
            if (obs == "field") {
                data_vector["field"].push_back(Phys.field());
            }
            if (obs == "C_2kf") {
                data_scalar["C_2kf"].push_back(Phys.C_2kf());
            }
            if (obs == "N_space") {
                data_scalar["N_space"].push_back(Phys.N_space());
            }
            if (obs == "N_time") {
                data_scalar["N_time"].push_back(Phys.N_time());
            }
            if (obs == "kappa") {
                data_scalar["kappa"].push_back(Phys.kappa());
            }
            if (obs == "rho_s") {
                data_scalar["rho_s"].push_back(Phys.rho_s());
            }
            if (obs == "Cx_phi") {
                data_vector["Cx_phi"].push_back(Phys.Cx_phi());
            }
            if (obs == "Ct_phi") {
                data_vector["Ct_phi"].push_back(Phys.Ct_phi());
            }
            if (obs == "algotime") {
                data_int["algotime"].push_back(Phys.algotime());
            }
            Theta_save.assign(L*B, 0);
        }
        number_sample++;
    }

    void store() {
        for (const auto& obs : OBSERVABLES) {
            std::ofstream outFile;
            char filename[200];
            sprintf_s(filename, sizeof(filename), "g%g_mu%g_K%g_L%d_B%d", g, mu, K, L, B);
            outFile.open(dir_root + "/" + filename + "__" + t_init + "__" + obs + ".data", std::ios_base::app);
            if (obs == "kappa" || obs == "rho_s" || obs == "C_2kf" || obs == "N_time" || obs == "N_space") {  // Scalar observables
                for (const auto& elem : data_scalar[obs]) {
                    outFile << elem << " ";
                }
                data_scalar[obs].clear();
            }
            if (obs == "Cx_phi" || obs == "Ct_phi") {  // Vector observables
                std::vector<std::vector<double>> averages = block_avg(data_vector[obs]);
                for (const auto& vec : averages) {
                    for (size_t i = 0; i < vec.size(); ++i) {
                        outFile << vec[i] << " ";
                    }
                    outFile << std::endl; // New line after each inner vector
                }
                data_vector[obs].clear();
            }
            if (obs == "field") {
                for (const auto& vec : data_vector[obs]) {
                    for (size_t i = 0; i < vec.size(); ++i) {
                        outFile << vec[i] << " ";
                    }
                    outFile << std::endl; // New line after each inner vector
                }
                data_vector[obs].clear();
            }
            if (obs == "algotime") {       // (Large) Int observables
                for (const auto& elem : data_int[obs]) {
                    outFile << elem << " ";
                }
                data_int[obs].clear();
            }
            if (obs == "C_theta") {
                for (const auto& elem : Theta_save) {
                    outFile << elem << " ";
                }
                Theta_save.assign(L*B, 0);
            }
            outFile.close();
        }
    }

    std::vector<std::vector<double>> block_avg(const std::vector<std::vector<double>>& input) {
        std::vector<std::vector<double>> averages;
        int totalSize = input.size();
        for (int i = 0; i < totalSize; i += data_block_avg) {
            size_t end = std::min(i + data_block_avg, totalSize);
            size_t vectorSize = input[0].size();
            std::vector<double> blockAverage(vectorSize, 0.0);
            for (size_t j = i; j < end; ++j) {
                for (size_t k = 0; k < vectorSize; ++k) {
                    blockAverage[k] += input[j][k];
                }
            }
            for (size_t k = 0; k < vectorSize; ++k) {
                blockAverage[k] /= (end - i);
            }
            averages.push_back(blockAverage);
        }
        return averages;
    }
};

