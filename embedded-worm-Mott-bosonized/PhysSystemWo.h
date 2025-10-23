#pragma once

#include "Imports.h"

class PhysSystem {
private:
    inline std::pair<int, int> coordinates(int k) {
        return std::make_pair(k % L, k / L);
    }

    int index(int i, int j) {
        if (i < 0) { i += L; }
        if (j < 0) { j += B; }
        return i % L + (j % B) * L;
    }

    inline double Relu(double x) {
        return std::max(0.0, x);
    }

    double time_quad_link(double d_phi) {
        double b;
        b = dir_vec * d_phi;
        double nu = UniformDouble(gen);
        return b + std::sqrt(Relu(-b) * Relu(-b) - std::log(nu) / inv_K_full);
    }

    std::pair<double, int> time_quad() {
        int i, j;
        std::tie(i, j) = coordinates(lift_label);
        std::vector<double> t_quad;
        std::vector<int> neighbors = { index(i,j + 1), index(i,j - 1), index(i - 1,j), index(i + 1,j) };
        std::vector<int> bs = { Jx[index(i,j)], -Jx[index(i,j - 1)], Jt[index(i - 1,j)], -Jt[index(i,j)] };
        for (size_t k = 0; k < 4; ++k) {
            t_quad.push_back(time_quad_link(bs[k] + f[neighbors[k]] - f[index(i, j)]));
        }
        double t_min = *std::min_element(t_quad.begin(), t_quad.end());
        int neighbor_label = neighbors[std::min_element(t_quad.begin(), t_quad.end()) - t_quad.begin()];
        return std::make_pair(t_min, neighbor_label);
    }

    double time_cos(double t_quad) {
        if (g == 0) {
            return t_quad + 1;
        }
        double nu = UniformDouble(gen);
        double A;
        A = dir_vec * f[lift_label];
        double sig = std::fmod(std::floor(2 * A), 2);
        if (sig < 0) {   //to fix the modulo in C++
            sig += 2.0;
        }
        double z = 2 * A - std::floor(2 * A);
        double B2 = -std::log(nu) / (2 * g) + sig + (1 - sig) / 2 * (1 - std::cos(Pi * z));
        double t_cos = std::floor(B2) - A + std::floor(A) + 1 / (2 * Pi) * std::acos(1 - 2 * (B2 - std::floor(B2)));
        return t_cos;
    }

    double time_refresh(double t_quad) {
        double nu = UniformDouble(gen);
        return -std::log(nu) / lambda_r;
    }

    double time_worm() {
        double nu = UniformDouble(gen);
        return -std::log(nu) / lambda_w;
    }

    double dE_b(int i, int j, int dx, int dt) { //energy difference when the worm goes from (i,j) to (i+dx, j + dt)       
        if (dt == 0) {
            return 2 * J * dx * (Jx[index(i + 1, j)] + dx - f[index(i + 1, j)] + f[index(i + 1, j + 1)]);
        }
        else if (dx == 0) {
            return 2 * J * dt * (Jt[index(i, j + 1)] + dt - f[index(i + 1, j + 1)] + f[index(i, j + 1)] - 2 * K * mu / Pi);
        }
    }

    void metropolis_J_move() {
        if (head_label == tail_label && UniformDouble(gen) < 0.5) {
            head_label = UniformLabel(gen);
            tail_label = head_label;
        }
        else {
            auto [ih, jh] = coordinates(head_label);
            std::vector<std::pair<int, int>> deltas = { {0,-1}, {0,1}, {-1,0}, {1,0} };
            auto [dx, dt] = deltas[UniformChoice(gen)];
            int neighbor_label = index(ih + dx, jh + dt);
            //r accounts for the a priori probability to do a head move
            double r = 1;
            if (head_label == tail_label) { r = 2; }
            else if (neighbor_label == tail_label) { r = 0.5; }
            //compute if metropolis accepted
            auto [i, j] = coordinates(index(ih - Relu(-dx), jh - Relu(-dt))); //index of the moving bond
            if (UniformDouble(gen) < r * std::exp(-dE_b(i, j, dx, dt))) {
                head_label = neighbor_label;
                Jx[index(i + 1, j)] += 2 * dx;
                Jt[index(i, j + 1)] += 2 * dt;
            }
        }
    }

public:
    std::vector<int> Jx;
    std::vector<int> Jt;
    std::vector<double> f;
    std::vector<double> varphi;
    int L;
    int B;
    double K;
    double g;
    double mu;
    double inv_K_full;
    double J;
    double lambda_r;
    double lambda_w;
    int dir_vec = 1;
    int lift_label = 0;
    int head_label = 0;
    int tail_label = 0;
    int64_t algo_time = 0;

    //Random number generator
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> UniformDouble;
    std::uniform_int_distribution<> UniformUnit;
    std::uniform_int_distribution<> UniformLabel;
    std::uniform_int_distribution<> UniformChoice;
    std::uniform_int_distribution<> UniformChoice12;
    std::normal_distribution<> Gaussian;

    // Constructor
    PhysSystem(int L, int B, double K, double g, double mu, double lambda_r, double lambda_w)
        : inv_K_full{ Pi / (8 * K) },
        J{ Pi / (4 * K) },
        L{ L },
        B{ B },
        K{ K },
        g{ g },
        mu{ mu },
        lambda_r{ lambda_r },
        lambda_w{ lambda_w },
        gen(rd()), UniformDouble(0.0, 1.0), UniformUnit(0.0, 1.0), UniformLabel(0.0, L* B - 1), Gaussian(0, 1), UniformChoice(0.0, 3), UniformChoice12(0.0, 11)
    {}

    void init_config() {
        f.assign(L * B, 0.0);
        varphi.assign(L * B, 0.0);
        Jx.assign(L * B, 0);
        Jt.assign(L * B, 0);
        head_label = UniformLabel(gen);
        tail_label = head_label;
        for (size_t i = 0; i < f.size(); ++i) {
            f[i] = Gaussian(gen);
        }
        update_varphi();
    }

    void redistribute_J_f() {
        int N2x = 0;
        int N2t = 0;
        for (size_t j = 0; j < B; ++j) {
            N2t -= Jx[index(0, j)];
        }
        for (size_t i = 0; i < L; ++i) {
            N2x += Jt[index(i, 0)];
        }
        std::vector<int> n(L * B);
        for (size_t i = 0; i < L; ++i) {
            if (i > 0) {
                n[index(i, 0)] = n[index(i - 1, 0)] - Jt[index(i - 1, 0)];
            }
            for (size_t j = 1; j < B; ++j) {
                n[index(i, j)] = n[index(i, j - 1)] + Jx[index(i, j - 1)];
            }
        }
        for (size_t k = 0; k < n.size(); ++k) {
            double n_and_f = n[k] + f[k];
            n[k] = std::lrint(n_and_f);
            f[k] = n_and_f - n[k];
        }
        for (size_t i = 0; i < L; ++i) {
            for (size_t j = 0; j < B - 1; ++j) {
                Jx[index(i, j)] = n[index(i, j + 1)] - n[index(i, j)];
            }
            Jx[index(i, B - 1)] = n[index(i, 0)] - n[index(i, B - 1)] - N2t;
        }
        for (size_t j = 0; j < B; ++j) {
            for (size_t i = 0; i < L - 1; ++i) {
                Jt[index(i, j)] = n[index(i, j)] - n[index(i + 1, j)];
            }
            Jt[index(L - 1, j)] = n[index(L - 1, j)] - n[index(0, j)] + N2x;
        }
    }

    std::tuple<double, std::string, int> get_next_event() {
        auto [t_quad, n_quad] = time_quad();
        std::vector<std::string> event_types = { "quad", "cos", "worm", "refresh" };
        std::vector<double> event_times = { t_quad, time_cos(t_quad), time_worm(), time_refresh(t_quad) };
        std::vector<int> event_neighbors = { n_quad, 0, 0, 0 };
        int min_index = std::min_element(event_times.begin(), event_times.end()) - event_times.begin();
        algo_time += 7;
        return std::make_tuple(event_times[min_index], event_types[min_index], event_neighbors[min_index]);
    }

    void update_lift_var(std::string& event_type, int partner_label) {
        if (event_type == "quad") {
            lift_label = partner_label;
        }
        else if (event_type == "cos") {
            dir_vec *= -1;
        }
        else if (event_type == "refresh") {
            dir_vec = 2 * UniformUnit(gen) - 1;
            lift_label = UniformLabel(gen);
        }
        else if (event_type == "worm") {
            metropolis_J_move();
            algo_time += 1;
        }   
    }

    void do_lift(double distance) {
        algo_time += 1;
        f[lift_label] = f[lift_label] + distance * dir_vec;
        if (head_label == tail_label && std::abs(f[lift_label]) > 1000) {//to avoid annoying precision issues
            redistribute_J_f();
        }
    }

    void update_varphi() {
        int nx = 0;
        for (size_t i = 0; i < L; ++i) {
            nx += Jt[index(i, 0)];
        }
        int nt = 0;
        for (size_t j = 0; j < B; ++j) {
            nt -= Jx[index(0, j)];
        }
        std::vector<int> n(L * B);
        for (size_t i = 0; i < L; ++i) {
            if (i > 0) {
                n[index(i, 0)] = n[index(i - 1, 0)] - Jt[index(i - 1, 0)];
            }
            for (size_t j = 1; j < B; ++j) {
                n[index(i, j)] = n[index(i, j - 1)] + Jx[index(i, j - 1)];
            }
        }
        for (size_t k = 0; k < n.size(); ++k) {
            auto [i, j] = coordinates(k);
            varphi[k] = Pi / 2 * (n[k] + f[k]) + Pi * nx * i / (2 * L) + Pi * nt * j / (2 * B);
        }
    }

    double N_space() {
        int n = 0;
        for (size_t i = 0; i < L; ++i) {
            n += Jt[index(i, 0)];
        }
        return n / 2;
    }

    double N_time() {
        int n = 0;
        for (size_t j = 0; j < B; ++j) {
            n -= Jx[index(0, j)];
        }
        return n / 2;
    }

    double C_2kf() {
        double com_phi = 0.0;
        for (size_t k = 0; k < L * B; ++k) {
            com_phi += varphi[k];
        }
        com_phi /= L * B;
        double C = 0.0;
        for (size_t k = 0; k < L * B; ++k) {
            C += std::cos(2 * (varphi[k] - com_phi));
        }
        C /= L * B;
        return C;
    }

    double kappa() {
        std::complex<double> phi_tilde_q(0.0, 0.0);
        for (size_t k = 0; k < L * B; ++k) {
            auto [i, j] = coordinates(k);
            const std::complex<double> phase(0.0, 2 * Pi / L * i);
            phi_tilde_q += varphi[k] * std::exp(phase);
        }
        double rho_q = 2 / Pi * std::sin(Pi / L) * std::abs(phi_tilde_q);
        return rho_q * rho_q / (B * L);
    }

    double rho_s() {
        std::complex<double> phi_tilde_w(0.0, 0.0);
        for (size_t k = 0; k < L * B; ++k) {
            auto [i, j] = coordinates(k);
            const std::complex<double> phase(0.0, 2 * Pi / B * j);
            phi_tilde_w += varphi[k] * std::exp(phase);
        }
        double rho_w = 2 / Pi * std::sin(Pi / B) * std::abs(phi_tilde_w);
        return rho_w * rho_w / (B * L);
    }

    std::vector<double> Cx_phi() {
        std::vector<double> corr_phi;
        corr_phi.assign(L, 0.0);
        for (int i0 = 0; i0 < L; ++i0) {
            for (int j = 0; j < B; ++j) {
                for (int i = 0; i < L; ++i) {
                    corr_phi[i0] += std::cos(varphi[index(i, j)] - varphi[index(i + i0, j)]);
                }
            }
            corr_phi[i0] /= L * B;
        }
        return corr_phi;
    }

    std::vector<double> Ct_phi() {
        std::vector<double> corr_phi;
        corr_phi.assign(B, 0.0);
        for (int j0 = 0; j0 < B; ++j0) {
            for (int j = 0; j < B; ++j) {
                for (int i = 0; i < L; ++i) {
                    corr_phi[j0] += std::cos(varphi[index(i, j)] - varphi[index(i, j + j0)]);
                }
            }
            corr_phi[j0] /= L * B;
        }
        return corr_phi;
    }

    int worm_size() {
        auto [ih, jh] = coordinates(head_label);
        auto [it, jt] = coordinates(tail_label);
        return index(ih - it, jh - jt);
    }

    int64_t algotime() {
        return algo_time;
    }

    std::vector<double> field() {
        std::vector<int> n(L * B);
        for (size_t i = 0; i < L; ++i) {
            if (i > 0) {
                n[index(i, 0)] = n[index(i - 1, 0)] - Jt[index(i - 1, 0)];
            }
            for (size_t j = 1; j < B; ++j) {
                n[index(i, j)] = n[index(i, j - 1)] + Jx[index(i, j - 1)];
            }
        }
        std::vector<double> phi(L * B);
        for (size_t k = 0; k < n.size(); ++k) {
            phi[k] = Pi / 2 * (n[k] + f[k]);
        }
        return phi;
    }
};