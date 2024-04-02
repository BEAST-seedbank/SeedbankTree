#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <libconfig.h++>
#include <unistd.h>
#include <vector>
#include <sstream>

std::string newick(int id, std::vector<int> island, std::vector<double> time_vector, 
        std::vector<int> anc, std::vector<int> dec_1, std::vector<int> dec_2) {
    std::string type;
    if (island[id] == 0) {
        type = "active";
    } else {
        type = "dormant";
    }

    double length;
    if (anc[id] == -1) {
        length = 0.0;
    } else {
        length = time_vector[anc[id]] - time_vector[id];
    }

    std::ostringstream newickStr;
    if (dec_1[id] != -1) {
        newickStr << "(" << newick(dec_1[id], island, time_vector, anc, dec_1, dec_2);
        if (dec_2[id] != -1) {
            newickStr << "," << newick(dec_2[id], island, time_vector, anc, dec_1, dec_2);
        }
        newickStr << ")";
    }
    newickStr << id << "[&type=\"" << type << "\"]:" << length;

    return newickStr.str();
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Call " << argv[0] << " <path to config file>" << std::endl;
        return 1;
    }
    int n, m, model;
    double K, theta_1, theta_2, c, beta;
    libconfig::Config cfg;
    cfg.readFile(argv[1]);
    cfg.lookupValue("model", model);
    cfg.lookupValue("sample_size_1", n);
    cfg.lookupValue("sample_size_2", m);
    cfg.lookupValue("island_2_relative_size", K);
    cfg.lookupValue("migration_rate", c);
    cfg.lookupValue("mutation_rate_1", theta_1);
    cfg.lookupValue("mutation_rate_2", theta_2);
    cfg.lookupValue("weak_seed_bank_mean_delay", beta);
    beta = 1.0 / beta;
    gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gen, time(NULL) * getpid());
    int next_parent = n + m;
    std::vector<int> anc(next_parent, -1);
    std::vector<int> dec_1(next_parent, -1);
    std::vector<int> dec_2(next_parent, -1);
    std::vector<int> island(next_parent, 0);
    std::vector<double> time_vector(next_parent, 0.0);
    double sim_time = 0.0;
    std::vector<int> active_1(n, -1), active_2(m, -1);
    for (int i = 0; i < n; i++) {
        active_1[i] = i;
    }
    for (int i = n; i < n + m; i++) {
        active_2[i - n] = i;
        island[i] = 1;
    }
    double total_rate, c_rate_1, c_rate_2, m_rate_1, m_rate_2, coin;
    int child_1, child_2;
    while (active_1.size() + active_2.size() > 1) {
        switch (model) {
            case 0: c_rate_1 = (double)(active_1.size() 
                        * (active_1.size() - 1)) * beta * beta;
                    c_rate_2 = 0.0;
                    m_rate_1 = 0.0;
                    m_rate_2 = 0.0;
                    break;
            case 1: c_rate_1 = (double)(active_1.size() 
                        * (active_1.size() - 1));
                    c_rate_2 = (double)(active_2.size() 
                        * (active_2.size() - 1)) * K;
                    m_rate_1 = (double)active_1.size() * c;
                    m_rate_2 = (double)active_2.size() * c * K;
                    break;
            case 2: c_rate_1 = (double)(active_1.size() 
                        * (active_1.size() - 1))/2;
                    c_rate_2 = 0.0;
                    m_rate_1 = (double)active_1.size() * c;
                    m_rate_2 = (double)active_2.size() * c * K;
                    break;
            default: std::cout << "Unrecognised model specification." << std::endl;
                    std::cout << "See dev.cfg for implemented models." << std::endl;
                    return 1;
        }
        total_rate = c_rate_1 + c_rate_2 + m_rate_1 + m_rate_2;
        sim_time += gsl_ran_exponential(gen, 1.0 / total_rate);
        coin = gsl_rng_uniform(gen);
        if (coin < c_rate_1 / total_rate) {
            child_1 = floor(active_1.size() * gsl_rng_uniform(gen));
            child_2 = floor(active_1.size() * gsl_rng_uniform(gen));
            while (child_2 == child_1) {
                child_2 = floor(active_1.size() * gsl_rng_uniform(gen));
            }
            time_vector.push_back(sim_time);
            anc.push_back(-1);
            island.push_back(0);
            anc[active_1[child_1]] = next_parent;
            anc[active_1[child_2]] = next_parent;
            dec_1.push_back(-1);
            dec_2.push_back(-1);
            dec_1[next_parent] = active_1[child_1];
            dec_2[next_parent] = active_1[child_2];
            active_1[std::min(child_1, child_2)] = next_parent;
            next_parent++;
            active_1.erase(active_1.begin() + std::max(child_1, child_2));
        } else if (coin < (c_rate_1 + c_rate_2) / total_rate) {
            child_1 = floor(active_2.size() * gsl_rng_uniform(gen));
            child_2 = floor(active_2.size() * gsl_rng_uniform(gen));
            while (child_2 == child_1) {
                child_2 = floor(active_2.size() * gsl_rng_uniform(gen));
            }
            time_vector.push_back(sim_time);
            anc.push_back(-1);
            island.push_back(1);
            anc[active_2[child_1]] = next_parent;
            anc[active_2[child_2]] = next_parent;
            active_2[std::min(child_1, child_2)] = next_parent;
            next_parent++;
            active_2.erase(active_2.begin() + std::max(child_1, child_2));
        } else if (coin < (total_rate - m_rate_2) / total_rate) {
            child_1 = floor(active_1.size() * gsl_rng_uniform(gen));
            anc.push_back(-1);
            island.push_back(1);
            time_vector.push_back(sim_time);
            anc[active_1[child_1]] = next_parent;
            dec_1.push_back(-1);
            dec_2.push_back(-1);
            dec_1[next_parent] = active_1[child_1];
            active_2.push_back(next_parent);
            active_1.erase(active_1.begin() + child_1);
            next_parent++;
        } else {
            child_1 = floor(active_2.size() * gsl_rng_uniform(gen));
            anc.push_back(-1);
            island.push_back(0);
            time_vector.push_back(sim_time);
            anc[active_2[child_1]] = next_parent;
            dec_1.push_back(-1);
            dec_2.push_back(-1);
            dec_1[next_parent] = active_2[child_1];
            active_1.push_back(next_parent);
            active_2.erase(active_2.begin() + child_1);
            next_parent++;
        }
    }
    // std::cout << "Decendant array 1";
    // for (int i = 0; i < next_parent; i++) {
    //     std::cout << "[" << dec_1[i] << "]";
    // }
    // std::cout << "\nDecendant array 2";
    // for (int i = 0; i < next_parent; i++) {
    //     std::cout << "[" << dec_2[i] << "]";
    // }
    for (int i = 0; i < next_parent; i++) {
        std::cout << i << " | " << dec_1[i] << " | " << dec_2[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << newick(next_parent-1, island, time_vector, anc, dec_1, dec_2) << std::endl;
    std::vector<std::vector<int> > hap(anc.size());
    std::vector<int> tmp(1, 0);
    hap.back() = tmp;
    next_parent = 1;
    int mut_count;
    for (int i = (int)anc.size() - 2; i > -1; i--) {
        tmp = hap[anc[i]];
        if (island[i] == 1) {
            mut_count = gsl_ran_poisson(gen, (time_vector[anc[i]] - time_vector[i]) * theta_2);
        } else {
            mut_count = gsl_ran_poisson(gen, (time_vector[anc[i]] - time_vector[i]) * theta_1);
        }
        for (int j = 0; j < mut_count; j++) {
            tmp.push_back(next_parent);
            next_parent++;
        }
        hap[i] = tmp;
    }
    std::vector<std::vector<int> > obs_hap(1);
    obs_hap[0] = hap[0];
    std::vector<int> obs_counts_1(1, 1);
    std::vector<int> obs_counts_2(1, 0);
    int check, j;
    for (int i = 1; i < n; i++) {
        check = 0;
        j = 0;
        while (check == 0 && j < (int)obs_hap.size()) {
            if (hap[i] == obs_hap[j]) {
                obs_counts_1[j]++;
                check = 1;
            }
            j++;
        }
        if (check == 0) {
            obs_hap.push_back(hap[i]);
            obs_counts_1.push_back(1);
            obs_counts_2.push_back(0);
        }
    }
    for (int i = n; i < n + m; i++) {
        check = 0;
        j = 0;
        while (check == 0 && j < (int)obs_hap.size()) {
            if (hap[i] == obs_hap[j]) {
                obs_counts_2[j]++;
                check = 1;
            }
            j++;
        }
        if (check == 0) {
            obs_hap.push_back(hap[i]);
            obs_counts_1.push_back(0);
            obs_counts_2.push_back(1);
        }
    }
    for (unsigned int i = 0; i < obs_hap.size(); i++) {
        for (unsigned int j = 0; j < obs_hap[i].size(); j++) {
            std::cout << obs_hap[i][j] << " ";
        }
        std::cout << obs_counts_1[i] << " " << obs_counts_2[i] << " -1" << std::endl;
        // -1 is an end-of-line character for the accompanying inference algorithm
    }
    gsl_rng_free(gen);
    return 1;
}
