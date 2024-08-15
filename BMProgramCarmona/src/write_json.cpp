#include "write_json.h"

#include "logger.h"

using json = nlohmann::json;
using namespace std;


void write_json(const std::string &filename,
                const Rede &bm,
                
                const std::vector<double> &av_s,
                const std::vector<double> &av_ss,
                const std::vector<double> &av_sss,
                
                const VecDoub &bm_av_s,
                const VecDoub &bm_av_ss,
                const VecDoub_IO &bm_av_sss)
{
    auto logger = get_logger();
    int n=bm.n;
    int npairs=bm.nbonds;
    int ntriplets=bm_av_sss.size();

    vector<int> s(n,0);
    vector<double> bmavs(n,0), bmavss(npairs,0), bmavsss(ntriplets,0);
    vector<double> h(n,0.0), J(npairs,0.0);

    // h_i and First Moment
    for (int i = 0; i < n; ++i){
        s[i] = bm.s[i];
        h[i] = bm.h[i];
        bmavs[i] = bm_av_s[i];
    }
    
    // Second moment and J_ij
    for (int i = 0; i < npairs; ++i){
        bmavss[i] = bm_av_ss[i];
        J[i] = bm.J[i];
    }
    
    // Third moment
    for (int i = 0; i < ntriplets; ++i)
        bmavsss[i] = bm_av_sss[i];

   

    json j;
    // number of spins
    j["n_spins"] = bm.n;
    j["n_bonds"] = bm.nbonds;
    // h_i
    j["h"] = h;
    // Jij
    j["J"] = J;
    
    j["σ"] = s;
    // First moment exp (magnetization)
    j["S_exp"] = av_s;
    // Second moment exp
    j["SS_exp"] = av_ss;
    // Thir moment exp
    j["SSS_exp"] = av_sss;
    
    // First moment ising (magnetization)
    j["S_ising"] = bmavs;
    // Second moment Ising
    j["SS_ising"] = bmavss;
    // Third momento ising
    j["SSS_ising"] = bmavsss;

    std::ofstream file;
    try
    {
        file.open(filename);
        file << j.dump(4);  // The '4' here is for pretty-printing with an indent of 4 spaces
        file.close();
        logger->info("wrote to {}",filename);
    }
    catch (const std::ofstream::failure& e) {
        logger->error("Exception occurred: {}", e.what());
        exit(1);
    }
}