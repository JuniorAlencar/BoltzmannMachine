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

    for (int i = 0; i < n; ++i){
        s[i] = bm.s[i];
        h[i] = bm.h[i];
        bmavs[i] = bm_av_s[i];
    }

    for (int i = 0; i < npairs; ++i){
        bmavss[i] = bm_av_ss[i];
        J[i] = bm.J[i];
    }
    
    for (int i = 0; i < ntriplets; ++i)
        bmavsss[i] = bm_av_sss[i];

    json j;
    j["n_spins"] = bm.n;
    j["n_bonds"] = bm.nbonds;
    j["h"] = h;
    j["J"] = J;

    j["σ"] = s;
    j["S_obs"] = av_s;
    j["SS_obs"] = av_ss;
    j["SSS_obs"] = av_sss;
    j["S"] = bmavs;
    j["SS"] = bmavss;
    j["SSS"] = bmavsss;

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