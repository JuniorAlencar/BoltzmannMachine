#include "read_json.h"
#include "logger.h"

using json = nlohmann::json;
using namespace std;


Rede read_json(const std::string filename) {
        auto logger = get_logger();
    try {

        // Open and parse the JSON file
        std::ifstream file(filename.c_str());
        json j;
        file >> j;
        file.close();
        int nspins = j["n_spins"];

	
        Rede bm(nspins, 0, 0, 0, 0, 0);


        bm.n = j["n_spins"];
        bm.nbonds = j["n_bonds"];
        vector<double> h = j["h"];
        vector<double> J = j["J"];
        vector<int> s = j["σ"];

        bm.h = h;
        bm.J = J;
        bm.s = s;

        logger->info("read_json: read {} spins and {} bonds from {}", bm.n, bm.nbonds, filename);

        // cout << "σs" << endl;
        // for (int i =0; i< nspins; ++i)
        //     cout << bm.s[i] << endl;

        // cout << "hs" << endl;
        // for (int i =0; i< nspins; ++i)
        //     cout << bm.h[i] << endl;
        
        // cout << "Js" << endl;
        // for (int i =0; i< bm.nbonds; ++i)
        //     cout << bm.J[i] << endl;
        
        return bm;

    } catch (const json::parse_error& e) {
        logger->error("read_json: JSON parsing error: {}", e.what());
        exit(1);
    } catch (const std::ifstream::failure& e) {
        logger->error("read_json: File I/O error: {}", e.what());
        exit(1);
    } catch (...) {
        logger->error("read_json: An unknown error occurred");
        exit(1);
    }

}
