#ifndef READ_INPUT_JSON_H
#define READ_INPUT_JSON_H

#include "logger.h"
#include <string>
#include <iostream>
#include <nlohmann/json.hpp>

struct Params
{
    /* data */
    std::string run_name;
    std::string input_data;
    std::string input_init_guess;
    std::string output_data;
    std::string output_err;
    bool use_exact;
    int iter_display;
    int iter_max;
    double eta_J;
    double eta_h;
    double min_err_S;
    double min_err_SS;
    int n_rep;
    int t_eq;
    int t_meas;

    Params(std::string rn = "run",
           std::string idata = "",
           std::string iguess = "",
           std::string oerr = "run_err",
           std::string odata = "run_data",
           bool exact = false,
           int idisplay = 100,
           int imax = 10000,
           double etaJ = 0.5,
           double etah = 0.3,
           double minerrS = 1e-3,
           double minerrSS = 1e-5,
           int nrep = 40,
           int teq = 3000,
           int tmeas = 120000) : run_name(rn),
                                 input_data(idata),
                                 input_init_guess(iguess),
                                 output_data(odata),
                                 output_err(oerr),
                                 use_exact(exact),
                                 iter_display(idisplay),
                                 iter_max(imax),
                                 eta_J(etaJ),
                                 eta_h(etah),
                                 min_err_S(minerrS),
                                 min_err_SS(minerrSS),
                                 n_rep(nrep),
                                 t_eq(teq),
                                 t_meas(tmeas){};
    void print() const
    {
        std::cout << "run_name : " << run_name << std::endl;
        std::cout << "input_data : " << input_data << std::endl;
        std::cout << "input_init_guess : " << input_init_guess << std::endl;
        std::cout << "output_data : " << output_data << std::endl;
        std::cout << "output_err : " << output_err << std::endl;
        std::cout << "use_exact : " << use_exact << std::endl;
        std::cout << "iter_display : " << iter_display << std::endl;
        std::cout << "iter_max : " << iter_max << std::endl;
        std::cout << "eta_J : " << eta_J << std::endl;
        std::cout << "eta_h : " << eta_h << std::endl;
        std::cout << "min_err_S : " << min_err_S << std::endl;
        std::cout << "min_err_SS : " << min_err_SS << std::endl;
        std::cout << "n_rep : " << n_rep << std::endl;
        std::cout << "t_eq : " << t_eq << std::endl;
        std::cout << "t_meas : " << t_meas << std::endl;
    }

    void log_info() const
    {
        auto logger = get_logger();

        logger->info("Params: run_name          {}", run_name);
        logger->info("Params: input_data        {}", input_data);
        logger->info("Params: input_init_guess  {}", input_init_guess);
        logger->info("Params: output_data       {}", output_data);
        logger->info("Params: output_err        {}", output_err);
        logger->info("Params: use_exact         {}", use_exact);
        logger->info("Params: iter_display      {}", iter_display);
        logger->info("Params: iter_max          {}", iter_max);
        logger->info("Params: eta_J             {}", eta_J);
        logger->info("Params: eta_h             {}", eta_h);
        logger->info("Params: min_err_S         {}", min_err_S);
        logger->info("Params: min_err_SS        {}", min_err_SS);
        logger->info("Params: n_rep             {}", n_rep);
        logger->info("Params: t_eq              {}", t_eq);
        logger->info("Params: t_meas            {}", t_meas);

        logger->flush();
    }
};

void read_input_json(const std::string &filename, Params &p);

#endif // READ_INPUT_JSON_H