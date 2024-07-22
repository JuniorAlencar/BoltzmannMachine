#ifndef ENERGY_H
#define ENERGY_H

#include "network.h"
double Energy (Rede &r);
double delta_E (Rede &r, const int s_flip);
double Magnetization (Rede &r);

#endif //ENERGY_H