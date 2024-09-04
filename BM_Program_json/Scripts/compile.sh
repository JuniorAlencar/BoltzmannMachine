#!/bin/bash

# flags to run create folders with libboost
g++ ProcessingData.cpp -o ProcessingData -lboost_filesystem -lboost_system -O3
g++ -O3 BMfinal.cpp -o BMfinal

g++ SpecificHeat.cpp -o SpecificHeat