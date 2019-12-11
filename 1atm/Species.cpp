#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include "Species.h"
#include <map>
#include <cmath>

Species::Species(std::string &line1, std::string &line2, std::string &line3, std::string &line4)
{
    std::stringstream ss(line1);
    std::vector <std::string> vline1;
    std::string item;

    while(getline(ss,item,' '))
    {
        if (!item.empty())
        {
            vline1.push_back(item);
        }
    }
    
    name = vline1[0];
    int nC = atoi(vline1[2].substr(0,vline1[2].size()-1).c_str());
    int nH = atoi(vline1[3].substr(0,vline1[3].size()-1).c_str());
    int nN = atoi(vline1[4].substr(0,vline1[4].size()-1).c_str());
    int nO = atoi(vline1[5].substr(0,vline1[5].size()-1).c_str());

    composition.insert(std::make_pair("C",nC));
    composition.insert(std::make_pair("H",nH));
    composition.insert(std::make_pair("N",nN));
    composition.insert(std::make_pair("O",nO));
    
    molecular_weight = 12.0107*nC + 1.0079*nH + 14.0067*nN + 15.9994*nO;

    T_low = atof(vline1[6].c_str());
    T_high = atof(vline1[7].c_str());
    T_mid = atof(vline1[8].c_str());   
   
    coeffs_high[0] = atof(line2.substr(0,15).c_str());
    coeffs_high[1] = atof(line2.substr(15,15).c_str());
    coeffs_high[2] = atof(line2.substr(30,15).c_str());
    coeffs_high[3] = atof(line2.substr(45,15).c_str());
    coeffs_high[4] = atof(line2.substr(60,15).c_str());
    coeffs_high[5] = atof(line3.substr(0,15).c_str());
    coeffs_high[6] = atof(line3.substr(15,15).c_str());
    
    coeffs_low[0] = atof(line3.substr(30,15).c_str());
    coeffs_low[1] = atof(line3.substr(45,15).c_str());
    coeffs_low[2] = atof(line3.substr(60,15).c_str());
    coeffs_low[3] = atof(line4.substr(0,15).c_str());
    coeffs_low[4] = atof(line4.substr(15,15).c_str());
    coeffs_low[5] = atof(line4.substr(30,15).c_str());
    coeffs_low[6] = atof(line4.substr(45,15).c_str());
   
}

double Species::Cp(double T) {
    double R = 1.9872;
    if (T < T_mid) {
        return R*(coeffs_low[0] + coeffs_low[1]*T + coeffs_low[2]*std::pow(T,2) + coeffs_low[3]*std::pow(T,3) + coeffs_low[4]*std::pow(T,4));
    }
    else {
        return R*(coeffs_high[0] + coeffs_high[1]*T + coeffs_high[2]*std::pow(T,2) + coeffs_high[3]*std::pow(T,3) + coeffs_high[4]*std::pow(T,4));
    }       
}

double Species::enthalpy(double T) {
    double R = 1.9872;
    if (T < T_mid) {
        return R*T*(coeffs_low[0] + coeffs_low[1]*T/2.0 + coeffs_low[2]*std::pow(T,2)/3.0 + coeffs_low[3]*std::pow(T,3)/4.0 + coeffs_low[4]*std::pow(T,4)/5.0 + coeffs_low[5]/T);
    }
    else {
        return R*T*(coeffs_high[0] + coeffs_high[1]*T/2.0 + coeffs_high[2]*std::pow(T,2)/3.0 + coeffs_high[3]*std::pow(T,3)/4.0 + coeffs_high[4]*std::pow(T,4)/5.0 + coeffs_high[5]/T);
    }       
}


