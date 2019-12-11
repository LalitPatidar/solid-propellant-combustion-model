#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include "Reactions.h"

Reactions::Reactions(std::string &line) {
    std::stringstream ss(line);
    std::vector <std::string> x;
    std::string item;

    while(getline(ss,item,' '))
    {
        if (!item.empty())
        {
            x.push_back(item);
        }
    }

    name = x[0];
    for_sym = atof(x[1].c_str());
    for_conv = atof(x[2].c_str());
    back_sym = atof(x[3].c_str());
    back_conv = atof(x[4].c_str());
    img_freq = atof(x[5].c_str());
    dHf = std::max(0.0,atof(x[8].c_str()));
    dHb = std::max(0.0,atof(x[9].c_str()));
    dGf = std::max(0.0,atof(x[10].c_str()));
    dGb = std::max(0.0,atof(x[11].c_str()));
    
    ff = std::max(0.0,atof(x[14].c_str()));
    fb = std::max(0.0,atof(x[15].c_str()));
    
    std::string reactants;
    std::string delimiter = "=";
    size_t pos = 0;
    std::string token;
    while ((pos = x[0].find(delimiter)) != std::string::npos) {
        token = x[0].substr(0, pos);
        reactants = token;
        x[0].erase(0, pos + delimiter.length());
    }
    
    std::string products = x[0];
    
    delimiter = "+";
    pos = 0;
    while ((pos = reactants.find(delimiter)) != std::string::npos) {
        token = reactants.substr(0, pos);
        reactants_list.push_back(token);
        reactants.erase(0, pos + delimiter.length());
    }
    reactants_list.push_back(reactants);
    
    while ((pos = products.find(delimiter)) != std::string::npos) {
        token = products.substr(0, pos);
        products_list.push_back(token);
        products.erase(0, pos + delimiter.length());
    }
    products_list.push_back(products);
    
    if (int(reactants_list.size()) > 2) {std::cout << "Warning: More than 2 reactants in the following reaction:" << name << std::endl;}
    if (int(products_list.size()) > 3) {std::cout << "Warning: More than 3 products in the following reaction:" << name << std::endl;}
}

