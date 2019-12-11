#ifndef REACTIONS_H
#define REACTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>

class Reactions {
    public:
        Reactions(std::string &line);
        std::string name;
        double for_sym, for_conv, back_sym, back_conv, img_freq, dHf, dHb, dGf, dGb, ff, fb;
        std::vector <std::string> reactants_list;
        std::vector <std::string> products_list;      
};

#endif
