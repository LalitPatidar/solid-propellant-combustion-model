#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <stdexcept>
#include "Phase.h"
#define _USE_MATH_DEFINES

Phase::Phase(std::string &file1, std::string &file2) {
    //-------------------------------------------------------------------------
    // Default initialization, CGS unit is used
    // Should be properly set after creating the Phase object using 
    // respective setters
    //-------------------------------------------------------------------------
    
    T = 298.15;
    P = 1013250.0;
    density = 1.0;
    viscosity = 1.0;
   
    //-------------------------------------------------------------------------
    // Read elements, species and reactions
    //-------------------------------------------------------------------------
    std::ifstream infile1,infile2;
    infile1.open(file1.c_str());

    std::string line;
    
    int flag_elements = 0;
    int flag_species = 0;
    int flag_reactions = 0;
    int line_number = 1;
    std::string line1, line2, line3, line4;
    
    while (getline(infile1,line)){
        if (line=="ELEMENTS") 
        {
            flag_elements = 1;
        }
        if (line=="THERMO") 
        {
            flag_species = 1;
        }
        if (line=="REACTIONS") 
        {
            flag_reactions = 1;
        }
        if (flag_elements==1 && line!= "ELEMENTS" && line!="END")
        {
            elements.push_back(line);
        }
        if (flag_species==1 && line!= "THERMO" && line!="END")
        {
            if (line_number==1){
                line1 = line;
            }
            if (line_number==2){
                line2 = line;
            }
            if (line_number==3){
                line3 = line;
            }
            if (line_number==4){
                line4 = line;
                species.push_back(Species(line1,line2,line3,line4));    
            }
            line_number = line_number + 1;
            if (line_number==5){line_number = 1;}
        }
        if (flag_reactions==1 && line!= "REACTIONS" && line!="END")
        {   
            char line0 = line[0];
            if (line0 != '!') {
                reactions.push_back(Reactions(line));
            }
        }
        if (line=="END")
        {
            flag_elements = 0;
            flag_species = 0;
            flag_reactions = 0;
        }
    }
    //-------------------------------------------------------------------------
    // Check for duplicate species and reactions (very crude as of now)
    // Check unused species or missing species
    //-------------------------------------------------------------------------
    for (int i=0; i< int(species.size());i++) {
        for (int j=0; j< int(species.size());j++) {
            if (j!=i && species[j].name==species[i].name) { 
                std::string errmsg = "Duplicate species found: " + species[j].name;
                throw std::invalid_argument(errmsg);
            }     
        }
    }
    for (int i=0; i< int(reactions.size());i++) {
        for (int j=0; j< int(reactions.size());j++) {
            if (j!=i && reactions[j].name==reactions[i].name) { 
                std::string errmsg = "Duplicate reaction found: " + reactions[j].name;
                throw std::invalid_argument(errmsg);
            }     
        }
    }
    int sfound;
    for (int i=0; i< int(species.size());i++) {
        sfound = 0;
        for (int j=0; j< int(reactions.size());j++) {
            for (int k=0; k< int(reactions[j].reactants_list.size());k++) {
                if (species[i].name == reactions[j].reactants_list[k]) {
                    sfound = 1;
                    break;
                }
                else {;}
            }
            for (int k=0; k< int(reactions[j].products_list.size());k++) {
                if (species[i].name == reactions[j].products_list[k]) {
                    sfound = 1;
                    break;
                }
                else {;}
            }
        }
        if (sfound==0 && species[i].name != "N2") {
            std::string errmsg = "Redundant/unused species found: " + species[i].name;
            throw std::invalid_argument(errmsg);
        }
    }
    
    //-------------------------------------------------------------------------
    // Create a vector of molecular weights
    //-------------------------------------------------------------------------
    for (int i=0; i< int(species.size());i++) {
        molecular_weights.push_back(species[i].molecular_weight);
    }
    //-------------------------------------------------------------------------
    // Create vectors for conversion factors, symmetry factors and 
    // imgainary frequencies and activation barriers of all reactions
    //------------------------------------------------------------------------- 
    
    for (int i=0; i< int(reactions.size());i++) {
        All_img_freq.push_back(reactions[i].img_freq);
        
        All_for_sym.push_back(reactions[i].for_sym);
        All_for_conv.push_back(reactions[i].for_conv);
        
        All_back_sym.push_back(reactions[i].back_sym);
        All_back_conv.push_back(reactions[i].back_conv);
        
        double kcal = 1000.0;
        // Use dHf and dHb for bond breaking reactions
        // which are identified using (dGf < 1e-10 or dGb < 1e-10) and (frequency of TS < 15 cm-1)
        if ((reactions[i].dGf < 1e-10 || reactions[i].dGb < 1e-10) && (reactions[i].img_freq < 15.0)) {
            if (int(reactions[i].reactants_list.size()) == 1 || int(reactions[i].products_list.size()) == 1) {
                Eaf.push_back(reactions[i].dHf + reactions[i].ff*kcal);
                Eab.push_back(reactions[i].dHb + reactions[i].fb*kcal);
            }
            else {
                if (reactions[i].dGf < 1e-10) { 
                    Eaf.push_back(0.0 + reactions[i].ff*kcal); 
                    Eab.push_back(reactions[i].dGb + reactions[i].fb*kcal);
                }
                if (reactions[i].dGb < 1e-10) { 
                    Eaf.push_back(reactions[i].dGf + reactions[i].ff*kcal); 
                    Eab.push_back(0.0 + reactions[i].fb*kcal);
                }
            }
        }
        else {
            Eaf.push_back(reactions[i].dGf + reactions[i].ff*kcal);
            Eab.push_back(reactions[i].dGb + reactions[i].fb*kcal);
        }
    }
    //-------------------------------------------------------------------------
    // Create vectors for reactants and products stoichiometric coeffs 
    //-------------------------------------------------------------------------
    reactants_stoic_coeffs.resize(reactions.size());
    products_stoic_coeffs.resize(reactions.size());
    
    for (int i=0; i < int(reactions.size()); i++){
        reactants_stoic_coeffs[i].resize(species.size());
        for (int j=0; j < int(reactions[i].reactants_list.size()); j++) {
            reactants_stoic_coeffs[i][species_index(reactions[i].reactants_list[j])] = reactants_stoic_coeffs[i][species_index(reactions[i].reactants_list[j])] + 1;
        }
        
        products_stoic_coeffs[i].resize(species.size());
        for (int j=0; j < int(reactions[i].products_list.size()); j++) {
            products_stoic_coeffs[i][species_index(reactions[i].products_list[j])] = products_stoic_coeffs[i][species_index(reactions[i].products_list[j])] + 1;
        }
    }
    //-------------------------------------------------------------------------
    // Read species radius 
    // Read in as angstrom and converted to cm
    //-------------------------------------------------------------------------
    infile2.open(file2.c_str());
    std::string name = " ";
    std::string rad;
    double rad1;
    std::vector<std::string> temp_names;
    std::vector<double> temp_radius;
    
    line_number = 1;
    
    while (getline(infile2,line)){
        if (line_number > 2) {
            name = line.substr(0,line.find(' '));
            
            rad = line.substr(line.find("mol)"),line.find("angstrom")-line.find("mol)"));
            rad1 = atof(rad.substr(4).c_str());   
        }
        line_number = line_number + 1;
        temp_names.push_back(name);
        temp_radius.push_back(rad1);        
    }
    
    for (int i=0; i < int(species.size()); i++) {
        for (int j=0; j < int(temp_names.size()); j++) {
            if (species[i].name == temp_names[j]) {
                species[i].radius = 1.0e-8*temp_radius[j];
            }
        }
    }   
    //-------------------------------------------------------------------------
    // create matrix for diffusion rate constants coefficients 
    // each element of kmat is kdiff divided by T/viscosity. The T/viscosity
    // factor is multiplied to kmat in calculating kdiff 
    //-------------------------------------------------------------------------
    kmat_reactant.resize(reactions.size());
    kmat_product.resize(reactions.size());
    
    //Boltzmann constant kB in CGS units and Avagadro number NA
    double kB = 1.38064852e-16;
    double NA = 6.022140857e23;
    double Rsum;
    
    int max_reactants = 0;
    int max_products = 0;
    
    for (int i=0; i <int(reactions.size()); i++) {
        if (int(reactions[i].reactants_list.size()) > max_reactants) { max_reactants = reactions[i].reactants_list.size();}
        if (int(reactions[i].products_list.size()) > max_products) { max_products = reactions[i].products_list.size();}
    }

    if (max_reactants > 2 || max_products > 3) {
        std::string errmsg = "Maximum number of reactants >2 or products>3: Check Warnings";
        throw std::invalid_argument(errmsg);
    }  

    for (int i=0; i <int(reactions.size()); i++) {
        kmat_reactant[i].resize(max_reactants-1);
        kmat_product[i].resize(max_products-1);
        if (reactions[i].reactants_list.size() > 1) {
            Rsum = 0.0;
            for (int j=0; j < int(reactions[i].reactants_list.size()-1); j++) {
                kmat_reactant[i][j].resize(2);
                Rsum = Rsum + species[species_index(reactions[i].reactants_list[j])].radius;
                kmat_reactant[i][j][0] = (2.0*NA*kB/3.0)*(Rsum + species[species_index(reactions[i].reactants_list[j+1])].radius)*(1.0/Rsum + 1.0/species[species_index(reactions[i].reactants_list[j+1])].radius);
                kmat_reactant[i][j][1] = (kB/2.0/M_PI)*(1.0/Rsum + 1.0/species[species_index(reactions[i].reactants_list[j+1])].radius)/(std::pow((Rsum + species[species_index(reactions[i].reactants_list[j+1])].radius),2.0));
            }
        }

        if (reactions[i].products_list.size() > 1) {
            Rsum = 0.0;
            for (int j=0; j < int(reactions[i].products_list.size()-1); j++) {
                kmat_product[i][j].resize(2);
                Rsum = Rsum + species[species_index(reactions[i].products_list[j])].radius;
                kmat_product[i][j][1] = (2.0*NA*kB/3.0)*(Rsum + species[species_index(reactions[i].products_list[j+1])].radius)*(1.0/Rsum + 1.0/species[species_index(reactions[i].products_list[j+1])].radius);
                kmat_product[i][j][0] = (kB/2.0/M_PI)*(1.0/Rsum + 1.0/species[species_index(reactions[i].products_list[j+1])].radius)/(std::pow((Rsum + species[species_index(reactions[i].products_list[j+1])].radius),2.0));
            }
        }
    }
    /*
    //-------------------------------------------------------------------------
    // Create vectors for diffusion rate constants - temporary 
    //-------------------------------------------------------------------------
    kdiff_f.resize(reactions.size());
    for (int i=0;i<int(reactions.size()); i++) {
        if (int(reactions[i].reactants_list.size())==1) { kdiff_f[i]=1e30;}
        if (int(reactions[i].reactants_list.size())==2) { kdiff_f[i]=5e8;}
        if (int(reactions[i].reactants_list.size())==3) { kdiff_f[i]=5e11;}
    }
    
    kdiff_b.resize(reactions.size());
    for (int i=0;i<int(reactions.size()); i++) {
        if (int(reactions[i].products_list.size())==1) { kdiff_b[i]=1e30;}
        if (int(reactions[i].products_list.size())==2) { kdiff_b[i]=5e8;}
        if (int(reactions[i].products_list.size())==3) { kdiff_b[i]=5e11;}
        if (int(reactions[i].products_list.size())==4) { kdiff_b[i]=5e14;}
    }*/
    //-------------------------------------------------------------------------
    // Modifications of activation barriers
    // Hardcoded section, comment out the whole secction to use the default
    // reaction rate parameters as provided in mechanism file  
    //-------------------------------------------------------------------------
    double kcal = 1000.0;
    std::string rxn;
    
    for (int i=0;i<int(reactions.size()); i++) {
        Eaf[i] += 0.0*kcal;
        Eab[i] += 0.0*kcal;
    }
            
}

//-----------------------------------------------------------------------------
// Getter setter for T, P, density and viscosity
//-----------------------------------------------------------------------------
double Phase::getT() {
    return T;
}

void Phase::setT(double value) {
    T = value;
}

double Phase::getP() {
    return P;
}

void Phase::setP(double value) {
    P = value;
}

double Phase::get_density() {
    return density;
}

void Phase::set_density(double value) {
    density = value;
}

double Phase::get_viscosity() {
    return viscosity;
}

void Phase::set_viscosity(double value) {
    viscosity = value;
}

double Phase::get_conductivity() {
    return conductivity;
}

void Phase::set_conductivity(double value) {
    conductivity = value;
}

//-----------------------------------------------------------------------------
// Getter for number of species and reactions, 
// and names and index of species, species molecular weights
// get reactants and products list for ith reaction
//-----------------------------------------------------------------------------
int Phase::get_n_species() {
    return species.size();
}

int Phase::get_n_reactions() {
    return reactions.size();
}

std::string Phase::species_name(int k) {
    return species[k].name;
}

std::string Phase::reaction_name(int k) {
    return reactions[k].name;
}

int Phase::species_index(std::string &input_name) {
    for(int indx=0; indx < int(species.size()); indx++) {
        if(species[indx].name == input_name) { return indx; }   
    }
    std::string errmsg = "Species missing: " + input_name;
    throw std::invalid_argument(errmsg);
    return -1;
}

int Phase::reaction_index(std::string &input_name) {
    for(int indx=0; indx < int(reactions.size()); indx++) {
        if(reactions[indx].name == input_name) { return indx; }   
    }
    std::string errmsg = "Reaction missing: " + input_name;
    throw std::invalid_argument(errmsg);
    return -1;
}

std::vector <double> Phase::get_molecular_weights() {
    return molecular_weights;
}

std::vector <std::string> Phase::get_reactants_list(int k) {
    return reactions[k].reactants_list;
}

std::vector <std::string> Phase::get_products_list(int k) {
    return reactions[k].products_list;
}

std::vector <double> Phase::get_Eaf() {
    return Eaf;
}

std::vector <double> Phase::get_Eab() {
    return Eab;
}
//-----------------------------------------------------------------------------
// Getter setter for X,Y, 
// and getters for concentrations and mean molecular weight
//-----------------------------------------------------------------------------
void Phase::setX(std::vector <double> values) {
    X = values;
    
    mean_molecular_weight = 0.0;
    for (int i = 0; i < int(species.size()); i++) {
        mean_molecular_weight = mean_molecular_weight + values[i]*molecular_weights[i];
    }
    
    Y.clear();
    for (int i = 0; i < int(species.size()); i++) {
        Y.push_back(values[i]*molecular_weights[i]/mean_molecular_weight);
    }
    
    concentrations.clear();
    for (int i = 0; i < int(species.size()); i++) {
        concentrations.push_back(X[i]*density/mean_molecular_weight);
    }
}

std::vector <double> Phase::getX() {
    return X;
}

void Phase::setY(std::vector <double> values) {
    Y = values;
    
    double sum = 0.0;
    for (int i = 0; i < int(species.size()); i++) {
        sum = sum + values[i]/molecular_weights[i];
    }
    mean_molecular_weight = 1.0/sum;
    
    X.clear();
    for (int i = 0; i < int(species.size()); i++) {
        X.push_back(values[i]*mean_molecular_weight/molecular_weights[i]);
    }
    
    concentrations.clear();
    for (int i = 0; i < int(species.size()); i++) {
        concentrations.push_back(X[i]*density/mean_molecular_weight);
    }    
}

std::vector <double> Phase::getY() {
    return Y;
}


double Phase::get_mean_molecular_weight() {
    return mean_molecular_weight;
}

std::vector <double> Phase::get_concentrations() {
    return concentrations;
}

//-----------------------------------------------------------------------------
// Getter setter for temperature dependent data members
//-----------------------------------------------------------------------------
std::vector <double> Phase::get_forward_rate_constants() {
    forward_rate_constants.clear();
    
    double h_planck = 6.626e-27;
    double kB = 1.38064852e-16;
    double R = 1.9872036;
    
    double wigner;
    
    for (int i = 0; i < int(reactions.size()); i++) {
        wigner = 1.0 + (1.0/24.0)*std::pow((h_planck*All_img_freq[i]*2.99792458e10/kB/T),2.0);
        forward_rate_constants.push_back((All_for_conv[i]*All_for_sym[i]*wigner*kB*T/h_planck)*(exp(-Eaf[i]/(R*T))));
    }
    return forward_rate_constants;
}

std::vector <double> Phase::get_backward_rate_constants() {
    backward_rate_constants.clear();
    
    double h_planck = 6.626e-27;
    double kB = 1.38064852e-16;
    double R = 1.9872036;
    
    double wigner;
    
    for (int i = 0; i < int(reactions.size()); i++) {
        wigner = 1.0 + (1.0/24.0)*std::pow((h_planck*All_img_freq[i]*2.99792458e10/kB/T),2.0);
        backward_rate_constants.push_back((All_back_conv[i]*All_back_sym[i]*wigner*kB*T/h_planck)*(exp(-Eab[i]/(R*T))));
    }
    return backward_rate_constants;
}

std::vector <double> Phase::get_kdiff_f() {
    return kdiff_f;
}

std::vector <double> Phase::get_kdiff_b() {
    return kdiff_b;
}

std::vector <std::vector<double> > Phase::get_kdiff() {
    // Only handles 2 reactants and 3 products
    // Need modifications to handle more reactants and products
    // make sure mechanism file has maximum of 2 reactants and 3 products
    kdiff.resize(reactions.size());
    
    std::vector <std::vector<double> > k_prod;
    k_prod.resize(reactions.size());
    std::vector <std::vector<double> > k_reac;
    k_reac.resize(reactions.size());
    
    double conc_E;
    
    for (int i=0; i < int(reactions.size()); i++) {
        k_prod[i].resize(2);
        if (int(reactions[i].products_list.size()) == 1) {
            k_prod[i][0] = 1e30;
            k_prod[i][1] = 1e30;
        }
        else if (int(reactions[i].products_list.size()) == 2) {
            k_prod[i][0] = kmat_product[i][0][0]*T/viscosity;
            k_prod[i][1] = kmat_product[i][0][1]*T/viscosity;
        }
        else {
            conc_E = concentrations[species_index(reactions[i].products_list[2])];
            k_prod[i][0] = (T/viscosity)*kmat_product[i][0][0]*kmat_product[i][1][0]/(kmat_product[i][1][0] + kmat_product[i][0][1]*conc_E);
            k_prod[i][1] = (T/viscosity)*kmat_product[i][0][1]*kmat_product[i][1][1]/(kmat_product[i][1][0] + kmat_product[i][0][1]*conc_E);
        }
    }
    
    for (int i=0; i < int(reactions.size()); i++) {
        k_reac[i].resize(2);
        if (int(reactions[i].reactants_list.size()) == 1) {
            k_reac[i][0] = 1e30;
            k_reac[i][1] = 1e30;
        }
        else {
            k_reac[i][0] = kmat_reactant[i][0][0]*T/viscosity;
            k_reac[i][1] = kmat_reactant[i][0][1]*T/viscosity;
        }
    }
    
    for (int i=0; i < int(reactions.size()); i++) {
        kdiff[i].resize(2);
        if (int(reactions[i].reactants_list.size()) == 1 && int(reactions[i].products_list.size()) == 1) {
            kdiff[i][0] = 1e30;
            kdiff[i][1] = 1e30;
        }
        else if (int(reactions[i].reactants_list.size()) == 1 && int(reactions[i].products_list.size()) > 1) {
            kdiff[i][0] = k_prod[i][0];
            kdiff[i][1] = k_prod[i][1];
        }
        else if (int(reactions[i].reactants_list.size()) > 1 && int(reactions[i].products_list.size()) == 1) {   
            kdiff[i][0] = k_reac[i][0];
            kdiff[i][1] = k_reac[i][1];
        }
        else {
            kdiff[i][0] = k_reac[i][0]*k_prod[i][0]/(k_prod[i][0]+k_reac[i][1]);
            kdiff[i][1] = k_reac[i][1]*k_prod[i][1]/(k_prod[i][0]+k_reac[i][1]);
        }
    }
    return kdiff;
}

std::vector <double> Phase::get_net_forward_rate_constants() {
    net_forward_rate_constants.clear();
    std::vector <double> frk = get_forward_rate_constants();
    std::vector <std::vector<double> > kd = get_kdiff();
    
    double kf_net;
    double limit_high = 1e18;
    double limit_low = 1e-18;  
     
    for (int i = 0; i < int(reactions.size()); i++) {
        //net_forward_rate_constants.push_back(1.0/((1.0/kd[i][0])+(1.0/frk[i])));
        kf_net = 1.0/((1.0/kd[i][0])+(1.0/frk[i]));
        if (kf_net > limit_high) { net_forward_rate_constants.push_back(limit_high); }
        else if (kf_net < limit_low) { net_forward_rate_constants.push_back(limit_low); }
        else { net_forward_rate_constants.push_back(kf_net); }
    }
    return net_forward_rate_constants;
}

std::vector <double> Phase::get_net_backward_rate_constants() {
    net_backward_rate_constants.clear();
    std::vector <double> brk = get_backward_rate_constants();
    std::vector <std::vector<double> > kd = get_kdiff();
    
    double kb_net;
    double limit_high = 1e18;
    double limit_low = 1e-18;
       
    for (int i = 0; i < int(reactions.size()); i++) {
        //net_backward_rate_constants.push_back(1.0/((1.0/kd[i][1])+(1.0/brk[i])));
        kb_net = 1.0/((1.0/kd[i][1])+(1.0/brk[i]));
        if (kb_net > limit_high) { net_backward_rate_constants.push_back(limit_high); }
        else if (kb_net < limit_low) { net_backward_rate_constants.push_back(limit_low); }
        else { net_backward_rate_constants.push_back(kb_net); }
    }
    return net_backward_rate_constants;
}

std::vector <double> Phase::get_net_rates_of_progress() {
    net_rates_of_progress.clear();
    double prodp, prodr;
    
    std::vector <double> nfrk = get_net_forward_rate_constants();
    std::vector <double> nbrk = get_net_backward_rate_constants();
    
    for (int i=0; i < int(reactions.size()); i++) {
        prodr = 1.0;
        prodp = 1.0;
        for (int j=0; j < int(species.size()); j++) {
            prodr = prodr*std::pow(concentrations[j],reactants_stoic_coeffs[i][j]);
            prodp = prodp*std::pow(concentrations[j],products_stoic_coeffs[i][j]);
        }
        net_rates_of_progress.push_back(nfrk[i]*prodr - nbrk[i]*prodp);     
    }
    return net_rates_of_progress;
}

std::vector <double> Phase::get_net_rates_of_production() {
    net_rates_of_production.clear();
    std::vector <std::vector<int> > vpr;
    vpr.resize(reactions.size());
    
    for (int i=0; i<int(reactions.size());i++) {
        vpr[i].resize(species.size());
        for (int j=0; j<int(species.size());j++) {
            vpr[i][j] = (products_stoic_coeffs[i][j] - reactants_stoic_coeffs[i][j]);
        }
    }
    std::vector <double> q = get_net_rates_of_progress();
    double sum;
    
    for (int j=0; j < int(species.size()); j++) {
        sum = 0.0;
        for (int i=0; i < int(reactions.size()); i++) {
            sum = sum + vpr[i][j]*q[i];    
        }
        net_rates_of_production.push_back(sum);
    }
    return net_rates_of_production;    
}

std::vector <double> Phase::get_Cps() {
    Cps.clear();
    for (int i=0; i< int(species.size()); i++) {
        Cps.push_back(species[i].Cp(T));
    }
    return Cps;
}

std::vector <double> Phase::get_enthalpies() {
    enthalpies.clear();
    for (int i=0; i< int(species.size()); i++) {
        enthalpies.push_back(species[i].enthalpy(T));
    }
    return enthalpies;
}

double Phase::get_Cpbar() {
    std::vector <double> Cps = get_Cps();
    double sum = 0.0;
    for (int i=0; i < int(species.size()); i++) {
        sum = sum + Cps[i]*Y[i];
    }
    return sum;
}
