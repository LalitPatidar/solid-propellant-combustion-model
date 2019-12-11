#ifndef PHASE_H
#define PHASE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include "Species.h"
#include "Reactions.h"

class Phase {
    public:
        // Constructor
        Phase(std::string &file1, std::string &file2);
        
        // Getter setter for T, P, density and viscosity
        double getT();
        void setT(double);
        
        double getP();
        void setP(double);

        double get_density();
        void set_density(double);
        
        double get_viscosity();
        void set_viscosity(double);
        
        double get_conductivity();
        void set_conductivity(double);
        
        // Get number of species, reactions, names and index of species, species molecular weights
        // get reactants and products list for ith reaction
        int get_n_species();
        int get_n_reactions();
        std::string species_name(int);
        std::string reaction_name(int);
        int species_index(std::string &name);
        int reaction_index(std::string &name);
        std::vector <double> get_molecular_weights();
        std::vector <std::string> get_reactants_list(int);
        std::vector <std::string> get_products_list(int);
        
        std::vector <double> get_Eaf();
        std::vector <double> get_Eab();
        
        // Getter setter for X, Y, and getter for concentrations and mean molecular weight
        void setX(std::vector <double>);
        std::vector <double> getX();
        void setY(std::vector <double>);
        std::vector <double> getY();
    
        double get_mean_molecular_weight();
        std::vector <double> get_concentrations();
        
        // getter for temperature dependent data members
        std::vector <double> get_forward_rate_constants();
        std::vector <double> get_backward_rate_constants();
        
        std::vector <double> get_kdiff_f();
        std::vector <double> get_kdiff_b();
        
        std::vector <std::vector<double> > get_kdiff();
        
        std::vector <double> get_net_forward_rate_constants();
        std::vector <double> get_net_backward_rate_constants();
        
        std::vector <double> get_net_rates_of_progress();
        std::vector <double> get_net_rates_of_production();
        
        std::vector <double> get_Cps();
        std::vector <double> get_enthalpies();
        double get_Cpbar();
               
    private:
        // data members based on mechanism file
        // fixed by constructors and are never changed later         
        std::vector <std::string> elements;
        std::vector <Species> species;
        std::vector <Reactions> reactions;
        
        std::vector <double> molecular_weights;
        std::vector <double> All_img_freq;
        
        std::vector <double> All_for_sym;
        std::vector <double> All_for_conv;
        std::vector <double> Eaf;
        
        std::vector <double> All_back_sym;
        std::vector <double> All_back_conv;
        std::vector <double> Eab;
        
        std::vector <std::vector<int> > reactants_stoic_coeffs;
        std::vector <std::vector<int> > products_stoic_coeffs;
        
        // data members defined after creating the object
        // set by respective setters
        double T,P,density,viscosity,conductivity;
        
        // linked data members such as mole fractions X, mass fractions Y and mean molecular weight
        // are linked in setter itself
        // X sets X, mean molecular weight, Y and concentrations
        // Y sets Y, mean molecular weight, X and concentrations
        std::vector <double> X;
        std::vector <double> Y;
        double mean_molecular_weight;
        std::vector <double> concentrations;
        
        // temperature dependent data members
        std::vector <double> forward_rate_constants;
        std::vector <double> backward_rate_constants;
        std::vector <double> kdiff_f;
        std::vector <double> kdiff_b;
        
        std::vector <std::vector< std::vector<double> > > kmat_reactant;
        std::vector <std::vector< std::vector<double> > > kmat_product;
        
        std::vector <std::vector<double> > kdiff;
        
        std::vector <double> Cps;
        std::vector <double> enthalpies;
        
        // temeprature and concentrations dependent data members
        std::vector <double> net_forward_rate_constants;
        std::vector <double> net_backward_rate_constants;
        
        std::vector <double> net_rates_of_progress;
        std::vector <double> net_rates_of_production;
};

#endif
