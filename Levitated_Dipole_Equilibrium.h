#ifndef LD_SOLVER_H
#define LD_SOLVER_H


#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>


class Levitated_Dipole_Equilibrium{
    private:
        // Geometrical Parameters _-_-_-_-_-
        double coil_minor_radius; //minor radius of the levitated dipole. For now we assume the coil is circular.
        double coil_major_radius; //major radius of the leviated dipole. 
        double coil_height; //height of the coil above the floor. 
        double chamber_width; //width of the dipole containment chamber. 
        double chamber_height; //height of the dipole containment chamber. 

        //Simulation Parameters _-_-_-_-_-
        int NR; //number of horizontal grid points.
        int NZ; //number of vertical grid points.
        double relaxation; //The relaxation paramter in the SOR solver.
        double tolerance; //User specified tolerance at which point to stop the solver.
        double max_iter; //The maximum number of iterations the user is willing to tolerate. 
        double DR; //space between R grid points. 
        double DZ; //space between Z grid points. 

        //plasma parameters
        double psi_max;
        double mu_0;


        //Grids needed for the solver _-_-_-_-_-
        std::vector<std::vector<double> > current_psi_grid; //This grid stores the most recent value of the flux at each grid point in the domain.  
        std::vector<std::vector<double> > current_pressure_grid; //This grid stores the most recent pressure at each gridpoint, calculated via pressure(psi). 
        std::vector<std::vector<double> > previous_psi_grid; //This grid stores the previous value of the flux at each grid point in the domain.  
        std::vector<std::vector<double> > previous_pressure_grid; //This grid stores the previous pressure at each gridpoint. 
        std::vector<std::vector<int> > label_grid; //This grid stores the label for each grid point. "0" = interior, "1" = outer boundaries, "2" = inner boundary, "3" = coil, "4" = vacuum

    public: 
        //The constructor. 
        Levitated_Dipole_Equilibrium();
        //The destructor.
        ~Levitated_Dipole_Equilibrium();

        //The various setters: 
        void set_coil_minor_radius(double new_coil_minor_radius);
        void set_coil_major_radius(double new_coil_major_radius);
        void set_coil_height(double new_coil_height);
        void set_chamber_width(double new_chamber_width);
        void set_chamber_height(double new_chamber_height);
        void set_NR(int new_NR);
        void set_NZ(int new_NZ);
        void set_relaxation(double new_relaxation);
        void set_tolerance(double new_tolerance);

        //The various getters: 
        double get_coil_minor_radius();
        double get_coil_major_radius();
        double get_coil_height();
        double get_chamber_width();
        double get_chamber_height();
        double get_NR();
        double get_NZ();
        double get_relaxation();
        double get_tolerance();        

        //Defining the actual calculation functions. 
        // plasma property calculators. _-_-_-_-_-
        double pressure(); //contains the functional form for out pressure <-> psi relationship. 

        //Setting up the simulation domain _-_-_-_-_-
        void initialise_psi_grid(); //builds up the psi grid based on a vacuum dipole field for our coil. 
        void initialise_pressure_grid(); //builds the pressure grid based on out functional form for p(psi). 
        void initialise_label_grid(); //builds the label grid based on the geometry of the setup. 

        //Iterative solver tools _-_-_-_-_-
        std::vector<std::vector<double> > calculate_pressure(std::vector<std::vector<double> >& generic_psi_grid);
        void single_iteration(); //Runs a single round of the SOR iterator. 
        double tolerance_check(); //Samples the resulting solution to make sure that the next iteration is necessary.
        void solver(); //Runs the SOR solver. 

        //Data Output tools _-_-_-_-_-
        void output_to_txt(std::string file_name, std::vector<std::vector<double> >& generic_grid);
};
#endif