#include "Levitated_Dipole_Equilibrium.h"

//Constructor
Levitated_Dipole_Equilibrium::Levitated_Dipole_Equilibrium() {
    //Default parameters. 
        coil_minor_radius = 1; //minor radius of the levitated dipole. For now we assume the coil is circular.
        coil_major_radius = 2; //major radius of the leviated dipole. 
        coil_height = 3.0; //height of the coil above the floor of the chamber. 
        chamber_width = 4.0; //width of the dipole containment chamber. 
        chamber_height = 6.0; //height of the dipole containment chamber. 

        //Simulation Parameters _-_-_-_-_-
        NR = 100; //number of horizontal grid points.
        NZ = 100; //number of vertical grid points.
        relaxation = 0.3; //The relaxation paramter in the SOR solver.
        tolerance = 1.e-5; //User specified tolerance at which point to stop the solver.

        //Grids needed for the solver _-_-_-_-_-
/*         current_psi_grid = Levitated_Dipole_Equilibrium::current_psi_grid; //This grid stores the most recent value of the flux at each grid point in the domain.  
        current_pressure_grid = Levitated_Dipole_Equilibrium::current_pressure_grid; //This grid stores the most recent pressure at each gridpoint, calculated via pressure(psi). 
        previous_psi_grid = Levitated_Dipole_Equilibrium::previous_psi_grid; //This grid stores the previous value of the flux at each grid point in the domain.  
        previous_pressure_grid = Levitated_Dipole_Equilibrium::previous_pressure_grid; //This grid stores the previous pressure at each gridpoint. 
        label_grid = Levitated_Dipole_Equilibrium::label_grid; //This grid stores the label for each grid point. "0" = interior, "1" = outer boundaries, "2" = inner boundary, "3" = coil, "4" = vacuum
 */
        current_psi_grid.resize(NR, std::vector<double>(NZ, 0.0));
        current_pressure_grid.resize(NR, std::vector<double>(NZ, 0.0));
        previous_psi_grid.resize(NR, std::vector<double>(NZ, 0.0));
        previous_pressure_grid.resize(NR, std::vector<double>(NZ, 0.0));
        label_grid.resize(NR, std::vector<int>(NZ, 0));
};

//Destructor 
Levitated_Dipole_Equilibrium::~Levitated_Dipole_Equilibrium(){};


//Defining the setters:

//Defining the getters:


void Levitated_Dipole_Equilibrium::initialise_label_grid(){
    //This function uses the geometry of the dipole and dipole chamber to work out what kind of grid point each 
    //point should be labeled as. 0 = interior, 1 = outer boundaries, 2 = inner boundary, 3 = coil, 4 = vacuum
    //interior points are those that get looped over in the SOR.
    //outer boundaries are the top, bottom, and outer edge of the chamber, and get set to some reference. 
    //coil is the outer edge of the coil, and is the minimum of flux. 
    //vacuum points get ignored, and are just the inside of the dipole coil. 
    
    //I couldn't think of a neater way of doing this than to just run through all the points in the grid and check what point they are.
    //This definitely could be optimised, but for now it's as fast as I can make it. 

    //loop to deal with the external point.
    double R = 0;
    double Z = 0;
    double DR = chamber_width/NR;
    double DZ = chamber_height/NZ; 

    for(int i; i < NR; i++){
        label_grid[i][0] = 1;
        label_grid[i][NZ-1] = 1;
        label_grid[NR-1][i] = 1;
        label_grid[0][i] = 0;
    };

    for(int i = 1; i < NR; i++){
        R = i*DZ;
        for (int j = 1; j < NZ; j++){
            Z = j*DZ;
            //labeling the coild interior. 
            if ((((R-coil_major_radius)*(R-coil_major_radius)) + ((Z-coil_height)*(Z-coil_height))) <= (coil_minor_radius*coil_minor_radius)){
                label_grid[i][j] = 4;
            } 
            
            //labeling coil boundary. 
            if (abs(((R-coil_major_radius)*(R-coil_major_radius)) + ((Z-coil_height)*(Z-coil_height)) - (coil_minor_radius*coil_minor_radius)) <= ((DR*DR) + (DZ*DZ))){
                label_grid[i][j] = 3;
            }
        }
    }
    for(int i = 0; i < NR; i++){
        R = i*DZ;
        for (int j = 0; j < NZ; j++){
            std::cout << label_grid[i][j];
        }
        std::cout << "\n";
    }
};  

void Levitated_Dipole_Equilibrium::initialise_psi_grid(){
    //This function sets an initial value of psi that we can then iterate on. A suitable choice seems to be a dipole field, so for 
    //now that's just what I'll use. It sets both the previous and current psi grid to be this dipole field. 

};

void Levitated_Dipole_Equilibrium::calculate_pressure(std::vector<std::vector<double> >& generic_psi_grid){
    //This runs through and calculates the pressure at each point in a psi grid that is given to it. It requires that we specify a 
    //a functional form for the dependence of pressure on psi. This can be tailored by the user.

};

void Levitated_Dipole_Equilibrium::initialise_pressure_grid(){
    //This function uses the calculate_pressure function to assign a pressure to each point in the previous and current pressure grids. 

};

void Levitated_Dipole_Equilibrium::single_iteration(){

};

void Levitated_Dipole_Equilibrium::tolerance_check(){

};

void Levitated_Dipole_Equilibrium::solver(){

};


int main(){
    std::cout << "did this shit work?" << "\n";
    Levitated_Dipole_Equilibrium equilibrium;
    std::cout << "did this shit work?" << "\n";
    equilibrium.initialise_label_grid();  
};