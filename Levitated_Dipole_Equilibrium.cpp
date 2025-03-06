#include "Levitated_Dipole_Equilibrium.h"

//Constructor
Levitated_Dipole_Equilibrium::Levitated_Dipole_Equilibrium() {
    //Default parameters. 
        coil_minor_radius = 0.01; //minor radius of the levitated dipole. For now we assume the coil is circular.
        coil_major_radius = 1.0; //major radius of the leviated dipole. 
        coil_height = 3.0; //height of the coil above the floor of the chamber. 
        chamber_width = 4.0; //width of the dipole containment chamber. 
        chamber_height = 6.0; //height of the dipole containment chamber. 

        //Simulation Parameters _-_-_-_-_-
        NR = 200; //number of horizontal grid points.
        NZ = 200; //number of vertical grid points.
        relaxation = 0.2; //The relaxation paramter in the SOR solver.
        tolerance = 0.0001; //User specified tolerance at which point to stop the solver.
        max_iter = 1000;

        DZ = chamber_height/NZ;
        DR = chamber_width/NR;
        //plasma parameters
        psi_max = 1.0;
        mu_0 = 4*M_PI*1.e-7;

        //Grids needed for the solver _-_-_-_-_-
/*         current_psi_grid = Levitated_Dipole_Equilibrium::current_psi_grid; //This grid stores the most recent value of the flux at each grid point in the domain.  
        current_pressure_grid = Levitated_Dipole_Equilibrium::current_pressure_grid; //This grid stores the most recent pressure at each gridpoint, calculated via pressure(psi). 
        previous_psi_grid = Levitated_Dipole_Equilibrium::previous_psi_grid; //This grid stores the previous value of the flux at each grid point in the domain.  
        previous_pressure_grid = Levitated_Dipole_Equilibrium::previous_pressure_grid; //This grid stores the previous pressure at each gridpoint. 
        label_grid = Levitated_Dipole_Equilibrium::label_grid; //This grid stores the label for each grid point. 0 = interior, 1 = outer boundary, 2 = inner boundary, 3 = upper boundary, 4 = lower boundary, 5 = coil, 6 = vacuum
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
    //point should be labeled as. 0 = interior, 1 = outer boundary, 2 = inner boundary, 3 = upper boundary, 4 = lower boundary, 5 = coil, 6 = vacuum
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
        label_grid[0][i] = 2;
        label_grid[i][0] = 4;
        label_grid[i][NZ-1] = 3;
        label_grid[NR-1][i] = 1;
    };
    label_grid[0][0] = 2; //fixing the issue where the corner doesn't get labeled. 
    for(int i = 1; i < NR; i++){
        R = i*DZ;
        for (int j = 1; j < NZ; j++){
            Z = j*DZ;
            //labeling the coil interior and a ring around it. 
            if ((((R-coil_major_radius)*(R-coil_major_radius)) + ((Z-coil_height)*(Z-coil_height))) <= (coil_minor_radius*coil_minor_radius)){
                label_grid[i+1][j] = 5;
                label_grid[i-1][j] = 5;
                label_grid[i][j+1] = 5;
                label_grid[i][j-1] = 5;
            } 
        } 
    }

    for(int i = 1; i < NR; i++){
        R = i*DZ;
        for (int j = 1; j < NZ; j++){
            Z = j*DZ;
            if ((((R-coil_major_radius)*(R-coil_major_radius)) + ((Z-coil_height)*(Z-coil_height))) <= (coil_minor_radius*coil_minor_radius)){
                label_grid[i][j] = 6;
            }
        }
    }

    label_grid[0][0] = 6;
    label_grid[0][NZ-1] = 6;
    label_grid[NZ-1][0] = 6;
    label_grid[NZ-1][NZ-1] = 6;


    for(int i = 0; i < NR; i++){
        for (int j = 0; j < NZ; j++){
            std::cout << label_grid[i][j];
        }
        std::cout << "\n";
    }
};  

void Levitated_Dipole_Equilibrium::initialise_psi_grid(){
    //This function sets an initial value of psi that we can then iterate on. A suitable choice seems to be a dipole field, so for 
    //now that's just what I'll use. It sets both the previous and current psi grid to be this dipole field. 
    double R = 0;
    double Z = 0;
    double r = 0; //The distance from the center of the dipole. 
    double cos_theta = 0; //Spherical coordinates definition of theta. 
    double dipole_psi = 0;
    //calculating the dipole field at interior points. 
    for(int i = 0; i < NR; i++){
        R = i*DZ;
        for (int j = 0; j < NZ; j++){
            if((label_grid[i][j] == 0)||(label_grid[i][j] == 2)){
                Z = j*DZ; //Why is there a divergence????
                r = sqrt((R*R) + ((Z-coil_height)*(Z-coil_height)));
                dipole_psi = R/(sqrt(r*r*r)+0.000001); //treating dipole and 4pi prefactor as = 1 for now. TODO: FIX THIS.
                current_psi_grid[i][j] = dipole_psi;
                previous_psi_grid[i][j] = dipole_psi;
            }
        }               
    }
    //Dealing with outer boundary points, where the normal derivative of psi has to go to zero. 
    for(int i; i < NR; i++){
        current_psi_grid[i][0] = current_psi_grid[i][1];
        previous_psi_grid[i][0] = previous_psi_grid[i][1];

        current_psi_grid[i][NZ-1] = current_psi_grid[i][NZ-2];
        previous_psi_grid[i][NZ-1] = previous_psi_grid[i][NZ-2];

        label_grid[NR-1][i] = 1;
        current_psi_grid[NR-1][i] = current_psi_grid[NR-2][i];
        previous_psi_grid[NR-1][i] = previous_psi_grid[NR-2][i];
    };
/*     for(int i = 0; i < NR; i++){
        R = i*DZ;
        for (int j = 0; j < NZ; j++){
            std::cout << current_psi_grid[i][j] << " ";
        }
        std::cout << "\n";
    } */


};

std::vector<std::vector<double> > Levitated_Dipole_Equilibrium::calculate_pressure(std::vector<std::vector<double> >& generic_psi_grid){
    //This runs through and calculates the pressure at each point in a psi grid that is given to it. It requires that we specify a 
    //a functional form for the dependence of pressure on psi. This can be tailored by the user.
    double psi = 0; //A temporary psi value for the computation.
    double pressure = 0; //A temporaty pressure value for the computation.
    std::vector<std::vector<double> > pressure_grid;
    pressure_grid.resize(NR, std::vector<double>(NZ, 0.0));
    for(int i = 0; i < NR; i++){
        for(int j = 0; j < NZ; j++){
            if(label_grid[i][j] == 0){
                psi = generic_psi_grid[i][j];
                pressure = (1 - cos((psi/psi_max)*M_PI)); //just a place holder that vanishes at the edge. If you change
                //this, make sure to change the derivative in the single_iteration function. 
                pressure_grid[i][j] = pressure;
                /* pressure_grid[i][j] = 0; //to test vacuum. */ 

            }
        }
    }
    return(pressure_grid);
};

void Levitated_Dipole_Equilibrium::initialise_pressure_grid(){
    //This function uses the calculate_pressure function to assign a pressure to each point in the previous and current pressure grids. 
    current_pressure_grid = calculate_pressure(current_psi_grid);
    previous_pressure_grid = calculate_pressure(previous_pressure_grid);

/*     for(int i = 0; i < NR; i++){
        for (int j = 0; j < NZ; j++){
            std::cout << current_pressure_grid[i][j] << " ";
        }
        std::cout << "\n";
    } */

};

void Levitated_Dipole_Equilibrium::single_iteration(){
    double gamma = 0;
    double temp_psi = 0;

    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            previous_psi_grid[i][j] = current_psi_grid[i][j];
        }
    for(int i = 1; i < NR; i++){
        for(int j =0; j < NZ; j++){
            if(label_grid[i][j] == 0){
                gamma = -mu_0*(i*DR*i*DR)*((M_PI/psi_max)*sin((previous_psi_grid[i][j]/psi_max))); //pressure derivative term 
  /*               gamma = 0; //to test vacuum.  */
                temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[i+1][j]+previous_psi_grid[i-1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[i][j+1]+previous_psi_grid[i][j-1])) - ((1/(2*i*DR*DR))*(previous_psi_grid[i+1][j] - previous_psi_grid[i-1][j])) - gamma);
                current_psi_grid[i][j] = (relaxation*previous_psi_grid[i][j]) + ((1-relaxation)*(temp_psi)); 
            }
        }
    }
    for(int j = 1; j < NZ-1; j++){
        gamma = -mu_0*(i*DR*i*DR)*((M_PI/psi_max)*sin((previous_psi_grid[i][j]/psi_max))); //pressure derivative term.
/*         gamma = 0; //to test vacuum.  */
        temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[1][j]+previous_psi_grid[1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[0][j+1]+previous_psi_grid[0][j-1])) - gamma);
         current_psi_grid[0][j] = ((1-relaxation)*previous_psi_grid[0][j]) + ((relaxation)*(temp_psi));
    }
    }    
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            previous_pressure_grid[i][j] = current_pressure_grid[i][j];
        }
    }


    current_pressure_grid = calculate_pressure(current_psi_grid);

/* 
    for(int i = 0; i < NR; i++){
        for (int j = 0; j < NZ; j++){
            std::cout << current_pressure_grid[i][j] << " ";
        }
        std::cout << "\n";
    } */
};

double Levitated_Dipole_Equilibrium::tolerance_check(){
    double epsilon = 0;
    double avg_prev_psi = 0; 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            epsilon += abs(current_psi_grid[i][j] - previous_psi_grid[i][j])/(NR*NZ);
            avg_prev_psi += previous_psi_grid[i][j]/(NR*NZ);
        }
    }
    return(epsilon/avg_prev_psi);
};
void Levitated_Dipole_Equilibrium::output_to_txt(std::string file_name, std::vector<std::vector<double> >& generic_grid){
    std::ofstream myfile(file_name); 
    for(int i = 0; i < NR; i++){
        for (int j = 0; j < NZ; j++){
            myfile << generic_grid[i][j] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
};

void Levitated_Dipole_Equilibrium::solver(){
    double epsilon = 1;
    int iteration_number = 0;
    while((epsilon > tolerance) && (iteration_number < max_iter)){
        single_iteration();
        epsilon = tolerance_check();
        iteration_number += 1;
        std::cout << iteration_number << "\n";
        std::cout << epsilon << "\n";
    }
    output_to_txt("pressure.txt", current_pressure_grid);
    output_to_txt("psi.txt", current_psi_grid);
};



int main(){
Levitated_Dipole_Equilibrium equilibrium;
    equilibrium.initialise_label_grid();
    equilibrium.initialise_psi_grid();
    equilibrium.initialise_pressure_grid();
    equilibrium.solver();
};