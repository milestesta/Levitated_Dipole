#include "Levitated_Dipole_Equilibrium.h"

//source: Downloads/Levitated_Dipole-main/no_pre_refresh
//compile with g++ -std=c++14 -I/opt/homebrew/include Levitated_Dipole_Equilibrium.cpp -o test
//Constructor
Levitated_Dipole_Equilibrium::Levitated_Dipole_Equilibrium(){
    //Default parameters. 
        coil_minor_radius = 0.1; //minor radius of the levitated dipole. For now we assume the coil is circular.
        coil_major_radius = 0.5; //major radius of the leviated dipole. 
        coil_height = 1.5; //height of the coil above the floor of the chamber. 
        chamber_width = 3.0; //width of the dipole containment chamber. 
        chamber_height = 3.0; //height of the dipole containment chamber. 

        //Simulation Parameters _-_-_-_-_-
        NR = 100; //number of horizontal grid points.
        NZ = 100; //number of vertical grid points.
        relaxation = 0.99; //The relaxation paramter in the SOR solver. Multiplies the contribution from the old psi values. 
        tolerance = 1e-6; //User specified percentage tolerance at which point to stop the solver.
        max_iter = 10000;

        DZ = chamber_height/NZ;
        DR = chamber_width/NR;

        //plasma parameters
        psi_max = 0.3; //set to be similar to the psi from the coil alone. 
        mu_0 = 4*M_PI*1e-7;
        pressure_max = 689476.0; //Pascals
        //pressure_max=1000000.0;
        //pressure_max = 0;
        wall_psi = 0;

        //coil parameters:
        coil_current = 1.0e6; //amperes
        coil_points = 100;
        current_element = (coil_current)/(coil_points); //current per unit poloidal length of the ring. 

    

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
        label_grid.resize(NR, std::vector<double>(NZ, 0));
        GS_grid.resize(NR, std::vector<double>(NZ, 0));
        source_grid.resize(NR, std::vector<double>(NZ, 0));
        coil_grid.resize(NR, std::vector<double>(NZ, 0));
        total_flux_grid.resize(NR, std::vector<double>(NZ, 0));
        R_grid.resize(NR);
        Z_grid.resize(NZ);

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
    //outputting coordinate vectors before I forget. 
    for(int i; i < NR; i++){
        R_grid[i] = i*DR;
    }

    for(int j= 0; j < NZ; j++){
        Z_grid[j] = j*DZ;
    }


    std::ofstream Rfile("R_grid.txt"); 
    for(int i = 0; i < NR; i++){
        Rfile << R_grid[i] << " ";
        }
    Rfile << "\n";
    Rfile.close();
    
    std::ofstream Zfile("Z_grid.txt"); 
    for(int j = 0; j < NZ; j++){
        Zfile << Z_grid[j] << " ";
        }
    Zfile << "\n";
    Zfile.close();

    for(int i; i < NR; i++){
        label_grid[0][i] = 2;
        label_grid[i][0] = 4;
        label_grid[i][NZ-1] = 3;
        label_grid[NR-1][i] = 1;
    };
    label_grid[0][0] = 2; //fixing the issue where the corner doesn't get labeled. 
    for(int i = 1; i < NR; i++){
        R = i*DR;
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
        R = i*DR;
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

};  

void Levitated_Dipole_Equilibrium::initialise_psi_grid(){
    //This function sets an initial value of psi that we can then iterate on. A suitable choice seems to be a dipole field, so for 
    //now that's just what I'll use. It sets both the previous and current psi grid to be this dipole field. 
    double R = 0;
    double Z = 0;
    double r = 0; //The distance from the center of the dipole. 
    double cos_theta = 0; //Spherical coordinates definition of theta.
    double sin_theta = 0; 
    double dipole_psi = 0;
    //just breaking the evaluation of A_phi into parts for troubleshooting. 
    double A_part_1_1 = 0;
    double A_part_2_1 = 0;
    double A_part_3_1 = 0;

    double A_part_1_2 = 0;
    double A_part_2_2 = 0;
    double A_part_3_2 = 0;
    double A_tot = 0;
    double k_1 = 0; //needed for the elliptic integrals. 
    double k_2 = 0; //needed for the elliptic integrals.  
    double d_theta = (2*M_PI/coil_points);
    double x_c; //holders for the location on the current ring. 
    double y_c;
    double x_c_2; //a fictitous point for the calculation for -ve coil. 
    double calc_shift = 2*coil_major_radius; //needed to account for the other coil, since k > 0.  
    double R_shifted; 
    double x_c_shifted; 
    //calculating the dipole field at interior points. 
    for(int i = 1; i < NR; i++){
        R = i*DR; //grid point R point.
        R_shifted = R + calc_shift;      
        for (int j = 0; j < NZ; j++){
            Z = (j*DZ); //gridpoint Z location.
            if((label_grid[i][j] == 0)||(label_grid[i][j] == 2)){//calculating the coil contribution.
                for (int p = 0; p < coil_points; p++){ // summing over the coil points.
                    //coil in +ve half . 
                    x_c = coil_major_radius + (coil_minor_radius*cos(p*d_theta)); //R source point on coil. 
                    y_c = coil_height + (coil_minor_radius*sin(p*d_theta)); //Z source point on coil. 

                    //coil in -ve half of plane (since we have are in one large chamber).
                    x_c_2 = calc_shift - coil_major_radius + (coil_minor_radius*cos(p*d_theta)); //R to source point of -ve coil. 




                    // r = sqrt(((R-x)*(R-x))+((Z-y)*(Z-y))); //distance from coil point to grid point. 

                    k_1 = sqrt((4*x_c*R)/(((x_c+R)*(x_c+R)) + ((Z-y_c)*(Z-y_c))));

                    k_2 = sqrt((4*x_c_2*R_shifted)/(((x_c_2+R_shifted)*(x_c_2+R_shifted)) + ((Z-y_c)*(Z-y_c))));

                    //psi is just r*the toroidal magnetic potential.
                    
                    //first coil:
                    A_part_1_1 = (mu_0*R*current_element); //prefactor to convert potential to psi. 
                    A_part_2_1 = (sqrt(R*x_c))/(2*M_PI*k_1);
                    A_part_3_1 = (((2 - (k_1*k_1))*boost::math::ellint_1(k_1)) - (2*boost::math::ellint_2(k_1)));

                    //second coil. 
                    A_part_1_2 = ((-1.0)*mu_0*R*current_element); //prefactor to convert potential to psi. -1 for coil with opposit current. 
                    A_part_2_2 = (sqrt(x_c_2*R_shifted))/(2*M_PI*k_2);
                    A_part_3_2 = (((2 - (k_2*k_2))*boost::math::ellint_1(k_2)) - (2*boost::math::ellint_2(k_2)));


                    A_tot = (A_part_1_1*A_part_3_1*A_part_2_1) +(A_part_1_2*A_part_3_2*A_part_2_2); //at this point this should be the 
                    current_psi_grid[i][j] += A_tot;
                    previous_psi_grid[i][j] += A_tot;
                    coil_grid[i][j] += A_tot;
                }
                
                
            }
        }               
    }
    //Dealing with outer boundary points, where the normal derivative of psi has to go to zero. 
    //Only works for NR = NZ. FIX THAT, IT'S EASY!
    for(int i; i < NR; i++){
/*         current_psi_grid[i][0] = current_psi_grid[i][1];
        previous_psi_grid[i][0] = previous_psi_grid[i][1];

        current_psi_grid[i][NZ-1] = current_psi_grid[i][NZ-2];
        previous_psi_grid[i][NZ-1] = previous_psi_grid[i][NZ-2];

        current_psi_grid[NR-1][i] = current_psi_grid[NR-2][i];
        previous_psi_grid[NR-1][i] = previous_psi_grid[NR-2][i]; */

        current_psi_grid[i][0] = wall_psi;
        previous_psi_grid[i][0] = wall_psi;

        current_psi_grid[i][NZ-1] = wall_psi;
        previous_psi_grid[i][NZ-1] = wall_psi;

        current_psi_grid[NR-1][i] = wall_psi;
        previous_psi_grid[NR-1][i] = wall_psi;
    };
    std::cout << "finished psi init" << "\n";

    //Running iterations to relax /psi to what it will be before adding the plasma. 
    //initial_relaxation();
    output_to_txt("initial_grid.txt", current_psi_grid);


};


//Pressure section. 

double Levitated_Dipole_Equilibrium::pressure_function(double flux){
    return(pressure_max*(1 - cos((flux/psi_max)*M_PI)));
}

double Levitated_Dipole_Equilibrium::pressure_function_derivative(double flux){
    return(pressure_max*M_PI/psi_max)*(sin((flux)*M_PI/psi_max));
}
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
                pressure = pressure_function(psi); //just a place holder that vanishes at the edge. If you change
                // pressure = pressure_max*exp(-(0.1)*((psi-psi_max)*(psi-psi_max)));
                //this, make sure to change the derivative in the single_iteration function. 
                pressure_grid[i][j] = pressure;
            }
        }
    }
    return(pressure_grid);
};

void Levitated_Dipole_Equilibrium::initialise_pressure_grid(){
    //This function uses the calculate_pressure function to assign a pressure to each point in the previous and current pressure grids. 
    current_pressure_grid = calculate_pressure(current_psi_grid);
    previous_pressure_grid = calculate_pressure(previous_pressure_grid);
    std::cout << "finished pressure init" << "\n";    

};

void Levitated_Dipole_Equilibrium::single_iteration(){
    double gamma = 0;
    double temp_psi = 0;

//storing the old grid, abd building the total flux. 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            previous_psi_grid[i][j] = current_psi_grid[i][j];
            total_flux_grid[i][j] = coil_grid[i][j] + current_psi_grid[i][j];
        }

//Calculating psi at the interior points. 
    for(int i = 1; i < NR; i++){
        for(int j =0; j < NZ; j++){
            if(label_grid[i][j] == 0){
                // gamma = -mu_0*(i*DR*i*DR)*((1/M_PI)*sin((current_psi_grid[i][j])/(M_PI*psi_max))); //pressure derivative term 
                gamma = -mu_0*(i*DR*i*DR)*(pressure_function_derivative(total_flux_grid[i][j]));
                temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[i+1][j]+previous_psi_grid[i-1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[i][j+1]+previous_psi_grid[i][j-1])) - ((1/(2*i*DR*DR))*(previous_psi_grid[i+1][j] - previous_psi_grid[i-1][j])) - gamma);
                current_psi_grid[i][j] = (relaxation*previous_psi_grid[i][j]) + ((1-relaxation)*(temp_psi)); 
                source_grid[i][j] = gamma;
            }
        }
    }

//calculating psi along the boundaries. 

//At the upper and lower boundaries. 
//I want the normal of psi to be 0, since the magnetic field shouldn't enter the chamber walls. Therefore I set psi(walls) = 0.

//At the far end boundary. 
    for(int j = 0; j < NZ-1; j++){
        //current_psi_grid[NR-1][j] = current_psi_grid[NR-2][j];
        current_psi_grid[NR-1][j] = wall_psi;
    }
//Along the floor and ceiling:
    for(int i = 0; i < NR-1; i++){
        //current_psi_grid[i][0] = current_psi_grid[i][1]; //floor
        //current_psi_grid[i][NZ-1] = current_psi_grid[i][NZ-2]; //ceiling. 

        current_psi_grid[i][0] = wall_psi; //floor
        current_psi_grid[i][NZ-1] = wall_psi; //ceiling
    }


//at the inner boundary. Using the left right symmetry of the dipole chamber. 
    for(int j = 1; j < NZ-1; j++){
        gamma = -mu_0*(i*DR*i*DR)*(pressure_function_derivative(total_flux_grid[i][j])); //pressure derivative term.
        temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[1][j]+previous_psi_grid[1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[0][j+1]+previous_psi_grid[0][j-1])) - gamma);
        current_psi_grid[0][j] = ((relaxation)*previous_psi_grid[0][j]) + ((1-relaxation)*(temp_psi));
        source_grid[i][j] = gamma;
    }
    }    
//

//updating the previous pressure grid to the new grid. 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            previous_pressure_grid[i][j] = current_pressure_grid[i][j];
        }
    }
//updating the pressure grid based on the new flux. 
    current_pressure_grid = calculate_pressure(total_flux_grid);
};

double Levitated_Dipole_Equilibrium::tolerance_check(){
    double epsilon = 0;
    double max_diff = 0;
    double diff = 0;
    double mean = 0;
    double N_points = 0; 
    double offset = 1e-8; //offset to avoid division by 0 in the tolerance. 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
                diff += ((abs(current_psi_grid[i][j] - previous_psi_grid[i][j])/(abs(current_psi_grid[i][j])+offset))*100); //percentage change 
                N_points += 1; 
            }
        }
    return(diff/N_points);
};

void Levitated_Dipole_Equilibrium::no_plasma_single_iterations(){
    double temp_psi = 0;

//storing the old grid. 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            previous_psi_grid[i][j] = current_psi_grid[i][j];
        }

//Calculating psi at the interior points. 
    for(int i = 1; i < NR; i++){
        for(int j =0; j < NZ; j++){
            if(label_grid[i][j] == 0){
                temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[i+1][j]+previous_psi_grid[i-1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[i][j+1]+previous_psi_grid[i][j-1])) - ((1/(2*i*DR*DR))*(previous_psi_grid[i+1][j] - previous_psi_grid[i-1][j])));
                current_psi_grid[i][j] = (relaxation*previous_psi_grid[i][j]) + ((1-relaxation)*(temp_psi)); 
            }
        }
    }

//calculating psi along the boundaries. 

//At the upper and lower boundaries. 
//I want the normal of psi to be 0, since the magnetic field shouldn't enter the chamber walls. Therefore I set psi(walls) = 0.

//At the far end boundary. 
    for(int j = 0; j < NZ-1; j++){
        //current_psi_grid[NR-1][j] = current_psi_grid[NR-2][j];
        current_psi_grid[NR-1][j] = wall_psi;
    }
//Along the floor and ceiling:
    for(int i = 0; i < NR-1; i++){
        //current_psi_grid[i][0] = current_psi_grid[i][1]; //floor
        //current_psi_grid[i][NZ-1] = current_psi_grid[i][NZ-2]; //ceiling. 

        current_psi_grid[i][0] = wall_psi;
        current_psi_grid[i][NZ-1] = wall_psi;
    }


//at the inner boundary. Using the left right symmetry of the dipole chamber. 
    for(int j = 1; j < NZ-1; j++){
        temp_psi = (1/((2/(DR*DR))+(2/(DZ*DZ))))*(((1/(DR*DR))*(previous_psi_grid[1][j]+previous_psi_grid[1][j])) + ((1/(DZ*DZ))*(previous_psi_grid[0][j+1]+previous_psi_grid[0][j-1])));
         current_psi_grid[0][j] = ((relaxation)*previous_psi_grid[0][j]) + ((1-relaxation)*(temp_psi));
    }
    }    
//
};

void Levitated_Dipole_Equilibrium::initial_relaxation(){
    double epsilon = 1;
    int iteration_number = 2000;
    while((epsilon > tolerance) && (iteration_number < max_iter)){
        no_plasma_single_iterations();
        epsilon = tolerance_check();
        iteration_number += 1;
        std::cout << "Pre-plasma Relaxation Iteration Number: " << iteration_number << "\n";
        std::cout << "Current tolerance: " << epsilon << "\n";
    }
};

void Levitated_Dipole_Equilibrium::GS_check(){
    //breaking the evaluation of the GS into three terms for clarity. 
    double term_1;
    double term_2;
    double term_3;
    double gamma; //pressure source term. 
    for(int i = 0; i < NR; i++){
        for(int j =0; j < NZ; j++){
            if(label_grid[i][j] == 0){ 
                // gamma = -mu_0*(i*DR*i*DR)*(pressure_max*(M_PI/psi_max)*sin((current_psi_grid[i][j]*M_PI/psi_max)))
                gamma = -mu_0*(i*DR*i*DR)*(pressure_function_derivative(total_flux_grid[i][j]));
                term_1 = (current_psi_grid[i+1][j] + current_psi_grid[i-1][j] - (2.0*current_psi_grid[i][j]))/(DR*DR);
                term_2 = (current_psi_grid[i][j+1] + current_psi_grid[i][j-1] - (2.0*current_psi_grid[i][j]))/(DZ*DZ);
                term_3 = (current_psi_grid[i+1][j] - current_psi_grid[i-1][j])/(2.0*DR*i*DR);
                GS_grid[i][j] = term_1 + term_2 - term_3; // - gamma; //should return 0 if the solution is valid.  
            }
        }
    }
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
        std::cout << "Iteration Number: " << iteration_number << "\n";
        std::cout << "Current tolerance: " << epsilon << "\n";
    }
    GS_check();
    output_to_txt("pressure.txt", current_pressure_grid);
    output_to_txt("psi.txt", total_flux_grid);
    output_to_txt("labels.txt", label_grid);
    output_to_txt("GS_check.txt", GS_grid);
    output_to_txt("source_grid.txt", source_grid);
};



int main(){
Levitated_Dipole_Equilibrium equilibrium;
    equilibrium.initialise_label_grid();
    equilibrium.initialise_psi_grid();
    equilibrium.initialise_pressure_grid();
    equilibrium.solver();
};