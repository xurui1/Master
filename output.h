/****************Here I output the concentration profile***********************/
void output(double dr, double **phi){
    
    int i;
    
    ofstream outputFile1("./results/phi.dat");
    
    for(i=0;i<Nr;i++){
        outputFile1 <<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outputFile1.close();
    
}

/*******************Here I output concentration profile for various radii***********************/
void outputphi(double **phi, double dr){
    
    ofstream outphi;
    string filename3;
    filename3="./results/phi/phi_r" + DoubleToStr(r_0)+ ".dat";
    outphi.open(filename3.c_str());
    
    for (int i=0;i<Nr;i++){
        outphi<<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outphi.close();
    
}

/*************************Here I output quadratic and quartic fit parameters*******************************/
void outputkappa(double *a1, double *a2, double *a3, double *a4, double *a5, double *a6, int nfa, double *chi){
    
    ofstream outkappa;
    string filename;
    
    if (Coord==2){
        if (poly==0){
            filename="./results/ABAcyl"+DoubleToStr(chi[0])+".dat";
            outkappa.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/ABcyl"+DoubleToStr(chi[0])+".dat";
            outkappa.open(filename.c_str());
        }
    }
    else if (Coord==3){
        if (poly==0){
            filename="./results/ABAsph"+DoubleToStr(chi[0])+".dat";
            outkappa.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/ABsph"+DoubleToStr(chi[0])+".dat";
            outkappa.open(filename.c_str());
        }
    }
    
    for (int i=0;i<nfa;i++){
        outkappa<<0.3+(double)i*0.02<<"  "<<a1[i]<<" "<<a2[i]<<" "<<a3[i]<<" "<<a4[i]<<" "<<a5[i]<<" "<<a6[i]<<endl;
    }
    
}