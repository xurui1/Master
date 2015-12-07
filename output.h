void output(double dr, double **phi){
    
    int i;
    
    ofstream outputFile1("./results/phi.dat");
    
    for(i=0;i<Nr;i++){
        outputFile1 <<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outputFile1.close();
    
}

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

void outputkappa(double *kappaM){
    
    ofstream outkappa;
    outkappa.open("./results/kappa.dat");
    
    for (int i=0;i<20;i++){
        outkappa<<0.3+(double)i*0.02<<"  "<<kappaM[i]<<endl;
    }
    
}