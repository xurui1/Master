double calc_excess(double **phi,double dr, double volume){
    
    double Omega=0.0;
    
    //need to calculate bulk concentration of block copolymers
    
    for (int i=0;i<Nr;i++){
        Omega+=(phi[0][i]+phi[1][i]+phi[2][i]+phi[3][i]+phi[4][i]- phi_bulk)*dV(i,dr);
    }
    
    Omega/=volume;
    Omega*=dr*(double)Nr;
    
    return Omega;
}