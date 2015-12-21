void set_radius(double *r_0vector, int nradii, double fA, double A, double B){
    
    double curvature=0.64;
    double dcurv = curvature/(double)nradii;
    
    for (int i=0;i<nradii;i++){
        r_0vector[i]=(4.3/curvature)-A-B*fA;
        
        curvature-=dcurv;
        
    }
    
    r_0vector[nradii]=100.0; //garbage entry
    
    
    
}