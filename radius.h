void set_radius(double *r_0vector, int nradii, double fA){
    
    double curvature=0.64;
    double dcurv = curvature/(double)nradii;
    
    for (int i=0;i<nradii;i++){
        if (poly==0){
            r_0vector[i]=(4.3/curvature)-7.05+2.36*fA;
        }
        else if(poly==1){
            r_0vector[i]=(4.3/curvature)-6.97+2.28*fA;
        }
        
        curvature-=dcurv;
        
    }
    
    r_0vector[nradii]=100.0; //garbage entry
    
    
    
}