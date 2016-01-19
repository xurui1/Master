//Here I set the inner radius of the computation cell
void set_radius(double *r_0vector, int nradii, double fA){
    
    double curvature=0.5;
    double dcurv = curvature/(double)nradii;
    
    for (int i=0;i<nradii;i++){
        r_0vector[i]=(4.3/curvature)-6.0;
        
        curvature-=dcurv;
        
    }
    
    r_0vector[nradii]=100.0; //garbage entry
    
    
    
}