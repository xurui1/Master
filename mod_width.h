void mod_width(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr, int nfa){
    
    
    double fE_hom;
    double volume;
    double displacer;
    
    int tempcoord = Coord;
    int tempPoly=poly;
    
    
    //Set coordinate system to cartesian
    Coord=1;
    
    double *fA=create_1d_double_array(nfa,"fA");
    double *AB_width=create_1d_double_array(nfa,"test_rad");
    double *ABA_width=create_1d_double_array(nfa,"ABA_width");
    
    ofstream outputwidth_fa;
    outputwidth_fa.open("./results/outputwidth_fa.dat");
    
    for (int i=0;i<2;i++){
    
        poly=i;
        parameters(chi,f,&ds,Ns,&dr,mu);
        
        int counter=0;
        
        for (int dds=0; dds<=80;dds+=4){
            //Set parameters
            updateparameters(f,Ns,dds);
            fA[counter]=f[0];
            fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
            omega(w);                                       //Initiate omega field
            secant(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,f,2*Nr/5);  //Find tensionless mmb
            volume=vol(dr);                                 //calculate volume
            r_0=1.0;
        
            volume=vol(dr);
            omega(w);
            
            displacer=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,2*Nr/5,0);
            int imax = mmbcentre(phi);
            int imax_left=mmbleft(phi,imax);
            int imax_right=mmbright(phi,imax);
            double width =  dr*((double)imax_right-(double)imax_left);
        
            if (poly==0){
                ABA_width[counter]=width;
            }
            else if (poly==1){
                AB_width[counter]=width;
            }
        
        
            counter++;
        
        }
        
    }
    
    for (int i=0;i<nfa;i++){
        outputwidth_fa<<fA[i]<<" "<<AB_width[i]<<" "<<ABA_width[i]<<endl;
    }
    
    outputwidth_fa.close();
    
    //reset parameters
    Coord=tempcoord;
    poly=tempPoly;
    
}