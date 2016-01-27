void mod_main(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double ds,double *chi,double dr,int nfa,double *A, double *B, double *C,double *D, double *E, double *F, int nradii, double *dFE, double *mu_vec){
    
    double *r_0vector=create_1d_double_array(nradii+1, "r_0vector");
    double *Rad=create_1d_double_array(nradii, "Rad");                  //Radius for fitting
    double *Rad2=create_1d_double_array(nradii, "Rad2");                  //Radius for fitting

    double *diameter=create_1d_double_array(nradii, "diameter");          //membrane diameter
    double *Curv=create_1d_double_array(nradii, "Curvsq");            //Curv sq. for fitting
    double *Curvsq=create_1d_double_array(nradii, "Curvsq");            //Curv sq. for fitting
    
    double *Curv2=create_1d_double_array(nradii, "Curvsq2");            //Curv sq. for fitting
    double *Curvsq2=create_1d_double_array(nradii, "Curvsq2");            //Curv sq. for fitting
    
    double *a1=create_1d_double_array(nfa, "a1");            //bending moduli const term
    double *a2=create_1d_double_array(nfa, "a2");            //bending moduli linear term
    double *a3=create_1d_double_array(nfa, "a3");            //bending moduli quad term
    double *a4=create_1d_double_array(nfa, "a4");            //bending moduli const term
    double *a5=create_1d_double_array(nfa, "a5");            //bending moduli quad term
    double *a6=create_1d_double_array(nfa, "a6");            //bending moduli quart term
    
    double *b1=create_1d_double_array(nfa, "a1");            //bending moduli const term
    double *b2=create_1d_double_array(nfa, "a2");            //bending moduli linear term
    double *b3=create_1d_double_array(nfa, "a3");            //bending moduli quad term
    double *b4=create_1d_double_array(nfa, "a4");            //bending moduli const term
    double *b5=create_1d_double_array(nfa, "a5");            //bending moduli quad term
    double *b6=create_1d_double_array(nfa, "a6");            //bending moduli quart term
    
    
    //open main output file
    ofstream outFile2;
    string filename2;
    filename2="./results/fA_test.dat";
    outFile2.open(filename2.c_str());
    
    //radius ouput
    ofstream radiout;
    radiout.open("./results/main_radius.dat");
    
    
    //output quadratic bending modulus
    ofstream bendingout;
    bendingout.open("./results/bending_mod.dat");
    
    //output average diameter
    ofstream diamout;
    diamout.open("./results/diameter_output.dat");
    
    
    int counter=0;
    double volume;
    double OP;
    double fE_hom;
    
    int ds_increment = 0.4*(double)Ns[2]/((double)nfa-1.0);
    
    for (int dds=0; f[0]<0.7 ;dds+= ds_increment){
        counter+=1;
        //Set parameters s
        updateparameters(f,Ns,dds);
        mu[2] = mu_vec[counter-1];                        //don't want to calc mu again
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        
        double pin_location=10.8-D[0]-(E[0]*f[0])-(F[0]*f[0]*f[0]);
        int pin = pin_location/dr;
        
        volume=vol(dr);                                 //calculate volume
        OP = calcOP(phi,dr,volume);                     //calculate order parameter
        //Set radius vector
        set_radius(r_0vector,nradii,f[0]);
        
        r_0=r_0vector[0];                                        //reset radius
        double avgradius=0.0;                                  //reset avgradius
        double avgmiddle=0.0;
        double avgdiameter=0.0;
        
        for (int radius=0;radius<nradii;radius++){
            volume=vol(dr);
            
            //initialize omega field
            omega(w);
            
            //calculate free energy minus homogeneneous free energy
            dFE[radius]=FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,volume,f,pin,1);
            OP = calcOP(phi,dr,volume);                    //calculate order parameter
            diameter[radius] = calc_excess(phi,dr,volume); //calculate copolymer excess
            avgdiameter+=diameter[radius];
            
            int imax=mmbcentre(phi);                       //membrane center (max phib)
            int ihalf=mmb_half(phi,imax,pin);                      //membrane middle (1/2 phib = phiA)
            
            avgradius+=(double)imax*dr;                    //avg membrane center
            avgmiddle+=(double)ihalf*dr;
            
            Rad[radius]=r_0;                               //set radius vectors
            Rad2[radius]=r_0;
            
            //output free energy data
            outFile2 <<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            cout<<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            
            //output concentration profile
            outputphi(phi,dr);
            
            //set new radius
            r_0=r_0vector[radius+1];
            
        }
        //output average radius
        avgradius/= (double)nradii;
        avgmiddle/= (double)nradii;
        radiout<<f[0]<<" "<<avgradius<<" "<<avgmiddle<<endl;
        
        //calculate average diameter
        avgdiameter/=nradii;
        diamout<<f[0]<<" "<<avgdiameter<<endl;
        
        
        //build curvature and curvature squared vectors
        for (int radius=0;radius<nradii;radius++){
            //Rad[radius]+=avgradius;   //membrane should be centered at 6, but just in case
            //Rad2[radius]+=avgmiddle;
            
            //Rad[radius]+=(A[0]+f[0]*B[0]+f[0]*f[0]*C[0]);
            //Rad2[radius]+=(D[0]+f[0]*E[0]+f[0]*f[0]*F[0]);
            
            Rad[radius]+=6.0;
            Rad2[radius]+=6.0;

            
            Curv[radius] =(4.3/Rad[radius]);
            Curvsq[radius] = pow(Curv[radius],2.0);
            
            Curv2[radius] =(4.3/Rad2[radius]);
            Curvsq2[radius] = pow(Curv2[radius],2.0);
        }
        
        //output files of free energy as a function of radius
        outputfE_FA(f[0],Curv,dFE,nradii);
        
        //quadratic curve fit with max hydrophobic
        curvefit(Curv,dFE,nradii,counter,a1,a2,a3);
        
        //quartic curve fit with max hydrophobic
        curvefit(Curvsq,dFE,nradii,counter,a4,a5,a6);
        
        //quadratic curve fit with max hydrophobic
        curvefit(Curv2,dFE,nradii,counter,b1,b2,b3);
        
        //quartic curve fit with max hydrophobic
        curvefit(Curvsq2,dFE,nradii,counter,b4,b5,b6);
        
        
        if (Coord==2){
            bendingout<<f[0]<<" "<<a5[counter-1]<<" "<<2.0*a5[counter-1]*200.0/(sqrt(chi[0])*avgdiameter*avgdiameter)<<endl;
        }
        else if (Coord==3){
            bendingout<<f[0]<<" "<<a5[counter-1]<<" "<<a5[counter-1]*200.0/(sqrt(chi[0])*avgdiameter*avgdiameter)<<endl;
        }
        
    }
    outFile2.close();
    radiout.close();
    bendingout.close();
    diamout.close();
    
    //output all curvefit results
    outputkappa(a1,a2,a3,a4,a5,a6,nfa,chi,0);
    outputkappa(b1,b2,b3,b4,b5,b6,nfa,chi,1);
    
    
}