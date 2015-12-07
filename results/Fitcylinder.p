f1(x) = (k/2.0)*(4.3/(x+5.84))**2           # define the function to be fit
k = 0.01; 
fit f1(x) 'fe(r)_OP_-1.dat' using 2:4 via k
