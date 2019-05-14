This is a repository for creating Data Assimilation codes for the Stechmann & Hottovy 2017 (GRL) stochastic model. 

Created By: Scott Hottovy, 14 MAY 2019, hottovy@usna.edu

List of Folders: 
                  GlobalFuns                -Contains m-files for global functions such as creating parabolic cylinder functions and the model code. 
                  
                  
List of Files: 
                  Hovmoller_uPrecip_DAMethod.m      -A MATLAB script to run a test code the Stechmann & Hottovy 2017 model in 
                                                     a format conducive to the language of the Majda and Harlim 2012 textbook. That
                                                     is a model given by a linear function F and covariance matrix R. 
                                                     
                  SH2017_LinearModel_DA_F.m         -A MATLAB function to give the model F_k for each wavenumber k (integer). 
                  Covariance_SH2017.m               -A MATLAB function to give the model covariance R_k for each wavenumber k (integer). 
                  RealizationFourierQMode_wF.m      -A MATLAB funciton to run the Stechmann & Hottovy 2017 model, used in the GRL paper. 
                                                     This method is much faster than the DA method. It is used in the script to generate 
                                                     the statistical steady state as the initial condition for the Hovmoller test. 
