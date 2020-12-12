Goal-oriented Coagulation Management (GCM) Algorithm
====================

Personalized Modulation of Coagulation Factors Using a Thrombin Dynamics Model to Treat Trauma-induced Coagulopathy

---------------------

This is the repository for the algorithm developed and proposed in the above mentioned paper. 

## System Requirements
# Hardware requirement 
This algorithm only requires a standard computer with enough RAM to support in-memory operations. 

#Software requirements
Running the GCM algorithm code requires a recent version of MATLAB. R2018a and newer is recommended.
The algorithm is tested on the following MATLAB versions: R2018a, R2018b, R2019b, R2020a

Install MATLAB according to MathWorks installation guide (https://www.mathworks.com/help/install/install-products.html)


## Initializing and User Instructions: 
1. Open MATLAB.
2. Boundedline function used to mark the normal region can be retrieved from the following: 
Kelly Kearney (2020). boundedline.m (https://github.com/kakearney/boundedline-pkg), GitHub. Retrieved December 12, 2020.
3. Extract boundedline-pkg to MATLAB directory. 
4. Add DetermineDelayTherapy.m function file to MATLAB directory. 
5. Add GCM_Algorithm_Parameter_Data.mat file to MATLAB directory. 
6. Open and run the algorithm file, GCM_Algorithm_for_TIC.m.
7. When prompted to input factor concentration values, use the trauma patient coagulation factor concentration measurements provided as a sample data set in TraumaSampleFactorConcentration_ForDemo.xlsx.
8. The algorithm will output the recommended concentrations, as well as CAT with current factors and CAT with normal factors. 
9. Estimated run time is less than a minute.

## Important Dates
15 December 2020 – Initial Submission <br /> 