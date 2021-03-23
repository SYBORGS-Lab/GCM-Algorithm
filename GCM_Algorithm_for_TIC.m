% Goaal-oriented Coagulation Management (GCM) Algorithm 
% Personalized Modulation of Coagulation Factors USing a Thrombin Dynamics
% Model to Treat Trauma-Induced Coagulopathy 
% Damon E. Ghetmiri, Mitchell J. Cohen, Amor A. Menezes 
% ©Copyright 2021 University of Florida Research Foundation, Inc. All Commercial Rights Reserved.

%% Trainining the Dynamic thrombin generation model parameter estimation, using matching pursuit greedy algorithm,  

% Load data that includes trained thrombin prediction model paramters. 

load('GCM_Algorithm_Parameter_Data.mat');

%% Forming the normal region: Identifying what is the desired normal range for thrombin generation trajector (CAT) 
% Mean of Parameters / Max and Min of Fit to Data

sys_est_Min_Nor=tf(Min_ModelParameters_Nor(4),[1 Min_ModelParameters_Nor(3) Min_ModelParameters_Nor(2) Min_ModelParameters_Nor(1)],'InputDelay',Min_ModelParameters_Nor(5));
sys_est_Max_Nor=tf(Max_ModelParameters_Nor(4),[1 Max_ModelParameters_Nor(3) Max_ModelParameters_Nor(2) Max_ModelParameters_Nor(1)],'InputDelay',Max_ModelParameters_Nor(5));
sys_est_Mean_Nor=tf(Mean_ModelParameters_Nor(4),[1 Mean_ModelParameters_Nor(3) Mean_ModelParameters_Nor(2) Mean_ModelParameters_Nor(1)],'InputDelay',Mean_ModelParameters_Nor(5));

T2 = linspace(0,42,124)';
Y_est_Min_Nor = 5*impulse(sys_est_Min_Nor,T2);
Y_est_Max_Nor = 5*impulse(sys_est_Max_Nor,T2);
Y_est_Mean_Nor = 5*impulse(sys_est_Mean_Nor,T2);

%Identify normal range 
[Range_peak_Min,i_m]=max(Y_est_Min_Nor);
Range_T_peak_Min= T2(i_m) ;
Range_Area_Min= trapz(T2,Y_est_Min_Nor ) ;
Range_Delay_Min=DetermineDelayTherapy(T2,Y_est_Min_Nor);

[Range_peak_Max,i_m]=max(Y_est_Max_Nor);
Range_T_peak_Max= T2(i_m) ;
Range_Area_Max= trapz(T2,Y_est_Max_Nor ) ;
Range_Delay_Max=DetermineDelayTherapy(T2,Y_est_Max_Nor);

[Range_peak_Mean,i_m]=max(Y_est_Mean_Nor);
Range_T_peak_Mean= T2(i_m) ;
Range_Area_Mean= trapz(T2,Y_est_Mean_Nor ) ;
Range_Delay_Mean=DetermineDelayTherapy(T2,Y_est_Mean_Nor);

Range_Peak=[Range_peak_Min Range_peak_Mean Range_peak_Max];
Range_T_Peak=[Range_T_peak_Max Range_T_peak_Mean Range_T_peak_Min];
Range_Area=[Range_Area_Min Range_Area_Mean Range_Area_Max];
Range_Delay=[Range_Delay_Max Range_Delay_Mean Range_Delay_Min];


%% Trauma Data: Import new trauma patient data that has not been used in the training 
%Order of factor concentrations: Factors II, V, VII, IX, X, VIII, ATIII, PC

factor_concen_Tra_sample_initial(1,1)=input('What is factor II concentration? ') ;
factor_concen_Tra_sample_initial(1,2)=input('What is factor V concentration? ') ;
factor_concen_Tra_sample_initial(1,3)=input('What is factor VII concentration? ') ;
factor_concen_Tra_sample_initial(1,4)=input('What is factor XI concentration? ') ;
factor_concen_Tra_sample_initial(1,5)=input('What is factor X concentration? ') ;
factor_concen_Tra_sample_initial(1,6)=input('What is factor VIII concentration? ') ;
factor_concen_Tra_sample_initial(1,7)=input('What is factor ATIII concentration? ') ;
factor_concen_Tra_sample_initial(1,8)=input('What is Protein C concentration? ') ;
 
%% CAT Variation with factor recommendation adjustments Sample by sample 
T3 = linspace(0,42,42001)';
FactorConcentration_History_TraSample=[];
Factor_tag={'II', 'V', 'VII' , 'IX', 'X', 'VIII', 'ATIII', 'PC'};
FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; factor_concen_Tra_sample_initial];

%STEP 1: Correct factor V to normal range 
if FactorConcentration_History_TraSample(1,6)>400
    FactorConcentration_History_TraSample(1,6)=FactorConcentration_History_TraSample(1,6)/2;
end
j=2; %factor V
recommend_factor_update=min(max(60, FactorConcentration_History_TraSample(end,j)),140);
FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

%STEP 2: Correct factor V to normal range 
j=3; %factor VII
recommend_factor_update=min(max(60, FactorConcentration_History_TraSample(end,j)),140);
FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

%STEP 3: adjusting factor II to correct peak 

j=1; %Factor II
Factor_Tra_changes_treat={};

factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 

%Estimate the CAT model using the update coagulation factor set 
Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));

Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);

Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];
% Forming the mappifng function from factor concentration to graph parameter

    for i=1:10 
        
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter= ((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
        
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
        
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];
        

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,2),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from peak to factor II
recommend_factor_update=min(max(60, polyval(factor_param_fit, Range_Peak(2))),140);

FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

%STEP 4: adding factor X to correct peak and peak time 

j=5;
Factor_Tra_changes_treat={};

factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 

Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));

Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);

Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    for i=1:10
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter= ((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
        
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
        
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,2),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from peak to factor X
recommend_factor_update=max(min(140, polyval(factor_param_fit, Range_Peak(3))),60);

FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;


%STEP 5: adjusting IX for time delay 
j=4;  % factor IX
Factor_Tra_changes_treat={};

factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 

Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));

Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);

Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    for i=1:10
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter=((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
              
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
        
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
        
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,5),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from time delay to PC
recommend_factor_update=max(min(140, polyval(factor_param_fit, Range_Delay(3))),60);

FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;


%STEP 6: adjusting VIII  based on peak time 

j=6;
Factor_Tra_changes_treat={};

factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 

Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));

Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);

Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    for i=1:10
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter= ((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        
        if Trau_1_coeff(4)<0
            change_parameter=-change_parameter;
            FactorSelected=FactorSelected+change_parameter ;
            factor_concen_Tra_sample(j)=FactorSelected;
            Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        end
        
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
        
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
        
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,3),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from peak time to factor VIII
recommend_factor_update=max(min(140, polyval(factor_param_fit, mean(Range_T_Peak(1:2)))),60);

FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

j=1;
Factor_Tra_changes_treat={};
factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 
Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    for i=1:10
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter= ((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);      
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,2),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from peak  to factor II
recommend_factor_update=max(min(140, polyval(factor_param_fit, Range_Peak(2))),60);
FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

% STEP 7 (IF NECESSARY)

Y_est_Min_Nor_Comparison = 5*impulse(sys_est_Min_Nor,T3);
Y_est_Max_Nor_Comparison = 5*impulse(sys_est_Max_Nor,T3);
Y_est_Mean_Nor_Comparison = 5*impulse(sys_est_Mean_Nor,T3);

%Identify normal range with T3 for comparison to understand where is the
%location of CAT 
Range_peak_Min_Comparison=max(Y_est_Min_Nor_Comparison);
Range_peak_Max_Comparison=max(Y_est_Max_Nor_Comparison);
Range_peak_Mean_Comparison=max(Y_est_Mean_Nor_Comparison);
Range_Peak=[Range_peak_Min_Comparison Range_peak_Mean_Comparison Range_peak_Max_Comparison];

%Current CAT after Step 5 adjustment 
Trau_1_coeff=FactorConcentration_History_TraSample(end,:)*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
Y_est_current_Tra= 5*impulse(sys_est_change_Tra1,T3);
Peak_current_Tra=max(Y_est_current_Tra);

% Identify if peak is in the range 
if Peak_current_Tra>Range_Peak(3)
    Current_CAT_Peak_Indicator=1;
elseif Range_Peak(1)>Peak_current_Tra
    Current_CAT_Peak_Indicator=-1;
else
    Current_CAT_Peak_Indicator=0;
end
% Identify if tail (between 10 and 20 minutes) is in the range 
if mean(Y_est_current_Tra(12500:22500)-Y_est_Max_Nor_Comparison(12500:22500))>0
    Current_CAT_Tail_Indicator=1;
elseif mean(Y_est_current_Tra(8000:18000)-Y_est_Min_Nor_Comparison(8000:18000))<0
    Current_CAT_Tail_Indicator=-1;
else
    Current_CAT_Tail_Indicator=0;
end


%CHoose what factor to change 
if Current_CAT_Peak_Indicator==1 && Current_CAT_Tail_Indicator~=-1
    j=8;
    Target_Peak_Value=mean(Range_Peak(2:3));
    NeedForStep6=1;
elseif Current_CAT_Peak_Indicator==1 && Current_CAT_Tail_Indicator==-1
    if FactorConcentration_History_TraSample(end,5)~=140
        j=5;
    else
        j=7;
    end
    Target_Peak_Value=mean(Range_Peak(2:3));
    NeedForStep6=1;
elseif Current_CAT_Peak_Indicator==-1 && Current_CAT_Tail_Indicator==-1
    j=7;
    Target_Peak_Value=mean(Range_Peak(2));
    NeedForStep6=1;
elseif Current_CAT_Peak_Indicator==0 && Current_CAT_Tail_Indicator==+1
    j=7;
    Target_Peak_Value=Range_Peak(2);
    NeedForStep6=1;
elseif Current_CAT_Peak_Indicator==0 && Current_CAT_Tail_Indicator==-1
    j=7;
    Target_Peak_Value=Range_Peak(2);
    NeedForStep6=1;
elseif Current_CAT_Peak_Indicator==-1 && Current_CAT_Tail_Indicator==0
    FactorConcentration_History_TraSample(end,6)=(FactorConcentration_History_TraSample(5,6)+FactorConcentration_History_TraSample(4,6))/9*5;
    NeedForStep6=0;
else
    NeedForStep6=0;
end


if NeedForStep6==1

Factor_Tra_changes_treat={};

factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
Delta_parameter=factor_concen_Tra_sample(j)/10 ; 

Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));

Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);
[Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
T_est_change_Tra_peak= T3(i_m) ;
AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);

Factor_Tra_changes_treat{1,j}=[factor_concen_Tra_sample(j) Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    for i=1:10
        factor_concen_Tra_sample=FactorConcentration_History_TraSample(end,:);
        FactorSelected=factor_concen_Tra_sample(j);
        change_parameter= ((-1)^i)*Delta_parameter*(ceil(i/2)) ;
        FactorSelected=FactorSelected+change_parameter ;
        factor_concen_Tra_sample(j)=FactorSelected;
        Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        
        if Trau_1_coeff(4)<0
            change_parameter=-change_parameter;
            FactorSelected=FactorSelected+change_parameter ;
            factor_concen_Tra_sample(j)=FactorSelected;
            Trau_1_coeff=factor_concen_Tra_sample*factor_coeff_All_Train+const_All_Train;
        end
        
        sys_est_change_Tra1=tf(Trau_1_coeff(4),[1 Trau_1_coeff(3) Trau_1_coeff(2) Trau_1_coeff(1)],'InputDelay',Trau_1_coeff(5));
        Y_est_change_Tra= 5*impulse(sys_est_change_Tra1,T3);       
        [Y_est_change_Tra_peak,i_m]=max(Y_est_change_Tra);
        T_est_change_Tra_peak= T3(i_m) ;
        AreaUnderCurve= trapz(T3,Y_est_change_Tra ) ;
        DelayCurve=DetermineDelayTherapy(T3,Y_est_change_Tra);     
        Factor_Tra_changes_treat{i+1,j}=[FactorSelected Y_est_change_Tra_peak T_est_change_Tra_peak AreaUnderCurve DelayCurve];

    end

Factor_Tra_changes_treat_Mat=cell2mat(Factor_Tra_changes_treat);
factor_param_fit=polyfit(Factor_Tra_changes_treat_Mat(:,2),Factor_Tra_changes_treat_Mat(:,1),2); %a second order polynomial from peak to required factor 
recommend_factor_update=max(min(140, polyval(factor_param_fit, Target_Peak_Value)),60);
FactorConcentration_History_TraSample=[FactorConcentration_History_TraSample; FactorConcentration_History_TraSample(end,:)];
FactorConcentration_History_TraSample(end,j)=recommend_factor_update;

end

Recommended_Coagulation_Factor_Set_S=FactorConcentration_History_TraSample(end,:);

%Display the recommended coagulation concentration set: 
fprintf('============================================ \n')
for factS=1:8
    fprintf('Recommended factor %s concnetration is: %4.1f \n',Factor_tag{factS},Recommended_Coagulation_Factor_Set_S(factS))
end

%Figure Visualization 
figure(1)
clf;
[hl, hp]=boundedline(T2, Y_est_Mean_Nor, [Y_est_Mean_Nor-Y_est_Min_Nor Y_est_Max_Nor-Y_est_Mean_Nor ], 'alpha');
outlinebounds(hl, hp)
box on
xlim([0 42])
ax = gca;
ax.FontSize = 20; 
xlabel('Time [min]')
ylabel('CAT')
title('Recommendation')

%Current CAT
Trau_coeff=factor_concen_Tra_sample_initial*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_coeff(4),[1 Trau_coeff(3) Trau_coeff(2) Trau_coeff(1)],'InputDelay',Trau_coeff(5));
Y_est_Tra= 5*impulse(sys_est_change_Tra1,T3);
figure(1)
hold on
p_c=plot(T3,Y_est_Tra,'r','LineWidth',3);

%Recommended CAT
Trau_coeff=FactorConcentration_History_TraSample(end,:)*factor_coeff_All_Train+const_All_Train;
sys_est_change_Tra1=tf(Trau_coeff(4),[1 Trau_coeff(3) Trau_coeff(2) Trau_coeff(1)],'InputDelay',Trau_coeff(5));
Y_est_Tra= 5*impulse(sys_est_change_Tra1,T3);
figure(1)
hold on
p_r=plot(T3,Y_est_Tra,'k','LineWidth',3);

legend([hp p_c p_r],{'Normal Region','Current CAT','Recommended CAT'})

