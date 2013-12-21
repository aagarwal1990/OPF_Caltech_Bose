clear all
close all
clc

ID = fopen('data_sum.dat','wt');
case_list = [9];
% case_list = [9,14,30,39,57];

fprintf(ID, '%10s %10s %12s %10s %15s \n','case_num', 'rel_type', 'obj_val', 'time', 'max_eig_ratio');
for ii = 1:length(case_list)
    case_num = strcat('case',int2str(case_list(ii)));
    
    solveOPF_chordal_relaxation_w
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','CH');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');
    
    solveOPF_SDP
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','SDP');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');
    
    solveOPF_SOCP_less_variables
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','SOCP');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');


end

fclose(ID);