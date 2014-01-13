clear all
close all
clc

ID = fopen('data_sum.dat','wt');
ID_2 = fopen('volt.dat','wt');
case_list = [9,14,30,39,57];

fprintf(ID, '%10s %10s %12s %10s %15s \n','case_num', 'rel_type', 'obj_val', 'time', 'max_eig_ratio');
fprintf(ID_2, '%5s %5s\n','case_num', 'rel_type');

for ii = 1:length(case_list)
    fprintf('\n');
    case_num = strcat('case',int2str(case_list(ii)));
    
    solveOPF_chordal_relaxation_w
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','CH');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');
    
    fprintf('\n');
    fprintf(ID_2,'%s ',case_num);
    fprintf(ID_2,'%10s ','CH');
    fprintf(ID_2,'\n');
    fprintf(ID_2,'Node\n');
    for jj = 1:n
     fprintf(ID_2,'%5f + %5fi\n', real(volt(jj)), imag(volt(jj)));
    end
    fprintf('\n');
    
    solveOPF_SDP
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','SDP');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');
    
    fprintf('\n');
    fprintf(ID_2,'%s ',case_num);
    fprintf(ID_2,'%8s ','SDP');
    fprintf(ID_2,'\n');
    fprintf(ID_2,'Node\n');
    for jj = 1:n
     fprintf(ID_2,'%5f + %5fi\n', real(volt(jj)), imag(volt(jj)));
    end
    fprintf('\n');
     
    solveOPF_SOCP_sparse
    fprintf(ID,'%10s ',case_num);
    fprintf(ID,'%10s ','SOCP');
    fprintf(ID,'%12f ',objective_value);
    fprintf(ID,'%10f ',elapsed_time);
    fprintf(ID,'%10e ',maxEigRatio);
    fprintf(ID,' \n');
    
    fprintf(ID_2,'\n');
    fprintf(ID_2,'%s ',case_num);
    fprintf(ID_2,'%10s ','SOCP');
    fprintf(ID_2,'\n');
    fprintf(ID_2,'Node\n');
    for jj = 1:n
     fprintf(ID_2,'%5f + %5fi\n', real(volt(jj)), imag(volt(jj)));
    end
    fprintf('\n');


end

fclose(ID);