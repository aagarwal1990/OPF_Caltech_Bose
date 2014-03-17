clear all
close all
clc

case_list = [9,14,30,39,57, 118];

for ii = 1:length(case_list)
    case_num = strcat('case',int2str(case_list(ii)));
    
    solveOPF_SDP_nonRank1
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
    



end

fclose(ID);