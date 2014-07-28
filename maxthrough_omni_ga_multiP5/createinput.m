% createinput.m is to create the input file(input.in) to be read by the
% nsga2.c code. 
function [] = createinput(M,N_f)
popsize = 1000;
gen = 2000;
nobj = 2;
pcross_real = 0.9;
pmut_real = 0.05;
eta_c = 15;
eta_m = 70; 
nbin = 0; % #binary variable
choice = 0; %choice for plotting
help_string = ''; %This string contains the order of values at teh end of input file
fid = fopen('input.in','w');
fprintf(fid,'%d\n',popsize);
help_string = sprintf('%s\n%s',help_string,'popsize');
fprintf(fid,'%d\n',gen);
help_string = sprintf('%s\n%s',help_string,'gen');
fprintf(fid,'%d\n',nobj);
help_string = sprintf('%s\n%s',help_string,'#obj');
fprintf(fid,'%d\n',M*N_f+M+N_f);
help_string = sprintf('%s\n%s',help_string,'#cons');
fprintf(fid,'%d\n',M*N_f);
help_string = sprintf('%s\n%s',help_string,'#realvar');
help_string = sprintf('%s\n%s',help_string,'#realvar lower-upper limit');
for i = 1:M*N_f
    fprintf(fid,'%f\n',-0.6);
    fprintf(fid,'%f\n',0.6);
end
help_string = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s',help_string,'real_crossover','real_mutation','eta_c','eta_m','#bin_variable','choice_for_plotting');
fprintf(fid,'%f\n%f\n',pcross_real,pmut_real);
fprintf(fid,'%d\n%d\n',eta_c,eta_m);
fprintf(fid,'%d\n%d\n',nbin,choice);
fprintf(fid,'%s',help_string);
end