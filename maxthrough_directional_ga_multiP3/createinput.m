function [] = createinput(M,N_f)
popsize = 1000;
gen = 2000;
nobj = 2;
pcross_real = 0.9;
pmut_real = 0.05;
eta_c = 15;
eta_m = 70;
nbin = 0;
choice = 0;

fid = fopen('input.in','w');
fprintf(fid,'%d\n',popsize);
fprintf(fid,'%d\n',gen);
fprintf(fid,'%d\n',nobj);
fprintf(fid,'%d\n',M*N_f+M+N_f);
fprintf(fid,'%d\n',M*N_f+M);
for i = 1:M*N_f
    fprintf(fid,'%f\n',-0.6);
    fprintf(fid,'%f\n',0.6);
end

for i = 1:M
    fprintf(fid,'%f\n',0);
    fprintf(fid,'%f\n',360);
end
fprintf(fid,'%f\n%f\n',pcross_real,pmut_real);
fprintf(fid,'%d\n%d\n',eta_c,eta_m);
fprintf(fid,'%d\n%d\n',nbin,choice);
end