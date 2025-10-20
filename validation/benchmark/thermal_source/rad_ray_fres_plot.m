%% Read in data
clear
filedir='./';
SOS_filename='vSOS_ray_fres';
rt3_output_filename='test_mw_fres.dat';

nheader1=8;
nquad_rt3=32;
fid = fopen(strcat(filedir,rt3_output_filename),'r');
header=textscan(fid,'%s',nheader1,'delimiter', '\n');
readdata=textscan(fid,'%f %f %f %f %f',1,'delimiter', '\n');
rt3_irrad_up_toa=readdata{1,4};
readdata=textscan(fid,'%f %f %f %f %f',1,'delimiter', '\n');
rt3_irrad_down_toa=readdata{1,4};
readdata=textscan(fid,'%f %f %f %f %f',nquad_rt3,'delimiter', '\n');
rt3mu=flipud(abs(readdata{1,3}));
rt3I=flipud(readdata{1,4});
rt3Q=flipud(readdata{1,5});
fclose(fid);


nheader1=16;

fid = fopen(strcat(filedir,SOS_filename),'r');
header=textscan(fid,'%s',nheader1,'delimiter', '\n');
SOS_irrad_up_toa=textscan(fid,'UPWELLING IRRADIANCE= %f',1,'delimiter', '\n');
header=textscan(fid,'%s',1,'delimiter', '\n');
SOSrad=textscan(fid,'%f %f %f %f %f %f',nquad_rt3,'delimiter', '\n');
fclose(fid);



%Plank_Epsilon=0.5*4.8147945192665320;
sosmu=cosd(SOSrad{1,1});
sosphi=SOSrad{1,2};
sosI=SOSrad{1,3}; %/Plank_Epsilon;
sosQ=SOSrad{1,4}; %/Plank_Epsilon;
delta_sos=sosQ./sosI;


%%
rt3I=interp1(rt3mu,rt3I,sosmu);
rt3Q=interp1(rt3mu,rt3Q,sosmu);
delta_rt3=rt3Q./rt3I;

%%
figure(31)
%subplot(2,2,1)
yyaxis left;
plot(sosmu,rt3I,'-ok',sosmu,sosI,'-+m','MarkerSize',10,'linewidth',2)
%legend('I RTSOS','I rt3')
xlabel('cos(\theta)','FontSize',18)
ylabel('Radiance','FontSize',18)
set(gca,'FontSize',18)
%%
figure(32)
%subplot(2,2,1)
yyaxis right;
plot(sosmu,rt3Q./rt3I,':sk',sosmu,sosQ./sosI,':xm','MarkerSize',10,'linewidth',2)
%legend('Q/I RTSOS','Q/I rt3')
xlabel('cos(\theta)','FontSize',18)
ylabel('Q/I','FontSize',18)
xlim([0.05 1])
legend('I RTSOS','I rt3','Q/I RTSOS','Q/I rt3','location','south')
set(gca,'FontSize',18)
exportgraphics(figure(31), 'rt3_comparison.pdf');
%%
figure(33)
%subplot(2,2,1)
plot(sosmu,100*(sosI-rt3I)./rt3I,'-+b',sosmu,100*(delta_sos-delta_rt3),'-xm','MarkerSize',10)
legend('I % diff','100 Q/I diff')
xlabel('cos(\theta)','FontSize',18)
ylabel('diff','FontSize',18)
set(gca,'FontSize',18)
