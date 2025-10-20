clear;

% modify the filedir to your directory
filedir='./';
kokfiledir='./Kok_data/SCIATRAN_BENCHMARK_RESULTS/';
szen=60;
ntangle=90;
refflag=1;
dataread=1;

nheader=19;
nheader2=5;
SOS_filename='vSOS_current';

if(dataread)
  fid = fopen(strcat(filedir,SOS_filename),'r');
  header=textscan(fid,'%s',nheader,'delimiter', '\n');
  SOSradref=textscan(fid,'%f %f %f %f %f %f',540,'delimiter', '\n');
  header2=textscan(fid,'%s',nheader2,'delimiter', '\n');
  SOSradtran=textscan(fid,'%f %f %f %f %f %f',540,'delimiter', '\n');
  fclose(fid);

%   fid = fopen(strcat(filedir,SOS_filename),'r');
%   header=textscan(fid,'%s',20,'delimiter', '\n');
%   SOSradref=textscan(fid,'%f %f %f %f %f %f',540,'delimiter', '\n');
%   header2=textscan(fid,'%s',4,'delimiter', '\n');
%   SOSradtran=textscan(fid,'%f %f %f %f %f %f',540,'delimiter', '\n');
%   fclose(fid);

  
  Kok_filename='aerosol_refl_N_240.dat';
  fid = fopen(strcat(kokfiledir,Kok_filename),'r');
  format = repmat('%f ', 1, 13);
  Kokradref=textscan(fid,format,ntangle,'delimiter', '\n');
  fclose(fid);

  Kok_filename='aerosol_trans_N_240.txt';
  fid = fopen(strcat(kokfiledir,Kok_filename),'r');
  format = repmat('%f ', 1, 13);
  Kokradtran=textscan(fid,format,ntangle,'delimiter', '\n');
  
  fclose(fid);
  
  mc_filename='MC_AEROSOL_TOA';
  fid = fopen(mc_filename,'r');
  header5=textscan(fid,'%s',17,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  mcradAref=textscan(fid,format,ntangle*6,'delimiter', '\n');
  fclose(fid);

  mc_filename='MC_AEROSOL_BOA';
  fid = fopen(mc_filename,'r');
  header6=textscan(fid,'%s',17,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  mcradAtran=textscan(fid,format,ntangle*6,'delimiter', '\n');
  fclose(fid);
  
  save('Kok_val_data');
else
  load 'Kok_val_data'
end



sosmuref=SOSradref{1,1};
sosphiref=SOSradref{1,2};
sosIref=SOSradref{1,3}./cosd(szen);
sosQref=-SOSradref{1,4}./SOSradref{1,3};
sosUref=SOSradref{1,5}./SOSradref{1,3};
sosVref=-SOSradref{1,6}./SOSradref{1,3};

sosmuref=reshape(sosmuref,180,3);
sosphiref=reshape(sosphiref,180,3);
sosIref=reshape(sosIref,180,3);
sosQref=reshape(sosQref,180,3);
sosUref=reshape(sosUref,180,3);
sosVref=reshape(sosVref,180,3);

sosmutran=SOSradtran{1,1};
sosphitran=SOSradtran{1,2};
sosItran=SOSradtran{1,3}./cosd(szen);
sosQtran=-SOSradtran{1,4}./SOSradtran{1,3};
sosUtran=SOSradtran{1,5}./SOSradtran{1,3};
sosVtran=-SOSradtran{1,6}./SOSradtran{1,3};

sosmutran=flipud(reshape(sosmutran,180,3));
sosphitran=flipud(reshape(sosphitran,180,3));
sosItran=flipud(reshape(sosItran,180,3));
sosQtran=flipud(reshape(sosQtran,180,3));
sosUtran=flipud(reshape(sosUtran,180,3));
sosVtran=flipud(reshape(sosVtran,180,3));

Kokangref=Kokradref{1,1};
KokIref0=Kokradref{1,2};
KokQref0=Kokradref{1,3}./Kokradref{1,2};
KokUref0=Kokradref{1,4}./Kokradref{1,2};
KokVref0=Kokradref{1,5}./Kokradref{1,2};

KokIref90=Kokradref{1,6};
KokQref90=Kokradref{1,7}./Kokradref{1,6};
KokUref90=Kokradref{1,8}./Kokradref{1,6};
KokVref90=Kokradref{1,9}./Kokradref{1,6};

KokIref180=Kokradref{1,10};
KokQref180=Kokradref{1,11}./Kokradref{1,10};
KokUref180=Kokradref{1,12}./Kokradref{1,10};
KokVref180=Kokradref{1,13}./Kokradref{1,10};

Kokangtran=Kokradtran{1,1};
KokItran0=Kokradtran{1,2};
KokQtran0=Kokradtran{1,3}./Kokradtran{1,2};
KokUtran0=Kokradtran{1,4}./Kokradtran{1,2};
KokVtran0=Kokradtran{1,5}./Kokradtran{1,2};

KokItran90=Kokradtran{1,6};
KokQtran90=Kokradtran{1,7}./Kokradtran{1,6};
KokUtran90=Kokradtran{1,8}./Kokradtran{1,6};
KokVtran90=Kokradtran{1,9}./Kokradtran{1,6};

KokItran180=Kokradtran{1,10};
KokQtran180=Kokradtran{1,11}./Kokradtran{1,10};
KokUtran180=Kokradtran{1,12}./Kokradtran{1,10};
KokVtran180=Kokradtran{1,13}./Kokradtran{1,10};

mcangref=mcradAref{1,1};
mcphiref=mcradAref{1,2};

mcIAref=mcradAref{1,3}*pi;
mcQAref=mcradAref{1,7}./mcradAref{1,3};
mcUAref=mcradAref{1,11}./mcradAref{1,3};
mcVAref=mcradAref{1,15}./mcradAref{1,3};
mcangref=reshape(mcangref,180,3);
mcIAref=reshape(mcIAref,180,3);
mcQAref=-reshape(mcQAref,180,3);
mcUAref=-reshape(mcUAref,180,3);
mcVAref=reshape(mcVAref,180,3);

mcangtran=mcradAtran{1,1};
mcphitran=mcradAtran{1,2};
mcIAtran=mcradAtran{1,3}*pi;
mcQAtran=-mcradAtran{1,7}./mcradAtran{1,3};
mcUAtran=-mcradAtran{1,11}./mcradAtran{1,3};
mcVAtran=mcradAtran{1,15}./mcradAtran{1,3};

mcangtran=flipud(reshape(mcangtran,180,3));
mcIAtran=flipud(reshape(mcIAtran,180,3));
mcQAtran=flipud(reshape(mcQAtran,180,3));
mcUAtran=flipud(reshape(mcUAtran,180,3));
mcVAtran=flipud(reshape(mcVAtran,180,3));

figure(1)
semilogy(Kokangref,sosIref(1:90,1),'+b',Kokangref,KokIref0,'ob',...
    Kokangref,sosIref(1:90,2),'+k',Kokangref,KokIref90,'ok',...
    Kokangref,sosIref(1:90,3),'+m',Kokangref,KokIref180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Radiance at TOA')
legend('SOS I \phi=0','Kok benchmark \phi=0',...
    'SOS I \phi=90','Kok benchmark \phi=90',...
    'SOS I \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(2)
plot(Kokangref,sosQref(1:90,1),'+b',Kokangref,KokQref0,'ob',...
    Kokangref,sosQref(1:90,2),'+k',Kokangref,KokQref90,'ok',...
    Kokangref,sosQref(1:90,3),'+m',Kokangref,KokQref180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Q at TOA')
legend('SOS Q \phi=0','Kok benchmark \phi=0',...
    'SOS Q \phi=90','Kok benchmark \phi=90',...
    'SOS Q \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(3)
plot(Kokangref,sosUref(1:90,1),'+b',Kokangref,KokUref0,'ob',...
    Kokangref,sosUref(1:90,2),'+k',Kokangref,KokUref90,'ok',...
    Kokangref,sosUref(1:90,3),'+m',Kokangref,KokUref180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('U at TOA')
legend('SOS U \phi=0','Kok benchmark \phi=0',...
    'SOS U \phi=90','Kok benchmark \phi=90',...
    'SOS U \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(4)
plot(Kokangref,sosVref(1:90,1),'+b',Kokangref,KokVref0,'ob',...
    Kokangref,sosVref(1:90,2),'+k',Kokangref,KokVref90,'ok',...
    Kokangref,sosVref(1:90,3),'+m',Kokangref,KokVref180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('V at TOA')
legend('SOS V \phi=0','Kok benchmark \phi=0',...
    'SOS V \phi=90','Kok benchmark \phi=90',...
    'SOS V \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(5)

sosIreferr0=100*(sosIref(1:90,1)-KokIref0)./KokIref0;
sosIreferr90=100*(sosIref(1:90,2)-KokIref90)./KokIref90;
sosIreferr180=100*(sosIref(1:90,3)-KokIref180)./KokIref180;

plot(Kokangref,sosIreferr0,'-+b',...
    Kokangref,sosIreferr90,'-+k',...
    Kokangref,sosIreferr180,'-+m','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at TOA')
legend('SOS I \phi=0',...
    'SOS I \phi=90',...
    'SOS I \phi=180')
ylim([-5 5])
set(gca,'FontSize',16)


figure(6)
semilogy(Kokangtran,sosItran(1:90,1),'+b',Kokangtran,KokItran0,'ob',...
    Kokangtran,sosItran(1:90,2),'+k',Kokangtran,KokItran90,'ok',...
    Kokangtran,sosItran(1:90,3),'+m',Kokangtran,KokItran180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Radiance at BOA')
legend('SOS I \phi=0','Kok benchmark \phi=0',...
    'SOS I \phi=90','Kok benchmark \phi=90',...
    'SOS I \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)


figure(7)
plot(Kokangtran,sosQtran(1:90,1),'+b',Kokangtran,KokQtran0,'ob',...
    Kokangtran,sosQtran(1:90,2),'+k',Kokangtran,KokQtran90,'ok',...
    Kokangtran,sosQtran(1:90,3),'+m',Kokangtran,KokQtran180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Q at BOA')
legend('SOS Q \phi=0','Kok benchmark \phi=0',...
    'SOS Q \phi=90','Kok benchmark \phi=90',...
    'SOS Q \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(8)
plot(Kokangtran,sosUtran(1:90,1),'+b',Kokangtran,KokUtran0,'ob',...
    Kokangtran,sosUtran(1:90,2),'+k',Kokangtran,KokUtran90,'ok',...
    Kokangtran,sosUtran(1:90,3),'+m',Kokangtran,KokUtran180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('U at BOA')
legend('SOS U \phi=0','Kok benchmark \phi=0',...
    'SOS U \phi=90','Kok benchmark \phi=90',...
    'SOS U \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(9)
plot(Kokangtran,sosVtran(1:90,1),'+b',Kokangtran,KokVtran0,'ob',...
    Kokangtran,sosVtran(1:90,2),'+k',Kokangtran,KokVtran90,'ok',...
    Kokangtran,sosVtran(1:90,3),'+m',Kokangtran,KokVtran180,'om','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('V at BOA')
legend('SOS V \phi=0','Kok benchmark \phi=0',...
    'SOS V \phi=90','Kok benchmark \phi=90',...
    'SOS V \phi=180','Kok benchmark \phi=180')
set(gca,'FontSize',16)

figure(10)

sosItranerr0=100*(sosItran(1:90,1)-KokItran0)./KokItran0;
sosItranerr90=100*(sosItran(1:90,2)-KokItran90)./KokItran90;
sosItranerr180=100*(sosItran(1:90,3)-KokItran180)./KokItran180;

plot(Kokangtran,sosItranerr0,'-+b',...
    Kokangtran,sosItranerr90,'-+k',...
    Kokangtran,sosItranerr180,'-+m','LineWidth',2)
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at BOA')
legend('SOS I \phi=0',...
    'SOS I \phi=90',...
    'SOS I \phi=180')
%ylim([-5 5])
set(gca,'FontSize',16)

% figure(11)
% semilogy(Kokangref,mcIAref(1:90,1),'+b',Kokangref,KokIref0,'ob',...
%     Kokangref,mcIAref(1:90,2),'+k',Kokangref,KokIref90,'ok',...
%     Kokangref,mcIAref(1:90,3),'+m',Kokangref,KokIref180,'om')
% xlabel('Viewing angle Deg')
% ylabel('Radiance at TOA')
% legend('mc I \phi=0','Kok benchmark \phi=0',...
%     'mc I \phi=90','Kok benchmark \phi=90',...
%     'mc I \phi=180','Kok benchmark \phi=180')
% 
% 
% figure(12)
% plot(Kokangref,mcQAref(1:90,1),'+b',Kokangref,KokQref0,'ob',...
%     Kokangref,mcQAref(1:90,2),'+k',Kokangref,KokQref90,'ok',...
%     Kokangref,mcQAref(1:90,3),'+m',Kokangref,KokQref180,'om')
% xlabel('Viewing angle Deg')
% ylabel('Q at TOA')
% legend('mc Q \phi=0','Kok benchmark \phi=0',...
%     'mc Q \phi=90','Kok benchmark \phi=90',...
%     'mc Q \phi=180','Kok benchmark \phi=180')
% 
% figure(13)
% plot(Kokangref,mcUAref(1:90,1),'+b',Kokangref,KokUref0,'ob',...
%     Kokangref,mcUAref(1:90,2),'+k',Kokangref,KokUref90,'ok',...
%     Kokangref,mcUAref(1:90,3),'+m',Kokangref,KokUref180,'om')
% xlabel('Viewing angle Deg')
% ylabel('U at TOA')
% legend('mc U \phi=0','Kok benchmark \phi=0',...
%     'mc U \phi=90','Kok benchmark \phi=90',...
%     'mc U \phi=180','Kok benchmark \phi=180')
% 
% figure(14)
% plot(Kokangref,mcVAref(1:90,1),'+b',Kokangref,KokVref0,'ob',...
%     Kokangref,mcVAref(1:90,2),'+k',Kokangref,KokVref90,'ok',...
%     Kokangref,mcVAref(1:90,3),'+m',Kokangref,KokVref180,'om')
% xlabel('Viewing angle Deg')
% ylabel('V at TOA')
% legend('mc V \phi=0','Kok benchmark \phi=0',...
%     'mc V \phi=90','Kok benchmark \phi=90',...
%     'mc V \phi=180','Kok benchmark \phi=180')
% 
% figure(15)
% 
% mcIreferr0=100*(mcIAref(1:90,1)-KokIref0)./KokIref0;
% mcIreferr90=100*(mcIAref(1:90,2)-KokIref90)./KokIref90;
% mcIreferr180=100*(mcIAref(1:90,3)-KokIref180)./KokIref180;
% 
% plot(Kokangref,mcIreferr0,'+b',...
%     Kokangref,mcIreferr90,'+k',...
%     Kokangref,mcIreferr180,'+m')
% xlabel('Viewing angle Deg')
% ylabel('Radiance % Error at TOA')
% legend('mc I \phi=0',...
%     'mc I \phi=90',...
%     'mc I \phi=180')
% 
% 
% 
% figure(16)
% semilogy(Kokangtran,mcIAtran(1:90,1),'+b',Kokangtran,KokItran0,'ob',...
%     Kokangtran,mcIAtran(1:90,2),'+k',Kokangtran,KokItran90,'ok',...
%     Kokangtran,mcIAtran(1:90,3),'+m',Kokangtran,KokItran180,'om')
% xlabel('Viewing angle Deg')
% ylabel('Radiance at BOA')
% legend('mc I \phi=0','Kok benchmark \phi=0',...
%     'mc I \phi=90','Kok benchmark \phi=90',...
%     'mc I \phi=180','Kok benchmark \phi=180')
% 
% 
% figure(17)
% plot(Kokangtran,mcQAtran(1:90,1),'+b',Kokangtran,KokQtran0,'ob',...
%     Kokangtran,mcQAtran(1:90,2),'+k',Kokangtran,KokQtran90,'ok',...
%     Kokangtran,mcQAtran(1:90,3),'+m',Kokangtran,KokQtran180,'om')
% xlabel('Viewing angle Deg')
% ylabel('Q at BOA')
% legend('mc Q \phi=0','Kok benchmark \phi=0',...
%     'mc Q \phi=90','Kok benchmark \phi=90',...
%     'mc Q \phi=180','Kok benchmark \phi=180')
% 
% figure(18)
% plot(Kokangtran,mcUAtran(1:90,1),'+b',Kokangtran,KokUtran0,'ob',...
%     Kokangtran,mcUAtran(1:90,2),'+k',Kokangtran,KokUtran90,'ok',...
%     Kokangtran,mcUAtran(1:90,3),'+m',Kokangtran,KokUtran180,'om')
% xlabel('Viewing angle Deg')
% ylabel('U at BOA')
% legend('mc U \phi=0','Kok benchmark \phi=0',...
%     'mc U \phi=90','Kok benchmark \phi=90',...
%     'mc U \phi=180','Kok benchmark \phi=180')
% 
% figure(19)
% plot(Kokangtran,mcVAtran(1:90,1),'+b',Kokangtran,KokVtran0,'ob',...
%     Kokangtran,mcVAtran(1:90,2),'+k',Kokangtran,KokVtran90,'ok',...
%     Kokangtran,mcVAtran(1:90,3),'+m',Kokangtran,KokVtran180,'om')
% xlabel('Viewing angle Deg')
% ylabel('V at BOA')
% legend('mc V \phi=0','Kok benchmark \phi=0',...
%     'mc V \phi=90','Kok benchmark \phi=90',...
%     'mc V \phi=180','Kok benchmark \phi=180')
% 
% figure(20)
% 
% mcItranerr0=100*(mcIAtran(1:90,1)-KokItran0)./KokItran0;
% mcItranerr90=100*(mcIAtran(1:90,2)-KokItran90)./KokItran90;
% mcItranerr180=100*(mcIAtran(1:90,3)-KokItran180)./KokItran180;
% 
% plot(Kokangtran,mcItranerr0,'+b',...
%     Kokangtran,mcItranerr90,'+k',...
%     Kokangtran,mcItranerr180,'+m')
% xlabel('Viewing angle Deg')
% ylabel('Radiance % Error at BOA')
% legend('mc I \phi=0',...
%     'mc I \phi=90',...
%     'mc I \phi=180')
