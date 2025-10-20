clear;

% modify the filedir to your directory

ntangle=20*3;
refflag=1;
dataread=1;
szen=78.4630;
taua=0.5;
tauo=0.5;

SOS_testfileme='vSOS_10202025';
%SOS_testfileme='sos_N100';

MCfilenameTOA='MCRAY_W0_wO_TOA';
MCfilenameBOA='MCRAY_W0_wO_BOA';
MCfilenameTOO='MCRAY_W0_wO_TOO';
MCfilenameBOO='MCRAY_W0_wO_BOO';

if(dataread)
  fid = fopen(SOS_testfileme,'r');
  header=textscan(fid,'%s',23,'delimiter', '\n');
  SOSradTOA=textscan(fid,'%f %f %f %f %f %f',60,'delimiter', '\n');
  header2=textscan(fid,'%s',5,'delimiter', '\n');
  SOSradBOA=textscan(fid,'%f %f %f %f %f %f',60,'delimiter', '\n');
  header3=textscan(fid,'%s',5,'delimiter', '\n');
  SOSradTOO=textscan(fid,'%f %f %f %f %f %f',60,'delimiter', '\n');
  header4=textscan(fid,'%s',5,'delimiter', '\n');
  SOSradBOO=textscan(fid,'%f %f %f %f %f %f',60,'delimiter', '\n');
  fclose(fid);

 
  fid = fopen(MCfilenameTOA,'r');
  header5=textscan(fid,'%s',20,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  MCradTOA=textscan(fid,format,ntangle,'delimiter', '\n');
  fclose(fid);

  fid = fopen(MCfilenameBOA,'r');
  header6=textscan(fid,'%s',20,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  MCradBOA=textscan(fid,format,ntangle,'delimiter', '\n');
  fclose(fid);
  
  fid = fopen(MCfilenameTOO,'r');
  header7=textscan(fid,'%s',20,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  MCradTOO=textscan(fid,format,ntangle,'delimiter', '\n');
  fclose(fid);
  
  fid = fopen(MCfilenameBOO,'r');
  header7=textscan(fid,'%s',20,'delimiter', '\n');
  format = repmat('%f ', 1, 18);
  MCradBOO=textscan(fid,format,ntangle,'delimiter', '\n');
  fclose(fid);
  
  save('L60_val_data');
else
%  load 'L60_val_data'
end

sosmuTOA=SOSradTOA{1,1};
sosphiTOA=SOSradTOA{1,2};
sosITOA=SOSradTOA{1,3};
sosQTOA=SOSradTOA{1,4}./SOSradTOA{1,3};
sosUTOA=-SOSradTOA{1,5}./SOSradTOA{1,3};
sosVTOA=-SOSradTOA{1,6}./SOSradTOA{1,3};

sosmuTOA=reshape(sosmuTOA,20,3);
sosphiTOA=reshape(sosphiTOA,20,3);
sosITOA=reshape(sosITOA,20,3);
sosQTOA=reshape(sosQTOA,20,3);
sosUTOA=reshape(sosUTOA,20,3);
sosVTOA=reshape(sosVTOA,20,3);

sosmuBOA=SOSradBOA{1,1};
sosphiBOA=SOSradBOA{1,2};
sosIBOA=SOSradBOA{1,3};
sosQBOA=SOSradBOA{1,4}./SOSradBOA{1,3};
sosUBOA=-SOSradBOA{1,5}./SOSradBOA{1,3};
sosVBOA=-SOSradBOA{1,6}./SOSradBOA{1,3};

sosmuBOA=reshape(sosmuBOA,20,3);
sosphiBOA=reshape(sosphiBOA,20,3);
sosIBOA=reshape(sosIBOA,20,3);
sosQBOA=reshape(sosQBOA,20,3);
sosUBOA=reshape(sosUBOA,20,3);
sosVBOA=reshape(sosVBOA,20,3);

sosmuTOO=SOSradTOO{1,1};
sosphiTOO=SOSradTOO{1,2};
sosITOO=SOSradTOO{1,3};
sosQTOO=SOSradTOO{1,4}./SOSradTOO{1,3};
sosUTOO=-SOSradTOO{1,5}./SOSradTOO{1,3};
sosVTOO=-SOSradTOO{1,6}./SOSradTOO{1,3};

sosmuTOO=reshape(sosmuTOO,20,3);
sosphiTOO=reshape(sosphiTOO,20,3);
sosITOO=reshape(sosITOO,20,3);
sosQTOO=reshape(sosQTOO,20,3);
sosUTOO=reshape(sosUTOO,20,3);
sosVTOO=reshape(sosVTOO,20,3);

sosmuBOO=SOSradBOO{1,1};
sosphiBOO=SOSradBOO{1,2};
sosIBOO=SOSradBOO{1,3};
sosQBOO=SOSradBOO{1,4}./SOSradBOO{1,3};
sosUBOO=-SOSradBOO{1,5}./SOSradBOO{1,3};
sosVBOO=-SOSradBOO{1,6}./SOSradBOO{1,3};

sosmuBOO=reshape(sosmuBOO,20,3);
sosphiBOO=reshape(sosphiBOO,20,3);
sosIBOO=reshape(sosIBOO,20,3);
sosQBOO=reshape(sosQBOO,20,3);
sosUBOO=reshape(sosUBOO,20,3);
sosVBOO=reshape(sosVBOO,20,3);

mcangTOA=MCradTOA{1,1};
mcphiTOA=MCradTOA{1,2};

mcITOA=MCradTOA{1,3}.*cosd(szen)*pi;
mcQTOA=MCradTOA{1,7}./MCradTOA{1,3};
mcUTOA=MCradTOA{1,11}./MCradTOA{1,3};
mcVTOA=MCradTOA{1,15}./MCradTOA{1,3};
mcangTOA=reshape(mcangTOA,20,3);
mcITOA=reshape(mcITOA,20,3);
mcQTOA=reshape(mcQTOA,20,3);
mcUTOA=reshape(mcUTOA,20,3);
mcVTOA=reshape(mcVTOA,20,3);

mcangBOA=MCradBOA{1,1};
mcphiBOA=MCradBOA{1,2};
mcIBOA=MCradBOA{1,3}.*cosd(szen)*pi;
mcQBOA=MCradBOA{1,7}./MCradBOA{1,3};
mcUBOA=MCradBOA{1,11}./MCradBOA{1,3};
mcVBOA=MCradBOA{1,15}./MCradBOA{1,3};

mcangBOA=reshape(mcangBOA,20,3);
mcIBOA=reshape(mcIBOA,20,3);
mcQBOA=reshape(mcQBOA,20,3);
mcUBOA=reshape(mcUBOA,20,3);
mcVBOA=reshape(mcVBOA,20,3);

mcangTOO=MCradTOO{1,1};
mcphiTOO=MCradTOO{1,2};
mcITOO=MCradTOO{1,3}.*cosd(szen)*pi;
mcQTOO=MCradTOO{1,7}./MCradTOO{1,3};
mcUTOO=MCradTOO{1,11}./MCradTOO{1,3};
mcVTOO=MCradTOO{1,15}./MCradTOO{1,3};

mcangTOO=reshape(mcangTOO,20,3);
mcITOO=reshape(mcITOO,20,3);
mcQTOO=reshape(mcQTOO,20,3);
mcUTOO=reshape(mcUTOO,20,3);
mcVTOO=reshape(mcVTOO,20,3);

mcangBOO=MCradBOO{1,1};
mcphiBOO=MCradBOO{1,2};
mcIBOO=MCradBOO{1,3}.*cosd(szen)*pi;
mcQBOO=MCradBOO{1,7}./MCradBOO{1,3};
mcUBOO=MCradBOO{1,11}./MCradBOO{1,3};
mcVBOO=MCradBOO{1,15}./MCradBOO{1,3};

mcangBOO=reshape(mcangBOO,20,3);
mcIBOO=reshape(mcIBOO,20,3);
mcQBOO=reshape(mcQBOO,20,3);
mcUBOO=reshape(mcUBOO,20,3);
mcVBOO=reshape(mcVBOO,20,3);

figure(1)
subplot(2,2,1)
semilogy(sosmuTOA(1:10,1),sosITOA(1:10,1),'+b',mcangTOA(1:10,1),mcITOA(1:10,1),'ob',...
    sosmuTOA(1:10,2),sosITOA(1:10,2),'+k',mcangTOA(1:10,2),mcITOA(1:10,2),'ok',...
    sosmuTOA(1:10,3),sosITOA(1:10,3),'+m',mcangTOA(1:10,3),mcITOA(1:10,3),'om')
xlabel('Viewing angle Deg')
ylabel('Radiance at TOA')
legend('SOS \phi=0','MC \phi=0',...
          'SOS \phi=90','MC \phi=90',...
          'SOS \phi=180','MC \phi=180','Location','NorthWest')

%figure(2)
subplot(2,2,2)
plot(sosmuTOA(1:10,1),sosQTOA(1:10,1),'+b',mcangTOA(1:10,1),mcQTOA(1:10,1),'ob',...
    sosmuTOA(1:10,2),sosQTOA(1:10,2),'+k',mcangTOA(1:10,2),mcQTOA(1:10,2),'ok',...
    sosmuTOA(1:10,3),sosQTOA(1:10,3),'+m',mcangTOA(1:10,3),mcQTOA(1:10,3),'om')
xlabel('Viewing angle Deg')
ylabel('Q/I at TOA')
% legend('SOS Q \phi=0','MC Q \phi=0',...
%           'SOS Q \phi=90','MC Q\phi=90',...
%           'SOS Q \phi=180','MC Q \phi=180')

%figure(3)
subplot(2,2,3)

plot(sosmuTOA(1:10,1),sosUTOA(1:10,1),'+b',mcangTOA(1:10,1),mcUTOA(1:10,1),'ob',...
    sosmuTOA(1:10,2),sosUTOA(1:10,2),'+k',mcangTOA(1:10,2),mcUTOA(1:10,2),'ok',...
    sosmuTOA(1:10,3),sosUTOA(1:10,3),'+m',mcangTOA(1:10,3),mcUTOA(1:10,3),'om')
xlabel('Viewing angle Deg')
ylabel('U/I at TOA')
% legend('SOS U \phi=0','MC U \phi=0',...
%           'SOS U \phi=90','MC U\phi=90',...
%           'SOS U \phi=180','MC U \phi=180')

%figure(4)
subplot(2,2,4)

plot(sosmuTOA(1:10,1),sosVTOA(1:10,1),'+b',mcangTOA(1:10,1),mcVTOA(1:10,1),'ob',...
    sosmuTOA(1:10,2),sosVTOA(1:10,2),'+k',mcangTOA(1:10,2),mcVTOA(1:10,2),'ok',...
    sosmuTOA(1:10,3),sosVTOA(1:10,3),'+m',mcangTOA(1:10,3),mcVTOA(1:10,3),'om')
xlabel('Viewing angle Deg')
ylabel('V/I at TOA')
% legend('SOS V \phi=0','MC V \phi=0',...
%           'SOS V \phi=90','MC V\phi=90',...
%           'SOS V \phi=180','MC V \phi=180')

figure(2)

sosITOAerr0=100*(sosITOA(1:10,1)-mcITOA(1:10,1))./mcITOA(1:10,1);
sosITOAerr90=100*(sosITOA(1:10,2)-mcITOA(1:10,2))./mcITOA(1:10,2);
sosITOAerr180=100*(sosITOA(1:10,3)-mcITOA(1:10,3))./mcITOA(1:10,2);

plot(sosmuTOA(1:10,1),sosITOAerr0,'+-b',...
    sosmuTOA(1:10,2),sosITOAerr90,'+-k',...
    sosmuTOA(1:10,3),sosITOAerr180,'+-m')
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at TOA')
 legend('SOS \phi=0',...
     'SOS \phi=90',...
     'SOS \phi=180')
%ylim([-1 1])


figure(3)
subplot(2,2,1)
semilogy(sosmuBOA(1:20,1),sosIBOA(1:20,1),'+b',mcangBOA(1:20,1),mcIBOA(1:20,1),'ob',...
    sosmuBOA(1:20,2),sosIBOA(1:20,2),'+k',mcangBOA(1:20,2),mcIBOA(1:20,2),'ok',...
    sosmuBOA(1:20,3),sosIBOA(1:20,3),'+m',mcangBOA(1:20,3),mcIBOA(1:20,3),'om')
xlabel('Viewing angle Deg')
ylabel('Radiance at BOA')
legend('SOS \phi=0','MC \phi=0',...
          'SOS \phi=90','MC \phi=90',...
          'SOS \phi=180','MC \phi=180','Location','NorthWest')


%figure(7)
subplot(2,2,2)

plot(sosmuBOA(1:20,1),sosQBOA(1:20,1),'+b',mcangBOA(1:20,1),mcQBOA(1:20,1),'ob',...
    sosmuBOA(1:20,2),sosQBOA(1:20,2),'+k',mcangBOA(1:20,2),mcQBOA(1:20,2),'ok',...
    sosmuBOA(1:20,3),sosQBOA(1:20,3),'+m',mcangBOA(1:20,3),mcQBOA(1:20,3),'om')
xlabel('Viewing angle Deg')
ylabel('Q/I at BOA')
% legend('SOS Q \phi=0','MC Q \phi=0',...
%           'SOS Q \phi=90','MC Q\phi=90',...
%           'SOS Q \phi=180','MC Q \phi=180')

%figure(8)
subplot(2,2,3)

plot(sosmuBOA(1:20,1),sosUBOA(1:20,1),'+b',mcangBOA(1:20,1),mcUBOA(1:20,1),'ob',...
    sosmuBOA(1:20,2),sosUBOA(1:20,2),'+k',mcangBOA(1:20,2),mcUBOA(1:20,2),'ok',...
    sosmuBOA(1:20,3),sosUBOA(1:20,3),'+m',mcangBOA(1:20,3),mcUBOA(1:20,3),'om')
xlabel('Viewing angle Deg')
ylabel('U/I at BOA')
% legend('SOS U \phi=0','MC U \phi=0',...
%           'SOS U \phi=90','MC U\phi=90',...
%           'SOS U \phi=180','MC U \phi=180')

%figure(9)
subplot(2,2,4)

plot(sosmuBOA(1:20,1),sosVBOA(1:20,1),'+b',mcangBOA(1:20,1),mcVBOA(1:20,1),'ob',...
    sosmuBOA(1:20,2),sosVBOA(1:20,2),'+k',mcangBOA(1:20,2),mcVBOA(1:20,2),'ok',...
    sosmuBOA(1:20,3),sosVBOA(1:20,3),'+m',mcangBOA(1:20,3),mcVBOA(1:20,3),'om')
xlabel('Viewing angle Deg')
ylabel('V/I at BOA')
% legend('SOS V \phi=0','MC V \phi=0',...
%           'SOS V \phi=90','MC V\phi=90',...
%           'SOS V \phi=180','MC V \phi=180')

figure(4)

sosIBOAerr0=100*(sosIBOA(1:20,1)-mcIBOA(1:20,1))./mcIBOA(1:20,1);
sosIBOAerr90=100*(sosIBOA(1:20,2)-mcIBOA(1:20,2))./mcIBOA(1:20,2);
sosIBOAerr180=100*(sosIBOA(1:20,3)-mcIBOA(1:20,3))./mcIBOA(1:20,2);

plot(sosmuBOA(1:20,1),sosIBOAerr0,'+-b',...
    sosmuBOA(1:20,2),sosIBOAerr90,'+-k',...
    sosmuBOA(1:20,3),sosIBOAerr180,'+-m')
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at BOA')
legend('SOS \phi=0',...
    'SOS \phi=90',...
    'SOS \phi=180')
%ylim([-1 1])

figure(10)
sosITOO(sosITOO==0)=NaN;
mcITOO(mcITOO==0)=NaN;
subplot(2,2,1)
plot(sosmuTOO(1:20,1),sosITOO(1:20,1),'+b',mcangTOO(1:20,1),mcITOO(1:20,1),'ob',...
    sosmuTOO(1:20,2),sosITOO(1:20,2),'+k',mcangTOO(1:20,2),mcITOO(1:20,2),'ok',...
    sosmuTOO(1:20,3),sosITOO(1:20,3),'+m',mcangTOO(1:20,3),mcITOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Radiance at TOO')
legend('SOS \phi=0','MC \phi=0',...
          'SOS \phi=90','MC \phi=90',...
          'SOS \phi=180','MC \phi=180','Location','NorthWest')


%figure(2)
subplot(2,2,2)
plot(sosmuTOO(1:20,1),sosQTOO(1:20,1),'+b',mcangTOO(1:20,1),mcQTOO(1:20,1),'ob',...
    sosmuTOO(1:20,2),sosQTOO(1:20,2),'+k',mcangTOO(1:20,2),mcQTOO(1:20,2),'ok',...
    sosmuTOO(1:20,3),sosQTOO(1:20,3),'+m',mcangTOO(1:20,3),mcQTOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Q at TOO')
% legend('SOS Q \phi=0','MC Q \phi=0',...
%           'SOS Q \phi=90','MC Q\phi=90',...
%           'SOS Q \phi=180','MC Q \phi=180')

%figure(3)
subplot(2,2,3)

plot(sosmuTOO(1:20,1),sosUTOO(1:20,1),'+b',mcangTOO(1:20,1),mcUTOO(1:20,1),'ob',...
    sosmuTOO(1:20,2),sosUTOO(1:20,2),'+k',mcangTOO(1:20,2),mcUTOO(1:20,2),'ok',...
    sosmuTOO(1:20,3),sosUTOO(1:20,3),'+m',mcangTOO(1:20,3),mcUTOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('U at TOO')
% legend('SOS U \phi=0','MC U \phi=0',...
%           'SOS U \phi=90','MC U\phi=90',...
%           'SOS U \phi=180','MC U \phi=180')

%figure(4)
subplot(2,2,4)

plot(sosmuTOO(1:20,1),sosVTOO(1:20,1),'+b',mcangTOO(1:20,1),mcVTOO(1:20,1),'ob',...
    sosmuTOO(1:20,2),sosVTOO(1:20,2),'+k',mcangTOO(1:20,2),mcVTOO(1:20,2),'ok',...
    sosmuTOO(1:20,3),sosVTOO(1:20,3),'+m',mcangTOO(1:20,3),mcVTOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('V at TOO')
% legend('SOS V \phi=0','MC V \phi=0',...
%           'SOS V \phi=90','MC V\phi=90',...
%           'SOS V \phi=180','MC V \phi=180')

figure(6)

sosITOOerr0=100*(sosITOO(1:20,1)-mcITOO(1:20,1))./mcITOO(1:20,1);
sosITOOerr90=100*(sosITOO(1:20,2)-mcITOO(1:20,2))./mcITOO(1:20,2);
sosITOOerr180=100*(sosITOO(1:20,3)-mcITOO(1:20,3))./mcITOO(1:20,2);

plot(sosmuTOO(1:20,1),sosITOOerr0,'+-b',...
    sosmuTOO(1:20,2),sosITOOerr90,'+-k',...
    sosmuTOO(1:20,3),sosITOOerr180,'+-m')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at TOO')
 legend('SOS \phi=0',...
     'SOS \phi=90',...
     'SOS \phi=180')
%ylim([-1 1])


figure(7)
sosIBOO(sosIBOO==0)=NaN;
mcIBOO(mcIBOO==0)=NaN;
subplot(2,2,1)
plot(sosmuBOO(1:20,1),sosIBOO(1:20,1),'+b',mcangBOO(1:20,1),mcIBOO(1:20,1),'ob',...
    sosmuBOO(1:20,2),sosIBOO(1:20,2),'+k',mcangBOO(1:20,2),mcIBOO(1:20,2),'ok',...
    sosmuBOO(1:20,3),sosIBOO(1:20,3),'+m',mcangBOO(1:20,3),mcIBOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Radiance at BOO')
legend('SOS \phi=0','MC \phi=0',...
          'SOS \phi=90','MC \phi=90',...
          'SOS \phi=180','MC \phi=180','Location','NorthWest')


%figure(2)
subplot(2,2,2)
plot(sosmuBOO(1:20,1),sosQBOO(1:20,1),'+b',mcangBOO(1:20,1),mcQBOO(1:20,1),'ob',...
    sosmuBOO(1:20,2),sosQBOO(1:20,2),'+k',mcangBOO(1:20,2),mcQBOO(1:20,2),'ok',...
    sosmuBOO(1:20,3),sosQBOO(1:20,3),'+m',mcangBOO(1:20,3),mcQBOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Q at BOO')
% legend('SOS Q \phi=0','MC Q \phi=0',...
%           'SOS Q \phi=90','MC Q\phi=90',...
%           'SOS Q \phi=180','MC Q \phi=180')

%figure(3)
subplot(2,2,3)

plot(sosmuBOO(1:20,1),sosUBOO(1:20,1),'+b',mcangBOO(1:20,1),mcUBOO(1:20,1),'ob',...
    sosmuBOO(1:20,2),sosUBOO(1:20,2),'+k',mcangBOO(1:20,2),mcUBOO(1:20,2),'ok',...
    sosmuBOO(1:20,3),sosUBOO(1:20,3),'+m',mcangBOO(1:20,3),mcUBOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('U at BOO')
% legend('SOS U \phi=0','MC U \phi=0',...
%           'SOS U \phi=90','MC U\phi=90',...
%           'SOS U \phi=180','MC U \phi=180')

%figure(4)
subplot(2,2,4)

plot(sosmuBOO(1:20,1),sosVBOO(1:20,1),'+b',mcangBOO(1:20,1),mcVBOO(1:20,1),'ob',...
    sosmuBOO(1:20,2),sosVBOO(1:20,2),'+k',mcangBOO(1:20,2),mcVBOO(1:20,2),'ok',...
    sosmuBOO(1:20,3),sosVBOO(1:20,3),'+m',mcangBOO(1:20,3),mcVBOO(1:20,3),'om')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('V at BOO')
% legend('SOS V \phi=0','MC V \phi=0',...
%           'SOS V \phi=90','MC V\phi=90',...
%           'SOS V \phi=180','MC V \phi=180')

figure(8)

sosIBOOerr0=100*(sosIBOO(1:20,1)-mcIBOO(1:20,1))./mcIBOO(1:20,1);
sosIBOOerr90=100*(sosIBOO(1:20,2)-mcIBOO(1:20,2))./mcIBOO(1:20,2);
sosIBOOerr180=100*(sosIBOO(1:20,3)-mcIBOO(1:20,3))./mcIBOO(1:20,2);

plot(sosmuBOO(1:20,1),sosIBOOerr0,'+-b',...
    sosmuBOO(1:20,2),sosIBOOerr90,'+-k',...
    sosmuBOO(1:20,3),sosIBOOerr180,'+-m')
xlim([0 180]);
xlabel('Viewing angle Deg')
ylabel('Radiance % Error at BOO')
 legend('SOS \phi=0',...
     'SOS \phi=90',...
     'SOS \phi=180')
%ylim([-1 1])