%% Start of script

close all;                          % close all figures
% clear;                              % clear all variables
% clc;                                % clear the command terminal
maindir=fileparts(fileparts(which('TestMethodsFromScript.m')));
%%
CompareMethods=[2:6]; % 1=YoungSooSuh, 2=madg, 3=valenti, 4=PRCF, 5=MadgwickAHRSclanek, 6=modif2
set=8;   %set1=1 set2=2 set12=3 set123=4 ALSdataset=5 synthetic2=6 synthetic3=7

optimAttempts=0;
RMS=0;
optimFreezeParam=0;
IMU=1; %usemag or not

if(set<5)
    load('Justa_dataset.mat')
elseif(set==5)
    load('ALS_dataset.mat')
elseif(set==6)
    load('Synthetic2.mat')
elseif(set==7)
    load('Synthetic3.mat')
elseif(set==8)
   data= readmatrix(fullfile(maindir,'Datasets','move01.csv'),'NumHeaderLines',5);
   Accelerometer=data(:,18:20)/1000;
   Gyroscope=data(:,21:23)*pi/180;
   Magnetometer=data(:,24:26)*1e5; %need to check conversion
   GlobalA=data(:,27:29)*pi/180;

   GlobalMarkers=readmatrix('framekinematics.txt','NumHeaderLines',5);
    for c=2:4
        nanx = isnan(GlobalMarkers(:,c));
        t    = 1:numel(GlobalMarkers(:,c));
        GlobalMarkers(nanx,c)= interp1(t(~nanx), GlobalMarkers(~nanx,c), t(nanx));
    end
    GlobalMarkers(:,1)=[];
    GlobalMarkers = unwrap(GlobalMarkers*pi/180,[],1)

   GlobalMarkers= interp1(t,GlobalMarkers,linspace(1,length(t),length(GlobalA))); 
   qViconReference=eul2quat(GlobalMarkers, 'XYZ');
   time=0:1/225:(length(t))/225;


%load('Wand.mat')
end

quaternionCountJ = zeros(length(time), 4);
test=zeros(length(time), 4);
test2=zeros(length(time), 4);

switch(set)
    case 1
        start=1;
        myEnd=2800-1;
    case 2
        start=2800;
        myEnd=5300-1;
    case 3
        start=1;
        myEnd=5300-1;
    otherwise
        start=1;
        myEnd=length(time)-1;
end


%%
AHRSs={};
if(any(CompareMethods==1))
    AHRS = YoungSooSuh_AHRS();
    AHRS.rg= 0.0370;
    AHRS.ra= 0.000031;
    AHRS.rm= 0.0000235;    

    
%     AHRS.rg= 0.0005631158677499999;
%     AHRS.ra= 0.0000018351604003246236;
%     AHRS.rm= 0.0001794678715498521;  
    
    AHRSs = cat(2,AHRSs,{AHRS});
end
%%
if(any(CompareMethods==2))
%     switch(set)
%         case 1
%             AHRS = MadgwickAHRS3('Beta',0.011421);
%         case 2
%             AHRS = MadgwickAHRS3('Beta',0.03226);
%         case 3
%             AHRS = MadgwickAHRS3('Beta',0.016630);
%         case 4
%             AHRS = MadgwickAHRS3('Beta',0.02113);
%         case 6
%             AHRS = MadgwickAHRS3('Beta',0.0328);
%     end
    AHRS = MadgwickAHRS3('Beta',0.011);
    AHRSs = cat(2,AHRSs,{AHRS});
end
%%
if(any(CompareMethods==3))
    AHRS = Valenti_AHRS();
    AHRS.wAcc=0.00116;%
    AHRS.wMag=1.012e-7;%
    switch(set)
        case 4
            AHRS.wAcc=0.00608;%
            AHRS.wMag=1.2759e-04;%
    end
    AHRSs = cat(2,AHRSs,{AHRS});
end
%%
if(any(CompareMethods==4))
    AHRS = JustaAHRSPureFast();
%     AHRS.wAcc=0.001;%
%     AHRS.wMag=0.01;%
    AHRS.wAcc=0.002906;%
    AHRS.wMag=5.219;%
    AHRS.gain=.0254;
    AHRSs = cat(2,AHRSs,{AHRS});
end
%%
if(any(CompareMethods==5))
    switch(set)
        case 1
            AHRS = MadgwickAHRSclanek('Beta',0.00166);
        otherwise
            AHRS = MadgwickAHRSclanek('Beta',0.01);
    end
    AHRSs = cat(2,AHRSs,{AHRS});
end
%%
if(any(CompareMethods==6))
    AHRS = JustaAHRSPureFastConstantCorr();
    AHRS.wAcc=0.0017;%
    AHRS.wMag=0.0001;%
    AHRSs = cat(2,AHRSs,{AHRS});
end

%%
scalErr=zeros(5,length(AHRSs));
names=[];
errorAngles=[];
for method = 1:length(AHRSs)
    names = [names,string(class(AHRSs{method}))];
    AHRSs{method}.Quaternion = qViconReference(start,:);
    quaternionCountJ = zeros(myEnd, 4);
    

tic
while toc<10

                   [ AHRSs{method},e,param,optimAttempts,optimFreezeParam]=...
                       optimFilterParams( AHRSs{method},1,set,optimAttempts,optimFreezeParam,IMU);
    end

    for t = start:myEnd
        %AHRS.iter=AHRS.iter+1;
        if(t==start)
            dt=time(1);
        else
            dt=time(t)-time(t-1);
        end
        AHRSs{method}.SamplePeriod=dt;
        
        AHRSs{method}.Update(Gyroscope(t,:), Accelerometer(t,:),Magnetometer(t,:));
        quaternionCountJ(t, :) = AHRSs{method}.Quaternion;
        
        %         if(CompareMethods==4 || CompareMethods==1 || CompareMethods==6)
                     test(t,:)=AHRS.test;
        %             test2(t,:)=AHRS.test2;
        %         end
    end
    qErr=quaternProd(qViconReference((start:myEnd),:),quaternConj(quaternionCountJ(start:myEnd,:)));
    qErr(qErr(:,1)<0,:)=-qErr(qErr(:,1)<0,:);
    qErr=real(qErr);
    uhel=abs(2*atan2(sqrt(sum(qErr(:,2:4).^2,2)),qErr(:,1))*180/pi).^2;
    errorAngles=[errorAngles uhel];
    %scalErr(4,method)=0;
    scalErr(5,method)=mean(uhel);
end

figure('DefaultAxesFontSize',18,'Position',[50,50,floor(1900*1/1),600]);
plot(time(1:end-1),quaternionCountJ,'LineWidth',1.0);
% xline(70.92,'LineWidth',1.0);
% xline(39.3,'LineWidth',1.0);
title('Estimated quaternion FSCF')
ylabel('[-]')
xlabel('Time (s)')
legend('w','x','y','z');

[minim2,ind]= min(scalErr(5,:));
disp(['Relative to the best:', num2str(scalErr(5,:)/minim2)]);
figure
plot(time(start:myEnd),movmean(errorAngles,30));
ylabel('MAE Error (deg)')
xlabel('Time (s)')
legend(names);
%% End of script