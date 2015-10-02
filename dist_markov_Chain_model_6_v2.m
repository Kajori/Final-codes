% 3 phase fault  -no
% Gaussian Noise - YES
% malacious node - if present node 5
% there will be a central guy.. 
MAX_DEV = 1e-6;

clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Defining the constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO_AREA=5;
for ATTACK_REGION=1:NO_AREA;
%for Attack_per=-5:0.5:5
%Attack=0.01*Attack_per; %0.01=1 % attack
NO_OF_EIGEN_VALUES=20;
deg=2*NO_OF_EIGEN_VALUES; 
for DV=-5:0.5:5

    DESIRED_VAL=ones(deg)*(DV*(exp(-02)));

SCALE = 0; %0.001=0.1 % noise
MAX_DEV = 1e-6;


PMU_PER_AREA=3;
SIZE_OF_voltage_data=1000;


TRUE=1;
FALSE=0;
area_line_no=[53,60,61;
              30,48,63;
              62,66,67;
              64,65,68;
              54,56,57];

Y=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);
Y_org=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);

m=60; % height of the Hankel matrix
H_pmu=zeros(NO_AREA,PMU_PER_AREA,m,deg);
C_pmu=zeros(NO_AREA,PMU_PER_AREA,m);

H_pdc=zeros(NO_AREA,PMU_PER_AREA*m,deg);
C_pdc=zeros(NO_AREA,PMU_PER_AREA*m);
%display('Gaussian Noise no 3 phase fault');

rho=10^-9;
disp('-----------------------------------------------');
tic;

TT=PMU_PER_AREA*m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Loading data from the input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
for row=1:NO_AREA
     for col=1:PMU_PER_AREA
           file_name=sprintf('no_3_phase_fault_area_%d_line_%d.txt',row,area_line_no(row,col));
           disp(file_name);
           Y_org(row,col,:)=importdata(file_name,'\n');
      end %end of for col
end %end of for row
%disp(sprintf( ' size(Y_org) =  %d ',size(Y_org)));
disp(' Loaded data from the input files ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Write statistics to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('stat.txt', 'a+');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %         adding Gausian Noise  to the data
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row=1:NO_AREA 
    for col=1:PMU_PER_AREA
             noise = randn(1); % noise with mean=0 and std=1;
             Y(row,col,:)= ( 1+  noise* SCALE )* Y_org(row,col,:);
    end %end of for col
end %end of for row        
%disp(sprintf( ' size(Y) = %d',size(Y)));

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Hankel Matrix Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row=1:NO_AREA
    for col=1:PMU_PER_AREA
               temp_1=reshape(Y(row,col,:),SIZE_OF_voltage_data,1);
               temp_2=fliplr(Hankel(temp_1,m,deg));
               H_pmu(row,col,:,:)=temp_2;
               temp_3=reshape(Y(row,col,1:m),1,m);
               C_pmu(row,col,:)=temp_3.';      
    end %end of for col
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Combining Data from PMU to PDC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for col=1:PMU_PER_AREA
          H_pdc(row,(col-1)*m+1: col*m,:)= reshape(H_pmu(row,col,:,:),m,deg);
          C_pdc(row,(col-1)*m+1: col*m) = reshape(C_pmu(row,col,:),m,1);
     end %end of for col
     
     %disp(text);
     %disp(size(H_pdc));
     %disp(size(C_pdc));
 end %end of for row



N=20000;
filename=sprintf('Model_6_clever_attacker_noise_%s_per_desired_val_%se-02_ATTACK_REGION_%s_N_%d.txt',num2str(SCALE*100),num2str(DV),num2str(ATTACK_REGION),N);
disp(filename);
fid = fopen(filename,'w');
fprintf(fid,sprintf(' Always area %s under attack \n Default number of iterations = %d',num2str(ATTACK_REGION),N));

for rn=-1:5 %rn indicates the removed node
w=zeros(NO_AREA*deg,1);
x=ones(NO_AREA*deg,N); %Everyone has equal initial x
z=zeros(deg,N);
%prev=0;%alternate attack



for j=1:N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         if we have not reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %text=sprintf(' -------------- %d --------------',j);
    %disp('.');
    
    for row=1:NO_AREA
        if(row~=rn)
            H_r=reshape(H_pdc(row,1:TT,:),TT,deg);
            C_r=reshape(C_pdc(row,1:TT,:),TT,1);
            temp=(H_r'*C_r)-w((row-1)*deg+1:row*deg,j)+rho*z(:,j);
            x((row-1)*deg+1:row*deg,j+1)=(H_r'*H_r+rho*eye(size(H_r,2)))\temp;
            if(row==ATTACK_REGION && rn~=-1)
                %x((row-1)*deg+1:row*deg,j+1)=(1+Attack)*x((row-1)*deg+1:row*deg,j+1);
%               
               % random_num = randi(2)-1; %X = randi(imax) returns a pseudorandom scalar integer between 1 and imax.
               % x((row-1)*deg+1:row*deg,j+1)=(1+random_num*Attack)*x((row-1)*deg+1:row*deg,j+1);
                for k=1:deg
                    diff=DESIRED_VAL(k)-x((row-1)*deg+k,j+1);
                    %disp(diff);
                    x((row-1)*deg+1:row*deg,j+1)=(1+0.5*diff)*x((row-1)*deg+k,j+1);
                end
                 %x((row-1)*deg+1:row*deg,j+1)=(1+Attack)*x((row-1)*deg+1:row*deg,j+1);
                 %x((row-1)*deg+1:row*deg,j+1)=(1+prev*Attack)*x((row-1)*deg+1:row*deg,j+1);
                 %prev=mod(prev+1,2);
                %disp(random_num);
            end
        end
        
   end %end row
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %       Remove one node from the calculation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   temp=zeros(deg,1);
   count=0;
   for row=1:NO_AREA
       if(row~=rn)
            temp=temp+x((row-1)*deg+1:row*deg,j+1);
            count=count+1;
       end
   end
   temp=temp/count;
   z(:,j+1)=transpose(temp');
   
   for row=1:NO_AREA
       if(row~=rn)
           w((row-1)*deg+1:row*deg,j+1)= w((row-1)*deg+1:row*deg,j)+rho*(x((row-1)*deg+1:row*deg,j+1)-z(:,j+1)); 
       end
   end %end row
    
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %         Check if we have reached convergence
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   convergence_flag=TRUE;
   for row=1:NO_AREA*deg
        if(row~=rn)
            if ((abs(x(row,j+1)-x(row,j)))>0.00001)
                convergence_flag=FALSE;
            end
        end
   end
   if( convergence_flag==TRUE)
         display(sprintf('Convergence Reached at %d rn =%d',j,rn));
         fprintf(fid,sprintf('\n Convergence Reached at %d rn =%d',j,rn));
        break;
   end 
end %end of j

if (rn==-1)    no_attack=x(:,j);
elseif (rn==0) x_0=x(:,j);
elseif (rn==1) x_1=x(:,j);
elseif (rn==2) x_2=x(:,j);
elseif (rn==3) x_3=x(:,j);
elseif (rn==4) x_4=x(:,j);
elseif (rn==5) x_5=x(:,j);
end    

end %for rn=1:NO_AREA


disp(' \n \n\n ');

display(' \n \n %%%%%%%%%%%%%% AREA 1 values  %%%%%%%%%%%%%%%%');
fprintf(fid,' \n \n %%%%%%%%%%%%%% AREA 1 values %%%%%%%%%%%%%%%%');
fprintf(fid,'\n Row    NoAttack               NoRemoval               RemoveArea1          RemoveArea2            RemoveArea3            RemoveArea4             RemoveArea5');


for row=1:deg
     
     fprintf(fid,'\n%3d)  %10d   %20d   %20d   %20d   %20d   %20d   %20d  ', row, no_attack(row),x_0(row),x_1(row),x_2(row),x_3(row),x_4(row),x_5(row));
end


display(' \n \n %%%%%%%%%%%%%% AREA 2 values  %%%%%%%%%%%%%%%%');
fprintf(fid,' \n \n %%%%%%%%%%%%%% AREA 2 values %%%%%%%%%%%%%%%%');
fprintf(fid,'\n Row    NoAttack               NoRemoval               RemoveArea1          RemoveArea2            RemoveArea3            RemoveArea4             RemoveArea5');



for row=deg+1:2*deg
     
     fprintf(fid,'\n%3d)  %10d   %20d   %20d   %20d   %20d   %20d   %20d  ', row, no_attack(row),x_0(row),x_1(row),x_2(row),x_3(row),x_4(row),x_5(row));
end

display(' \n \n %%%%%%%%%%%%%% AREA 3 values  %%%%%%%%%%%%%%%%');
fprintf(fid,' \n \n %%%%%%%%%%%%%% AREA 3 values %%%%%%%%%%%%%%%%');
fprintf(fid,'\n Row    NoAttack               NoRemoval               RemoveArea1          RemoveArea2            RemoveArea3            RemoveArea4             RemoveArea5');

for row=2*deg+1:3*deg
     
    fprintf(fid,'\n%3d)  %10d   %20d   %20d   %20d   %20d   %20d   %20d  ', row, no_attack(row),x_0(row),x_1(row),x_2(row),x_3(row),x_4(row),x_5(row));
end

display(' \n \n %%%%%%%%%%%%%% AREA 4 values  %%%%%%%%%%%%%%%%');
fprintf(fid,' \n \n %%%%%%%%%%%%%% AREA 4 values %%%%%%%%%%%%%%%%');
fprintf(fid,'\n Row    NoAttack               NoRemoval               RemoveArea1          RemoveArea2            RemoveArea3            RemoveArea4             RemoveArea5');


for row=3*deg+1:4*deg
     
     fprintf(fid,'\n%3d)  %10d   %20d   %20d   %20d   %20d   %20d   %20d ', row, no_attack(row),x_0(row),x_1(row),x_2(row),x_3(row),x_4(row),x_5(row));
end

display(' \n \n %%%%%%%%%%%%%% AREA 5 values  %%%%%%%%%%%%%%%%');
fprintf(fid,' \n \n %%%%%%%%%%%%%% AREA 5 values %%%%%%%%%%%%%%%%');
fprintf(fid,'\n Row    NoAttack               NoRemoval               RemoveArea1          RemoveArea2            RemoveArea3            RemoveArea4             RemoveArea5');


for row=4*deg+1:5*deg
     
     fprintf(fid,'\n%3d)  %10d   %20d   %20d   %20d   %20d   %20d   %20d ', row, no_attack(row),x_0(row),x_1(row),x_2(row),x_3(row),x_4(row),x_5(row));
end

%RMS Error
fprintf(fid,sprintf('\n\n Reporting RMS Error '));
min_RMS_Region=99;
min_RMS_Val=99999;
min_RMS_Vec=ones(deg)*999;
for rn=-1:NO_AREA
    if (rn==-1)  
        Vec=no_attack;
        fprintf(fid,sprintf('\n\n ## No attack scenario'));
    elseif (rn==0) 
        Vec=x_0;
        fprintf(fid,sprintf('\n\n No area is removed from the computation'));
    elseif (rn==1) 
        Vec=x_1;
        fprintf(fid,sprintf('\n\n Area 1 is removed from the computation'));
    elseif (rn==2) 
        Vec=x_2;
        fprintf(fid,sprintf('\n\n Area 2 is removed from the computation'));
    elseif (rn==3) 
        Vec=x_3;
        fprintf(fid,sprintf('\n\n Area 3 is removed from the computation'));
    elseif (rn==4) 
        Vec=x_4;
        fprintf(fid,sprintf('\n\n Area 4 is removed from the computation'));
    elseif(rn==5)
        Vec=x_5;
        fprintf(fid,sprintf('\n\n Area 5 is removed from the computation'));
    end 
    
    
    V_diff=zeros(deg);
    for j=1:deg
        V_max=intmin('int64');
        V_min=intmax('int64');
        for i=1:NO_AREA
            if(i~=rn)
                if (Vec((i-1)*deg+j)>V_max)
                    V_max=Vec((i-1)*deg+j);
                end
                if (Vec((i-1)*deg+j)<V_min)
                    V_min=Vec((i-1)*deg+j);
                end
            end
        end
        V_diff(j)=(V_max-V_min); % Error
    end
    RMSE= rms(V_diff);
    fprintf(fid,sprintf('\n rms = %d ',RMSE(1)));
    if(rn>0)
        V_diff=zeros((NO_AREA-1)*deg);
    else
        V_diff=zeros(NO_AREA*deg);
    end
    k=1;
    for i=1:NO_AREA
        if(i~=rn)
            for j=1:deg
                V_diff(k)=no_attack((i-1)*deg+j)-Vec((i-1)*deg+j); % Error
                k=k+1;
            end
        end
    end
    RMSE= rms(V_diff);
    fprintf(fid,sprintf('\n Actual RMSE area %d is %d ',rn,RMSE(1)));
    
end%rn
    fclose(fid);
%ACTUAL RMS ERROR
end % Attack amount
end %Attack Region




%  RMSE_sum=0;
%     for i=1:NO_AREA-1
%         if(i~=rn)
%             for j=i+1:NO_AREA
%                 if(j~=rn)
%                     V_1=Vec((i-1)*deg+1:i*deg);
%                     V_2=Vec((j-1)*deg+1:j*deg);
%                     diff=zeros(deg);
%                     for k=1:deg
%                         diff(k)=(V_1(k) - V_2(k)); % Error
%                     end
%                     RMSE= rms(diff);
%                     RMSE_sum=RMSE+RMSE_sum;
%                     fprintf(fid,sprintf('\n i=%d j= %d rms = %d ',i,j,RMSE(1)));
%                     %fprintf(fid,sprintf('\n %d ',RMSE));
%                 end
%             end
%         end
%     end
% if (rn>-1 && min_RMS_Val>RMSE(1))
%        min_RMS_Region=rn;
%        min_RMS_Vec=Vec;
%        min_RMS_Val=RMSE(1);
%     end
% 	
