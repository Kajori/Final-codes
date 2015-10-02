% Gaussian Noise - No
% Individual - Yes
% Algorithm 1 'cal_outlier_individual_same_val' : all nodes get its neighbours values, and also the vlues of the
% other nodes. If it finds that the neighbours values are among the
% outliers then it replaces the neighbours values with its own value. It takes indivudual values
%
% Algorithm 2 'cal_outlier_individual_next_val' : all nodes get its neighbours values, and also the values of the
% other nodes. Then it replaces the neighbours values if they are in the
% outliers
clear all; clc; close all;
N=20000;
Attack_per=5;

SCALE = 0; %0.001 = 0.1 % noise
NO_AREA=5;
for ATTACK_REGION=5:NO_AREA
    for Attack_per=5:5
        Attack=0.01*Attack_per; %0.01=1 % attack
        %load n16 % PST model
        %load IEEE68_dyn

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         Defining the constants
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        PMU_PER_AREA=3;
        SIZE_OF_voltage_data=1000;
        NO_OF_EIGEN_VALUES=20;

        PRECISION=100;
        TRUE=1;
        FALSE=0;
        area_line_no=[53,60,61;
                    30,48,63;
                    62,66,67;
                    64,65,68;
                    54,56,57];

        Y=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);
        Y_org=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);
        deg=2*NO_OF_EIGEN_VALUES;
        m=60; % height of the Hankel matrix
        H_pmu=zeros(NO_AREA,PMU_PER_AREA,m,deg);
        C_pmu=zeros(NO_AREA,PMU_PER_AREA,m);

        H_pdc=zeros(NO_AREA,PMU_PER_AREA*m,deg);
        C_pdc=zeros(NO_AREA,PMU_PER_AREA*m);
        display('Gaussian Noise no 3 phase fault');

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
                Y_org(row,col,:)=importdata(file_name,'\n');
            end %end of for col
        end %end of for row
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
        end %end of for row

        
        N=20000;
        filename=sprintf('Model_7_constant_attacker_noise_%s per_attack_%sper_ATTACK_REGION_%s_N_%d.txt',num2str(SCALE*100),num2str(Attack*100),num2str(ATTACK_REGION),N);
        disp(filename);
        fid = fopen(filename,'w');
        fprintf(fid,sprintf(' Always area %s under attack \n Default number of iterations = %d',num2str(ATTACK_REGION),N));
        for attack_mode=-1:1
            % -1 No attack + No outlier elimination  ,
            % 0 = attack + No outlier elimination ,
            % 1 = attack + Outlier elimination
           w=zeros(deg,1); %Everyone has equal initial w_ij
        w12_1(:,1)=w; w12_2(:,1)=w; w15_1(:,1)=w;  w15_5(:,1)=w;
        w23_2(:,1)=w; w23_3(:,1)=w; w34_4(:,1)=w;  w34_3(:,1)=w;
        w45_4(:,1)=w; w45_5(:,1)=w;

        x(:,1)=ones(NO_AREA*deg,1);
        convergence_flag=FALSE;
         for j=1:N
               


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         if we have not reached convergence
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp(sprintf(' -------------- %d --------------',j));

                outliers=zeros(NO_AREA,2);
                for row=1:NO_AREA

                    [left,right]=cal_left_right(row);
                    disp(sprintf(' area = %d left =%d right =%d',row,left,right));
                    if(attack_mode<1)
                        left_val= x((left-1)*deg+1:left*deg,j);
                        right_val= x((right-1)*deg+1:right*deg,j);
                    else
                        left_val= eliminate_outlier_individual(x,j,NO_AREA,deg,left);
                        right_val= eliminate_outlier_individual(x,j,NO_AREA,deg,right);
                     end

                    Pre_succ=2;
                    dummy=reshape(H_pdc(row,1:TT,:),TT,deg);
                    dummy_0=dummy.';
                    dummy_01=dummy_0*dummy;
                    dummy_1=dummy_01+rho*Pre_succ*eye(size(dummy,2));
                    dummy_3=reshape(C_pdc(row,1:TT,:),TT,1);
                    dummy_4=dummy'*dummy_3;
                    disp(size(left_val));
                    disp(size(right_val));
                    if(row==1)
                        dummy_5=dummy_4+(w12_1(:,j) + w15_1(:,j))+rho*(left_val+right_val);
                    elseif(row==2)
                        dummy_5=dummy_4+(w23_2(:,j) - w12_2(:,j))+rho*(left_val+right_val);
                    elseif(row==3)
                        dummy_5=dummy_4+(w34_3(:,j) - w23_3(:,j))+rho*(left_val+right_val);
                    elseif(row==4)
                        dummy_5=dummy_4+(w45_4(:,j) - w34_4(:,j))+rho*(left_val+right_val);
                    elseif ( row==5)
                        dummy_5=dummy_4+((-1)*w15_5(:,j) - w45_5(:,j))+rho*(left_val+right_val);
                    end %end of if(row==1)

                    if(ATTACK_REGION==row && attack_mode>-1)
                        disp(sprintf(' I am here ATTACK_REGION =%d',ATTACK_REGION));
                        %disp(size(dummy_5));
                        %disp(size(dummy_1));
                        temp=dummy_1\dummy_5;
                        for k=1:deg
                            x((row-1)*deg+k,j+1)=(1+33*Attack)*temp(k);
                        end
                        disp(x((row-1)*deg+1:row*deg,j+1));
                    else 
                        x((row-1)*deg+1:row*deg,j+1) =dummy_1\dummy_5;
                    end

                    if(row==2)
                        disp(x(deg+1:2*deg,j+1));
                        w12_2(:,j+1)=w12_2(:,j)-rho*(left_val-x(deg+1:2*deg,j+1));   % w12 update as 2 is the successor of 1
                    elseif(row==3)
                        w23_3(:,j+1)=w23_3(:,j)-rho*(left_val-x(2*deg+1:3*deg,j+1));   % w23 update as 2 is the predecessor of 3
                    elseif(row==4)
                        w34_4(:,j+1)=w34_4(:,j)-rho*(left_val-x(3*deg+1:4*deg,j+1));
                    elseif ( row==5)
                        w15_5(:,j+1)=w15_5(:,j)-rho*(left_val-x(4*deg+1:5*deg,j+1));  % w14 update as 1 is the predecessor of 4
                        w45_5(:,j+1)=w45_5(:,j)-rho*(right_val-x(4*deg+1:5*deg,j+1));   % w14 update as 1 is the predecessor of 4
                    end %end of if
                end %end of for row

                convergence_flag=TRUE;
                for row=1:NO_AREA*deg
                    if ((abs(x(row,j+1)-x(row,j)))>0.00001)
                        convergence_flag=FALSE;
                        %if(j>N-3) disp(sprintf('abs(x(%d,%d)-x(%d,%d) = %d - %d =%d  ',row,j+1,row,j,x(row,j+1),x(row,j),abs(x(row,j+1)-x(row,j)))); end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         Check if we have reached convergence
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if( convergence_flag==TRUE)
                    display(sprintf('Convergence Reached at %d ',j));
                    break;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         after receiving from all successor
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for row=1:NO_AREA

                    [left,right]=cal_left_right(row);
                    if(attack_mode<1)
                        left_val= x((left-1)*deg+1:left*deg,j+1);
                        right_val= x((right-1)*deg+1:right*deg,j+1);
                    else
                        left_val= eliminate_outlier_individual(x,j+1,NO_AREA,deg,left);
                        right_val= eliminate_outlier_individual(x,j+1,NO_AREA,deg,right);
                    end
                    if(row==1)
                        w15_1(:,j+1)=w15_1(:,j)-rho*(x(1:deg,j+1)-left_val);
                        w12_1(:,j+1)=w12_1(:,j)-rho*(x(1:deg,j+1)-right_val);
                    elseif(row==2)
                        w23_2(:,j+1)=w23_2(:,j)-rho*(x(deg+1:2*deg,j+1)-right_val);
                    elseif(row==3)
                        w34_3(:,j+1)=w34_3(:,j)-rho*(x(2*deg+1:3*deg,j+1)-right_val);
                    elseif(row==4)
                        w45_4(:,j+1)=w45_4(:,j)-rho*(x(3*deg+1:4*deg,j+1)-right_val);
                    end  %  end if(row==1)
                end%for row after receiving from all successor
            end%for j

            if(attack_mode==-1)
                no_attack_no_elimination=x(:,j+1);
            elseif(attack_mode==0)
                attack_no_elimination=x(:,j+1);
            else
                attack_elimination=x(:,j+1);
            end
        end %attack_mode
disp(' \n \n\n ');

        for row=1:NO_AREA*deg

            if(mod(row,40)==1)
                fprintf(fid,sprintf(' \n \n %%%%%%%%%%%%%% AREA %d values %%%%%%%%%%%%%%%%',(mod(row,40)+1)));
                fprintf(fid,'\n Row    NoAttack_Noelimination Attack_NoElimination Attack_Elimination');
            end
               fprintf(fid,'\n%3d)  %10d   %20d %20d', row, no_attack_no_elimination(row),attack_no_elimination(row),attack_elimination(row));
               %fprintf(fid,'\n%3d)   %20d', row, attack_no_elimination(row));
        end
% 

%         disp(' \n \n\n ');
% 
%         for row=1:NO_AREA*deg
% 
%             if(mod(row,40)==1)
%                 fprintf(fid,sprintf(' \n \n %%%%%%%%%%%%%% AREA %d values %%%%%%%%%%%%%%%%',(mod(row,40)+1)));
%                 fprintf(fid,'\n Row    NoAttack_Noelimination Attack_NoElimination Attack_Elimination');
%             end
%                 fprintf(fid,'\n%3d)  %10d   %20d %20d', row, no_attack_no_elimination(row),attack_no_elimination(row),attack_elimination(row));
%         end
% 
%         %RMS Error
%         fprintf(fid,sprintf('\n\n Reporting RMS Error '));
%         for attack_mode=-1:1
%             if (attack_mode==-1)
%                 Vec=no_attack_no_elimination;
%                 fprintf(fid,sprintf('\n\n ## No attack No elimination scenario'));
%             elseif(attack_mode==0)
%                 Vec=attack_no_elimination;
%                 fprintf(fid,sprintf('\n\n ## Attack  but No elimination scenario'));
%             elseif(attack_mode==1)
%                 Vec=attack_elimination;
%                 fprintf(fid,sprintf('\n\n ## Attack and  elimination scenario'));
%             end
% 
%             V_diff=zeros(deg);
%             for j=1:deg
%                 V_max=intmin('int64');
%                 V_min=intmax('int64');
%                 for i=1:NO_AREA
%                     if (Vec((i-1)*deg+j)>V_max)
%                         V_max=Vec((i-1)*deg+j);
%                     end
%                     if (Vec((i-1)*deg+j)<V_min)
%                         V_min=Vec((i-1)*deg+j);
%                     end
%                 end
%                 V_diff(j)=(V_max-V_min); % Error
%             end %for j=1:deg
%             RMSE= rms(V_diff);
%             fprintf(fid,sprintf('\n rms = %d ',RMSE(1)));
%         end    %attack_mode
% 
%         %ACTUAL RMS ERROR
%         V_diff=zeros(NO_AREA*deg);
%         for j=1:NO_AREA*deg
%             V_diff(j)=no_attack_no_elimination(j)-attack_no_elimination(j); % Error
%         end
%         RMSE= rms(V_diff);
%         fprintf(fid,sprintf('\n Attack but no elimination scenario  = %d ',RMSE(1)));
% 
% 
%         V_diff=zeros(NO_AREA*deg);
%         for j=1:NO_AREA*deg
%             V_diff(j)=no_attack_no_elimination(j)-attack_elimination(j); % Error
%         end
%         RMSE= rms(V_diff);
%         fprintf(fid,sprintf('\n Attack and  elimination scenario = %d ',RMSE(1)));
% 
%         fclose(fid);

    end % Attack amount
end %Attack Region
