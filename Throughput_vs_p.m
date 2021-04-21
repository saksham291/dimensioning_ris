%SimRIS_v18(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz)
clearvars

tic
Environment=1;                     % 1 Indoor (InH - Indoor Office) / 2 Outdoor (UMi - Street Canyon)
Scenario=1;                        % 1 (RIS in xz plane - left side wall) or 2 (RIS in yz plane - opposite wall)
Frequency=28;                      % Operating Frequency (in GHz)
N=256;                              % Number of RIS  Elements (sqrt(N) should be an integer)
Nt=1;                              % should be an integer power of 2'
Nr=1;                              % should be an integer power of 2

ArrayType=2;                       % Uniform Linear Array=1 or Uniform Planar Array=2
Tx_xyz=[1,1,1];                  % Tx coordinates (in meters)
Rx_xyz=[51,1,1];                % Rx coordinates (in meters)
%RIS_xyz=[5,2,2];       % RIS coordinates (in meters)
Nsym = 30000;
shifts = 49;
hei = 3;         %keep hei=2 for finding spread
I = eye(Nr);
H = zeros(N,Nt,Nsym);
G = zeros(Nr,N,Nsym);
D = zeros(Nr,Nt,Nsym);

H_LOS = zeros(N,Nt,Nsym);
G_LOS = zeros(Nr,N,Nsym);
D_LOS = zeros(Nr,Nt,Nsym);

R_ris = zeros(Nsym,shifts);
R_no_ris = zeros(Nsym,shifts);
M_ris = zeros(shifts,hei);
M_no_ris = zeros(shifts,hei);

R_ris_LOS = zeros(Nsym,shifts);
R_no_ris_LOS = zeros(Nsym,shifts);
M_ris_LOS = zeros(shifts,hei);
M_no_ris_LOS = zeros(shifts,hei);

%comp_phi = zeros(N,N,Nsym);
%theta = complex(comp_phi,0);

for y = 1:hei
    for x = 1:shifts
        RIS_xyz=[x+1,y+1,1];                
        parfor iter = 1:Nsym  
                 %INDOORS
            [H(:,:,iter),G1,D(:,:,iter)] = SimRIS_v18_1(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz);
            G(:,:,iter) = transpose(G1);
            %phi_N = angle(H(:,:,iter));
            %theta(:,:,iter) = diag(exp(-i*phi_N));
            Ch = G(:,:,iter)*H(:,:,iter) + D(:,:,iter); 
            R_ris(iter,x) = log2(I+det((Ch'*Ch)*1e13))             %1e13 to be replaced with Pt/Pn
            R_no_ris(iter,x) = log2(I+det((D(:,:,iter)'*D(:,:,iter))*1e13))

        end
        M_ris(x,y) = mean(abs(R_ris(:,x)))                               %Calculation of mean of x values
        M_no_ris(x,y) = mean(abs(R_no_ris(:,x)))
        if Environment==2  
            M_ris_LOS(x,y) = mean(abs(R_ris_LOS(:,x)))                               %Calculation of mean of x values
            M_no_ris_LOS(x,y) = mean(abs(R_no_ris_LOS(:,x)))
        end
    end
end

%Save variables to file
%save('sim_data.mat','Environment','Scenario','Frequency','ArrayType','N','Nt','Nr','Tx_xyz','Rx_xyz','RIS_xyz','shifts', 'M_ris','M_no_ris','D','H','G','R_ris','R_no_ris')
toc
x = linspace(1,shifts,shifts)  
%p = linspace(2,22,6)
Tadd = M_ris-M_no_ris
% for j = 1:hei
%     M(j) = sum(Tadd(:,j))
% end
%plot(p,M)
plot(x,M_ris(:,1))
hold on 
plot(x,M_ris(:,2))
plot(x,M_ris(:,3))
plot(x,M_no_ris(:,1))
hold off
%xlim([1 26])
legend({'RIS, z=1 m','RIS, z=2 m','RIS, z=3 m','No RIS'},'Location','north','FontSize',13)
xlabel('y (m)','FontSize',12) 
ylabel('Achievable Rate [bit/s/Hz]','FontSize',12) 
grid on
print -depsc epsFig2
save('Data_Throughput_p.mat')
 


% Coordinates (with respect to reference (0,0,0) point (Scenario 1 - Indoors)
%Tx_xyz=[0,25,2];       % Tx coordinates (in meters)
%Rx_xyz=[38,48,1];     % Rx coordinates (in meters)
%RIS_xyz=[40,50,2];       % RIS coordinates (in meters)

% for test (Scenario 1 - Indoors)
%Tx_xyz=[0,25,3];       % Tx coordinates (in meters)
%Rx_xyz=[70,45,0.75];     % Rx coordinates (in meters)
%RIS_xyz=[65,50,0.5];       % RIS coordinates (in meters)

% for test (Scenario 2 - Indoors)
%Tx_xyz=[0,25,2];       % Tx coordinates (in meters)
%Rx_xyz=[70,30,1];     % Rx coordinates (in meters)
%RIS_xyz=[75,25,2];       % RIS coordinates (in meters)

% for test (Scenario 1 - Outdoors)
%Tx_xyz=[0,25,30];       % Tx coordinates (in meters)
%Rx_xyz=[140,120,1];     % Rx coordinates (in meters)
%RIS_xyz=[120,150,10];       % RIS coordinates (in meters)

% % for test (Scenario 2 - Outdoors)
%Tx_xyz=[0,25,30];       % Tx coordinates (in meters)
%Rx_xyz=[150,100,1];     % Rx coordinates (in meters)
%RIS_xyz=[180,40,10];       % RIS coordinates (in meters)
