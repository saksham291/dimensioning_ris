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
% Rx_xyz=[51,1,1];                % Rx coordinates (in meters)
%RIS_xyz=[5,2,2];       % RIS coordinates (in meters)
Nsym = 40000;
shifts = 49;
hei = 2;         %keep hei=2 for finding spread
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
    rx_pos = 50/y+1;
    Rx_xyz=[rx_pos,1,1];
    for x = 1:shifts
        RIS_xyz=[((x+1)/y),4-y,1];                
        parfor iter = 1:Nsym  
              %INDOORS
            [H(:,:,iter),G1,D(:,:,iter)] = SimRIS_v18_1(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz);
            G(:,:,iter) = transpose(G1);
            %phi_N = angle(H(:,:,iter));
            %theta(:,:,iter) = diag(exp(-i*phi_N));
            Ch = G(:,:,iter)*H(:,:,iter) + D(:,:,iter); 
            R_ris(iter,x) = det((Ch'*Ch)*1e13)              %1e13 to be replaced with Pt/Pn
            R_no_ris(iter,x) = det((D(:,:,iter)'*D(:,:,iter))*1e13)
        end
        M_ris(x,y) = mean(abs(R_ris(:,x)))                               %Calculation of mean of x values
        M_no_ris(x,y) = mean(abs(R_no_ris(:,x)))
    end
end


%Save variables to file
%save('sim_data.mat','Environment','Scenario','Frequency','ArrayType','N','Nt','Nr','Tx_xyz','Rx_xyz','RIS_xyz','shifts', 'M_ris','M_no_ris','D','H','G','R_ris','R_no_ris')
toc

x = linspace(1,25,49);
y = linspace(1,49,49);
Padd_1 = M_ris(:,1)-M_no_ris(:,1)
Padd_2 = M_ris(:,2)-M_no_ris(:,2)
fac = max(Padd_2)/max(Padd_1)
plot(x./25,Padd_1)
hold on 
plot(y./50,Padd_2./fac)
hold off
legend({'RIS, r = 25m','RIS, r = 50m'},'Location','north','FontSize',13)
grid on
xlabel('p','FontSize',13) 
ylabel('Additional received power (P_{RIS})','FontSize',13)
print -depsc Additional_pow_equivalence
matlab2tikz('equivalence.tex');




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
%M_RIS = [7.46356127776573;7.27778319830253;9.00000914818980;8.99757768606307;8.68444880578030;8.52454951937622;8.23235147402740;8.02895209329259;7.96896759919037;7.87332452427982;7.67482650686835;7.55299517167672;7.56605730975148;7.56497768261837;7.60757442521537;7.53925720186601;7.54797987875802;7.24855877752752;7.29065506283447;7.39363484204512;7.29778566289320;7.22971846178936;7.39729923196148;7.38306227651305;7.26034078038203;7.36977870631079;7.29123581236765;7.24837886666015;7.32141232399802;7.38964756415677;7.46156174685431;7.48987583774850;7.25875365307961;7.37915454166918;7.55261191524338;7.52840830363829;7.61161239291015;7.57542404528963;7.82102856356991;7.89947396324445;7.96687297805145;8.05049818667457;8.25456142948732;8.45860371440874;8.66099316259031;8.94267989168175;8.77701419966748;7.17488595499654;7.11509224078134]
%M_NO_RIS = [6.73237176115224;6.56182846145994;6.52163957767655;6.56682719410527;6.40602118532028;6.54987009606421;6.46841098112398;6.45283654099653;6.57535224597201;6.61365037432031;6.47207919093562;6.39802397701397;6.42243831632479;6.35753885895349;6.54807936642373;6.49619804091775;6.56115757383270;6.30314811881750;6.42340800074278;6.42253891276012;6.30949390726293;6.38536224141056;6.51083448547461;6.50739256920484;6.48801548944346;6.35748311934655;6.33946329926568;6.35570951058994;6.43689593900774;6.48491541832041;6.51528088771299;6.45429368972614;6.34330075704473;6.33698650594920;6.51059227599846;6.46115716887643;6.31249801622814;6.48482869087168;6.49849360688567;6.61761098666348;6.42804036791032;6.48249124921072;6.29814906867709;6.48754139915296;6.38284012234981;6.33380167937887;6.44901345239657;6.49127042439024;6.49507741318783]
%M_RIS_2 = [10.5804344949954;10.3276006648989;11.5643517912155;11.6139770751351;11.4019993006667;11.2004660194734;11.1713526926444;11.2016717869078;10.8848645168174;10.9624362562303;10.9065812819977;10.7076452087282;10.8389268839635;10.9127084900573;11.0392774525815;10.8510137900620;11.0947543179534;11.0797261023563;11.3197770833119;11.2964041303007;11.5475854702915;11.4578178639744;10.2626005029469;10.2636080512668]
%M_NO_RIS_2 = [10.1902082956104;10.0026020647632;10.1182323596494;10.1705692585415;10.0780534440315;10.1278152767466;10.1111392996746;10.0775836253963;9.90113683577106;9.85959878444276;9.95100126452210;9.87074771751728;9.98063283216849;10.0985824421786;10.0630943341440;9.96077264918525;10.0138221860139;9.89452439069472;10.0629869253358;9.98987533593480;10.0473887847911;9.90333868463048;9.87767471481777;9.89746562196095]

