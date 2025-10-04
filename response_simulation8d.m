function dof8_resp_nonoise = response_simulation8d(param, DOF)
num_rel=1; % Set number of realization
deltat = 0.01; %seconds, also sampling period
inputduration = 1e5 ; % number of data point, in this case 1000 seconds

%build mass, stiffnes and damping matrices
n = 8; % to test eigenvalue
%number of DOF = n
ninput = 8;
noutput = n;
Vecsensor = [1,2,3,4,5,6,7,8];

m0 = 1;
m_values=[m0,m0,m0,m0,m0,m0,m0,m0]; % mass
k_list=param(1:DOF); % stiffness
xivalue = param(DOF+1:2*DOF); % damping ratios

% set percentage noise
percRMS = 0.1;

loc_history=zeros(num_rel,1);
dof8_resp_nonoise=zeros(inputduration,n,num_rel);
dof8_resp_10per=zeros(inputduration,n,num_rel);

loc=1:8;

actuatorlocation = [1,2,3,4,5,6,7,8];
vecsen = zeros(1,n);
vecact = zeros(1,n);
for i=1:n;
    if Vecsensor(i)> 0;
        vecsen(1,i) = 1;
    else
        vecsen(1,i) = 0;
    end
end
vecact(1, actuatorlocation) = 1;
forcemagnitude = 1;

multiplot = 1;
pigreco = 4.0 * atan(1);
Mass = zeros(n,n);
Stiffness = zeros(n,n);
Damping = zeros(n,n);
 for i=1:n-1;
     Mass(i,i) = m_values(i);
 end
%Mass(n,n) = mvalue/2;
Mass(n,n) = m_values(n);
 for i=1:n;
     if i==n
         Stiffness(i,i)=k_list(i);
     else
         Stiffness(i,i) = k_list(i)+k_list(i+1);
         Stiffness(i, i+1) = - k_list(i+1);
         Stiffness(i+1, i) = - k_list(i+1);
     end
 end

[Modesh, Fre2] = eig(Stiffness,Mass);
Fre = sqrtm(Fre2);
% Modaldampingmatrix = 2 * xivalue * Fre;
Modaldampingmatrix = 2 * (xivalue .* eye(n)) .* Fre;
Damping = inv(Modesh') * Modaldampingmatrix * inv(Modesh);
%
%   build the force location matrix
%
Flocat = zeros(n,ninput);
Flocat(loc,ninput)=1; % forces applied at either 1st or 4th mass
%
%   Building the sensor location matrix
%
Senslocat=zeros(noutput,n);
for i=1:noutput;
    Senslocat(i,Vecsensor(i)) = 1;
end
%disp(' Building the system matrices: End') 
%disp(' Transforming the system in discrete time: Start') 
%
%   Transform 2nd order model to 1st order model
%
Acont = zeros(2*n,2*n);
Bcont = zeros(2*n,ninput);
Ccont = zeros(noutput,2*n);
Ccontwork = zeros(n,2*n);
Dcontwork = zeros(n,ninput);
Dcont = zeros(noutput,ninput);% sdof is 1x2
MinvK = zeros(n,n);
MinvC = zeros(n,n);
MinvF = zeros(n,ninput);
MinvK = inv(Mass) * Stiffness;
MinvC = inv(Mass) * Damping;
MinvF = inv(Mass) * Flocat;
for i=1:n
    Acont(i,i+n)=1;
end
Acont(n+1:2*n,1:n)= - MinvK;
Acont(n+1:2*n,n+1:2*n)= -MinvC;
% Bcont(n+1:2*n,ninput) = MinvF;
Bcont(n+1:2*n,:) = MinvF;
Ccontwork(1:n,1:n)= - MinvK;
Ccontwork(1:n,n+1:2*n) = - MinvC;
Ccont = Senslocat * Ccontwork;
Dcontwork = MinvF;
Dcont = Senslocat * Dcontwork;
[Adisc,Bdisc] = c2d(Acont,Bcont,deltat);
Cdisc = Ccont;
Ddisc = Dcont;
Observcont = obsv(Acont, Ccont);
Controcont = ctrb(Acont, Bcont);
rankObservcont = rank(Observcont);
rankControcont = rank(Controcont);
Observdisc = obsv(Adisc, Cdisc);
Controdisc = ctrb(Adisc, Bdisc);
rankObservdisc = rank(Observdisc);
rankControdisc = rank(Controdisc);
%disp(' Transforming the system in discrete time: End') 
%
% Generate input force
%
%disp(' Generating Input Force: Start') 
%randgenvect = wgn(inputduration,1,1,1,0);
randgenvect1 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect2 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect3 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect4 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect5 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect6 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect7 = idinput(inputduration,'RGS',[0 1],[-1  1]);
randgenvect8 = idinput(inputduration,'RGS',[0 1],[-1  1]);

inputforce1 = forcemagnitude * randgenvect1;
inputforce2 = forcemagnitude * randgenvect2;
inputforce3 = forcemagnitude * randgenvect3;
inputforce4 = forcemagnitude * randgenvect4;
inputforce5 = forcemagnitude * randgenvect5;
inputforce6 = forcemagnitude * randgenvect6;
inputforce7 = forcemagnitude * randgenvect7;
inputforce8 = forcemagnitude * randgenvect8;
inputforce = [inputforce1, inputforce2, inputforce3, inputforce4,...
    inputforce5, inputforce6, inputforce7, inputforce8];
%disp(' Generating Input Force: End') 
%
%   Simulation of system response
%
%disp(' Simulating System Response: Start') 
respoutput=dlsim(Adisc, Bdisc, Cdisc, Ddisc, inputforce);
%disp(' Simulating System Response: End') 
dof8_resp_nonoise(:,:,1)=respoutput;
%
%   Plot input and output

%
%   Modify response and input to account for noise
%
%disp(' Adding Measurement Noise: Start') 
noisesequence = idinput(inputduration,'RGS',[0 1],[-1 1]);
%noisesequence = wgn(inputduration,1,1,1,0);
inputnoise = percRMS * noisesequence/std(noisesequence)* std(inputforce);
inputforcenoise = inputforce + inputnoise;
respoutputno=dlsim(Adisc, Bdisc, Cdisc, Ddisc, inputforcenoise);
for i=1:noutput;
%    noiseoutput(:,i) = wgn(inputduration,1,1,1,i);
    noiseoutput(:,i) = idinput(inputduration,'RGS',[0 1],[-1 1]);
    responsenoise(:,i) = percRMS * noiseoutput(:,i)/std(noiseoutput(:,i)) * std(respoutput(:,i));
end
respoutputnoise = respoutputno + responsenoise; % input noise + measure noise
%disp(' Adding Measurement Noise: End')

dof8_resp_10per(:,:,1)=respoutputnoise;

end