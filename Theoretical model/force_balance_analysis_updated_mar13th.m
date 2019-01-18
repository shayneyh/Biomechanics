%Model for joint reaction forces when carrying a load
%Assumptions: 
%1. PCSA remain constant when muscle exerts forces
%2. f_PCSA max is 20 for all muscle fibres
%3. f_PCSA is only dependent on length
%4. Palm stays flat
%5. Wrist range from 0-90 while elbow ranges from 0-180

%function [JRF_diff, JRFr] = force_balance_analysis_updated_mar13th(angle_elbow, angle_wrist, adjust_f_pcsa)

function [force JRFx JRFy JRFr moment_elbow] = force_balance_analysis_updated_mar13th(angle_elbow_deg, angle_wrist_deg, adjust_f_pcsa)

% angle conversion
angle_elbow = angle_elbow_deg*pi/180;
angle_wrist = angle_wrist_deg*pi/180;

%define constants (cm)
L_H_Biceps = 31;
L_H_Brachialis = 10;
L_H_Brachioradialis = 8;
L_H_Extensor = 3;
L_F_Biceps = 8;
L_F_Brachialis = 5;
L_F_Brachioradialis = 24;
L_F_Extensor = 25;
%maximum force per PCSA
Max_Force = 20; %N/cm^2
%resting lengths
L_R_Biceps = 31.3; % (avg of 27.7 and 34.9cm)
%http://journals1.scholarsportal.info.myaccess.library.utoronto.ca/pdf/03635023/v36i0005/881_sauama.xml   
L_R_Brachialis = 15; % +-0.8cm http://journals1.scholarsportal.info.myaccess.library.utoronto.ca/pdf/03635023/v36i0005/881_sauama.xml
L_R_Brachioradialis = 19.9;
L_R_Extensor = 15.5; %http://journals1.scholarsportal.info.myaccess.library.utoronto.ca/pdf/00219290/v28i0007/791_tbpehwms.xml

%Physiological cross sectional area (cm^2)
PCSA_Biceps = 12.3;
PCSA_Brachialis = 13.0;
PCSA_Brachioradialis = 2.9;
PCSA_Extensor = 3.6;

%weight (assume 10kg)
weight = 27.44;

L_Forearm = 25;
%length from wrist to where weight is loaded
L_Hand = 6;

%muscle lengths (resting length, Max F produced)
l_m_biceps = sqrt(L_H_Biceps^2 + L_F_Biceps^2 - 2*L_H_Biceps*L_F_Biceps*cos(angle_elbow));
l_m_brachialis = sqrt(L_H_Brachialis^2 + L_F_Brachialis^2 - 2*L_H_Brachialis*L_F_Brachialis*cos(angle_elbow));
l_m_brachioradialis = sqrt(L_H_Brachioradialis^2 + L_F_Brachioradialis^2 - 2*L_H_Brachioradialis*L_F_Brachioradialis*cos(angle_elbow));
l_m_extensor = sqrt(L_H_Extensor^2 + L_F_Extensor^2 - 2*L_H_Extensor*L_F_Extensor*cos(angle_elbow));

%angle between muscle and forearm: 1 biceps 2 brachialis 3 brachioradialis
%4extensor
angle_m_f = zeros(1,4);
angle_m_f(1) = acos((L_F_Biceps^2 + l_m_biceps^2 - L_H_Biceps^2)/(2*L_F_Biceps*l_m_biceps));
angle_m_f(2) = acos((L_F_Brachialis^2 + l_m_brachialis^2 - L_H_Brachialis^2)/(2*L_F_Brachialis*l_m_brachialis));
angle_m_f(3) = acos((L_F_Brachioradialis^2 + l_m_brachioradialis^2 - L_H_Brachioradialis^2)/(2*L_F_Brachioradialis*l_m_brachioradialis));
angle_m_f(4) = acos((L_F_Extensor^2 + l_m_extensor^2 - L_H_Extensor^2)/(2*L_F_Extensor*l_m_extensor));


%f(L) = F_max - k(Lmax-L)^2 (k = 0.025, Fmax = Max_Force, L = l_m)
f_per_pcsa = zeros(1,4);
    
if adjust_f_pcsa == 0
f_per_pcsa(1) = Max_Force;
f_per_pcsa(2) = Max_Force;
f_per_pcsa(3) = Max_Force;
f_per_pcsa(4) = Max_Force;
else
f_per_pcsa(1) = 20 - 0.025*(L_R_Biceps - l_m_biceps)^2;
f_per_pcsa(2) = 20 - 0.025*(L_R_Brachialis - l_m_brachialis)^2;
f_per_pcsa(3) = 20 - 0.025*(L_R_Brachioradialis - l_m_brachioradialis)^2;
f_per_pcsa(4) = 20 - 0.025*(L_R_Extensor - l_m_extensor)^2;
end

disp([l_m_biceps l_m_brachialis l_m_brachioradialis l_m_extensor])
disp(angle_m_f)
disp(f_per_pcsa)

%Physiological cross sectional area
pcsa = [PCSA_Biceps PCSA_Brachialis PCSA_Brachioradialis PCSA_Extensor];

%force exerted by each muscle
force = zeros(1,4);
for i = 1:4
    force(i) = f_per_pcsa(i)*pcsa(i);
end


JRFx = 0;
JRFy = 0;

for i = 1:4
    JRFx = JRFx + f_per_pcsa(i)*pcsa(i)*cos(angle_m_f(i));
    JRFy = JRFy + f_per_pcsa(i)*pcsa(i)*sin(angle_m_f(i));
end
 
JRFx = JRFx + weight*cos(pi/2 - angle_wrist);
JRFy = JRFy - weight*sin(pi/2 - angle_wrist);

disp(JRFx)
disp(JRFy)

%resultant elbow joint reaction forces
JRFr = sqrt(JRFx^2 + JRFy^2);
% JRF_diff = abs(JRFx - JRFy);

moment_arm = zeros(1,4);
moment_arm(1) = L_F_Biceps;
moment_arm(2) = L_F_Brachialis;
moment_arm(3) = L_F_Brachioradialis;
moment_arm(4) = L_F_Extensor;


%+ve for CCW
moment_elbow = 0;
for i = 1:4
    moment_elbow = moment_elbow + f_per_pcsa(i)*pcsa(i)*moment_arm(i)/100*sin(angle_m_f(i));
end

l_weight_moment_arm = sqrt((L_Forearm/100)^2 + (L_Hand/100)^2 - 2*(L_Forearm/100)*(L_Hand/100)*cos(pi-angle_wrist));
angle_between_weight_and_forearm = acos((L_Forearm/100)*sin(pi-angle_wrist)/l_weight_moment_arm);
moment_elbow = moment_elbow - weight*l_weight_moment_arm*angle_between_weight_and_forearm;



end

