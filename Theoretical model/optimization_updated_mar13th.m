function [min_JRFr, angle_elbow, angle_wrist, angle_elbow_diff, angle_wrist_diff] = optimization_updated_mar13th()

[elbow, wrist] =  meshgrid(0:180/399:180,0:180/799:90);
JRFr = zeros(length(elbow(:,1)), length(wrist(:,1)));
% JRF_diff = zeros(length(elbow(:,1)), length(wrist(:,1)));
JRFr_adjusted = zeros(length(elbow(:,1)), length(wrist(:,1)));
JRFx = zeros(length(elbow(:,1)), length(wrist(:,1)));
JRFy = zeros(length(elbow(:,1)), length(wrist(:,1)));
moment = zeros(length(elbow(:,1)), length(wrist(:,1)));
force_biceps = zeros(length(elbow(:,1)), length(wrist(:,1)));
force_brachialis = zeros(length(elbow(:,1)), length(wrist(:,1)));
force_brachioradialis = zeros(length(elbow(:,1)), length(wrist(:,1)));
force_ECRL = zeros(length(elbow(:,1)), length(wrist(:,1)));

for i = 1:length(elbow)
    for j = 1:length(wrist)
%         [~, ~, ~, JRFr(i,j), ~] = force_balance_analysis_updated_mar13th(elbow(1,j),wrist(i,1), 0);
        [force, JRFx(i,j), JRFy(i,j) JRFr_adjusted(i,j) moment(i,j)] = force_balance_analysis_updated_mar13th(elbow(1,j),wrist(i,1), 1);
%         force_biceps(i,j) = force(1);
%         force_brachialis(i,j) = force(2);
%         force_brachioradialis(i,j) = force(3);
%         force_ECRL(i,j) = force(4);
    end
end

elbow_index = 1:0.5:180;
for i = 1:length(elbow_index)
    [force_2d, ~, ~, ~, ~] = force_balance_analysis_updated_mar13th(elbow_index(i),0, 1);
    force_biceps_2d(i) = force_2d(1);
    force_brachialis_2d(i) = force_2d(2);
    force_brachioradialis_2d(i) = force_2d(3);
    force_ECRL_2d(i) = force_2d(4);
end

%2-D graphs muscle generating capacity changes on biceps
figure
plot(elbow_index, force_biceps_2d,elbow_index, force_brachialis_2d,elbow_index, force_brachioradialis_2d,elbow_index, force_ECRL_2d)
xlabel('Elbow angle (degree)')
ylabel('Muscle generating capacity(N)')
legend('Biceps', 'Brachialis', 'Brachioradialis', 'ECRL')


%plot muscle generating capacity of each muscle at diff. elbow and wrist angles
figure
mesh(elbow, wrist, force_biceps);
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Bicep force (N)')

figure
mesh(elbow, wrist, force_brachialis);
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Brachialis force (N)')

figure
mesh(elbow, wrist, force_brachioradialis);
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Brachioradialis force (N)')

figure
mesh(elbow, wrist, force_ECRL);
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('ECRL force (N)')


%Plot JRFr before and after force adjustment due to change in muscle length
figure
hold on
%create color map
C1 =zeros(length(elbow(:,1)), length(wrist(:,1)), 3);
C1(:,:,1) = 0.0;
C1(:,:,2) = 0.5;
C1(:,:,3) = 1;
mesh(elbow, wrist, JRFr, C1);

C2 =zeros(length(elbow(:,1)), length(wrist(:,1)), 3);
C2(:,:,1) = 1;
C2(:,:,2) = 0.5;
C2(:,:,3) = 0;
mesh(elbow, wrist, JRFr_adjusted, C2);
hold off
title('Graph of elbow joint reaction forces at various angle positions')
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Elbow Joint Reaction Force (JRF)')
legend('Before muscle force adjustment', 'After muscle force adjustment');

%Plot fo JRFr
figure
mesh(elbow, wrist, JRFr_adjusted);
title('Graph of elbow joint reaction forces at various angle positions')
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Elbow Joint Reaction Force (JRF)')
text([pi pi pi/2], [0 pi/2 0], {'(180, 0, 299.7)', '(180, 90, 620.5)', '(90, 0, 308.1)'})

%plot JRFx and JRFy
% figure
% mesh(elbow, wrist, JRFx, gradient(JRFx));
% title('Graph of elbow JRF in axial direction')
% xlabel('Elbow angle (degree)')
% ylabel('Wrist angle (degree)')
% zlabel('Elbow JRFx')
% 
% figure
% mesh(elbow, wrist, JRFy, gradient(JRFy));
% title('Graph of elbow JRF in lateral direction')
% xlabel('Elbow angle (degree)')
% ylabel('Wrist angle (degree)')
% zlabel('Elbow JRFy')

%plot JRFx with a z=0 plane
figure
hold on
mesh(elbow, wrist, JRFx, gradient(JRFx));
mesh(elbow, wrist, zeros(length(elbow(:,1)), length(wrist(:,1))), 'CDataMapping','direct');
hold off
title('Graph of elbow JRF in axial direction')
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Elbow JRFx')

%plot JRFy with a z=0 plane
figure
hold on
mesh(elbow, wrist, JRFy, gradient(JRFy));
mesh(elbow, wrist, zeros(length(elbow(:,1)), length(wrist(:,1))), 'CDataMapping','direct');
hold off
title('Graph of elbow JRF in normal direction')
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Elbow JRFy')

%plot moment with a z=0 plane
figure
hold on
mesh(elbow, wrist, zeros(length(elbow(:,1)), length(wrist(:,1))));
mesh(elbow, wrist, moment, gradient(moment));
hold off
title('Graph of elbow moment')
xlabel('Elbow angle (degree)')
ylabel('Wrist angle (degree)')
zlabel('Elbow moment (CCW+)')


% min_JRFr = min(JRFr(:)); 
min_JRFr_adjusted = min(JRFr_adjusted(:)); 

% [pos_y pos_x] = find(JRFr == min_JRFr); 
[pos_y pos_x] = find(JRFr_adjusted == min_JRFr_adjusted); 

angle_wrist = wrist(pos_y,1); 
angle_elbow = elbow(1,pos_x);

% min_JRF_diff = min(JRF_diff(:)); 
% [pos_y_diff pos_x_diff] = find(JRF_diff == min_JRF_diff); 
% angle_wrist_diff = wrist(pos_y_diff,1); 
% angle_elbow_diff = elbow(1,pos_x_diff);

% disp(min_JRFr)
disp(min_JRFr_adjusted)
end
