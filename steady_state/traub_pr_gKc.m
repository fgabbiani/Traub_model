function traub_pr_gKc

%601 x 2 matrix contains cai, g_Kc
load gKc_conductance.txt

%set v to a value where the activation c = 1
v = 60;
c = 0:0.5:300;

%assume gmax is equal to 1
g = min(c/250,1);

figure;
plot(c,g);
hold on;
plot(gKc_conductance(:,1),gKc_conductance(:,2),'r--');
ylabel('g');
xlabel('calcium concentration (arbitrary)');
title('gKc conductance');
