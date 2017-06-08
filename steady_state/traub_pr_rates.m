%
%  Rates of the Pinsky-Rinzel model and the Traub model since they are 
% the same. Modified from prsolve_rk4.m in
% invariance_project/pinsky_rinzel_matlab
%
% References: [PR94] Pinsky PF, Rinzel J (1994) Intrinsic and Network
% Rhythmogenesis in a Reduced Traub Model for CA3 Neurons. J Comp Neurosci
% 1:39-60. Erratum in J Comp Neurosci 2:275, 1995. 
%
% Traub RC, et al. (1991) A model of a CA3 hippocampal pyramidal neuron 
% incorporating voltage clamp data on intrinsic conductances  

function traub_pr_rates(channel)

v = -20:0.1:130;
c = 0:0.5:550;

switch(channel)
    case 'na'
        plot_na(v);
    
    case 'dr'
        plot_dr(v);
        
    case 'ka'
        plot_ka(v);

    case 'ca'
        plot_ca(v);
        
    case 'kc'
        plot_kc(v);
        
    case 'kahp'
        plot_kahp(c);

end

return

function plot_kahp(c)

%1501 x 3 matrix contains v, q_inf, tau_q
load gKahp_ss.txt

%activation variables
tau_q = 1./(aq(c) + bq(c));
q_inf = aq(c).*tau_q;

figure;
subplot(2,1,1);
plot(c, q_inf);
hold on;
plot(gKahp_ss(:,1),gKahp_ss(:,2),'r--');
ylabel('activation');
title('gKahp activation');

subplot(2,1,2);
plot(c,tau_q);
hold on;
plot(gKahp_ss(:,1),gKahp_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('calcium concentration (arbitrary)');
ylabel('time constant (ms)');

function plot_kc(v)

%1501 x 3 matrix contains v, c_inf, tau_c
load gKc_ss.txt

%activation variables
tau_c = 1./(ac(v) + bc(v));
c_inf = ac(v).*tau_c;

figure;
subplot(2,1,1);
plot(v, c_inf);
hold on;
plot(gKc_ss(:,1),gKc_ss(:,2),'r--');
ylabel('activation');
title('gKc activation');

subplot(2,1,2);
plot(v,tau_c);
hold on;
plot(gKc_ss(:,1),gKc_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

function plot_ca(v)

%1501 x 5 matrix contains v, a_inf, tau_a, b_inf, tau_b
load gCa_ss.txt

%activation variables
tau_s = 1./(as(v) + bs(v));
s_inf = as(v).*tau_s;

figure;
subplot(2,1,1);
plot(v, s_inf);
hold on;
plot(gCa_ss(:,1),gCa_ss(:,2),'r--');
ylabel('activation');
title('gCa activation');

subplot(2,1,2);
plot(v,tau_s);
hold on;
plot(gCa_ss(:,1),gCa_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

%inactivation variables
tau_r = 1./(ar(v) + br(v));
r_inf = ar(v).*tau_r;

figure;
subplot(2,1,1);
plot(v, r_inf);
hold on;
plot(gCa_ss(:,1),gCa_ss(:,4),'r--');
ylabel('inactivation');
title('gCa inactivation');

subplot(2,1,2);
plot(v,tau_r);
hold on;
plot(gCa_ss(:,1),gCa_ss(:,5),'r--');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

function plot_ka(v)

%1501 x 5 matrix contains v, a_inf, tau_a, b_inf, tau_b
load gKa_ss.txt

%activation variables
tau_a = 1./(aa(v) + ba(v));
a_inf = aa(v).*tau_a;

figure;
subplot(2,1,1);
plot(v, a_inf);
hold on;
plot(gKa_ss(:,1),gKa_ss(:,2),'r--');
ylabel('activation');
title('gKa activation');

subplot(2,1,2);
plot(v,tau_a);
hold on;
plot(gKa_ss(:,1),gKa_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

%inactivation variables
tau_b = 1./(ab(v) + bb(v));
b_inf = ab(v).*tau_b;

figure;
subplot(2,1,1);
plot(v, b_inf);
hold on;
plot(gKa_ss(:,1),gKa_ss(:,4),'r--');
ylabel('inactivation');
title('gKa inactivation');

subplot(2,1,2);
plot(v,tau_b);
hold on;
plot(gKa_ss(:,1),gKa_ss(:,5),'r--');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

function plot_dr(v)

%1501 x 3 matrix contains v, n_inf, tau_n
load gKdr_ss.txt

%activation variables
tau_n = 1./(an(v) + bn(v));
n_inf = an(v).*tau_n;

figure;
subplot(2,1,1);
plot(v, n_inf);
hold on;
plot(gKdr_ss(:,1),gKdr_ss(:,2),'r--');
ylabel('activation');
title('gKdr activation');

subplot(2,1,2);
plot(v,tau_n);
hold on;
plot(gKdr_ss(:,1),gKdr_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

function plot_na(v)

%1501 x 5 matrix contains v, m_inf, tau_m, h_inf, tau_h
load gNa_ss.txt

%activation variables
tau_m = 1./(am(v) + bm(v));
m_inf = am(v).*tau_m;

figure;
subplot(2,1,1);
plot(v, m_inf);
hold on;
plot(gNa_ss(:,1),gNa_ss(:,2),'r--');
ylabel('activation');
title('gNa activation');

subplot(2,1,2);
plot(v,tau_m);
hold on;
plot(gNa_ss(:,1),gNa_ss(:,3),'r--');
legend('PR matlab code','NEURON code');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

%inactivation variables
tau_h = 1./(ah(v) + bh(v));
h_inf = ah(v).*tau_h;

figure;
subplot(2,1,1);
plot(v, h_inf);
hold on;
plot(gNa_ss(:,1),gNa_ss(:,4),'r--');
ylabel('inactivation');
title('gNa inactivation');

subplot(2,1,2);
plot(v,tau_h);
hold on;
plot(gNa_ss(:,1),gNa_ss(:,5),'r--');
xlabel('membrane potential (mV)');
ylabel('time constant (ms)');

%
%For following rate constants, see eq. 6 of [PR94] and erratum
%

%forward rate constant for fast sodium
function val = am(v)
val = 0.32*(13.1-v)./(exp((13.1-v)/4)-1);

%backward rate constant for fast sodium
function val = bm(v)
val = 0.28*(v-40.1)./(exp((v-40.1)/5)-1);

%forward rate constant for DR activation
function val = an(v)
val = 0.016*(35.1-v)./(exp((35.1-v)/5)-1);

%backward rate constant for DR activation
function val = bn(v)
val = 0.25*exp(0.5-0.025*v);

%forward rate constant for sodium inactivation
function val = ah(v)
val = 0.128*exp((17-v)/18);

%backward rate constant for sodium inactivation
function val = bh(v)
val = 4./(1+exp((40-v)/5));

%forward rate constant for Ca current activation
function val = as(v)
val = 1.6./(1+exp(-0.072*(v-65)));

%backward rate constant for Ca current activation
function val = bs(v)
val = 0.02*(v-51.1)./(exp((v-51.1)/5)-1);

%forward rate constant for Ca current inactivation
%note: the following will not work correctly in a vectorized case 
%if v <= 0
%    val = 0.005;
%else
%    val = exp(-v/20)/200;
%end

function val = ar(v)
val = (v <= 0)*0.005 + (v > 0).*exp(-v/20)/200;

%backward rate constant for Ca current inactivation
function val = br(v)
val = (v > 0).*(0.005 - ar(v));

%forward rate constant for fast Ca and V dependent K current
% if v <= 50
%    %see erratum 
%    val = exp((v-10)/11-(v-6.5)/27)/18.975;
% else
%    val = 2*exp((6.5-v)/27);
% end
% Here is a parallelized version of the code
function val = ac(v)
val = (v <= 50).*exp((v-10)/11-(v-6.5)/27)/18.975 + ...
    (v > 50).*2.*exp((6.5-v)/27);

%backward rate constant for fast Ca and V dependent K current
function val = bc(v)
val = (v<= 50).*(2*exp((6.5-v)/27) - exp((v-10)/11-(v-6.5)/27)/18.975);


%forward rate constant for slow Ca dependent K current
function val = aq(v)
val = min((0.00002)*v,0.01);

%backward rate constant for slow Ca dependent K current
function val = bq(v)
val = 0.001;

%forward rate constant for A-type current activation
%(variable a of Traub et al. 1991)
function val = aa(v)
val = 0.02*(13.1-v)./(exp( (13.1-v)/10 ) - 1);

%backward rate constant for A-type current activation 
%(variable a of Traub et al. 1991)
function val = ba(v)
val = 0.0175*(v-40.1)./(exp( (v-40.1)/10 ) - 1);

%forward rate constant for A-type current inactivation
%(variable b of Traub et al. 1991)
function val = ab(v)
val = 0.0016*exp( (-13 -v)/18 );

%backward rate constant for A-type current inactivation 
%(variable a of Traub et al. 1991)
function val = bb(v)
val = 0.05./( 1 + exp( (10.1-v)/5 ) );
