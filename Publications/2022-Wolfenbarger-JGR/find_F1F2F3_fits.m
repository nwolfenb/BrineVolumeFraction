% find_F1F2F3_fits.m
% Used to find the polynomial cubic spline defining the phase behavior
% functions for a solution composition and pressure defined by FREZCHEM
% version 13.3 or 15.1 over the temperature range from the freezing point
% to the eutectic point. Fits are published in the Supporting Information
% for
%
% Wolfenbarger, N. S., Fox-Powell, M. G., Buffo, J. J., Soderlund, K. M., &
% Blankenship, D. D. (2022). Compositional controls on the distribution of
% brine in Europa's ice shell. Journal of Geophysical Research: Planets,
% 127, e2022JE007305. https://doi.org/10.1029/2022JE007305

clear all; close all; clc
addpath('..\..\FREZCHEM')
addpath('..\..\Ice Ih Properties')
addpath('..\..\Water Properties')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FREZCHEM v13.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Earth
fn = '..\..\Freezing Simulations\FREZCHEM\v13.3\Seawater\Gitterman';
[T,F1,F2,F3] = F1F2F3_FREZCHEM(fn);
T = round(T,1);

figure
subplot(3,1,1)
plot(T,F1,'k','linewidth',1)

T1 = T(T>=-22.8);
F11 = F1(T>=-22.8);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F11;
C = H;
d = y;
T0 = T(1);
F10 = F1(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.8 & T>=-33.3);
F12 = F1(T<=-22.8 & T>=-33.3);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F12;
C = H;
d = y;
T0 = T2(1);
F10 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

T3 = T(T<-33.3);
F13 = F1(T<-33.3);
H3 = [T3.^3 T3.^2 T3 ones(size(T3))];
y3 = F13;
x3 = H3\y3;
f3 = H3*x3;

hold on
plot(T(T>=-22.8),f1,'r--','linewidth',1)
plot(T(T<=-22.8 & T>=-33.3),f2,'r--','linewidth',1)
plot(T(T<-33.3),f3,'r--','linewidth',1)

ax1 = gca;
ax1.XDir = 'reverse';
ax1.XTickLabel = [];
ax1.TickLabelInterpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.Title.Interpreter = 'latex';
axis tight
ylabel('$F_1(T)$ (kg/m$^{3}$)')
title('Earth Seawater')

fprintf('Earth Seawater\n')
fprintf('F1(T)\n')
fprintf('T>=-22.8\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.8 & T>=-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('T<-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',x3)

subplot(3,1,2)
plot(T,F2,'k','linewidth',1)
T1 = T(T>=-22.8);
F21 = F2(T>=-22.8);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F21;
C = H;
d = y;
T0 = T(1);
F20 = F2(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.8 & T>=-33.3);
F22 = F2(T<=-22.8 & T>=-33.3);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F22;
C = H;
d = y;
T0 = T2(1);
F20 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

T3 = T(T<-33.3);
F23 = F2(T<-33.3);
H3 = [T3.^3 T3.^2 T3 ones(size(T3))];
y3 = F23;
x3 = H3\y3;
f3 = H3*x3;

hold on
plot(T(T>=-22.8),f1,'r--','linewidth',1)
plot(T(T<=-22.8 & T>=-33.3),f2,'r--','linewidth',1)
plot(T(T<-33.3),f3,'r--','linewidth',1)
ax2 = gca;
ax2.XDir = 'reverse';
ax2.XTickLabel = [];
ax2.TickLabelInterpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.Title.Interpreter = 'latex';
axis tight
ylabel('$F_2(T)$ ')

fprintf('F2(T)\n')
fprintf('T>=-22.8\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.8 & T>=-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('T<-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',x3)
fprintf('\n')

subplot(3,1,3)
plot(T,F3,'k','linewidth',1)

T1 = T(T>=-22.8 & T<=-6.4);
F31 = F3(T>=-22.8 & T<=-6.4);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F31;
C = H;
d = y;
T0 = T1(1);
F30 = F31(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F30;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.8 & T>=-33.3);
F32 = F3(T<=-22.8 & T>=-33.3);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F32;
C = H;
d = y;
T0 = T2(1);
F20 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

T3 = T(T<-33.3);
F33 = F3(T<-33.3);
H3 = [T3.^3 T3.^2 T3 ones(size(T3))];
y3 = F33;
x3 = H3\y3;
f3 = H3*x3;

hold on
plot(T(T>=-22.8 & T<=-6.4),f1,'r--','linewidth',1)
plot(T(T<=-22.8 & T>=-33.3),f2,'r--','linewidth',1)
plot(T(T<-33.3),f3,'r--','linewidth',1)
ax3 = gca;
ax3.XDir = 'reverse';
ax3.TickLabelInterpreter = 'latex';
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';
ax3.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature, $T$ ($^{\circ}$C)')
ylabel('$F_3(T)$')

fprintf('F3(T)\n')
fprintf('T>=-22.8 & T<=-6.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.8 & T>=-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('T<-33.3\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',x3)
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FREZCHEM v15.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NaCl
fn = '..\..\Freezing Simulations\FREZCHEM\v15.1\NaCl\NaCl_1ppt';
[T,F1,F2,~] = F1F2F3_FREZCHEM(fn);
T = round(T,1);

figure
subplot(2,1,1)
plot(T,F1,'k','linewidth',1)
hold on
T1 = T;
F11 = F1;
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F11;
C = H;
d = y;
T0 = T(1);
F10 = F1(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);
plot(T,f1,'r--','linewidth',1)
ax1 = gca;
ax1.XDir = 'reverse';
ax1.TickLabelInterpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_1(T)$ (kg/m$^{3}$)')
title('NaCl')

fprintf('NaCl\n')
fprintf('F1(T)\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)

subplot(2,1,2)
plot(T,F2,'k','linewidth',1)
T1 = T;
F21 = F2;
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F21;
C = H;
d = y;
T0 = T(1);
F20 = F2(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);
hold on
plot(T,f1,'r--','linewidth',1)
ax2 = gca;
ax2.XDir = 'reverse';
ax2.TickLabelInterpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_2(T)$')

fprintf('F2(T)\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('\n')

%% MgSO4
fn = '..\..\Freezing Simulations\FREZCHEM\v15.1\MgSO4\MgSO4_1ppt';
[T,F1,F2,~] = F1F2F3_FREZCHEM(fn);
T = round(T,1);

figure
subplot(2,1,1)
plot(T,F1,'k','linewidth',1)
hold on
T1 = T;
F11 = F1;
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F11;
C = H;
d = y;
T0 = T(1);
F10 = F1(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);
plot(T,f1,'r--','linewidth',1)
ax1 = gca;
ax1.XDir = 'reverse';
ax1.TickLabelInterpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_1(T)$ (kg/m$^{3}$)')
title('MgSO$_4$')

fprintf('MgSO4\n')
fprintf('F1(T)\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)

subplot(2,1,2)
plot(T,F2,'k','linewidth',1)
T1 = T;
F21 = F2;
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F21;
C = H;
d = y;
T0 = T(1);
F20 = F2(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);
hold on
plot(T,f1,'r--','linewidth',1)
ax2 = gca;
ax2.XDir = 'reverse';
ax2.TickLabelInterpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_2(T)$')

fprintf('F2(T)\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('\n')

%% Europa (Cl-Dominated)
fn = '..\..\Freezing Simulations\FREZCHEM\v15.1\Cl_Na_Mg_SO4_seawater\Cl_Na_Mg_SO4_seawater';
[T,F1,F2,F3] = F1F2F3_FREZCHEM(fn);
T = round(T,1);

figure
subplot(3,1,1)
plot(T,F1,'k','linewidth',1)

T1 = T(T>=-22.4);
F11 = F1(T>=-22.4);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F11;
C = H;
d = y;
T0 = T(1);
F10 = F1(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.4);
F12 = F1(T<=-22.4);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F12;
C = H;
d = y;
T0 = T2(1);
F10 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);


hold on
plot(T(T>=-22.4),f1,'r--','linewidth',1)
plot(T(T<=-22.4),f2,'r--','linewidth',1)
ax1 = gca;
ax1.XDir = 'reverse';
ax1.XTickLabel = [];
ax1.TickLabelInterpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.Title.Interpreter = 'latex';
axis tight
ylabel('$F_1(T)$ (kg/m$^{3}$)')
title('Cl-Dominated')

fprintf('Cl-Dominated\n')
fprintf('F1(T)\n')
fprintf('T>=-22.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')

subplot(3,1,2)
plot(T,F2,'k','linewidth',1)
T1 = T(T>=-22.4);
F21 = F2(T>=-22.4);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F21;
C = H;
d = y;
T0 = T(1);
F20 = F2(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.4);
F22 = F2(T<=-22.4);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F22;
C = H;
d = y;
T0 = T2(1);
F20 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

hold on
plot(T(T>=-22.4),f1,'r--','linewidth',1)
plot(T(T<=-22.4),f2,'r--','linewidth',1)
ax2 = gca;
ax2.XDir = 'reverse';
ax2.XTickLabel = [];
ax2.TickLabelInterpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.Title.Interpreter = 'latex';
axis tight
ylabel('$F_2(T)$')

fprintf('F2(T)\n')
fprintf('T>=-22.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')

subplot(3,1,3)
plot(T,F3,'k','linewidth',1)
T1 = T(T>=-22.4 & T<=-6.1);
F31 = F3(T>=-22.4 & T<=-6.1);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F31;
C = H;
d = y;
T0 = T1(1);
F30 = F31(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F30;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-22.4);
F32 = F3(T<=-22.4);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F32;
C = H;
d = y;
T0 = T2(1);
F20 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);


hold on
plot(T(T>=-22.4 & T<=-6.1),f1,'r--','linewidth',1)
plot(T(T<=-22.4),f2,'r--','linewidth',1)
ax3 = gca;
ax3.XDir = 'reverse';
ax3.TickLabelInterpreter = 'latex';
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';
ax3.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_3(T)$')

fprintf('F3(T)\n')
fprintf('T>=-22.4 & T<=-6.1\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-22.4\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')

%% Europa (SO4-Dominated)
fn = '..\..\Freezing Simulations\FREZCHEM\v15.1\SO4_Mg_Na_Cl\SO4_Mg_Na_Cl';
[T,F1,F2,F3] = F1F2F3_FREZCHEM(fn);
T = round(T,1);

figure
subplot(3,1,1)
plot(T,F1,'k','linewidth',1)
T1 = T(T>=-5.6);
F11 = F1(T>=-5.6);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F11;
C = H;
d = y;
T0 = T(1);
F10 = F1(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-5.6);
F12 = F1(T<=-5.6);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F12;
C = H;
d = y;
T0 = T2(1);
F10 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F10;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

hold on
plot(T(T>=-5.6),f1,'r--','linewidth',1)
plot(T(T<=-5.6),f2,'r--','linewidth',1)
ax1 = gca;
ax1.XDir = 'reverse';
ax1.XTickLabel = [];
ax1.TickLabelInterpreter = 'latex';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax1.Title.Interpreter = 'latex';
axis tight
ylabel('$F_1(T)$ (kg/m$^{3}$)')
title('SO$_4$-Dominated')

fprintf('SO4-Dominated\n')
fprintf('F1(T)\n')
fprintf('T>=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')

subplot(3,1,2)
plot(T,F2,'k','linewidth',1)
T1 = T(T>=-5.6);
F21 = F2(T>=-5.6);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F21;
C = H;
d = y;
T0 = T(1);
F20 = F2(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-5.6);
F22 = F2(T<=-5.6);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F22;
C = H;
d = y;
T0 = T2(1);
F20 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F20;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

hold on
plot(T(T>=-5.6),f1,'r--','linewidth',1)
plot(T(T<=-5.6),f2,'r--','linewidth',1)
ax2 = gca;
ax2.XDir = 'reverse';
ax2.XTickLabel = [];
ax2.TickLabelInterpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.Title.Interpreter = 'latex';
axis tight
ylabel('$F_2(T)$')

fprintf('F2(T)\n')
fprintf('T>=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')

subplot(3,1,3)
plot(T,F3,'k','linewidth',1)
T1 = T(T>=-5.6 & T<=-3.4);
F31 = F3(T>=-5.6 & T<=-3.4);
H = [T1.^3 T1.^2 T1 ones(size(T1))];
y = F31;
C = H;
d = y;
T0 = T(1);
F30 = F31(1);
Aeq = [T0^3 T0^2 T0 1];
beq = F30;
p1 = lsqlin(C,d,[],[],Aeq,beq);
f1 = polyval(p1,T1);

T2 = T(T<=-5.6);
F32 = F3(T<=-5.6);
H = [T2.^3 T2.^2 T2 ones(size(T2))];
y = F32;
C = H;
d = y;
T0 = T2(1);
F30 = polyval(p1,T0);
Aeq = [T0^3 T0^2 T0 1];
beq = F30;
p2 = lsqlin(C,d,[],[],Aeq,beq);
f2 = polyval(p2,T2);

hold on
plot(T(T>=-5.6 & T<=-3.4),f1,'r--','linewidth',1)
plot(T(T<=-5.6),f2,'r--','linewidth',1)
ax3 = gca;
ax3.XDir = 'reverse';
ax3.TickLabelInterpreter = 'latex';
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';
ax3.Title.Interpreter = 'latex';
axis tight
xlabel('Temperature ($^{\circ}$C)')
ylabel('$F_3(T)$')

fprintf('F3(T)\n')
fprintf('T>=-3.4 & T<=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p1)
fprintf('T<=-5.6\n')
fprintf('%1.4eT^3 + %1.4eT^2 + %1.4eT + %1.4e\n',p2)
fprintf('\n')
