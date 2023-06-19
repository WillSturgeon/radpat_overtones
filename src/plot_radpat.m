clear all;close all;clc;

% Matlab script to plot Rayleigh wave radiation patterns
% Made by Will Sturgeon 19/06/2023. UCL.

n='0';
modes=[028,043,061,098,202,271];
periods=[275,200,151,100,50,37.6];periods=string(periods);

clrs=parula(length(modes));

for i=1:length(modes)
mode=sprintf('%03d',modes(i));

instr=strcat('202201281114A_prem_noocean_n',n,'_l',mode,'.txt');
data=textread(instr);

az=data(:,1);
rad_pat=data(:,2);

figure(1)
p = polarplot(deg2rad(az),rad_pat);hold on;
set(p,'markersize',20,'LineWidth',4)
p.Color=clrs(i,:)

ax = ancestor(p, 'polaraxes')
ax.ThetaZeroLocation = 'bottom';
ax.ThetaDir = 'counterclockwise';

title('PREM N=0. Spectral amplitude in mHz')
set(gca,'Fontsize',20);

end
hold off;
colstr=strcat('T=',periods,' s');
colorbar('TickLabels',colstr)