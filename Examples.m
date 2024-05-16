%% requirements:
% 1. gcc for compiling c codes. Load before starting Matlab:
%    module load gcc/6.3.0
% 2. Use Matlab 2020a
%% Setup the directory where the membrane object is located and add the directory to Matlab's function pool 
dir_mod = '/home2/s171152/codes/matlab/mine/git/memCompCourse/codes2024/Session2Codes';
addpath(genpath(dirMod));
%--------------------------------------------------------------------------
% create 'unit' u using the unit module, and 'membrane' m using the membrane module 
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); 
m=ModMembrane(2,'unit',u);
%% Plot the membrane 'm'. Note that Matlab autonatically recognize m is an 'object' and apply m's own plot function  
fig=figure;
subplot(1,2,1);
plot(m,'f',fig);
subplot(1,2,2);
col=rand(m.var.n_coord,3);
plot(m,'f',fig,'col',col,'col_min',0,'col_max',1,'colBar',true);
%--------------------------------------------------------------------------
%% Compute the area, volume, Helfrich bending energy using the provided functions under m
A=Area(m);
V=Volume(m);
r=mean(sqrt(sum(m.var.coord(:,1).^2+m.var.coord(:,2).^2+m.var.coord(:,3).^2,2)));
fprintf('A: %f, %f\n',sum(A),4*pi*r^2);
fprintf('V: %f, %f\n',sum(V),4/3*pi*r^3);
figure;
n=10;
H=zeros(n,1);
for i=1:n
    m.pm.k_c=i;
    H_temp=Helfrich(m);
    H(i)=sum(H_temp);
end
plot((1:n),H,'.'); hold on;
plot((1:n),(1:n)*8*pi,'-');
legend('computed H','line with 8\pi slope');
xlabel('k_c');ylabel('H');
%% Compute the internal force and potential that regulates the edge length of the triangular meshes
m.pm.Vdh.V0=0.02; %adjusting the internal force
[Fi] = Finternal(m,'plot_or_not',true);
%% Remeshing, blue: edge manipulated, green: replacement edges 
[~,remeshed,edg_add] = Remesh(m,0,1,[min(Fi.rg),max(Fi.rg)],'plot_or_not',true);
[~,remeshed,edg_add] = Remesh(m,1,1,[min(Fi.rg),max(Fi.rg)],'plot_or_not',true);