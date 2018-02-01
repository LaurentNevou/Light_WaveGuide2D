%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 31Jan2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Solving Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=1.55E-6;       %% Wavelength (meter)

neff_guess=2;         %% Guess of the effective index solutions (NOT WORKING)
neff_min=1.5;         %% filter the solutions where the effective index is superior than
neff_max=4;           %% filter the solutions where the effective index is inferior than
nmodes=10;            %% number of solutions asked 

AbsorbingBoundaryCondition=0;     %% 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AAbs=0;               %% Plot abs(E,H)
RReal=1;              %% Plot real(E,H)
IImag=0;              %% Plot imag(E,H)

plotEx=1;
plotEy=1;
plotEz=0;
plotHx=0;
plotHy=0;
plotHz=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grid definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=61;                %% number of points in the direction x
Ny=51;                %% number of points in the direction y

Dx=1E-6;              %% map
Dy=4E-6;              %% map

x = linspace(-Dx, Dx, Nx);
%y = linspace(-Dy/1.02, Dy, Ny);
y = linspace(-0.5e-6, 1.5e-6, Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);

[X,Y] = meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx=1e-6;
Ly=1e-6;

x0=0;y0=Ly/2;

idx  = 1 > abs((X-x0)/Lx*2);
idy  = 1 > abs((Y-y0)/Ly*2);
idXY = idx.*idy ;

n = (2.5)*idXY + (1)*(1-idXY);  %% ridge optical index

n(Y<=0)=1.5+0i;                %% substrate optical index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AbsorbingBoundaryCondition==1
    LOSS=1;
    n(:,1)         = n(:,1)         + LOSS*i;
    n(:,end)       = n(:,end)       + LOSS*i;
    n(1,2:end-1)   = n(1,2:end-1)   + LOSS*i;
    n(end,2:end-1) = n(end,2:end-1) + LOSS*i;
end

eps=n.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure('position',[100 500 700 500])
pcolor(x*1e6,y*1e6,abs(n))
%pcolor(x*1e6,y*1e6,abs(eps))
xlabel('x (um)')
ylabel('y (um)')
colormap(jet)
colorbar
%shading flat
%break

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hx=[];Hy=[];Hz=[];Ex=[];Ey=[];Ez=[];

[Ex,Ey,neff]=WG2D_SolveExEy_f(x,y,eps,lambda,nmodes,neff_guess,neff_min,neff_max);

display(strcat('neff = ',num2str(neff,'%.5f')   ));

if plotEz==1;
  for i=1:length(neff)
      Ez(:,:,i)=WG2D_GetEz_f(x,y,Ex(:,:,i),Ey(:,:,i),eps,neff(i),lambda);
  end
end

if plotHx==1 || plotHy==1 || plotHz==1;
  for i=1:length(neff)
      [Hx(:,:,i),Hy(:,:,i),Hz(:,:,i)]=WG2D_GetHxHyHz_f(x,y,Ex(:,:,i),Ey(:,:,i),eps,neff(i),lambda);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AAbs==1
  HHx=abs(Hx);HHy=abs(Hy);HHz=abs(Hz);
  EEx=abs(Ex);EEy=abs(Ey);EEz=abs(Ez);
end
if RReal==1
  HHx=real(Hx);HHy=real(Hy);HHz=real(Hz);
  EEx=real(Ex);EEy=real(Ey);EEz=real(Ez);  
end
if IImag==1
  HHx=imag(Hx);HHy=imag(Hy);HHz=imag(Hz);
  EEx=imag(Ex);EEy=imag(Ey);EEz=imag(Ez);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0fig=-1900; Y0fig=100;
Wfig=1800;Hfig=1000;

if plotEx==1

figure('Name','Ex','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)
    
    subplot(2,6,i)
    hold on
    A=EEx(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    plot([-1 1]*Dx*1e6,[0 0],'w','linewidth',2)
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Ex'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotEy==1
  
figure('Name','Ey','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)

    subplot(2,6,i)
    hold on
    A=EEy(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Ey'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotEz==1

figure('Name','Ez','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)

    subplot(2,6,i)
    hold on
    A=EEz(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Ez'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotHx==1

figure('Name','Hx','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)
    
    subplot(2,6,i)
    hold on
    A=HHx(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    plot([-1 1]*Dx*1e6,[0 0],'w','linewidth',2)
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Hx'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotHy==1
  
figure('Name','Hy','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)

    subplot(2,6,i)
    hold on
    A=HHy(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Hy'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotHz==1

figure('Name','Hz','position',[X0fig Y0fig Wfig Hfig])

for i=1:length(neff)

    subplot(2,6,i)
    hold on
    A=HHz(:,:,i);
    pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Hz'  )    )
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
