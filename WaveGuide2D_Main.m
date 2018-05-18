%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 26Feb2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%
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
nmodes=6;            %% number of solutions asked 

AbsorbingBoundaryCondition=0;     %% 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AAbs=0;               %% Plot abs(E,H)
RReal=1;              %% Plot real(E,H)
IImag=0;              %% Plot imag(E,H)

plotElectricalField=1;
plotMagneticalField=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grid definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=41;                %% number of points in the direction x
Ny=31;                %% number of points in the direction y

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
    LOSS=1e-5;
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

if plotElectricalField==1;
  for i=1:length(neff)
      Ez(:,:,i)=WG2D_GetEz_f(x,y,Ex(:,:,i),Ey(:,:,i),eps,neff(i),lambda);
  end
end

if plotMagneticalField==1
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

X0fig=50; Y0fig=100;
Wfig=1800;Hfig=1000;

if plotElectricalField==1

figure('Name','Electrical Field','position',[X0fig Y0fig Wfig Hfig])


for i=1:6%length(neff)
    
    subplot(3,6,i)
    hold on
    A=EEx(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    plot([-1 1]*Dx*1e6,[0 0],'w','linewidth',2)
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Ex)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Ex)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Ex)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:6%length(neff)

    subplot(3,6,i+6)
    hold on
    A=EEy(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Ey)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Ey)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Ey)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:6%length(neff)

    subplot(3,6,i+12)
    hold on
    A=EEz(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Ez)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Ez)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Ez)'  )    )
    end
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

if plotMagneticalField==1

figure('Name','Magnetical Field','position',[X0fig Y0fig Wfig Hfig])

for i=1:6%length(neff)
    
    subplot(3,6,i)
    hold on
    A=HHx(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    plot([-1 1]*Dx*1e6,[0 0],'w','linewidth',2)
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Hx)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Hx)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Hx)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for i=1:6%length(neff)

    subplot(3,6,i+6)
    hold on
    A=HHy(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Hy)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Hy)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Hy)'  )    )
    end
    shading flat
    xlabel('x (um)')
    ylabel('y (um)')
    colormap(jet)
    colorbar
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:6%length(neff)

    subplot(3,6,i+12)
    hold on
    A=HHz(:,:,i);
    if sum(A(:)==0)
      pcolor(x*1e6,y*1e6,A)
    else
      pcolor(x*1e6,y*1e6,A/max(abs(A(:))))
    end
    contour(x*1e6,y*1e6,abs(n),1,'linewidth',2,'linecolor','w')
    
    if RReal==1 || IImag==1
        caxis([-1 1])
    elseif AAbs==1
        caxis([0 1])
    end
    
    if AAbs==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Abs(Hz)'  )    )
    elseif RReal==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Re(Hz)'  )    )
    elseif IImag==1
      title(  strcat('neff=',  num2str( neff(i) ,'%.3f'),',  Im(Hz)'  )    )
    end
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
