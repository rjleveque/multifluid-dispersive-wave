function plotslice(m,j,k)
%
% read time and number of grids:
n1 = j+10000;
fname = ['fort.',num2str(n1)];
fname(6) = 't';
%
fid = fopen(fname);
t1 = fscanf(fid,'%g',1);        fscanf(fid,'%s',1);
meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
fclose(fid);
%
qmin = 1e6;
qmax = -1e6;
%
fname(6) = 'g';
fid = fopen(fname);
mx = fscanf(fid,'%d',1);  
my = fscanf(fid,'%d',1);     
grid = fscanf(fid,'%g %g',[2 inf]);
status = fclose(fid);
grid = grid';
%
x = reshape(grid(:,1),mx,my);
y = reshape(grid(:,2),mx,my);
x = x/1.0e3;
y = y/1.0e3;
%
% data set 1
fname(6) = 'q';
fid      = fopen(fname);
%
gridno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
level = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
mx = fscanf(fid,'%d',1);         fscanf(fid,'%s',1);
my = fscanf(fid,'%d',1);         fscanf(fid,'%s',1);   

xlow = fscanf(fid,'%g',1);       fscanf(fid,'%s',1);
ylow = fscanf(fid,'%g',1);       fscanf(fid,'%s',1);
dx = fscanf(fid,'%g',1);         fscanf(fid,'%s',1);
dy = fscanf(fid,'%g',1);         fscanf(fid,'%s',1);
%
data = fscanf(fid,'%g',[meqn,mx*my]);
data = data';
%
if m==1
   q1 = reshape(data(:,1),mx,my)+...
        reshape(data(:,2),mx,my);
elseif m==2
   q1 = reshape(data(:,3),mx,my);  
elseif m==3
   q1 = reshape(data(:,4),mx,my);  
elseif m==4
   q1 = reshape(data(:,5),mx,my);  
elseif m==5
   q1 = reshape(data(:,6),mx,my);  
end
%
status = fclose(fid);
%
vof1 = reshape(data(:,6),mx,my);  

xslice = x(k,1)
yslice = y(k,:);
qslice = q1(k,:);

clf
plot(yslice,qslice,'bx',...
     'LineWidth',1)

framest = [num2str(j)];
if m ==1
   title(['density at t=', num2str(t1),'s'],...
      'fontsize',20,'interpreter','latex')
%
   fname = ['crater_5eqns_transport_rho' framest];
elseif m == 2
   title(['velocity ($u$) at t=', num2str(t1),'s'],...
      'fontsize',20,'interpreter','latex')
%
   fname = ['crater_5eqns_transport_u' framest];
elseif m == 3
   title(['velocity ($v$) at t=', num2str(t1),'s'],...
      'fontsize',20,'interpreter','latex')
%
   fname = ['crater_5eqns_transport_v' framest];
elseif m == 4
   title(['pressure at t=', num2str(t1),'s'],...
      'fontsize',20,'interpreter','latex')
%
   fname = ['crater_5eqns_transport_p' framest];
elseif m == 5
   title(['air volume fraction at t=', num2str(t1),'s'],...
      'fontsize',20,'interpreter','latex')
   fname = ['crater_5eqns_transport_vof' framest];
end
%
legend(['$x$ =', num2str(xslice)],...
      'fontsize',20,'interpreter','latex',...
      'Location','NorthEast',...
      'Box','off')
%
xlabel('$y$(km)','fontsize',20,'interpreter','latex')
%ylabel('$y$(km)','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%axis off
%
%printpng(fname)
