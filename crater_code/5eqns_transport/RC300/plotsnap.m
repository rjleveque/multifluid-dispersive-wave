function plotsnap(m,j)  
n1 = j+10000;
fname = ['fort.',num2str(n1)];
fname(6) = 't';
fid  = fopen(fname);
t1   = fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
fclose(fid);
%
fname(6) = 'c';
fid    = fopen(fname);
data1  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
data1 = data1';
%
framest = [num2str(j)];
%
p0 = 0.123e6;
h0 = 4000.0;
km = 1000.0;
%
clf
if m==1
   plot(data1(:,1)/km,data1(:,3)-h0,'b-',...
        'LineWidth',1)
%
  eta_minmax = [min(min(data1(:,3))) max(max(data1(:,3)))]-h0
%
   title(['surface displacement at time t=', num2str(t1),'s'],...
             'fontsize',20,'interpreter','latex')
%
   pname = ['crater_5eqns_transport_RC300_eta' framest];
elseif m==2
   plot(data1(:,1)/km,data1(:,2),'b-',...
        'LineWidth',1)

   ur_minmax = [min(min(data1(:,2))) max(max(data1(:,2)))]

   title(['depth average radial velocity at time $t=0$ and t=', num2str(t1),'s'],...
             'fontsize',20,'interpreter','latex')
%
   pname = ['crater_5eqns_transport_ur' framest];
end
%
xlabel('$r(km)$','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
%xmin = min(min(data0(:,1)));
%xmax = max(max(data0(:,1)));
%xlim([xmin xmax])
%
grid on
%printeps(pname)      
end
