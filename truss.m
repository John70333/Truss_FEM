clc; clear;
ui1 = [0,0]; ui2 = [2.5,2.5*3^.5]; ui3 = [5,0]; ui4 = [7.5,2.5*3^.5]; ui5 = [10,0];
ui = [ui1,ui2,ui3,ui4,ui5]'; conn = {[2,3],[1,3,4],[1,2,4,5],[2,3,5],[3,4]}';
Ra = [0,-10,0,-10,0,-10,0]'*1e3/3; ai = 3:9; ub = [0,0,0]'; bi = [1,2,10];
A = zeros(length(conn)); E = A; nu = .3; ns = 100;
for i = 1:length(conn)
    for j = conn{i}
        A(i,j) = 1e-4;
        E(i,j) = 1e10;
    end
end
[R,K,u,uf,lci,lcf,eps] = tfem(ui,conn,A,E,nu,Ra,ai,ub,bi,ns);
draw(ui,uf(:,end),conn);
%%
function [R,K,u,uf,lci,lcf,eps] = tfem(ui,conn,A,E,nu,Ra,ai,ub,bi,ns)
n = length(conn); uc = ui; Ras = Ra./ns;
a = @(u)atan2(u(2),u(1)); l = @(u)sum(u.^2)^.5; p = @(i)2*i-1:2*i;
c = @(u1,u2)l(u2-u1)^-1*[cos(a(u2-u1)),sin(a(u2-u1))]'*[cos(a(u2-u1)),sin(a(u2-u1))];
R = zeros(2*n,1,ns); K = R; u = R; uf = R; lci = zeros(n,n,ns); lcf = lci; eps = lci;
for pg = 1:ns
    for i = 1:n
        s = zeros(2);
        for j = conn{i}
            if pg == 1
                lci(i,j) = l(uc(p(i))-uc(p(j))); pc = 1;
            else
                pc = (1-nu*eps(i,j,pg-1))^2;
            end
            k = E(i,j)*A(i,j)*pc*c(uc(p(i)),uc(p(j)));
            K(p(i),p(j),pg) = -k; s = s+k;
        end
        K(p(i),p(i),pg) = s;
    end
    Kaa = K(ai,ai,pg); Kab = K(ai,bi,pg); Kba = K(bi,ai,pg); Kbb = K(bi,bi,pg);
    ua = Kaa\(Ras-Kab*ub); Rb = Kba*ua+Kbb*ub;
    u(ai,:,pg) = ua; u(bi,:,pg) = ub; R(ai,:,pg) = Ras; R(bi,:,pg) = Rb;
    uf(:,:,pg) = uc+u(:,:,pg); uc = uf(:,:,pg);
    for i = 1:n
        for j = conn{i}
            lcf(i,j,pg) = l(uf(p(i),:,pg)-uf(p(j),:,pg));
            eps(i,j,pg) = lcf(i,j,pg)./lci(i,j)-1;
        end
    end
end
end
function draw(ui,uf,conn)
xi = ui(1:2:end-1); yi = ui(2:2:end);
xf = uf(1:2:end-1,:,end); yf = uf(2:2:end,:,end);
figure('windowstate','maximized'); hold on; grid on;
xlim([min([xi;xf])-1,max([xi;xf])+1]); 
ylim([min([yi;yf])-1,max([yi;yf])+1]); 
for i = 1:length(ui)/2
    for j = conn{i}
        if j < i
            plot([xi(i),xi(j)],[yi(i),yi(j)],'b-');
            plot([xf(i),xf(j)],[yf(i),yf(j)],'r-');
        end
    end
end
saveas(gcf,'truss.png');
end
