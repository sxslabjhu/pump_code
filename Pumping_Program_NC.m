% Program for: Kidney Epithelial Cells are Active Mechano-biological Fluid Pumps

% JactiveMode: 1: Jcactive = 0 or a constant
%              2: Jcactive = gammac*c0, gammac is constant
%              3: Jcactive = gammac*c0, gammac is linear in c0
%              4: Jcactive = -gammac (Delta muc - Delta mu)
%              5: Jcactive = gammac*c0*(1 + kp*(p0 - p0r))
%              6: Jcactive = gammac*((cin(1)-cin(N))/c0r*kp*(p0 - p0r));
%              7: Jcactive = -gammac*(Delta muc - Delta mu)*(1 + kp*(p0 - p0r))

clear
clc

RankCheck = 0;
JactiveMode = 6;

Iter = 21;
R0 = 5d-6;              % (m) radius of the cell
b = 10.d-6;             % (m) cell width
w = 10.d-6;             % (m) cell depth
L = 5.d-6;              % (m) cell length
h = 0.5d-6;             % (m) membrane thickness
mu = 1d-3;              % (Pa s) dynamic viscosity of the fluid
p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
p0r = 0*1d5;            % (Pa) referece pressure where gammac = 1
kp = 2d-2;              % (1/Pa) slope of pressure on gammac
c0f = 340;              % (mol/m^3 = mM) external ion concentration at the front
c0b = 340;              % (mol/m^3 = mM) external ion concentration at the back
c0r = 340;              % (mol/m^3 = mM) reference ion concentration where gammac = 1
c0c = 0;                % (mol/m^3 = mM) critical ion concentration where gammac = 0
deltaccf = 1;           % (mol/m^3 = mM) critical ion concentration difference at the front
deltaccb = 1;           % (mol/m^3 = mM) critical ion concentration difference at the back
cinb = 340.35;          % expected value, not prescribed
cinf = 340.87;          % expected value, not prescribed
Dcommon = 1d-3;         % (nm m/s) diffusion constant for ions
vc0 = 0;                % (nm/s) expected average cytosol velocity

R = 8.31451;            % (J/mol K) Ideal gas constant
T = 300;                % (K) absolute temperature
kB = 1.38d-23;          % (J/K) Boltzmann constant
NA = 6.02d23;           % (1/mol) Avogadios number

alphaf = 1d0*1.d-1;     % (nm/Pa/s) coefficient of water permeation
alphab = 1d0*1.d-1;     % (nm/Pa/s) coefficient of water permeation
gf = 1d-0*5d3;          % (nm/s) passive channel coefficient at the front
gb = 1d-0*5d3;          % (nm/s) passive channel coefficient at the back
gammacf = 5d-3*2d4;     % (nm/s) active channel coefficient at the front
gammacb = 5d-3*2d4;     % (nm/s) active channel coefficient at the back
if JactiveMode == 4
    alphaf = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
    alphab = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
    gf = 1d3;           % (nm/s) passive channel coefficient at the front
    gb = 1d3;           % (nm/s) passive channel coefficient at the back
    gammacf = 1d2;      % (nm/s) active channel coefficient at the front
    gammacb = 1d2;      % (nm/s) active channel coefficient at the back
end

if JactiveMode == 6
    gammacf = 2d4;      % (nm/s) active channel coefficient at the front
    gammacb = 2d4;      % (nm/s) active channel coefficient at the back
end
Jcactivef = 0d-5;
Jcactiveb = 0d-5;

fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell

N = 21;
dx = L/(N-1);
x = linspace(0,L,N);

D = Dcommon*ones(N,1);  % (nm m/s) diffusion constant for ions

N1 = 31;
N2 = 33;
C0F = linspace(200,340,N1);
Ym = C0F;
P0B = p0b + linspace(0,500,N2);
Xm = P0B;
[Xmesh,Ymesh] = meshgrid(Xm,Ym);

JWATERF = zeros(N1,N2);
PCF = zeros(N1,N2);
PCB = zeros(N1,N2);
CINB = zeros(N1,N2);
CINF = zeros(N1,N2);
POWER = zeros(N1,N2);

for loop2 = 1:N2
    
    p0b = P0B(loop2);
    
    for loop1 = 1:N1
        
        c0f = C0F(loop1);
        
        if JactiveMode == 2
            Jcactivef = gammacf*c0f;
            Jcactiveb = gammacb*c0b;
        elseif JactiveMode == 3
            Jcactivef = gammacf*(c0f-c0c)/(c0r-c0c)*c0f;
            Jcactiveb = gammacb*(c0b-c0c)/(c0r-c0c)*c0b;
        elseif JactiveMode == 5
            Jcactivef = gammacf*c0f*(kp*(p0f - p0r));
            Jcactiveb = gammacb*c0b*(kp*(p0b - p0r));
        end
        
        DF = zeros(N+2,N+2);
        Fn = zeros(N+2,1);
        
        % initial guess
        if loop1 == 1 && loop2 == 1
            cin = linspace(cinb,cinf,N)';
            pb = 1.5d3;
            pf = 1.5d3;
        elseif loop1 == 1 && loop2 > 1
            cin = temp_cin;
            pb = temp_pb;
            pf = temp_pf;
        end
        X = [cin; pb; pf];
        
        iter = 0;
        ITER = true;
        while ITER
            iter = iter + 1;
            
            if JactiveMode == 4
                Jcactivef = -gammacf*R*T*(cin(N)-c0f-deltaccf);
                Jcactiveb = -gammacb*R*T*(cin(1)-c0b-deltaccb);
            elseif JactiveMode == 6
                Jcactivef = gammacf*((cin(1)-cin(N))/c0r*kp*(p0f - p0r));
                Jcactiveb = gammacb*((cin(1)-cin(N))/c0r*kp*(p0b - p0r));
            elseif JactiveMode == 7
                Jcactivef = -gammacf*R*T*(cin(N)-c0f-deltaccf)*((cin(1)-cin(N))/c0r*kp*(p0f - p0r));
                Jcactiveb = -gammacb*R*T*(cin(1)-c0b-deltaccb)*((cin(1)-cin(N))/c0r*kp*(p0b - p0r));
            end
            
            Fn(1) = -gb*(cin(1)-c0b) + Jcactiveb +(D(1)+D(2))/2*(cin(2)-cin(1))/dx ...
                + R0^2/8/mu/L*(pf-pb)*(cin(1)+cin(2))/2;
            Fn(2:N-1) = -(D(2:N-1)+D(1:N-2))/2/dx.*(cin(2:N-1)-cin(1:N-2)) ...
                - R0^2/8/mu/L*(pf-pb)*(cin(2:N-1)+cin(1:N-2))/2 ...
                + (D(2:N-1)+D(3:N))/2/dx.*(cin(3:N)-cin(2:N-1)) ...
                + R0^2/8/mu/L*(pf-pb)*(cin(2:N-1)+cin(3:N))/2;
            Fn(N) = -gf*(cin(N)-c0f) + Jcactivef - (D(N)+D(N-1))/2*(cin(N)-cin(N-1))/dx ...
                - R0^2/8/mu/L*(pf-pb)*(cin(N)+cin(N-1))/2;
            Fn(N+1) = alphaf*(pf - p0f) - alphaf*R*T*(cin(N) - c0f) + R0^2/8/mu/L*(pf-pb);
            Fn(N+2) = -alphab*(pb - p0b) + alphab*R*T*(cin(1) - c0b) ...
                -alphaf*(pf - p0f) + alphaf*R*T*(cin(N) - c0f);
            
            DF(1,1) = -gb - 1/2/dx*(D(1)+D(2)) + R0^2/16/mu/L*(pf - pb);
            DF(1,2) = 1/2/dx*(D(1)+D(2)) + R0^2/16/mu/L*(pf - pb);
            DF(1,N+1:N+2) = [-1, 1]*R0^2/16/mu/L*(cin(1)+cin(2));
            DF(N,N-1) = 1/2/dx*(D(N)+D(N-1)) - R0^2/16/mu/L*(pf - pb);
            DF(N,N) = -gf - 1/2/dx*(D(N)+D(N-1)) - R0^2/16/mu/L*(pf - pb);
            DF(N,N+1:N+2) = [1,-1]*R0^2/16/mu/L*(cin(N)+cin(N-1));
            DF(N+1,N) = -alphaf*R*T;
            DF(N+1,N+1) = -R0^2/8/mu/L;
            DF(N+1,N+2) = alphaf + R0^2/8/mu/L;
            DF(N+2,[1,N]) = R*T*[alphab, alphaf];
            DF(N+2,N+1:N+2) = -[alphab, alphaf];
            
            if JactiveMode == 4
                DF(1,1) = DF(1,1) - gammacb*R*T;
                DF(N,N) = DF(N,N) - gammacf*R*T;
            elseif JactiveMode == 6
                DF(1,1) = DF(1,1) + gammacb/c0r*kp*(p0b-p0r);
                DF(1,N) = DF(1,N) - gammacb/c0r*kp*(p0b-p0r);
                DF(N,1) = DF(N,1) + gammacf/c0r*kp*(p0f-p0r);
                DF(N,N) = DF(N,N) - gammacf/c0r*kp*(p0f-p0r);
            elseif JactiveMode == 7
                DF(1,1) = DF(1,1) - gammacb*R*T*(cin(1)-c0b-deltaccb)/c0r*kp*(p0b-p0r) ...
                    - gammacb*R*T*((cin(1)-cin(N))/c0r*kp*(p0b - p0r));
                DF(1,N) = DF(1,N) + gammacb*R*T*(cin(1)-c0b-deltaccb)/c0r*kp*(p0b-p0r);
                DF(N,1) = DF(N,1) - gammacf*R*T*(cin(N)-c0f-deltaccf)/c0r*kp*(p0f-p0r);
                DF(N,N) = DF(N,N) + gammacf*R*T*(cin(N)-c0f-deltaccf)/c0r*kp*(p0f-p0r) ...
                    - gammacf*R*T*((cin(1)-cin(N))/c0r*kp*(p0f - p0r));
            end
            
            for i = 2:N-1
                DF(i,i-1) = 1/2/dx*(D(i)+D(i-1)) - R0^2/16/mu/L*(pf - pb);
                DF(i,i) = -1/2/dx*(D(i)+D(i-1)) - 1/2/dx*(D(i)+D(i+1));
                DF(i,i+1) = 1/2/dx*(D(i)+D(i+1)) + R0^2/16/mu/L*(pf - pb);
                DF(i,N+1:N+2) = [1,-1]*R0^2/16/mu/L*(cin(i-1)-cin(i+1));
            end
            
            DF = sparse(DF);
            X = X - DF\Fn;
            
            cin = X(1:N);
            pb = X(N+1);
            pf = X(N+2);
            
            if iter > 1
                error = abs((X-temp_X)./X);
                error = sum(error)/(N+2);
                if error < 1d-4 || iter == Iter
                    ITER = false;
                end
            end
            temp_X = X;
        end
        
        if loop1 == 1
            temp_cin = cin;
            temp_pb = pb;
            temp_pf = pf;
        end
        
        Jwaterf = -alphaf*(pf-p0f-R*T*cin(N)+R*T*c0f);
        Jwaterb = -alphab*(pb-p0b-R*T*cin(1)+R*T*c0b);
        power = Jwaterf*p0b*pi*R0^2;
        
        if RankCheck == 1
            DF = full(DF);
            RankCheck = (N+2) - rank(DF);
            if RankCheck ~= 0
                fprintf('RankCheck = %d\n',RankCheck);
            end
        end
        
        CINB(loop1,loop2) = cin(1);
        CINF(loop1,loop2) = cin(N);
        JWATERF(loop1,loop2) = Jwaterf;
        PCF(loop1,loop2)  = pf;
        PCB(loop1,loop2)  = pb;
        POWER(loop1,loop2)  = power;
        
    end
end

%%
PumpForce = (PCF-PCB)*pi*R0^2*1d9;  % nano-Newton

figure
for i = 1:3:N1
    plot(P0B,JWATERF(i,:)/max(max(JWATERF)),'-','linewidth',2,'Color',[1,1,1]*(1-i/N1)); hold on
    set(gca,'fontsize',13);
    xlabel('p_0^b (Pa)');
    ylabel('J_{water} (nm/s)')
end

