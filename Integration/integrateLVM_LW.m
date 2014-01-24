function q = integrateLVM_LW(F, p, c, psi0)

% m==0 is Kelvin Mode
for m=0:p.maxMermodes;
    %Make BC coefficients
    eref(m+1) = erefcof(m);
    wref(m+1) = wrefcof(m);

    %Make Phase Speed
    cn(m+1) = c/(2*(m) + 1);
    if m>=1
       cn(m+1) = -cn(m+1); % Rossby Wave
    end

   
end

%Make Damping
dampn = p.damp*(2.6/c).^2; % XXX --- Need to come back to this.

n = length(p.lons);                            % number of space gridpoints
dt = (p.time(2)-p.time(1)).*p.timetosec;       % time step
% tf = p.time(end).*p.timetosec;                 % final time

x = p.lons.*p.deg;                             % space in meters
h = x(2)-x(1);                                 %dx
r = cn(1).*dt/h; disp(sprintf('Courant number: %0.2f',r)) %Must be < 1 for stability
u = zeros(size(x));                            % IC
% u0 = u;

%Matrix
I = eye(n);
R = diag(ones(1,n-1),1);
Dxc = (R+diag([-1 zeros(1,n-2) 1])-R')/(2*h); % centered diff in x
Dxx = (R-diag([1 2*ones(1,n-2) 1])+R')/h^2;

uall = zeros([p.maxMermodes+1 length(u)]);          % All waves
uforc  = uall;                                      % Only Forced waves
urefl = uforc;                                      % Only reflected waves
urere = uforc;                                      % Only 2x relfected waves

%Preallocate
q = NaN * zeros([4 size(uall) length(p.time)]);      % Total Wave Field



for tn = 1:length(p.time)                   % TIME LOOP
    if (mod(tn, p.tdisps) == 0)
        display(['Integrating Timestep: ', num2str(tn), '/', num2str(length(p.time))]);
    end
    
    for m = 0:p.maxMermodes;                % MERIDIONAL MODE LOOP
        
        % Enforce BC
        if (m==0)
            % West BC
            wbcforced = 0;
            wbcall = 0;
            for mm = 0:p.maxMermodes
                wbcforced = wbcforced + wref(mm+1)*uforc(mm+1,1);
                wbcall = wbcall + wref(mm+1)*uall(mm+1, 1);
            end
            uall(m+1,1) = wbcall;   % Assign BC to Kelvin Mode
            uforc(m+1, 1) = 0;      % No reflection if forced only, note this is at the upwind boundary
            urefl(m+1, 1) = wbcforced; % Reflection only of directly forced
            urere(m+1, 1) = wbcall; % Reflected plus reflected
        else
            % East BC
            urefl(m+1, end) = eref(m+1).*uforc(1, end);
            urere(m+1, end) = eref(m+1).*uall(1,end); %XXX Look at this, maybe should come from urefl
            uforc(m+1, end) = 0;
            uall(m+1,end) = eref(m+1).*uall(1,end); % East BC, assigned to Rossby Mode

        end

%        % Lax Wendroff
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        % ut = u - c*dt*Dx(u) + c^2*dt^2/2*Dxx(u) + dt/2*(F(t+1) + F) -
%        %                                                c^2*dt^2/2 * Dx(F)
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        
%        
%        A = -cn(m+1).*Dxc+(cn(m+1).^2.*dt/2)*Dxx;  % A matrix
%        D = dampn.*I;                         % Damping Matrix
%        M = I+dt*(A - D);                          % Total Differential Equation
%                                                   % Make forcing matrix
%        if (tn~=length(p.time))
%            Fm = 1/2*(squeeze(F(m+1,:,tn+1)) + squeeze(F(m+1,:,tn)))'...
%                -cn(m+1).^2*dt/2*Dxx*(squeeze(F(m+1,:,tn))');
%        else
%            Fm = squeeze(F(m+1,:,tn))'...
%                -cn(m+1).^2*dt/2*Dxx*(squeeze(F(m+1,:,tn))');
%        end   
       
       % Upwind
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % XXX Testing for model comparison
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       if m==0
           R = diag(ones(1,n-1), -1);
       else
           R = diag(ones(1,n-1), 1);
       end
       A = abs(cn(m+1)).*(R - I)./h;
       D = dampn.*I;
       M = I + dt*(A - D);
       Fm = squeeze(F(m+1,:,tn))';
       

       
       Fn = Fm.*psi0./(p.rho.*p.H); % Forcing, dimensionalized
       
       uall(m+1,:) = M*squeeze(uall(m+1,:))'+dt.*Fn; % Time step forward.
       uforc(m+1,:) = M*squeeze(uforc(m+1,:))'+dt.*Fn; 
       urefl(m+1,:) = M*squeeze(urefl(m+1,:))'; % No forcing for reflected waves
       urere(m+1,:) = M*squeeze(urere(m+1,:))';
       
    end

       %Assign output
       q(1,:,:,tn) = uall;
       q(2,:,:,tn) = uforc;
       q(3,:,:,tn) = urefl;
       q(4,:,:,tn) = urere - urefl; %Twice reflected is Reflected - (Reflected once)
end


end


% Generates the reflection coefficients for reflection of Kelvin Waves
function eref = erefcof(m)
    if(m==1)
        eref=sqrt(3/2);
    elseif (mod(m, 2)==0)
        eref = 0;
    else
        mm = (m-1)/2;
        tmp1 = sqrt (factorial(2*mm-1)*(4*mm+3));
        tmp2 = sqrt ( 2^(2*mm) *factorial(mm-1)*factorial(mm+1));
        eref = real(tmp1/tmp2);
        
    end
end

%Generates the reflection coefficients for reflection of Rossby Waves
function wref = wrefcof(m)

   if (mod(m, 2)==0)
        wref = 0;
    else
        mm = (m-1)/2;
        tmp1 = sqrt (2*(mm + 1)*factorial(2*mm));
        tmp2 = sqrt ( 4*mm + 3)*(2^(mm+1))*factorial(mm+1);
        wref = real(tmp1/tmp2);
        
    end
end