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
    dampn = p.damp*(2.5/c).^2; % XXX --- Need to come back to this.

n = length(p.lons);                            % number of space gridpoints
dt = (p.time(2)-p.time(1)).*p.timetosec;       % time step
tf = p.time(end).*p.timetosec;                 % final time

x = p.lons.*p.deg;                             % space in meters
h = x(2)-x(1);                                 %dx
r = cn(1).*dt/h; disp(sprintf('Courant number: %0.2f',r)) %Must be < 1 for stability
u = zeros(size(x));                            % IC
u0 = u;

%Matrix
I = eye(n);
R = diag(ones(1,n-1),1);
Dxc = (R+diag([-1 zeros(1,n-2) 1])-R')/(2*h); % c*centered diff in x
Dxx = (R-diag([1 2*ones(1,n-2) 1])+R')/h^2;

ur = zeros([p.maxMermodes+1 length(u)]);       %Assign IC
q = NaN * zeros([size(ur) length(p.time)]);    %Preallocate

for tn = 1:length(p.time)                   % TIME LOOOP
    for m = 0:p.maxMermodes;                % MERIDIONAL MODE LOOP
        
        % Enforce BC
        if (m==0)
            % West BC
            wbc = 0;
            for mm = 0:p.maxMermodes
                wbc = wbc + wref(mm+1)*ur(mm+1,1);
            end
            ur(m+1,1) = p.wref.*wbc; % Assign BC to Kelvin Mode
        else
            ur(m+1,end) = p.eref.*eref(m+1).*ur(1,end); % East BC, assigned to Rossby Mode
        end

       % Lax Wendroff
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % ut = u - c*dt*Dx(u) + c^2*dt^2/2*Dxx(u) + dt/2*(F(t+1) + F) -
       %                                                c^2*dt^2/2 * Dx(F)
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
       A = -cn(m+1).*Dxc+(cn(m+1).^2.*dt/2)*Dxx;  % A matrix
       D = dampn.*I;                         % Damping Matrix
       M = I+dt*(A - D);                          % Total Differential Equation
                                                  % Make forcing matrix
       if (tn~=length(p.time))
           Fm = 1/2*(squeeze(F(m+1,:,tn+1)) + squeeze(F(m+1,:,tn)))'...
               -cn(m+1).^2*dt/2*Dxx*(squeeze(F(m+1,:,tn))');
       else
           Fm = squeeze(F(m+1,:,tn))'...
               -cn(m+1).^2*dt/2*Dxx*(squeeze(F(m+1,:,tn))');
       end
       
       Fn = Fm.*psi0./(p.rho.*p.H); % Forcing, dimensionalized
       
       ur(m+1,:) = M*squeeze(ur(m+1,:))'+dt.*Fn; % Time step forward.
       
    end
       %Plotting (temporary)
%        clf
%        subplot(3,1,1)
%        plot(x./p.deg, Fn, 'k');
%        ylim([-300 300])
%        xlim([x(1)./p.deg x(end)./p.deg]);
%        subplot(3,1,2:3)
%        plot(x./p.deg,u0,'b:',x./p.deg,ur(1,:),'r.-')
%        hold on
%               plot(x./p.deg,ur(2,:),'b.-')
%               plot(x./p.deg, ur(4,:), 'g');
%         hold off
%     %    axis([-1 1 -.1 1.1])
%        xlim([x(1)./p.deg x(end)./p.deg]);
%        ylim([-.5*1e8 .5*1e8]);
%        title(sprintf('%s , t=%0.2f','Lax-Wendroff', tn))
%        drawnow
%        
       %Assign output
       q(:,:,tn) = ur;
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