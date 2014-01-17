function [wmodes pnormal ce] = calcVertModes(ds, p)

% Vert modes

eqlon = 22;
eqlat = 11;
nsq = nanmean(squeeze(nanmean(squeeze(ds.argo.nsqr(eqlon,eqlat-5:eqlat+5,:,:)))),2);
nsq = [nsq; nsq(end)];
ts = nanmean(squeeze(nanmean(squeeze(ds.argo.temp(eqlon,eqlat-5:eqlat+5,:,:)))),2);
interpdepths = 2:10:3300;
nsqi = interp1([ds.argo.nsqrp; 3300], nsq, interpdepths);

% [wmodes, pmodes, ce] = dynmodes(nsq, [ds.argo.nsqrp; 4000], 10);
[wmodes, pmodes, ce] = dynmodes(nsqi, interpdepths, p.maxVermodes);

for i=1:p.maxVermodes
    pnorm(:,i) = trapz([0; interpdepths'], pmodes(:,i).^2);
%     pnormal(:,i) = pmodes(:,i)/(pnorm(:,i)); %All modes normalized to integral one
    pnormal(:,i) = pmodes(:,i)./std(pmodes(:,i));
    if (pnormal(1,i) <0)
        pnormal(:,i) = -pnormal(:,i);
    end
end
end