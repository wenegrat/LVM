function qvec = makeqvec(q, p, force, refl, rere)

%Make the q vector depending on options
qvec = zeros([p.maxVermodes p.maxMermodes+1 length(p.lons).*length(p.time)]);

if (force)
    % Reshape the q vectors           
    qvecforc = reshape(squeeze(q(2,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
    qvec = qvec + qvecforc;
end
if (refl)
    qvecrefl = reshape(squeeze(q(3,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
    qvecrefl(:,1,:) = qvecrefl(:, 1,:).*p.wref; % Apply western BC reflectivity to Kelvin mode
    qvecrefl(:,2:end,:) = qvecrefl(:,2:end,:).*p.eref; % Apply Eastern BC reflectivity to all Rossby modes;
    qvec = qvec + qvecrefl;
end
if (rere)
    qvecrere = reshape(squeeze(q(4,:,:,:,:)), p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
    qvecrere = qvecrere.*p.wref.*p.eref; % Apply both BCs to the reflect + reflect;
    qvec = qvec + qvecrere;
end

end