function [Uu] = compute_basis_field(fieldIndex, snapshotsDir, dirName, outDir, tol, format, varargin)

maxSnapshots = Inf;

if (nargin == 7)
    maxSnapshots = varargin{1};
end

% useNorm = 0;

count = 0;
snap = [];

while (exist([snapshotsDir,'/param',num2str(count)],'dir') && count < maxSnapshots)
    tentativeFile = [snapshotsDir,'/param',num2str(count),'/',dirName,'/field',num2str(fieldIndex),'.snap'];
    
    if (exist(tentativeFile,'file'))   
        isEmpty = 0;
        % check if it is empty 
        fid = fopen(tentativeFile);
        if all(fgetl(fid) == -1)
            isEmpty = 1;
        else
            fseek(fid,0,-1); % rewind it
        end

        if (~isEmpty)
            disp(tentativeFile)
            snap = [snap;csvread(tentativeFile)];
        end
    end
    count = count + 1;
end

% if (useNorm) 
%     normVelocity = spconvert(load('basis/tube_1x3_h0.08/norm0.m'));
%     H = chol(normVelocity);
% end

snap = snap';

disp(['Performing svd for ',dirName,' field ',num2str(fieldIndex),' ...'])
disp(['n snapshots = ',num2str(size(snap,2)),' ...'])
[U,S,~] = svd(snap,'econ');

% if (useNorm)
%     U = H \ U;
% end

totalenergy = sum(diag(S).^2);

partialSum = 0;
Nu = 1;
while (tol * tol < 1.0 - partialSum / totalenergy)
    partialSum = partialSum + S(Nu,Nu) * S(Nu,Nu);
    Nu = Nu + 1;
end
Nu = Nu - 1;

disp(['selected vectors = ',num2str(Nu),' ...'])
Uu = U(:,1:Nu);

Uu(abs(Uu) < 1e-15) = 0;
disp(['Dumping result for ',dirName,' field ',num2str(fieldIndex),' ...'])
dlmwrite([outDir,'/field',num2str(fieldIndex),'.basis'],Uu','delimiter', ',','precision', format);
dlmwrite([outDir,'/svd',num2str(fieldIndex),'.txt'],diag(S),'delimiter', ',','precision', '%.16g');