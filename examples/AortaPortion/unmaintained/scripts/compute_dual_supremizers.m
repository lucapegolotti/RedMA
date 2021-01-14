function compute_dual_supremizers(fieldIndex,basis,matricesDir,meshname,outDir, format)
normMatrix = spconvert(load([matricesDir,'/',meshname,'/norm',num2str(fieldIndex),'.m']));

suprdual = [];
count = 1;
while (exist([matricesDir,'/',meshname,'/dualConstraint',num2str(count),'.m'],'file'))
    curConstraint = spconvert(load([matricesDir,'/',meshname,'/dualConstraint',num2str(count),'.m']));
    % pad constraintMatrices
    curConstraint = [curConstraint; ...
                     sparse(size(basis,1)-size(curConstraint,1), ...
                     size(curConstraint,2))];
                 
    suprdual = [suprdual normMatrix \ curConstraint];

    count = count + 1;
end
size(suprdual)
for i = 1:size(suprdual,2)
    disp(['Normalizing supremizer ',num2str(i)])
    % normalize
    for j = 1:size(basis,2)
        suprdual(:,i) = suprdual(:,i) - (suprdual(:,i)' * basis(:,j))/norm(basis(:,j)) * basis(:,j);
    end

    for j = 1:(i-1)
        suprdual(:,i) = suprdual(:,i) - (suprdual(:,i)' * suprdual(:,j))/norm(suprdual(:,j)) * suprdual(:,j);
    end

    suprdual(:,i) = suprdual(:,i) / norm(suprdual(:,i));
end

suprdual(abs(suprdual) < 1e-15) = 0;

dlmwrite([outDir,'/dual_supremizers',num2str(fieldIndex),'.basis'],full(suprdual'),'delimiter', ',','precision', format);