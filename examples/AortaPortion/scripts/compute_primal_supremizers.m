function [basisAugmented] = compute_primal_supremizers(fieldIndex_augmented,fieldIndex_constraint,basisAugmented,basisConstraint,matricesDir,meshname,outDir,format)
constraintMatrix = spconvert(load([matricesDir,'/',meshname,'/primalConstraint.m']));
normMatrix = spconvert(load([matricesDir,'/',meshname,'/norm',num2str(fieldIndex_augmented),'.m']));

suprprimal = normMatrix \ (constraintMatrix * basisConstraint);

N0 = size(basisAugmented,2);
N1 = size(basisConstraint,2);

for i = 1:N1
    disp(['Normalizing supremizer ',num2str(i)])
    % normalize
    for j = 1:N0
        suprprimal(:,i) = suprprimal(:,i) - (suprprimal(:,i)' * basisAugmented(:,j))/norm(basisAugmented(:,j)) * basisAugmented(:,j);
    end

    for j = 1:(i-1)
        suprprimal(:,i) = suprprimal(:,i) - (suprprimal(:,i)' * suprprimal(:,j))/norm(suprprimal(:,j)) * suprprimal(:,j);
    end

    suprprimal(:,i) = suprprimal(:,i) / norm(suprprimal(:,i));
end

suprprimal(abs(suprprimal) < 1e-15) = 0;

basisAugmented = [basisAugmented suprprimal];

dlmwrite([outDir,'/primal_supremizers_',num2str(fieldIndex_augmented),'_',num2str(fieldIndex_constraint),'.basis'],suprprimal','delimiter', ',','precision', format);