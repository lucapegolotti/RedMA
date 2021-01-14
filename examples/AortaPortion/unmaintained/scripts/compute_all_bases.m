clear
clc

outDir = '../basis';
snapshotsDir = '../snapshots';
matricesDir = '../matricesForOffline';

format = '%0.16g';

mkdir(outDir);

subdirs = dir([snapshotsDir,'/param0']);

% param 0 must exist and be complete
for file = dir([snapshotsDir,'/param0/'])'
    if (file.isdir && ~strcmp(file.name,'.') && ~strcmp(file.name, '..'))
        curOutDir = [outDir,'/',file.name];
        mkdir(curOutDir)
        U0 = compute_basis_field(0,snapshotsDir,file.name,curOutDir,1e-5,format);
        U1 = compute_basis_field(1,snapshotsDir,file.name,curOutDir,1e-5,format);
        
        U0 = compute_primal_supremizers(0,1,U0,U1,matricesDir,file.name,curOutDir,format);
        compute_dual_supremizers(0,U0,matricesDir,file.name,curOutDir,format);
    end
end
