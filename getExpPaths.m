function expPaths = getExpPaths(subject)

% takes in a base directory and subject string and returns all relevant
% directories for that subject 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseDir = '/home/kelly/ShockAwe/data';

expPaths = struct();

expPaths.baseDir = baseDir;
expPaths.subj        = fullfile(baseDir, subject);  % subject directory

% subject directories
expPaths.afni        = fullfile(expPaths.subj, 'afni/');
expPaths.design_mats = fullfile(expPaths.subj, 'design_mats/');
expPaths.dti96trilin = fullfile(expPaths.subj, 'dti96trilin/');
expPaths.fibers      = fullfile(expPaths.subj, 'fibers/');
expPaths.raw         = fullfile(expPaths.subj, 'raw/');
expPaths.raw_func    = fullfile(expPaths.subj, 'raw_func/');
expPaths.raw_dti     = fullfile(expPaths.subj, 'raw_dti/');
expPaths.regs        = fullfile(expPaths.subj, 'regs/');
expPaths.results     = fullfile(expPaths.subj, 'results/');
expPaths.results_hab = fullfile(expPaths.subj, 'results_hab/');
expPaths.ROIs        = fullfile(expPaths.subj, 'ROIs/');
expPaths.stimtimes   = fullfile(expPaths.subj, 'stimtimes/');
expPaths.slicetimes  = fullfile(expPaths.subj, 'slicetimes/');
expPaths.funcROIs    = fullfile(expPaths.subj, 'funcROIs/');
expPaths.t1          = fullfile(expPaths.subj, 't1/');

end


