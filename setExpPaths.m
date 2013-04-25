function setExpPaths(subject)

% gets relevant subject directories as specified in getExpPaths and creates 
% them if they don't already exist. If directory already does exist, matlab 
% gives a warning to say so and doesn't overwrite it.  At least this is the 
% case on Rexy...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expPaths = getExpPaths(subject);

dirNames = fieldnames(expPaths);

for d = 2:numel(dirNames)   % 1st directory is base dir (already exists)
    mkdir(expPaths.(dirNames{d}));
end

fprintf(['\n\nmade exp directories for subject ',subject,'\n\n']);

end
