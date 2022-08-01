function GenerateDoc(cliArgs, Args)
    arguments
        cliArgs     string  = ""
        Args.root       string  = pwd
        Args.venv       string  = fullfile(pwd,"../../venv")
        Args.docPath    string
        Args.dataPath   string
    end
    if isfield(Args,'docPath')
        cliArgs = strcat('--doc ',Args.docPath,' ',cliArgs);
    end
    if isfield(Args,'dataPath')
        cliArgs = strcat('--data ',Args.dataPath,' ',cliArgs);
    end
    try
        if pyenv().Status ~= "Loaded"
            pyenv(Version=fullfile(Args.venv,'bin','python'));
        end
        pyrunfile(sprintf('GenerateDoc.py %s %s',Args.root,cliArgs));
    catch err
        rethrow(err)
    end
end