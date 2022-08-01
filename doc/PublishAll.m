function PublishAll(pkgRoot,docRoot,Args)
    arguments
        pkgRoot = fullfile(pwd, "../")
        docRoot = pwd
        Args.htmlPath = fullfile(docRoot, "html")
        Args.styleSheet = fullfile(docRoot, "StyleSheet.xsl")
    end
    % Execute `publish` in all subfolders of 

    addpath(pkgRoot);
    try
        [root, files, folders] = getFilesFolders(docRoot);
        [files, folders] = filterFilesFolders(files, folders);
        walk(@PublishFiles, root, files, folders);
        rmpath(pkgRoot);
    catch err
        rmpath(pkgRoot);
        rethrow(err);
    end
    function PublishFiles(root,files,~)
        % Execute `publish` on `files`

        relPath = erase(root,docRoot);
        curPath = pwd;
        cd(root);
        try
            files = filterMFiles(files);
            for file = files(:)'
                [~,fileName] = fileparts(file);
                publish(file, format="html",...
                    outputDir=fullfile(Args.htmlPath, relPath),...
                    stylesheet=Args.styleSheet);
                if file == "consturctor.m"
                    % TODO: Add symbolic link for convenience
                end
                fprintf("Generated %s from %s\n",...
                    fullfile(relPath, fileName + ".html"), file);
            end
            cd(curPath);
        catch err2
            cd(curPath);
            rethrow(err2);
        end
    end
end
function walk(func,root,files,folders)
    % Execute `func` recursively over `folders`

    func(root,files,folders);
    for nextFolder = folders(:)'
        [curRoot,curFiles,curFolders] = getFilesFolders(fullfile(root,nextFolder));
        walk(func,curRoot,curFiles,curFolders);
    end
end
function [root,files,folders] = getFilesFolders(path)
    % Extract files and folders in `path` similar to python's os.walk

    files = dir(path);
    % First object always exist and can get absolute path
    root = files(1).folder;
    % Filter current and parent folders
    files = files(~({files(:).name} == "." | {files(:).name} == ".."));
    % Split files and folders
    folders = files([files.isdir]);
    folders = string({folders.name});
    folders = folders(:);
    files = files(~[files.isdir]);
    files = string({files.name});
    files = files(:);
end
function [mFiles,otherFiles] = filterMFiles(files)
    % Filter the .m files

    [~,~,ext] = fileparts(files);
    mFiles = files(ext == ".m");
    otherFiles = files(ext ~= ".m");
end
function [files,folders] = filterFilesFolders(files,folders)
    % Filter the relevant documentation files

    files = filterMFiles(files);
    excludeFolders = ["html", "template"];
    folders = folders(~any(folders == excludeFolders, 2));
    excludeFiles = ["PublishAll.m"];
    files = files(~any(files == excludeFiles, 2));
end