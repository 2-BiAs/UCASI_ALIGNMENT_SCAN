clear
clc

printf('\nBegin Initialization\n');

printf('\nBuilding Libraries\n');

exec('bin\BuildLibs.sci');
BuildLibs();

printf('\nLoading Libraries\n');

exec('bin\LoadLibs.sci');
listLibs = LoadLibs();

//Setting Libraries to variables
for i=1:2:size(listLibs)
    sTemp = sprintf('%s = listLibs(%d)', listLibs(i+1), i);
    mprintf('\n\t%s\n' , listLibs(i+1));
    execstr(sTemp);
end

printf('\nInitialization Complete\n');
