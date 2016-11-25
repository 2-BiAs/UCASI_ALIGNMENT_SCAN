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
for i=1:2:size(listLibs) / 2
    execstr(listLibs(i+1) + ' = listLibs(i)');
end

printf('\nInitialization Complete\n');
