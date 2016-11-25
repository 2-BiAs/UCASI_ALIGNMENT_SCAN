function listOutput = LoadLibs()

    [output, bOK] = dos("FOR /D %I IN (BIN\*) DO @ECHO %~nxI");
    output = stripblanks(output);
    listLibs = list();
    for i=1:size(output,1)
        try
            listLibs($+1) = lib("BIN\" + output(i));
            listLibs($+1) = output(i);
        catch
            printf('WARNING: lib: %s library cannot be loaded\n', output(i));
        end
    end
    listOutput = listLibs;
endfunction
