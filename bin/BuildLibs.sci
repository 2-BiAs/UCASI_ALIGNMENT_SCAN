function BuildLibs()
    dos("\BIN\BUILD_NAMES.BAT");
    [output, bOK] = dos("FOR /D %I IN (BIN\*) DO @ECHO %~nxI");
    output = stripblanks(output);
    for i=1:size(output,1)
        genlib(output(i),path="BIN\" + output(i));
    end
endfunction
