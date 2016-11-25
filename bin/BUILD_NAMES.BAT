@ECHO OFF
Setlocal EnableDelayedExpansion
FOR /D %%I IN (*) DO (	
	IF EXIST %%I\names DEL %%I\names
	SET log=%%I\names
	FOR %%J IN (%%I\*.sci) DO (
		ECHO %%~nJ >> "!log!"
		@ECHO %%~nJ))