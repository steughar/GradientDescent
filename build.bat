@echo off

set CommonLinkerFlags= -incremental:yes -opt:ref
IF NOT EXIST ..\build mkdir ..\build
pushd ..\build
del *.pdb > NUL 2> NUL

cl -Z7 -nologo -Fmeval_circle.map ..\code\Source.c ..\code\eval_circle.c

popd
