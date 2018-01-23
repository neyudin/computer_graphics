#!/bin/bash

pictures='../pictures'
results='../results'

for x in `ls "$pictures"`
do
    if [ ${x:(-4):4} == ".bmp" ]
    then
        name=${x::-4};
### Align
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.bmp" --align;
### Postprocessing
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.gray.bmp" --gray-world
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.sharp.bmp" --unsharp
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.contr.bmp" --autocontrast 0.3
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.median.bmp" --median-linear 1
### Align and postprocessing
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.gray.bmp" --align --gray-world;
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.sharp.bmp" --align --unsharp;
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.contr.01.bmp" --align --autocontrast 0.1;
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.contr.03.bmp" --align --autocontrast 0.3;
### Mirroring
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.mirror.bmp" --align --mirror;
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.gray.mirror.bmp" --align --gray-world --mirror;
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.sharp.mirror.bmp" --align --unsharp --mirror;
        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.contr.01.mirror.bmp" --align --autocontrast 0.1 --mirror;
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.contr.03.mirror.bmp" --align --autocontrast 0.3 --mirror;
### Resize
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.resize.03.bmp" --resize 0.3;
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.resize.07.bmp" --resize 0.7;
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.resize.19.bmp" --resize 1.9;
#        ./build/bin/align "$results/$name.align.bmp" "$results/$name.resize.24.bmp" --resize 2.4;
### Subpixel
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.subpix.10.bmp" --align --subpixel 1.0;
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.subpix.20.bmp" --align --subpixel 2.0;
#        ./build/bin/align "$pictures/$name.bmp" "$results/$name.align.subpix.30.bmp" --align --subpixel 3.0;
    fi
done
