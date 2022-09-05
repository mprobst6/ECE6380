rm *.png
rm -rf outputfiles
rm -rf inputfiles

julia variablemesh.jl 
julia femtot.jl 
