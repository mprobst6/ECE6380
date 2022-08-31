#!/bin/bash
for resolution in 5 10 15 20 30 100
do
    julia mesh1D.jl $resolution
    julia fem.jl $resolution
done