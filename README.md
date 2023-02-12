# MasterChemicalMechanism

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://john-waczak.github.io/MasterChemicalMechanism.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://john-waczak.github.io/MasterChemicalMechanism.jl/dev/)
[![Build Status](https://github.com/john-waczak/MasterChemicalMechanism.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/john-waczak/MasterChemicalMechanism.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/john-waczak/MasterChemicalMechanism.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/john-waczak/MasterChemicalMechanism.jl)


A package for ingesting models from the [Master Chemical Mechanism](http://chmlin9.leeds.ac.uk/MCM/project.htt) and producing Coupled ODEs to simulate gas phase chemical kinetics. 

# Steps 

1. Define config parameters 
2. Parse mechanism file 
3. Generate R.H.S. reaction rate functions 
4. Construct model via `Catalyst.jl`
5. Convert to coupled ODEs and structural simplify 
6. Simulate 
7. Analyze Results 
8. Combine with data to perform *Data Assimilation*, i.e. fitting the model to real data.
