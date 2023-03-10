var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MasterChemicalMechanism","category":"page"},{"location":"#MasterChemicalMechanism","page":"Home","title":"MasterChemicalMechanism","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MasterChemicalMechanism.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MasterChemicalMechanism]","category":"page"},{"location":"#MasterChemicalMechanism.parse_rxns-Tuple{Any}","page":"Home","title":"MasterChemicalMechanism.parse_rxns","text":"function parse_rxns(rxns)\n\nGiven a list of equations rxns, parse equation returning a list of unique species in the entire mechanism and a list of reaction components, i.e. tuples of (reactants, reactants_stoich, products, products_stoich, rate_equation)\n\n\n\n\n\n","category":"method"},{"location":"#MasterChemicalMechanism.read_fac_file-Tuple{String}","page":"Home","title":"MasterChemicalMechanism.read_fac_file","text":"function read_fac_file(path::String)\n\nParse a mechanism file in the FACSIMILE format (ending with .fac). This assumes that the file has all necessary sections to define the master chemical mechanism including\n\nVOCs: a list of VOCs tracked in the simulation\ngeneric_rate_coefficients\ncomplex_rate_coefficients\nperoxy_radicals: the species being summed to form RO2\nreaction_definitions: The list of chemical reactions in our mechanism\n\nOutput\n\nReturns a dictionary containing the above sections as text.\n\n\n\n\n\n","category":"method"},{"location":"#MasterChemicalMechanism.species_and_stoich-Tuple{Any}","page":"Home","title":"MasterChemicalMechanism.species_and_stoich","text":"function species_and_stoich(rxn_half)\n\nGiven half of a reaction (either reactants or products), parse the equation and return a list of species involved and associated stoichiometric coefficients. If no species is involved (production or destruction), set the species to nothing.\n\n\n\n\n\n","category":"method"}]
}