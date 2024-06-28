# plotlyMol

A package to create molecular representations using the base Python Plotly. 

consider the ability to be able to pass multiple inputs of the same time... [smiles, smiles, etc] <- need to handle single inputs and lists of inputs 

For the future:

- 2D: give chemdraw like structures
- 3D: 
	- draw 3D structures 
		- add in the ability to scale the position of the half-bond by the van der waals radii of the atoms at either end. 
	- add orbitals and other surfaces
	- for smiles: need to be able to handle add structures
	- PROBLEM: the rdkit xyzblock --> mol object can fail, for things like NITRO.  Need to find a new way to handle xyz (maybe bable xyz --> smiles???)
	- Need to add in the ability to draw orbitals.  Code works, just needs to be implemented. 
- 4D: Animations
- GUI
	- ability to toggle aspects of the presentation/hover/etc
	- ability to dynamically enter mol information
