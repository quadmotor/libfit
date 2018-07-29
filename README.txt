
o888  o88   oooo       ooooooooooo o88    o8  
 888  oooo   888ooooo   888    88  oooo o888oo
 888   888   888    888 888ooo8     888  888  
 888   888   888    888 888         888  888  
o888o o888o o888ooo88  o888o       o888o  888o

Created: June 2017
Author: Cameron Mackintosh
Contact: cmackinto at posteo daht net

Description:
	The libFit library parameterically represents structured meshes 
	using Bezier surfaces, and generates arbitrary new structured meshes 
	from these parametric representations. 
	

File Hierarchy: 
	build/: object files created during compilation
	include/: module files created during compilation
	lib/: libfit library created during compilation
	src/: libfit source files
	examples/: example programs using libfit

Future Development:
	- Structured meshes are assumed to have rectilinear block connectivity. 
	  This should be changed to support arbitrary connectivity.
	- Interfacing blocks is done by solving for control points at the 
	  boundary then solving for control points on the blocks, subtracting 
	  out the known control points. When generating new meshes, first the 
	  blocks’ xyz values are calculated, then the boundaries’ xyz values are 
	  calculated, then the boundary values replace the block values at 
	  boundary locations. 
	- The process of calculating boundaries’ xyz values DOES NOT use arc 
	  length tables, merely using the uv values of one of the two adjacent 
	  blocks.
