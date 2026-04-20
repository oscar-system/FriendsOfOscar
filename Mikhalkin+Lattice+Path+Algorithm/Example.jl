include("RationalPlaneTropicalCurves.jl")

#Fix the degree d and then determine the valid paths using all_paths(), where we can also see their (non-zero) multiplicities.
d=3;
AP=all_paths(d);

#Build binary trees coming from Mikhalkin's lattice path algorithm using l_build_tree (for left turning) and r_build_tree (for right turnings).
#The initial path should be one of multiplicity >0. 
initial_path=AP[3]
root_tree_upper=l_build_tree(initial_path,d); 
root_tree_lower=r_build_tree(initial_path,d); 

#Collect paths in the above binary trees that start in initial_path and end in a leaf. Each pair (upper, lower) of those paths (of paths) produces a Newton subdivision (and a tropical curve) associated to the initial_path. Summing up their multiplicity (as subdivisions) we obtain the multiplicity of the initial path
allpaths_upper=collect_paths(root_tree_upper);
allpaths_lower=collect_paths(root_tree_lower);

#Prints the paths from the root to the leaves in the binary tree.
print_paths(allpaths_upper)
print_paths(allpaths_lower)

#Shows the size of allpaths_upper and the size of allpaths_lower, to determine which numbers can be plugged into in the visulisation functions (e.g. we can combine the partial subdivision in the first position in allpaths_upper with the second partial subdivision in allpaths_lower).
length(allpaths_upper) 
length(allpaths_lower)

#Note! The function visualize_curve might give you curves with cycles, but those are fake cycles coming from the parallelograms in the Newton subdivision. Note that the number of interior points minus the number of parallelograms in each corresponding Newton subdivision is zero, as expected.
#The first input integer corresponds to the position in allpaths_upper, the second input integer to the position in allpaths_lower.
#Now we can visualize the Newton Subdivision as well as an exemplary corresponding tropical curve.
#Note that the visualized curve does not necessarily pass through points in tropical general position.
visualize_subdivision(1,2,d)
visualize_curve(1,2,d)
