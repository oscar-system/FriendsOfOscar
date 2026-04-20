using Oscar
#This file is structured as follows: We start with all functions needed to implement the recursive part of the algorithm starting from a lambda-increasing lattice path to the two boundary paths D1 and D2. The recursive computation of the multiplicity follows.
#The next part of the code filters out the paths that contribute with positive multiplicity.
#In order to build the corresponding Newton subdivisions and associated tropical curves, the next part of this file will introduce a tree structure that captures the recursive process for each initial lambda-increasing lattice path.
#The final section of the code builds and visualizes the Newton subdivisions and associated tropical curves.

#####
#In this section of the code, we determine when a path turns to the left/right for the first time and cut and fold accordingly. This is required for the calculation of multiplicity.

#The input is a matrix, whose columns should be thought of as vectors, and an integer k. The function calculates the normalized area of the polygon that has as vertices the k-1,k,k+1 columns of the input matrix.
function area(M::ZZMatrix,k::Integer)::Integer
    E=zero_matrix(ZZ,2,2)
    if k<2 || k>size(M,1)-1
        print("Impossible calculation")
        return nothing
    else
        for j in 1:2
            E[1,j]=M[k,j]-M[k-1,j]
            E[2,j]=M[k+1,j]-M[k,j]
        end
    end 
    return E[1,1]*E[2,2]-E[1,2]*E[2,1]
end

#The input is a matrix representing a path in the plane. The function determines where the input path turns for the first time to the left.
function left_turn(M::ZZMatrix)::Integer
    l=0 
    for i in 2:size(M,1)-1 
        if area(M,i)>0
            l=i
            break
        end
    end
    return l
end

#The input is a matrix representing a path in the plane. The function determines where the input path turns for the first time to the right.
function right_turn(M::ZZMatrix)::Integer
    r=0 
    for i in 2:size(M,1)-1 
        if area(M,i)<0
            r=i
            break
        end
    end
    return r
end

#The input is a matrix representing a path in the plane. The function cuts the vertex corresponding to the first left turn of the path and returns the new path. If cutting is not possible, the function returns nothing.
function lcut(M::ZZMatrix)::Union{ZZMatrix,Nothing}
    N=zero_matrix(ZZ,size(M,1)-1,2)
    if left_turn(M)<1
        println("No turns to the left")
        return nothing
    end
    for i in 1:size(N,1)
        if i<left_turn(M) 
            N[i,1]=M[i,1]
            N[i,2]=M[i,2]
        else
            N[i,1]=M[i+1,1]
            N[i,2]=M[i+1,2]
        end
    end
    return N
end

#The input is a matrix representing a path in the plane. The function cuts the vertex corresponding to the first right turn of the path and returns the new path. If cutting is not possible, the function returns nothing.
function rcut(M::ZZMatrix)::Union{ZZMatrix,Nothing}
    M1=zero_matrix(ZZ,size(M,1)-1,2)
    if right_turn(M)<1
        println("No turns to the right")
        return nothing
    end
    for i in 1:size(M1,1)
        if i<right_turn(M) 
            M1[i,1]=M[i,1]
            M1[i,2]=M[i,2]
        else
            M1[i,1]=M[i+1,1]
            M1[i,2]=M[i+1,2]
        end
    end
    return M1
end

#The input is a matrix representing a path in the plane. The function finds the first left turn, say column l, which is a point in the plane. Then, it reflects that point along the line containing the points represented by the (l-1) and (l+1) colums of the input path. This procedure is called folding. The function returns the new path. If folding is not possible because the path would leave the given Newton polygon, the function returns nothing.
function lfold(M::ZZMatrix,d::Integer)::Union{ZZMatrix,Nothing}
    M2=zero_matrix(ZZ,size(M,1),2)
    for i in 1:size(M2,1)
        if i!=left_turn(M) 
            M2[i,1]=M[i,1]
            M2[i,2]=M[i,2]
        else
            M2[i,1]=M[i-1,1]+M[i+1,1]-M[i,1]
            M2[i,2]=M[i-1,2]+M[i+1,2]-M[i,2]
            if M2[i,1]<0 || M2[i,2]<0 || M2[i,2]>-M2[i,1]+d
                return nothing
            end
        end
    end
    return M2
end

#The input is a matrix representing a path in the plane. The function finds the first right turn, say column l, which is a point in the plane. Then, it reflects that point along the line containing the points represented by the (l-1) and (l+1) colums of the input path. This procedure is called folding. The function returns the new path. If folding is not possible because the path would leave the given Newton polygon, the function returns nothing.
function rfold(M::ZZMatrix,d::Integer)::Union{ZZMatrix,Nothing}
    M2=zero_matrix(ZZ,size(M,1),2)
    for i in 1:size(M2,1)
        if i!=right_turn(M) 
            M2[i,1]=M[i,1]
            M2[i,2]=M[i,2]
        else
            M2[i,1]=M[i-1,1]+M[i+1,1]-M[i,1]
            M2[i,2]=M[i-1,2]+M[i+1,2]-M[i,2]
            if M2[i,1]<0 || M2[i,2]<0 || M2[i,2]>-M2[i,1]+d
                return nothing
            end
        end
    end
    return M2
end

#######
#In this section of the code, we recursively calculate the multiplicity of a given path.

#The input is a matrix representing a path in the plane. The function determines whether the input path is one of the two extremal paths in the boundary of our triangle connecting the points where the linear map lambda obtains its maximal and minimal value.
function D1orD2(M::Union{ZZMatrix,Nothing},d::Int64)::Integer
    #We first build the boundary paths
    D1=zero_matrix(ZZ,2*d+1,2)
    D2=zero_matrix(ZZ,d+1,2)
    for i in 1:(d+1)
        D1[i,1]=0
        D1[i,2]=d+1-i
    end
    for i in (d+2):size(D1,1)
        D1[i,1]=i-(d+1)
        D1[i,2]=0
    end
    for i in 1:d+1
        D2[i,1]=i-1
        D2[i,2]=d-(i-1)
    end
    #then we compare
    if D1==M || D2==M
        return 1
    else 
        return 0
    end
end

#The input is a matrix representing a path in the plane. The function calculates recursively the multiplicity associated to left turns of the input path. This follows the algorithm described in the references.
function lmultiplicity(M::Union{ZZMatrix,Nothing},d::Integer)::Integer
    if D1orD2(M,d)==1
        return 1
    end
    if isnothing(M)
        return 0
    end
    for i in 1:size(M,1)
    if M[i,1]<0 || M[i,2]<0 || M[i,2]>-M[i,1]+d
       return 0
    end
    end
    if D1orD2(M,d)==0 && left_turn(M)==0 
        return 0
    else
        s=abs(area(M,left_turn(M)))*lmultiplicity(lcut(M),d) + lmultiplicity(lfold(M,d),d)
        return s
    end
end

#The input is a matrix representing a path in the plane. The function calculates recursively the multiplicity associated to right turns of the input path. This follows the algorithm described in the references.
function rmultiplicity(M::Union{ZZMatrix,Nothing},d::Integer)::Integer
    if D1orD2(M,d)==1
        return 1
    end
    if isnothing(M)
        return 0
    end
    for i in 1:size(M,1)
        if M[i,1]<0 || M[i,2]<0 || M[i,2]>-M[i,1]+d
            return 0
        end
    end
    if D1orD2(M,d)==0 && right_turn(M)==0 
        return 0
    else
        s=abs(area(M,right_turn(M)))*rmultiplicity(rcut(M),d) + rmultiplicity(rfold(M,d),d)
        return s
    end
end

#The input is a matrix representing a path in the plane. The function calculates recursively the total multiplicity of the input path. This follows the algorithm described in the references.
function multiplicity(M::Union{ZZMatrix,Nothing},d::Integer)::Integer
    return lmultiplicity(M,d)*rmultiplicity(M,d)
end

######
#In this section of the code, we construct the paths that contribute with positive multiplicity to the curve counting problem for a fixed degree d.

#This is an auxiliary function that makes the indexing in the function building_path() easier.
function helper_sum(i::Integer,d::Integer)::Integer
    s=0
    if i==0
        return 0
    end
    for j in 1:i
        s=s+d+1-(j-1)
    end
    return s
end

#The input is a vector of positive integers. The function constructs a path P by taking the vertices, described by the input vector, of a base path B. Namely, the i-th vertex of P is defined to be the v[i]-th vertex of B.
function building_path(v::Vector,d::Integer)::ZZMatrix
    P=zero_matrix(ZZ,3*d,2)
    B=zero_matrix(ZZ,Int((d+1)*(d+2)/2),2)
    l=d+1
    for i in 0:d
        for k in helper_sum(i,d)+1:helper_sum(i,d)+(d+1-i)
            l=k-helper_sum(i,d)
            B[k,1]=i
            B[k,2]=(d-i)-l+1
        end
    end
    for i in 1:3*d
        P[i,1]=B[v[i],1]
        P[i,2]=B[v[i],2]
    end
    return P
end

#The function generates all possible lambda-increasing paths of length 3d-1 that contribute positively to the associated curve counting problem.
function all_paths(d::Integer)::Vector{ZZMatrix}
    s=collect(2:Int((d+1)*(d+2)/2-1))
    v=subsets(s::Vector,3*d-2::Int)
    allpaths=Vector{ZZMatrix}()
    S=0
    l=0
    for i in 1:length(v)
        w=vec(hcat(1,v[i]...,Int((d+1)*(d+2)/2)))
        P=building_path(w,d)
        m=multiplicity(P,d)
        if m>0
            println("Path number ", i, " has multiplicity ", m)
            S=S+m
        end
        push!(allpaths,P)
    end
    println(" We have total contribution of ", S)
    return allpaths
end

#######
#In this section of the code, we construct the binary trees that are formed by following the recursive formula for the calculation of the left and right multiplicity. Furthermore, we save the paths from the root to the leaves of the binary trees.

#We want to save the binary tree formed in our recursion, so that we are able to recover paths from the root to the leaves, which later produce the Newton subdivision.
#Each entry of this structure contains the information of a path M, its multiplicity mult, and the information of its childern.
mutable struct TreeNode
    M::Any             
    mult::Float64      
    children::Vector{TreeNode}  
end

#The input is a matrix representing a path in the plane. The function calculates the multiplicity (from left turns) of the input path recursively. Furthermore, it saves the childern coming from the recursion.
function l_build_tree(M::Union{ZZMatrix,Nothing},d::Integer)::Union{TreeNode,Nothing}
    if isnothing(M)
        return nothing
    end
    # leaf case
    if D1orD2(M,d) == 1
        return TreeNode(M, 1, TreeNode[]) #Creates node with information M, multiplicity 1 and no childern.
    end
    # multipicity 0
    for i in 1:size(M,1)
        if M[i,1]<0 || M[i,2]<0 || M[i,2]>-M[i,1]+d
            return nothing
        end
    end
    # multipicity 0
    if D1orD2(M,d) == 0 && left_turn(M) == 0
        return nothing
    end
    # recursion
    child1=l_build_tree(lcut(M),d)
    child2=l_build_tree(lfold(M,d),d)
    mult1=0
    mult2=0
    if !isnothing(child1) 
        mult1=child1.mult
    end
    if !isnothing(child2)
        mult2=child2.mult
    end
    total_mult=abs(area(M, left_turn(M)))*mult1 + mult2
    if total_mult > 0 #We only care about non-trivial multiplicity
        children = TreeNode[]
        if !isnothing(child1)
            push!(children, child1)
        end
        if !isnothing(child2)
            push!(children, child2)
        end
        return TreeNode(M, total_mult, children)
    else
        return nothing
    end
end

#The input is a matrix representing a path in the plane. The function calculates the multiplicity (from right turns) of the input path recursively. Furthermore, it saves the childern coming from the recursion.
function r_build_tree(M::Union{ZZMatrix,Nothing},d)::Union{TreeNode,Nothing}
    if isnothing(M)
        return nothing
    end
    # leaf case
    if D1orD2(M,d) == 1
        return TreeNode(M, 1, TreeNode[]) #Creates node with information M, multiplicity 1 and no childern
    end
    # multipicity 0
    for i in 1:size(M,1)
        if M[i,1]<0 || M[i,2]<0 || M[i,2]>-M[i,1]+d
            return nothing
        end
    end
    # multipicity 0
    if D1orD2(M,d) == 0 && right_turn(M) == 0
        return nothing
    end
    # recursion
    child1=r_build_tree(rcut(M),d)
    child2=r_build_tree(rfold(M,d),d)
    mult1=0
    mult2=0
    if !isnothing(child1)
        mult1=child1.mult
    end
    if !isnothing(child2)
        mult2=child2.mult
    end
    total_mult=abs(area(M, right_turn(M)))*mult1 + mult2
    if total_mult > 0 #We only care about non-trivial multiplicity
        children = TreeNode[]
        if !isnothing(child1)
            push!(children, child1)
        end
        if !isnothing(child2)
            push!(children, child2)
        end
        return TreeNode(M, total_mult, children)
    else
        return nothing
    end
end

#The function collects all paths in the binary tree associated to the multiplicity calculation from the root, namely initial_path, to leaves (here we display also the childern).
function collect_paths(node::TreeNode, path::Vector{TreeNode}=Vector{TreeNode}())::Vector{Vector{TreeNode}}
    current_path = copy(path) 
    push!(current_path, node) 
    paths = Vector{Vector{TreeNode}}() 
    if isempty(node.children)
        push!(paths, current_path) #leaf node
    else
        for child in node.children
            append!(paths, collect_paths(child, current_path)) #we add children
        end
    end
    return paths
end

#The input is a vector, whose entries are paths in the binary tree from the root to a leaf. The function prints these paths in the binary tree.
function print_paths(paths::Vector{Vector{TreeNode}})
    for (i, path) in enumerate(paths)
        println("Path $i:")
        for node in path
            println("  Node multiplicity: ", node.mult)
            display(node.M)
        end
        println("------")
    end
end

#######
#In this section of the code, we construct Newton subdvisions of the triangle with vertices {(0,d),(0,0),(d,0)} for each of the paths in the binary tree from the root to a leaf. Furthermore, we implement functions to visualize the Newton subdivision and its corresponding tropical curve.

#The input is a vector, whose entries are paths in the plane. The desired Newton subdivision associated to a path (fron initial_path to a leaf) comes from triangles (formed when cutting) and parallelograms (formed when folding), which are the cells of the Newton subdivision. This function returns the described cells for a path in the binary tree coming from the multiplicity calculation.
function polyhedron_generation(k,path::Vector{TreeNode})::ZZMatrix
    if size(path[k].M,1)==size(path[k+1].M,1)
        parallelogram=zero_matrix(ZZ,4,3)
        for i in 1:size(path[k].M,1)
            if path[k].M[i,1] != path[k+1].M[i,1] || path[k].M[i,2] != path[k+1].M[i,2]
                parallelogram[1,1]=path[k].M[i-1,1]
                parallelogram[1,2]=path[k].M[i-1,2]
                parallelogram[2,1]=path[k].M[i,1]
                parallelogram[2,2]=path[k].M[i,2]
                parallelogram[3,1]=path[k].M[i+1,1]
                parallelogram[3,2]=path[k].M[i+1,2]
                parallelogram[4,1]=parallelogram[1,1]+parallelogram[3,1]-parallelogram[2,1]
                parallelogram[4,2]=parallelogram[1,2]+parallelogram[3,2]-parallelogram[2,2]
                return parallelogram
            end
        end
    else
        triangle=zero_matrix(ZZ,3,3)
        for i in 1:size(path[k].M,1)
            if path[k].M[i,1] != path[k+1].M[i,1] || path[k].M[i,2] != path[k+1].M[i,2]
                triangle[1,1]=path[k].M[i-1,1]
                triangle[1,2]=path[k].M[i-1,2]
                triangle[2,1]=path[k].M[i,1]
                triangle[2,2]=path[k].M[i,2]
                triangle[3,1]=path[k].M[i+1,1]
                triangle[3,2]=path[k].M[i+1,2]
                return triangle
            end
        end
    end
end

#The input are two vectors of paths in the plane, which correspond to paths in the binary trees coming from the left and right multiplicity calculations. The function returns the subdivision of points associated to upper and lower paths (from the initial_path to a leaf). This is used to construct the cells from the desired Newton subdivision.
function Newton_subdivision(path1::Vector{TreeNode},path2::Vector{TreeNode},d::Integer)::SubdivisionOfPoints{QQFieldElem}
    cells=Vector{Vector{Int64}}()
    base_set=zero_matrix(ZZ,Int((d+1)*(d+2)/2),3)
    l=d+1
    for i in 0:d
        for k in helper_sum(i,d)+1:helper_sum(i,d)+(d+1-i)
            l=k-helper_sum(i,d)
            base_set[k,1]=i
            base_set[k,2]=(d-i)-l+1
        end
    end
    BB=[collect(base_set[i,:]) for i in 1:nrows(base_set)]
    vec_dict=Dict(v=> i for (i,v) in enumerate(BB))
    for k in 1:size(path1,1)-1 
        P=polyhedron_generation(k,path1)
        PP=[collect(P[i,:]) for i in 1:nrows(P)]
        indices=[vec_dict[vec] for vec in PP]
        push!(cells,Int.(indices))
    end
    for k in 1:size(path2,1)-1 
        P=polyhedron_generation(k,path2)
        PP=[collect(P[i,:]) for i in 1:nrows(P)]
        indices=[vec_dict[vec] for vec in PP]
        push!(cells,Int.(indices))
    end
    return subdivision_of_points(base_set, cells)
end

#The input are two integers, that correspond to the entries of the vectors allpaths_upper and of allpaths_lower that we want to combine for the Newton subvision. The function displays the corresponding Newton subdivision.
function visualize_subdivision(m::Integer,n::Integer,d::Integer)
    @req m>0 && n>0 && m<=length(allpaths_upper) && n<=length(allpaths_lower) "The given entries are out of bound"
    Nsub=Newton_subdivision(allpaths_upper[m],allpaths_lower[n],d)
    complex=Nsub.pm_subdivision.POLYHEDRAL_COMPLEX
    Polymake.visual(complex)
end

#The input are two integers, that correspond to the entries of the vectors allpaths_upper and of allpaths_lower that we want to combine for the Newton subvision and for the corresponding tropical curve. The function displays the tropical curve.
function visualize_curve(m::Integer,n::Integer,d::Integer) 
    @req m>0 && n>0 && m<=length(allpaths_upper) && n<=length(allpaths_lower) "The given entries are out of bound"
    base=zero_matrix(ZZ,Int((d+1)*(d+2)/2),2)
    l=d+1
    for i in 0:d
        for k in helper_sum(i,d)+1:helper_sum(i,d)+(d+1-i)
            l=k-helper_sum(i,d)
            base[k,1]=i
            base[k,2]=(d-i)-l+1
        end
    end
    Nsub=Newton_subdivision(allpaths_upper[m],allpaths_lower[n],d)
    w=min_weights(Nsub)
    subdiv=subdivision_of_points(base,w)
    curve=tropical_hypersurface(subdiv)
    visualize(curve)
end
