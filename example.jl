include("vertex_enumeration.jl")
A = [3 3; 2 5; 0 6]
B = [3 2; 2 6; 3 1]
vertex_enumeration(A, transpose(B))