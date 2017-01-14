using LRSLib

function HtoV(A, b)
    """
    Change the H representation of a Polyhedron to
    V representation.

    Parameters
    ----------
    A : Array{Int64,2}
        The matrix of left-hand side of matrix inequality:
        Ax <= b

    b : Array{Int64,1}
        The vector of right-hand side of matrix inequality.

    Return
    ------
    V : Polyhedra.SimpleVRepresentation{2,Rational{BigInt}}
        The vertex matrix which is V representation of 
        Ax <= b 

    """
    linset = IntSet()
    H = Polyhedra.SimpleHRepresentation(A, b, linset)
    H = LRSMatrix(H)
    V = LRSGeneratorMatrix(H)
    V = Polyhedra.SimpleVRepresentation(V)

    return V
end

function vertex_enumeration(A, B)
    """
    Given payoff matrixs of two-player strategic form game,
    extend the matrixs with nonnegativity constraint of mixed actions,
    and do vertex_enumeration to find NEs. 

    Parameters
    ----------
    A : Array{Int64,2}
        The payoff matrix of Player0.

    B : Array{Int64,2}
        The payoff matrix of Player1.

    Return
    ------
    NE : Array{Polyhedra.SimpleVRepresentation{2,Rational{BigInt}}, 1}
        The NEs found by vertex enumeration, which are represented by 
        vertexs.

    """
    NE = []
    
    n1 = size(A)[2]
    n2 = size(B)[2]
    b1 = ones(size(A)[1])
    b2 = ones(size(B)[1])
    
    A = vcat(A, -eye(Int64, n1))
    b1 = vcat(b1, zeros(n1))
    B = vcat(-eye(Int64, n2), B)
    b2 = vcat(zeros(n2), b2)

    Q = HtoV(A, b1)
    P = HtoV(B, b2)

    for i in 1:size(A)[1]
        for j in 1:size(B)[1]
            v1 = Q.V[i, :]
            v2 = P.V[j, :]
        
            label1 = *(A, v1) - b1
            ind1 = label1 .!= 0
            ind2 = label1 .== 0
            label1[ind1] = 0
            label1[ind2] = 1
        
            label2 = *(B, v2) - b2
            ind3 = label2 .!= 0
            ind4 = label2 .== 0
            label2[ind3] = 0
            label2[ind4] = 1
        
            if label1 + label2 == [1., 1., 1., 1., 1.]
                push!(NE, (v1, v2))
            end
        end
    end

    return NE
end