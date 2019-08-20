using PyCall

# Block compressed row matrix
struct BSR_Matrix
    matrix::PyCall.PyObject
end

struct PrecomputedSolver
    alpha
    dcurve
    X
    Y
    Lgrid
    PUXStor
    interior
    solve
end

struct Solution
    # Store solution
    U1
    U2
    # Derivatives
    dU1dx
    dU1dy
    dU2dx
    dU2dy
    # Corresponding RHS
    RHS1
    RHS2
    # BCs
    ubc1
    ubc2
end
