import Base.*

# Overload matrix-vector multiplication for BSR matrix
function *(A::BSR_Matrix, x::AbstractVector{Float64})
    convert(Vector{Float64}, A.matrix*x)
end
