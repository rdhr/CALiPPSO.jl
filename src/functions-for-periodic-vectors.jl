using StaticArrays
include("PeriodicNumber.jl")


#= 
Definition of norm of vectors whose elements are PeriodicNumbers.
Note that in this way, the same function can be applied to usual 'Vector' arrays, but also to 'StaticVector' types as will be used afterwards
=#
import LinearAlgebra.norm
"norm function for a vector of `PeriodicNumber` elements"
norm(v::AbstractVector{<:PeriodicNumber}) = norm(value.(v))
norm(PeriodicNumber.(3.0*rand(5), 2.0))

#########################################################################################################
#########################################################################################################
### Functions to create static vectors whose components are of type 'PeriodicNumber'
#########################################################################################################
#########################################################################################################

"""
    PeriodicVector(vec::Vector{T}, L::T=1.0)

Convert 'vec' to a Static Vector of the same size, but with elements of 'PeriodicNumber{T}' type.


The periodicity of the numbers, 'L', is given as second argument and defaults to 1.0.
"""
function PeriodicVector(vec::Vector{T}, L::T=1.0) where {T<:Real}
    dim=length(vec)
    if typeof(L)==T
        SVector{dim}(PeriodicNumber.(vec, L))
    else
        SVector{dim}(PeriodicNumber.(vec, convert(T, L)))
    end    
end
PeriodicVector(rand(4), 2.0)
PeriodicVector([1,3,53,4], 3)


"""
    PeriodicVectors(mat::Matrix{T}, L::T=1.0) where {T<:Real}

Convert 'mat' to a Vector of elements of type SVector{d, PeriodicNumber{T}}.

mat should be of size d x N. The output is a 1-dimensional array of N elements, each of which
consists of StaticVector's of size d, and PeriodicNumber{T} as data.
The periodicity of the numbers, 'L', is given as second argument and defaults to 1.0.
"""
function PeriodicVectors(mat::Matrix{T}, L::T=1.0) where {T<:Real}
    dim, N = size(mat)
    @assert dim<N "d>N. N should be larger than the dimensionality of vectors"
    [PeriodicVector(mat[:,i], L) for i in 1:N]
end
PeriodicVectors(rand(4, 10), 2.0)
PeriodicVectors(rand(1:100, 4, 10), 2)

#### Format the output of periodic vectors
function Base.show(io::IO, PV::SVector{d, <:PeriodicNumber}) where {d}
    L = PV[1].L
    vals = value.(PV)
    if L==1.0
        print(io, "PV(", vals, ")")
    else
        print(io, "PV(", vals, "; L=", L, ")" )
    end
end

# vt1 = PeriodicVector(rand(4), 0.5)
# vt2 = PeriodicVector(rand(4), 0.5)

#########################################################################################################
#########################################################################################################
### Functions to compute the distance between two vectors with 'PeriodicNumber' components
#########################################################################################################
#########################################################################################################

"""
    trit_to_shift(x::Char, L::Real=1.0)

Associate a base-3 digit (as a character) to a shift.

The idea is to use a ternary representation (trit is for base 3, what bit is for base 2) to 
identify a system's images by shifts.
The association rule is the following:
    0: for no shift
    1: for a shift of size L; so the image to the right
    2: for a shift of size -L; so the image to the left

In Julia 0.7 and later transforming a number to a different base outputs a string, that's why
it is useful to assume that x is type 'Char'
"""
function trit_to_shift(x::Char, L::Float64=1.0)
    if x=='0'
        return 0.0
    elseif x=='1'
        return L
    elseif x=='2'
        return -L
    else
        return error("Wrong digit as input! Only '0', '1', '2' are possible.")
    end
end


"""
    trit_to_shift_vector(s::String, L::T=1.0)

Split 's' into characters and apply 'trit_to_shift' to each of them.

The output is thus a vector of the same number of elements as characters in 's'.
The value of the shift, 'L', defaults to 1.0.

See also: [`trit_to_shift`](@ref)
"""
function trit_to_shift_vector(s::String, L::Float64=1.0) 
    trit_to_shift.(collect(s), L)
end


"""
    generate_system_images(d::Int64, L::T=1.0)

Generate the 3^d-1 (virtual) images of a d-dimensional system of periodicity 'L' in each side.

The output is a Vector whose elements are of type SVector{d,T}. 
The value of 'L' defaults to 1.0.0
For a d-dimensional system, there are 3^d-1 virtual images. Each of them is identified with a
single number to which a shift vector is associated.

See also: [`trit_to_shift`](@ref), [`trit_to_shift_vector`](@ref)
"""
function generate_system_images(d::Int64, L::Float64=1.0)
    images = Vector{SVector{d, Float64}}(undef, 3^d-1)
    # Note that i does not begin in '0'; so we're excluding the case of no shift (the most common one for big systems)
    for i in 1:3^d-1::Int
        str = string(i, base=3, pad=d)
        images[i] = SVector{d}(trit_to_shift_vector(str, L) )
    end
    images
end
generate_system_images(2, 5.0)
# imgst = generate_system_images(4, 0.5)



"""
    MIC_distance(V1::SVector{d, PeriodicNumber{T}}, V2::SVector{d, PeriodicNumber{T}}, images::Vector{SVector{d, T}})

Compute the distance between 'V1' and 'V2' considering periodic boundary conditions and using 
the so called Minimum Image Convention (MIC).

To avoid unnecessary calls in other functions, the system's set of virtual images should be 
given as third argument, as an array of SVector's. 
It is assumed that both vectors are contained in a d-dimensional box of size L. 
The MIC guarantees that all the possible periodic shifts or 'virtual images' of the relative vectors 
are considered when calculating their distance.

The output is the index of the image of minimum distance (0 if no shift is needed) and the value of such distance.

See also: [`MIC_vector`](@ref).
"""
function MIC_distance(V1::SVector{d, PeriodicNumber{T}}, V2::SVector{d, PeriodicNumber{T}}, images::Vector{SVector{d, T}}) where {d, T<:Float64}
    L1 = V1[1].L; L2 = V2[1].L # periodicities of each vector
    #check that both vectors have the same periodicity
    (L1 != L2) && error("Comparing periodic vectors of different periodicity")

    X1 = value.(V1); X2 = value.(V2) # values of the components of the first and second vector
    diffs =  abs.(X1 .- X2) # absolute value of the difference of each of the vectors components

    if all(x -> x<=L1/2, diffs) # if all the components differ by less than half the box size, there is no need to consider the system images
        return 0, norm(diffs) # the output is the 0-th image (so the system itself), and the usual norm of the vector difference

    else # if at least one of the components is larger than L/2, then check all the distances with all the possible images
        #TODO: Try something smarter, so that not all the 3^d-1 shifts are evaluated, but only the relevant ones
        min_dist, close_image = findmin( [norm(X1 .- (X2.+shift) ) for shift in images])
        return close_image, min_dist # the output is the index of the image that produced the smallest distance, and the value of such distance
    end
end
MIC_distance(PeriodicVector(rand(4), 0.5), PeriodicVector(rand(4), 0.5), generate_system_images(4, 0.5))


"""
    MIC_vector(V1::SVector{d, PeriodicNumber{T}}, V2::SVector{d, PeriodicNumber{T}}, images::Vector{SVector{d, T}})

Compute the vector joining V2 to V1 (i.e. V1-V2) considering periodic boundary conditions
and using the so called Minimum Image Convention (MIC). 

To avoid unnecessary calls in other functions, the system's set of virtual images should be 
given as third argument, as an array of SVector's. 
It is assumed that both vectors are contained in a d-dimensional box of size L. 
The MIC guarantees that all the possible periodic shifts or 'virtual images' of the vectors 
are considered when calculating their distance.

The output is the MIC vector (as a `SVector` type) and the index of the image of minimum
 distance (0 if no such image is needed).

See also: [`MIC_distance`](@ref).
"""
function MIC_vector(V1::SVector{d, PeriodicNumber{T}}, V2::SVector{d, PeriodicNumber{T}}, images::Vector{SVector{d, T}}) where {d, T<:Float64}
    L1 = V1[1].L; L2 = V2[1].L # periodicities of each vector
    #check that both vectors have the same periodicity
    (L1 != L2) && error("Comparing periodic vectors of different periodicity")

    X1 = value.(V1); X2 = value.(V2) # values of the components of the first and second vector
    diffs =  abs.(X1 .- X2) # absolute value of the difference of each of the vectors components

    if all(x -> x<=L1/2, diffs) # if all the components differ by less than half the box size, there is no need to consider the system images
        return SVector{d}(X1 .- X2), 0   # the output is the usual vector difference and the 0-th image (so the system itself), and 
    else # if at least one of the components is larger than L/2, then check all the distances with all the possible images
        #TODO: Try something smarter, so that not all the 3^d-1 shifts are evaluated, but only the relevant ones
        min_dist, close_image = findmin( [norm(X1 .- (X2.+shift) ) for shift in images])
        
        # the output is the vector difference using a virtual image that produces the smallest distance, and the index of such image
        return SVector{d}(X1 .- (X2 .+ images[close_image])), close_image 
    end
end
MIC_vector(PeriodicVector(rand(4), 0.5), PeriodicVector(rand(4), 0.5), generate_system_images(4, 0.5))



#########################################################################################################
#########################################################################################################
### Functions to check distances between set of vectors
#########################################################################################################
#########################################################################################################
# Temporary values and arrays to perform first, simple calls on the following functions
# dt= 4; Lt=1.2
# Xst = PeriodicVectors(5 .*rand(dt, 20) .- 2, 1.2)

"""
    distances_between_centers(Xs::Vector{SVector{d, PeriodicNumber{T}}})
    
Compute the distance between all pair of SVector's.

If 'Xs' is an array of length N, the output is a NxN symmetric matrix. 
When this function is called each element of 'Xs' plays the role of the center of a particle.

We make use of the 'Symmetric' constructor, so only N*(N-1)/2 distances are calculated.

See also [`distances_between_particles`](@ref)
"""
function distances_between_centers(Xs::Vector{SVector{d, PeriodicNumber{T}}}) where {d, T<:Real}
    N = length(Xs)
    distances=zeros(N,N)
    for i in 1:N, j in i+1:N
        distances[i,j]=norm(Xs[i] .- Xs[j])
    end
    return Symmetric(distances)
end
# distances_between_centers(Xst)



"""
    distances_between_centers(Xs::Vector{SVector{d, PeriodicNumber{T}}},  images::Vector{SVector{d, T}})
    
Compute the distance between all pair of SVector's and the image of MIC distance between them.

If 'Xs' is an array of length N, the output is a NxN symmetric matrix containing the distances
and a NxN symmetric matrix containing the MIC images. 
When this function is called each element of 'Xs' plays the role of the center of a particle.

We make use of the 'Symmetric' constructor, so only N*(N-1)/2 distances are calculated.

See also: [`MIC_distance`](@ref).
"""
function distances_between_centers(Xs::Vector{SVector{d, PeriodicNumber{T}}}, images::Vector{SVector{d, T}}) where {d, T<:Real}
    N = length(Xs)
    distances=zeros(N,N)
    image_indices = zeros(Int64, N, N)
    for i in 1:N, j in i+1:N
        img, dist = MIC_distance(Xs[i], Xs[j], images)
        distances[i,j] = dist
        image_indices[i,j] = img
        
        if img != 0
            min_val, mirror_img = findmin(norm.( images .+ [images[img]] ))
            image_indices[j, i] = mirror_img
        else
            image_indices[j, i] = 0
        end
    end
    return Symmetric(distances), image_indices
end
# distances_between_centers(Xst, generate_system_images(dt, Lt))



"""
    check_for_overlaps(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, tolerance::T)

Check whether there is an overlap between all pairs of particles centred at 'Xs' and of radius
'R'. 

Each element of 'Xs' plays the role of the position of a particle's center, and all of them
are assumed to be of the same size, i.e. radius 'R'.
A tolerance to determine whether there is an overlap or not should  be passed as third argument.

# Output

- 'overlap': a boolean that is 'true' only when an overlap is present
- 'message': a string  that contains some information about the overlapping particles.
- 'particles': a tuple containing the indices of the overlapping particles; '(0, 0)' if no overlap
"""
function check_for_overlaps(Xs::Vector{SVector{d, PeriodicNumber{T}}}, R::T, tolerance::T) where {d, T<:Float64}
    σ = 2*R
    N = length(Xs)
    distances = distances_between_centers(Xs)
    overlap=false; message = "No overlap is present"; particles = (0, 0)
    for i in 1:N, j in i+1:N
        if distances[i,j]-σ < -tolerance
            gap = distances[i,j]-σ
            Xi = value.(Xs[i])
            Xj = value.(Xs[j])
            message = "Overlap between particles $i and $j;\n Centers at: $Xi; \t $Xj. \n\t\tradius:  $R; \t Overlap = $gap ."

            overlap = true
            particles = (i, j)
            break
        end
    end
    return overlap, message, particles
end
# check_for_overlaps(Xst, 0.0, 1e-10)
# check_for_overlaps(Xst, 1.2, 1e-10)
