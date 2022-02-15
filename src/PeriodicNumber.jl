# Import the functions for which we define methods to be called when using PeriodicNumbers
import Base.show
import Base: +, -, *, /, abs, abs2
using LinearAlgebra

#######################################################################################################################
#######################################################################################################################
### Definition of 'PeriodicNumber' with its respective methods for basic operations
#######################################################################################################################
#######################################################################################################################
"""
    PeriodicNumber{T}(value::T, L::T)

Define a number whose value is always contained between 0 and L; more precisely ∈ [0,L).


When it is initialized with value>L or value<0 the value is transformed to be its L-modulus.

"""
struct PeriodicNumber{T<:Real} <:Number
    value::T # Value to be converted if not ∈ [0,L)
    L::T # upper bound
    function PeriodicNumber{T}(value, L) where {T<:Real}
        zero_T = convert(T, 0)
        if L<=zero_T
            error("The period should be a positive number.")
        end
        if (value<zero_T)||(value>=L)
            value = mod(value, L)
            if value==L
                value=zero_T
            end
        end
    new(value, L)
    end
end

# Simple constructor for periodic number
"Constructor for 'PeriodicNumber', with L=1 by default."
PeriodicNumber(x::T, L::T=1.0) where {T<:Real} = PeriodicNumber{T}(x, L)

#### Format the output of periodic numbers
function Base.show(io::IO, PN::PeriodicNumber)
    L = PN.L
    if L==1.0
        print(io, "PN(", PN.value, ")" )
    else
        print(io, "PN(", PN.value, " L=", L, ")" )
    end
end
PeriodicNumber(1.2)


#################
### Definition of basic operations
#################

# Define sum of periodic values among themselves or with real numbers. The result is always a PeriodicNumber type
function Base.:+(PN1::PeriodicNumber, PN2::PeriodicNumber)::PeriodicNumber
    x1=PN1.value; x2=PN2.value
    L1=PN1.L; L2=PN2.L
    if (L1 !=L2)
        error("Comparing PeriodicNumbers of different periodicity")
    end

    return PeriodicNumber(x1+x2, L1)
end


function Base.:+(PN1::PeriodicNumber, x2::Real)::PeriodicNumber
    x1=PN1.value;  L=PN1.L
    return PeriodicNumber(x1+x2, L)
end

Base.:+(x1::Real, PN2::PeriodicNumber)::PeriodicNumber = PN2+x1


# Define product periodic numbers among themselves or with real numbers. The result is always a PeriodicNumber type
function Base.:*(PN1::PeriodicNumber, PN2::PeriodicNumber)::PeriodicNumber
    x1=PN1.value; x2=PN2.value
    L1=PN1.L; L2=PN2.L
    if (L1 !=L2)
        error("Comparing quantities of different periodicity")
    end
    return PeriodicNumber(x1*x2, L1)
end


function Base.:*(PN1::PeriodicNumber, x2::Real)::PeriodicNumber
    x1=PN1.value; L=PN1.L
    return PeriodicNumber(x1*x2, L)
end

Base.:*(x1::Real, PN2::PeriodicNumber)::PeriodicNumber = PN2*x1


#= Define difference of periodic values among themselves or with real numbers. The result is always a PeriodicNumber type, and hence positive.
The difference between periodic values is somewhat tricky, because we need to consider that two particles far apart in the real space might be very close in the periodic one.
For instance, if x1=0 and x2=L, the real distance is L, but in the periodic space it's 0
# NB: The following function is essentially computing a distance in the periodic space, since the value is always positive and
does not depend on the order, i.e. x1-x2 = x2-x1.
=#
"""
    -(PN1::PeriodicNumber, PN2::PeriodicNumber)

Return the difference of two periodic numbers, also as a 'PeriodicNumber' type.

Given that the output is a periodic number, and hence its value is positive, this 'difference
function' is actually commutative. 
Moreover, the output's value can never be larger than L/2.
"""
function Base.:-(PN1::PeriodicNumber, PN2::PeriodicNumber)::PeriodicNumber
    x1=PN1.value; x2=PN2.value
    L1=PN1.L; L2=PN2.L
    if (L1 !=L2)
        error("Comparing quantities of different periodicity")
    end
    arg=x1-x2
    arg=min(PeriodicNumber(arg, L1).value, PeriodicNumber(-arg, L1).value )
    return PeriodicNumber(arg, L1)
end


#= 
NB: This is also an UNsigned difference. 
Thus, it is made to work only with PeriodicNumber's whose lower limit is equal to 0 
=#
function Base.:-(PN1::PeriodicNumber, x2::Real)::PeriodicNumber
    x1=PN1.value; L=PN1.L
    dif=abs(x1-x2)
    while dif>0.5*L
        dif -= L
    end
    return PeriodicNumber(dif, L)
end


#NB: Note also that since the difference of a periodic number and any real number is UNsigned, the '-' operator is commutative
Base.:-(x2::Real, PN1::PeriodicNumber)::PeriodicNumber = PN1-x2

#Simple first calls for compilation purposes
PeriodicNumber(0.1, 3.0)+ PeriodicNumber(2.5, 3.0); PeriodicNumber(2.5, 3.0) - PeriodicNumber(0.1, 3.0)
PeriodicNumber(0.1, 3.0) +2.0; PeriodicNumber(0.1, 3.0) - 2.0
rand()*PeriodicNumber(2.5, 3.0); PeriodicNumber(2.5, 3.0)*rand()


# Define other useful methods for computing absolute values, norm, etc.
value(PN::PeriodicNumber{<:Real}) = PN.value
abs(PN::PeriodicNumber{<:Real}) = Base.abs(value(PN))
abs2(PN::PeriodicNumber{<:Real}) = Base.abs2(value(PN))
# value(pn1); abs(pn1); abs2(pn2) # call for compilation

# Define a 'distance' function to work with numbers with periodic BC
"""
    periodic_distance(PN1::PeriodicNumber, PN2::PeriodicNumber)

Compute the *un*-signed distance between a pair of periodic numbers.

It is implemented simply as the 'value' field of '-' method for periodic numbers.
Since such method always returns a non-negative value and is commutative, this function 
is a well defined distance.
"""
periodic_distance(P1::PeriodicNumber, P2::PeriodicNumber)::Float64=(P1-P2).value
periodic_distance(PeriodicNumber(0.1, 3.0), PeriodicNumber(2.5, 3.0)) # call for compilation
