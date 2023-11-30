module SUT
#=
A structure for the elements of a supply-use table (SUT). Additional (intensity) elements are calculated based on absolute SUT data. Changing individual elements after original setup will not recalculate dependent variables.

---

# usage as for example:

include("SUT_structure.jl")

### to get the make matrix
foo = <SUT.>structure(U,V,F,e,pii)
foo.V

# or

### to sum up the use matrix and commodity final demand
bar = <SUT.>structure(U,V,F,e,pii;t="physical")
sum(bar.U,dims=2) + bar.e
=#

    using LinearAlgebra

    export structure

    mutable struct structure

        U::Matrix{Float64} # Use matrix
        V::Matrix{Float64} # Make matrix (transpose of supply matrix)
        F::Matrix{Float64} # Direct factor requirements
        e::Matrix{Float64} # Commodity final demand
        pii::Matrix{Float64} # Factor prices
        g::Union{Matrix{Float64},Nothing} # Total industry output
        q::Matrix{Float64} # Total commodity output
        f::Matrix{Float64} # Total factor requirements
        ϕ::Union{Matrix{Float64},Nothing} # Factor limits (endowments)
        B::Union{Matrix{Float64},Nothing} # Use coefficients
        C::Union{Matrix{Float64},Nothing} # 
        D::Matrix{Float64} # Industry share
        S::Union{Matrix{Float64},Nothing} # Factor intensities
        t::Union{String,Nothing} # Unit of commodity-related SUT elements
        isActive

        function structure(dict, t; ϕ=nothing)
            #=
            if size(V) != size(U')
                @error "The make and use matrices have different dimensions."
            elseif size(V,1) != size(F,2)
                @error "The factor exchanges have a different industry dimension than the make and use matrices."
            end

            if !isequal(
                names(dict["V"][:,2:end]), dict["U"][:,1], dict["e"][:,1]
                )
                @error "The make, use, and final demand matrices cover different numbers of commodities."
            end
            if dict["V"][:,1] != names(dict["U"][:,2:end]) != names(dict["F"][:,2:end])
                @error "The make, use, and final demand matrices cover different numbers of industries."
            end
            #if dict["F"][:,1] == dict["pii"][:,1]
            =#

            
            this = new()
            try 
                #SUT type: monetary, physical, or mixed
                this.t = t
                #SUT system
                this.U = dict["U"]
                this.V = dict["V"]
                this.F = dict["F"]
                this.e = dict["e"]
                this.pii = dict["pii"]
                this.q = sum(transpose(dict["V"]),dims=2)
                this.f = sum(dict["F"],dims=2)
                isnothing(ϕ) && (this.ϕ = this.f)
                this.D = this.V*inv(diagm(vec(this.q)))
                
                # Calculate the following only when data is given in monetary units
                if this.t == :"monetary"
                    this.g = sum(dict["V"],dims=2)
                    this.B = dict["U"]*inv(diagm(vec(this.g)))
                    this.C = transpose(dict["V"])*inv(diagm(vec(this.g)))
                    this.S = dict["F"]*inv(diagm(vec(this.g)))

                    if sum(this.V',dims=2) != sum(this.U,dims=2) + this.e && sum(this.V,dims=2) != transpose(sum(this.U,dims=1) + this.pii'*this.F)
                        @warn("The SUT-system is not balanced.")
                    end
                else
                    this.g = nothing
                    this.B = nothing
                    this.C = nothing
                    this.S = nothing
                    @info "The SUT system is not given in monetary units such that no total industry output and dependent matrices may be calculated."
                    sum(this.V',dims=2) == sum(this.U,dims=2) + this.e || @warn("The SUT-system is not balanced.")
                end

            catch
                rethrow()
            end

            @info("A structure for supply-use elements was set up. Changing individual elements will not change others automatically.")
            return this
        end
    end
end
