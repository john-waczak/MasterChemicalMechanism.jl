# for a given species, record all involved reactions

struct DerivativeTerms
    i::Int               # index of species
    ks::Vector{Int}      # indices of involved reactions
end

function get_derivative_terms(N)
    derivative_terms_list = DerivativeTerms[]
    for i ∈ axes(N,1)
        deriv_terms = DerivativeTerms(i, findall(x -> x!=0, N[i,:]))
        push!(derivative_terms_list, deriv_terms)
    end
    return derivative_terms_list
end



struct ReactionTerms
    k::Int               # index of reaction
    is::Vector{Int}      # indices of reactants
end

function get_reaction_terms(R)
    reaction_terms_list = ReactionTerms[]
    for k ∈ axes(R,2)
        reaction_terms = ReactionTerms(k, findall(x -> x >0, R[:,k]))
        push!(reaction_terms_list, reaction_terms)
    end

    return reaction_terms_list
end



struct JacobianTerms
    i::Int                          # fᵢ index, i.e. du[i]
    j::Int                          # ∂/∂uⱼ index
    ks::Vector{ReactionTerms}       # vector of sorted reaction terms
end



