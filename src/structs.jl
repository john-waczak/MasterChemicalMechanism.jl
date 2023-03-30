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



function get_jac_prototype(N,R)
    Jprototype = zeros(Float64, size(R,1), size(R,1))

    In,Kn,Vn = findnz(N)

    for (i,k) ∈ zip(In, Kn)
        for j ∈ axes(R,1)
            if N[i,k] != 0 && R[j,k] > 0
                Jprototype[i,j] = 1.0
            end
        end
    end


    Jprototype = sparse(Jprototype)
    println("nonzero %: ", 100*length(nonzeros(Jprototype))/(size(Jprototype,1)^2))
    return Jprototype
end


function get_jac_terms(Jprototype, N, R)
    jac_terms = JacobianTerms[]
    Ijac,Jjac,Vjac = findnz(Jprototype)
    for (i,j) ∈ zip(Ijac, Jjac)
        rxn_terms = ReactionTerms[]
        for k ∈ axes(R,2)
            if N[i,k] != 0 && R[j,k] > 0
                reaction_indices = findall(x -> x > 0, R[:,k])
                indices_out = [j]
                for idx_rxn ∈ reaction_indices
                    if idx_rxn != j
                        push!(indices_out, idx_rxn)
                    end
                end
                push!(rxn_terms, ReactionTerms(k, indices_out))
            end
        end
        push!(jac_terms, JacobianTerms(i,j,rxn_terms))
    end
end



