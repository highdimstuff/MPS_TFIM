using ITensors 
using ITensorMPS
using Printf

function Ising_Hamiltonian(L::Int64,J::Float64,h::Float64;bc,noise_on,
    psi_init::Union{MPS,Nothing}=nothing,sites = nothing)
    @show bc # to show the boundary conditon
    # Condition for the boundary type 
    if !in(bc,("OBC","PBC"))
        error("Incorrect boundary condition type! Only OBC and PBC are allowed.")
    end 

    # define sites and its indices with the Pauli matrices
    if isnothing(sites)
        sites = siteinds("S=1/2",L)
    end 

    # OpSum for operator sum
    os = OpSum()
    for i in 1:L-1
        os += -J,"Z",i,"Z",i+1
        os += -h,"X",i#,"Id",i+1
    end 
    os += -h,"X",L

    # adding a small perturbation for classical ising result
    for i in 1:L
        os += -1e-6,"Z",i
    end 

    # Periodic Boundary Condition
    if bc=="PBC"
        os += -J,"Z",1,"Z",L
    end 

    # create the MPO for the Hamiltonian
    H = MPO(os,sites)

    # Take a psi with the above sites but with random values
    if isnothing(psi_init)
        state = ["Up" for i = 1:L]
        ψ = random_mps(sites,state)  
    else
        ψ = psi_init
    end 
    nsweeps = 20
    maxdim = [10, 20, 100, 200, 200]
    cutoff = 1.0e-14

    if noise_on
        # Adding noise in the drmg to converge into global minima 
        println("Noise activated!")
        noise = []
        for i in 1:nsweeps
            if i<round(Int,nsweeps/2)
                push!(noise,0.5*10.0^(-i))
            else 
                push!(noise,0)
            end 
        end 
        energy, ψ = dmrg(H,ψ; nsweeps,maxdim,cutoff,noise)
    else 
        # No noise in dmrg
        energy, ψ = dmrg(H,ψ; nsweeps,maxdim,cutoff)
    end
    return energy, ψ, sites
end 

# Calculation of entropy 
function entanglement_entropy(psi::MPS,l::Int64)
    psi_ortho = orthogonalize(psi,l) # shiftes the orthogonality center to cite l 
    U,S,V = svd(psi_ortho[l], (linkinds(psi_ortho,l-1)...,siteinds(psi_ortho,l)...))
    SvN = 0.0 # Von Neumann entropy S_{vN}
    for n=1:dim(S,1)
        p = S[n,n]^2
        SvN -= p*log(p)/log(2.0)
    end 
    return SvN
end 

# Calculation of expectation value
function expectation_m(psi,oper_type)
    # expectation_m calculates the magnetization <Mz>
    # It can also calucate for <Mx> and other operators 
    m_type = expect(psi,oper_type) 
    N = length(psi)
    # total mag 
    exp_m_type = sum(m_type)/float(N)
    return exp_m_type
end

# Fixed Parameters
J = 1.0; h_upper = 2.0;h_lower = 1e-6; bc = "OBC"; noise_on = true
L = 30; L_half = round(Int64,L/2)
# step size of h (order param)
del_h = 0.025
# Number of iterations
N_iter = round(Int64,h_upper/del_h)

println("Number of iterations:",N_iter)
open("ME_Incp_obc_del_h_0pt025_L_$L.dat","w") do file 
    @printf(file, "%10s %15s %15s %10s %15s %15s \n", "h1", "Mz_1","entropy_1", "h2","Mz_2","entropy_2")
    energy1,psi1,site1 = Ising_Hamiltonian(L,J,h_lower;bc,noise_on)
    energy2,psi2,site2 = Ising_Hamiltonian(L,J,h_upper;bc,noise_on)
    m1 = expectation_m(psi1,"Z")
    m2 = expectation_m(psi2,"Z")
    entropy1 = entanglement_entropy(psi1,L_half)
    entropy2 = entanglement_entropy(psi2,L_half)
    @printf(file,"%10.6f %15.8f %15.8f %10.6f %15.8f %15.8f\n",h_lower,m1,entropy1,h_upper,m2,entropy2)
    let 
        h1 = h_lower
        h2 = h_upper 
        for i in 1:N_iter
            #increasing h
            h1 += del_h
            energy1,psi1,site1 = Ising_Hamiltonian(L,J,h1;bc,noise_on,psi_init = psi1,sites = site1)
            m1 = expectation_m(psi1,"Z")
            entropy1 = entanglement_entropy(psi1,L_half)
            #decreasign h
            h2 -= del_h
            energy2,psi2,site2 = Ising_Hamiltonian(L,J,h2;bc,noise_on,psi_init = psi2,sites = site2)
            m2 = expectation_m(psi2,"Z")
            entropy2 = entanglement_entropy(psi2,L_half)
            @printf(file,"%10.6f %15.8f %15.8f %10.6f %15.8f %15.8f \n",h1,m1,entropy1,h2,m2,entropy2)
        end 
    end 
end 
