function Tableau(nqubits::Integer)
    tab = BitArray(undef, 2 * nqubits, 2 * nqubits + 1)
    for i in 1:(2 * nqubits)
        tab[i, i] = true
    end
    return tab
end

function nqubits(tableau::BitMatrix)
    return Int(size(tableau)[1] / 2)
end

"Does an inplace rowsum as described in Aaronson and Gottesman"
function rowsum!(tableau::BitMatrix, h::Integer, i::Integer)
    nq = nqubits(tableau)
    if h < 1
        throw("h=$(h) must be positive.")
    end
    if h > nq
        throw("h=$(h) must be less than or equal to the number of qubits ($(nq)).")
    end 
    if i < 1
        throw("i=$(i) must be positive.")
    end
    if i > nq
        throw("i=$(i) must be less than or equal to the number of qubits ($(nq)).")
    end

    # Get the sum of g.
    gsum = 0
    for j in 1:nq
        x1 = Int(tableau[i, j])
        z1 = Int(tableau[i, j + nq])
        x2 = Int(tableau[h, j])
        z2 = Int(tableau[h, j + nq])
        if x1 == 0 && z1 == 0
            gsum += 0
        elseif x1 == 1 && z1 == 1
            gsum += z2 - x2
        elseif x1 == 1 && z1 == 0
            gsum += z2 * (2 * x2 - 1)
        else
            gsum += x2 * (1 - 2 * z2)
        end
    end
    # Set the sign bit.
    ri = Int(tableau[i, 2 * nq + 1])
    rh = Int(tableau[h, 2 * nq + 1])
    if (2 * rh + 2 * ri + gsum) % 4 == 0
        # Set rh = 0
        tableau[h, 2 * nq + 1] = false
    else
        # Set rh = 1
        tableau[h, 2 * nq + 1] = true
    end
    # Set all x and z bits.
    for j in 1:nq
        # x_hj = x_ij XOR x_hj
        tableau[h, j] = xor(tableau[i, j], tableau[h, j])
        # z_hj = z_ij XOR x_hj
        tableau[h, nq + j] = xor(tableau[i, nq + j], tableau[h, nq + j])
    end
end

"Inplace CNOT gate controlled by qubit a and acting on qubit b."
function cnot!(tableau::BitMatrix, a::Integer, b::Integer)
    nq = nqubits(tableau)
    if a < 1
        throw("a=$(a) must be positive.")
    end
    if a > nq
        throw("a=$(a) must be less than or equal to the number of qubits ($(nq)).")
    end 
    if b < 1
        throw("b=$(b) must be positive.")
    end
    if b > nq
        throw("b=$(b) must be less than or equal to the number of qubits ($(nq)).")
    end 

    for i in 1:(2 * nq)
        # Set r_i = r_i XOR x_ia z_ia (x_ib XOR z_ia XOR 1)
        temp = xor(tableau[i, b], xor(tableau[i, a + nq], true))
        tableau[i, 2 * nq + 1] = xor(tableau[i, 2 * nq + 1], tableau[i, a] & tableau[i, b + nq] & temp)
        # Set x_ib = x_ib XOR x_ia
        tableau[i, b] = xor(tableau[i, b], tableau[i, a])
        # Set z_ia = z_ia XOR z_ib
        tableau[i, a + nq] = xor(tableau[i, a + nq], tableau[i, b + nq])
    end
end

"Inplace H gate on qubit a."
function hadamard!(tableau::BitMatrix, a::Integer)
    nq = nqubits(tableau)
    if a < 1
        throw("a=$(a) must be positive.")
    end
    if a > nq
        throw("a=$(a) must be less than or equal to the number of qubits ($(nq)).")
    end 

    for i in 1:(2 * nq)
        # Set r_i = r_i XOR x_ia z_ia.
        tableau[i, 2 * nq + 1] = xor(tableau[i, 2 * nq + 1], tableau[i, a] & tableau[i, a])
        # Swap x_ia and z_ia.
        temp = tableau[i, a] # Store x_ia
        tableau[i, a] = tableau[i, a + nq] # Set x_ia = z_ia
        tableau[i, a + nq] = temp
    end
end