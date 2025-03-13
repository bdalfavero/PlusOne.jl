"""
Build a tableau in the all-zero state. The tableau has one row for each
stabilizer/destabilizer (2n rows) as well as an extra row for measurement
scratch space. There are 2n columns for the (de)stabilizers, and a final
column for the sign bits.
"""
function Tableau(nqubits::Integer)
    tab = BitArray(undef, 2 * nqubits + 1, 2 * nqubits + 1)
    for i in 1:(2 * nqubits)
        tab[i, i] = true
    end
    return tab
end

"Returns the number of qubits for this tableau."
function nqubits(tableau::BitMatrix)
    return Int((size(tableau)[1] - 1) / 2)
end

"Does an inplace rowsum as described in Aaronson and Gottesman"
function rowsum!(tableau::BitMatrix, h::Integer, i::Integer)
    nq = nqubits(tableau)

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

"Inplace S gate on qubit a."
function phase!(tableau::BitMatrix, a::Integer)
    nq = nqubits(tableau)
    if a < 1
        throw("a=$(a) must be positive.")
    end
    if a > nq
        throw("a=$(a) must be less than or equal to the number of qubits ($(nq)).")
    end 

    for i in 1:(2 * nq)
        # Set r_i = r_i XOR x_ia z_ia
        tableau[i, 2 * nq + 1] = xor(tableau[i, 2 * nq + 1], tableau[i, a] & tableau[i, a + nq])
        # Set z_ia + z_ia XOR x_ia
        tableau[i, a + nq] = xor(tableau[i, a + nq], tableau[i, a])
    end
end

"Measure the a^th qubit, updating the tableau in place."
function measure!(tableau::BitMatrix, a::Integer)
    nq = nqubits(tableau)
    # Find a p in [n+1, ..., 2n] s.t. x_pa = 1.
    p_exists = false
    p = 0
    for pprime in (nq + 1):(2 * nq)
        if tableau[pprime, a]
            p_exists = true
            p = pprime
            break
        end
    end
    if p_exists
        # Measurement is random!
        # Call rowsum on the appropriate rows.
        for i in 1:(2 * nq)
            if (i != p) && tableau[i, a]
                rowsum!(tableau, i, p)
            end
        end
        # Set the (p-n)^th row to the p^th row.
        tableau[p - nq, :] = tableau[p, :]
        # Set the x and z bits p^th row to zero.
        for i in 1:(2 * nq)
            tableau[p, i] = false
        end
        # Set r_p to 0 or 1 with a coin flip.
        tableau[p, end] = rand([0, 1]) != 0
        # Set z_pa = 1.
        tableau[p, a + nq] = true
        return tableau[p, 2 * nq + 1]
    else
        # Measurement is deterministic.
        for j in 1:(2 * nq + 1)
            tableau[2 * nq + 1, j] = false
        end
        for i in 1:nq
            if tableau[i, a]
                rowsum!(tableau, 2 * nq + 1, i + nq)
            end
        end
        return tableau[2 * nq + 1, 2 * nq + 1]
    end
end