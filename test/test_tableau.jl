import PlusOne
import Test

function test_bell_parity()
    tab = PlusOne.Tableau(2)
    PlusOne.hadamard!(tab, 1)
    PlusOne.cnot!(tab, 1, 2)
    m1 = PlusOne.measure!(tab, 1)
    m2 = PlusOne.measure!(tab, 2)
    return m1 == m2
end

Test.@test test_bell_parity()