
@enum SolvingMode mtz scf mcf cec dcc

function solvingmode_from_string(x::AbstractString)
    i::Int = -1
    for mode in instances(SolvingMode)
        if string(mode) == x
            i = Int(mode)
            break
        end
    end
    return SolvingMode(i)
end
