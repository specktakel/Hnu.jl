#__precompile__(true)
module CosmoUnits

using Unitful

@unit pc "pc" pc 3.0857e16u"m" true

Unitful.register(@__MODULE__)
function __init__()
    return Unitful.register(@__MODULE__)
end

end
