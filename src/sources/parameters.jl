module Parameters

export LINEAR, LOG

const LINEAR = 0
const LOG = 1

abstract type Parameter end


struct LinearParameter<:Parameter
    name::String
    min::Float64
    max::Float64
    step::Float64
end


struct LogarithmicParameter<:Parameter
    name::String
    min::Float64
    max::Float64
    points::Int
    # fixed::Bool=false # TODO implement
end

function make_parameter_grid(p::LinearParameter)
    # step = (p.max - p.min) / p.points
    return Vector(range(p.min, p.max, step=p.step))
end

function make_parameter_grid(p::LogarithmicParameter)
    step = (log(p.max) - log(p.min)) / p.points
    return exp.(Vector(range(log(p.min), log(p.max), step)))
end

function make_linear_parameter_grid(
    name::String,
    min::Float64,
    max::Float64,
    points::Int
)
    param = LinearParameter(name, min, max, points)
    return make_parameter_grid(param)
end

end