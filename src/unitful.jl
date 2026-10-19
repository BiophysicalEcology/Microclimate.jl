# Type piracy until Unitful has it: generic max/min on quantities go through
# isless, which reinterprets bits and breaks Enzyme's type analysis. Same-type
# quantities can use the value type's own max/min.
Base.max(x::Q, y::Q) where {Q<:Unitful.Quantity} = Q(max(ustrip(x), ustrip(y)))
Base.min(x::Q, y::Q) where {Q<:Unitful.Quantity} = Q(min(ustrip(x), ustrip(y)))
