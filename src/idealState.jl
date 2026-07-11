# idealState.jl - Ideal Gas States

# Structure (type) definition
# ---------------------------

struct IdealState{ℙ <: FLOAT}
    G::IdealGas{ℙ}
    M::Quantity{ℙ, dimension(u"kg/kmol"), typeof(u"kg/kmol")}
    R::Quantity{ℙ, dimension(u"kJ/kg/K"), typeof(u"kJ/kg/K")}
    P::Quantity{ℙ, dimension(u"kPa"), typeof(u"kPa")}
    T::Quantity{ℙ, dimension(u"K"), typeof(u"K")}
    v::Quantity{ℙ, dimension(u"m^3/kg"), typeof(u"m^3/kg")}
    ρ::Quantity{ℙ, dimension(u"kg/m^3"), typeof(u"kg/m^3")}
    u::Quantity{ℙ, dimension(u"kJ/kg"), typeof(u"kJ/kg")}
    h::Quantity{ℙ, dimension(u"kJ/kg"), typeof(u"kJ/kg")}
    s::Quantity{ℙ, dimension(u"kJ/kg/K"), typeof(u"kJ/kg/K")}
    cp::Quantity{ℙ, dimension(u"kJ/kg/K"), typeof(u"kJ/kg/K")}
    cv::Quantity{ℙ, dimension(u"kJ/kg/K"), typeof(u"kJ/kg/K")}
    Pr::ℙ
    vr::ℙ
    # Internal, validating constructors
end

# "ﬆ" is U+FB06
