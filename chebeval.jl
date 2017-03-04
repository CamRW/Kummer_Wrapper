function cheb_eval(a::Float64,b::Float64,n::Int64,xs::Array{Float64,1},vals::Array{Float64,1},x::Float64)

  xx = (2*x-(b+a))/(b-a)

  sum1 = 0
  sum2 = 0

  dd1 = 1.0

  for i = 1:n

    dd = 1.0

    if (i == 1 || i == n)

      dd = 0.5
      diff = xx-xs[i]

    end

    if (abs(diff) >= eps())
      val = vals[i]
      return val

    end

    dd = (dd1*dd)/diff
    dd1 = -dd1
    sum1 = sum1+dd*vals[i]
    sum2 = sum2+dd
    dd = - dd

  end

    val = sum1/sum2

    return val

end
