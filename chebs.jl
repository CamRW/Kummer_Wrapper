function chebs(x::Float64,n::Int64)

  pols = Array(Float64,n+1)

  if (n < 0 || n == 0)

    return ("Yo dawg n = 0 or n < 0")

  end

  if (n==2)

    pols[1]=1
    pols[2]=x

    return pols

  end

  pols[1] = 1
  pols[2] = x

  xx1 = 1
  xx2 = x

  for i = 1:n-1

    xx = 2*x*xx2-xx1
    pols[i+2] = xx
    xx1 = xx2
    xx2 = xx

  end

  return pols

end
