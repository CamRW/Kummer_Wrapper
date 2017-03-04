
#function cheb_adap{T<:AbstractFloat,I<:Integer}(fun::Function,eps::T,k::I,xs::Array{T,1},u::Array{T,2},a::T,b::T)

function cheb_adap(fun::Function,eps::Float64,k::Int64,xs::Array{Float64,1},u::Array{Float64,2},a::Float64,b::Float64)

# Output Declarations
# [0,1] sqrt(x) test or 1/sqrt(x)

#  nints::Int64


#
  maxints::Int64
  nn::Int64
  nn = k/2
  maxints = 1000000

  ab0 = Array(Float64,2,maxints)
  ab1 = Array(Float64,2,maxints)
  vals = Array(Float64,k)
  coefs = Array(Float64,k)
  ts = Array(Float64,k)
  ab = Array(Float64,2,k)


  nints0 = 1

  ab0[1,1] = a
  ab0[2,1] = b

  nints1 = 0
  nints = 0
while (nints0 > 0)
    a0 = ab0[1,nints0]
    b0 = ab0[2,nints0]
    c0 = (a0+b0)/2

      if (b0-a0 == 0)

        ier = 1024
        return "Pew Pew"

      end

      nints0 = nints0 - 1

      ts = xs*(b0-a0)/2 + (b0+a0)/2

      for i = 1:k
        vals[i] = fun(ts[i])
      end

        coefs = *(u,vals)
#
# Measure the relative error in the trailing coefficients
#
  #    println(coefs)

      dd2 = maxabs(coefs)+1
      dd1 = maxabs(coefs[k-nn+1:k])
      dd = dd1/dd2

      if (dd > eps)

        if (nints0+2 > maxints)

          ier = 4
          return ier
        end

      c0 = (a0+b0)/2

      nints0 = nints0 + 1
      ab0[1,nints0] = c0
      ab0[2,nints0] = b0
      nints0 = nints0 + 1
      ab0[1,nints0] = a0
      ab0[2,nints0] = c0

      else

        if (nints1+1 > maxints)

          ier = 4
          return ier

        end

        nints1 = nints1 + 1
        ab1[1,nints1] = a0
        ab1[2,nints1] = b0

      end

      nints = nints1
      ab = ab1[:,1:nints1]

end


return nints,ab

end
