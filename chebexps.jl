# Written by Cameron Weigel 6/28/2016
#
# ** Might make this a module containing all the necessary Chebyshev code
# rewritten from the Fortran Chebyshev.f90 file, currently contains just
# the chebexps function with an include statement for the necessary chebs
# function located in chebs.jl **
#
#
#
# Construct the n-point Clenshaw-Cutris quadrature with the nodes
#
#   xⱼ = cos ( π * j/(n-1))
#
# on the interval [-1,1], as well as the matrices tkaing the values of an
# n-term Chebyshev expansion at the nodes xⱼ to its coefficients and vice-vera
# and two spectral integration matrices.
#
# Input parameters:
#
#   n - the length of the quadrature formula to construct
#
# Output parameters:
#
#   xs - an array containing the n quadrature nodes
#   whts - an array containing the n quadrature weights
#   u - the (n,n) matrix which takes the values of an n-term Chebyshev
#     expansion at the n quadrature nodes to the n expansion coefficients
#   v - the (n,n) matrix which takes the coefficients of an n-term Chebyshev
#     expansion to it values at the n quadrature nodes
#   aintl - the "left" spectral integration matrix which takes the values
#     of a function f(t) on the Chebyshev nodes to the value of the function
#     g(t) defined via the formula
#
#              t
#       g(t) = ∫ f(u) du
#              a
#
#   ainr - the "right" spectral integration matrix which takes the values
#     of a function f(t) on the Chebyshev nodes to the value of the function
#     g(t) defined via the formula
#
#              t
#       g(t) = ∫ f(u) du
#              a
#
#
#

#
# pick smooth functions [0,1] cos(10pi*x)/sqrt(x)
#

function cheb_exps(n::Int64)



# Variable/Array Definitions


 # Output parameters

  xs = Array(Float64,n)
  whts = Array(Float64,n)
  v = Array(Float64,n,n)
  u = Array(Float64,n,n)
  aintl = Array(Float64,n,n)
  aintr = Array(Float64,n,n)



  pols = Array(Float64,n+1)
  c = zeros(n,n)
  d = Array(Float64,n,n)

  xx = Array(Float64,n,n)


  h = pi/(n-1)

#

  for i = 1:n
    xs[n-i+1] = cos(h*(i-1))
  end

  for i = 1:n

    x = xs[i]
    pols = chebs(x,n-1)

    for j = 1:n

      u[j,i] = pols[j]
      v[i,j] = pols[j]

    end
  end

  u[1,:]=u[1,:]./2
  u[n,:]=u[n,:]./2
  u[:,1]=u[:,1]./2
  u[:,n]=u[:,n]./2

  u = u.*2.0./(n-1)

  # Construct the weights by multiplying u^2 on the left by the
  # integrals of the Chebyshev polynomials


  c[1,1] = 2.0
  for i = 2:2:n-1
    c[1,i+1] = 1.0/(i+1)-1.0/(i-1)
  end

  d = *(c,u)
  whts = d[1,1:n]

# Form the matrix which takes the values of the function f(t) to the values of
# the integral g(t) = f(u)du (from a to t)

# Left spectral integration matrix

  for i in 1:n

    pols = chebs(xs[i],n)

    xx[i,1] = xs[i]
    xx[i,2] = xs[i]^(2)/2

      for j in 3:n

        xx[i,j] = .5*(pols[j+1]/j-pols[j-1]/(j-2))

      end

    end

    for i in 2:n

      xx[i,:] = xx[i,:] - xx[1,:]

    end

    xx[1,:] = 0

    aintl = *(xx,u)

# Right spectral integration matrix

    for i in 1:n

      pols = chebs(xs[i],n)

      xx[i,1] = xs[i]
      xx[i,2] = xs[i]^(2)/2

        for j in 3:n

          xx[i,j] = .5*(pols[j+1]/j-pols[j-1]/(j-2))

        end

      end

      for i in 1:n

        xx[i,:] = xx[i,:] - xx[n,:]

      end

      xx[n,:] = 0

      aintr = *(xx,u)

      

return xs, whts, aintl, aintr, u, v

# more to go here

end
