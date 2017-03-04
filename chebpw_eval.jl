include("print0.jl")

function chebpw_eval(nints::Int64,ab::Array{Float64,2},k::Int64,xscheb::Array{Float64,1},vals::Array{Float64,2},x::Float64)

#
# Variable Declarations
#

#  ninters::Int64
#  intl::Int64
#  intr::Int64
#  a::Float64
#  b::Float64
#  xx::Float64
#  sum1::Float64
#  sum2::Float64
#  dd1::Float64
#  dd::Float64
#  diff::Float64
  val::Float64
  int = 0
  int0 = 0
  diff = 0.0
#
#
  intl::Int64
  ninters = 8
  intl = 1
  intr = nints

#  int = fld((intl+intr),2)

for iter = 1:ninters
int = div((intl+intr),2)
c   = ab[1,int]
if (x > c)
intl = int
else
if (int > 1)
intr = int-1
end
end
end

#
println("")
for int = intl:intr
  a = ab[1,int]
  b = ab[2,int]

  println(a,"   ",b,"    ",x)

  if (x <= b && x >= a)
    int0 = int
  end
end

if (int0 == 0)
println(int0)
print0("ab =",ab)
println(x)
exit(0)
end
  #

  a     = ab[1,int0]
  b     = ab[2,int0]

  #

  xx    = (2*x-(b+a))/(b-a)

  sum1  = 0.0
  sum2  = 0.0

  dd1   = 1.0

for i = 1:k

  dd     = 1.0

  if (i == 1) || (i == k)
    dd    = 0.5
    diff  = xx-xscheb[i]

    #

    if (abs(diff) < eps())
      val = vals[i,int0]
      return val
    end
  end

  dd    = (dd1*dd)/diff
  dd    = -dd1
  sum1  = sum1+dd*vals[i,int0]
  sum2  = sum2+dd
  dd    = -dd

end
  val = sum1/sum2
return val
end
