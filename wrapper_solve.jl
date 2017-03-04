#
#  Construct a  function q with chebadap
#
include("chebs.jl")
include("chebexps.jl")
include("chebadap.jl")
include("chebpw_eval.jl")
include("print0.jl")


function wrapper_solve(fun::Function,k::Int64,epsil::Float64,a::Float64,b::Float64)


(xs,whts,aintl,aintr,u,v) = cheb_exps(k)



(nintsin,abin)            = cheb_adap(fun,epsil,k,xs,u,a,b)
qvals                     = Array(Float64,k,nintsin);

for int=1:nintsin
a0 = abin[1,int]
b0 = abin[2,int]
for i=1:k
t = xs[i]*(b0-a0)/2 + (b0+a0)/2
qvals[i,int] = fun(t)
end
end

print0("qvals =",qvals)

#
#  Call the kummer phase function
#

println(abin)

ccall( (:__wrapper_MOD_phase_function,"./libkummer.so"), Ptr{Void},
  (Ptr{Int64},Ptr{Int64},Ptr{Float64},Ptr{Float64}),
  &k, &nintsin, abin, qvals)



nints = ccall( (:__wrapper_MOD_phase_intervals,"./libkummer.so"), Int64, ())

println("k=",k)
println("nints=",nints)



ab     = Array(Float64,2,nints)
alpha  = Array(Float64,k,nints)
alphap = Array(Float64,k,nints)
alphapp = Array(Float64,k,nints)




 ccall( (:__wrapper_MOD_phase_data,"./libkummer.so"), Void,
   (Ptr{Int64},Ptr{Int64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
   &nints,&k,ab,alpha,alphap,alphapp)

   print0("ab = ",ab)
   print0("alpha = ",alpha)
   print0("alphap = ",alphap)
   print0("alphapp =",alphapp)




#c_1 = y_init*sqrt(alphap[1,1])
#c_2 = y_prime_init/sqrt(alphap[1,1])+y_init*(alphapp[1,1]/(2*alphap[1,1]^(1.5)))

#a_1 = sqrt(c_1^2+c_2^2)
#a_2 = atan2(c_1,c_2)

#if (a_2 > pi)
#  a_1 = -a_1
#  a_2 = a_2 - pi
#end

#if (a_2 < 0)
#  a_1 = -a_1
#  a_2 = a_2 + pi
#end


#println("ts = ",ts)



#  a  = chebpw_eval(nints,ab,k,xs,alpha,ts)
#  ap = chebpw_eval(nints,ab,k,xs,alphap,ts)
#  val = a_1*sin(a+a_2)/sqrt(ap)

#  println("val =", val)

#  print0("ab =",ab)




  return nints,ab,k,xs,alpha,alphap,alphapp

end
