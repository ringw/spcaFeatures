# MIT License
# 
# Copyright (c) 2018 Lauren Berk
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

using LinearAlgebra


function projectMatrix(sig, xVal)
	I = 1.0*Matrix(LinearAlgebra.I,size(sig,1),size(sig,1))
	xx = (I - xVal*xVal')
	A3 = xx*sig*xx
    return A3

end


function multiOptimalSPCA(prob, k, ncomp, searchCap=500000)
	myprob = problem(copy(prob.data),copy(prob.Sigma))
	st = time()
	all_x = zeros(ncomp, size(myprob.Sigma,1))

	for rnd = 1:ncomp
		~, xVal, ~, ~, ~, ~, ~ = branchAndBound(myprob, k, outputFlag=0, searchCap=searchCap) ;
		all_x[rnd,:] = xVal

		myprob.Sigma = projectMatrix(myprob.Sigma, xVal)
	end

	return all_x, time()-st
end


function adjVarExplained(Sigma, xVal)

	depth = size(xVal,1)
	varexplained = zeros(depth)
	sig = copy(Sigma)

	for component = 1:depth
		mypc = xVal[component,:]
		mypc = mypc/norm(mypc)
		varexplained[component] = mypc'*sig*mypc
		sig = projectMatrix(sig, mypc)
	end

	return varexplained
end
