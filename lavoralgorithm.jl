using IntervalArithmetic, IntervalRootFinding, Distributions

function initialize(S::IntervalBox, f::Function)
    "Returns the array where the partitions of the search space S will be stored together with their lower bounds"
    S1, S2 = bisect(S)
    a = (S1,f(S1).lo)
    b = (S2,f(S2).lo)
    L1 = [a,b]
end

function pick_lowerbound!(A::Array)
    "Extract the box *min_box* from the array A with the less lower bound"
        min_bound = indmin([A[i][2] for i in 1:length(A)])
        min_box = A[min_bound][1]
        deleteat!(A,min_bound) 
        return min_box
end

function gradient_check(x::IntervalBox, gradient::Function)
    "Check that zero is present in each component of the gradient"
    g = gradient(x)
    if g[1].lo < 0. < g[1].hi && g[2].lo < 0. < g[2].hi
        return true
    else
        return false
    end
end

function local_optimization(x::IntervalBox, f::Function)
    "Choose a random point inside the box"
        u1 = Uniform(x[1].lo, x[1].hi)
        u2 = Uniform(x[2].lo, x[2].hi)
        m = [rand(u1),rand(u2)]
        fest = f(m)
end

function minimize(X::IntervalBox, f::Function, g::Function; tolerance = 0.05)
    L1 = initialize(X, f)
    L2 = typeof(L1)(0) #Array where boxes that contain the global minimum with size less than the given tolerance will be stored
    fest = Inf   ##First estimate of the global minumum
    while !isempty(L1)
        X = pick_lowerbound!(L1)
        if gradient_check(X,g)
            fupdate = local_optimization(X, f)
            fest = min(fest, fupdate)
            filter!(x->x[2] <= fest, L1)  ##Reject any box B with lb([f](B)) > fest
            if diam(X) < tolerance
                push!(L2, (X,f(X).lo))
            else
                X1,X2 = bisect(X) 
                push!(L1, (X1, f(X1).lo))
                push!(L1, (X2, f(X2).lo))
            end
        end
    end
    filter!(x->x[2] <= fest, L2) 
    return [L2[i][1] for i in 1:length(L2)], fest  #Boxes that may contain the global minimum (tolerance--dependent) and the global minimum
end
    



