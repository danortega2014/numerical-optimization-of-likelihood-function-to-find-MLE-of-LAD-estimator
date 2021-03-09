# numerical-optimization-of-likelihood-function-to-find-MLE-of-LAD-estimator

For econometrics, Maximum Likelihood Estimation is a very popular estimating technique, in which not all estimators have closed form solutions. Maximum likelihood estimates is the estimate which is the maximum of the likelihood of receiving a specific parameter estimate given the data.  Meaning the data is most probable in the statistical model given a specific point in the parameter space of the model. The data is also assumed follow a certain statistical model.
 The example I am doing is the MLE of the LAD estimator. Which can be represented by this parametic model:

![image](https://user-images.githubusercontent.com/64437206/110543875-7d710f80-80f0-11eb-94bb-79b6772f9449.png)


I will be maximizing the LAD log likelihood function to find the maximum likelihood estimate. 

![image](https://user-images.githubusercontent.com/64437206/110543912-88c43b00-80f0-11eb-9108-881533a16fe0.png)

This log of the likeihood function is derived from taking the log of the likelihood pdf of y. Which in this case is:

![image](https://user-images.githubusercontent.com/64437206/110544043-be692400-80f0-11eb-9273-5e08baf66069.png)



Since the log likelihood function is an absolute value function, it’s derivative cannot be taken hence no closed form solution. We’re trying to find the max likelihood of theta given x,y. I managed to construct an algorithm, in which I have created random samples of X and Y using the random and distributions packages.  
```
using Random, Distributions

x= rand(Normal(),100)
y= rand(Normal(),100)
```

The only parameter that must be maximized is beta. I created two functions, first the regular log likelihood function and the second the first order conditions for the log likelihood function. Again, since the likelihood function has abs value the derivative cannot be taken, but an approximation can be created using Newton’s difference quotient.
```
function loglike(Beta)
    llbeta  = -length(x)*log(2)-sum(abs.(y .-x*Beta))
    return -llbeta
end

#newton's difference quotient
function derivat(Beta)
    h = .00001
    g2 = (loglike(Beta + h) - loglike(Beta) ) / h
    return g2
end
```

With this function, I run it through the bisection function, and I run the log likelihood function through the optimize function, my results are in Figure.In which a maximum for theta is .0202 for both techniques.
