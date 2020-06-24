# Statistical-Functions


Random walk Metropolis-Hastings algorithm and Iterated weighted least squares functions (Statistical Computing MATH6166, University of Southampton):

**Random walk Metropolis-Hastings algorithm** - one of the most common Markov chain Monte Carlo algorithms  <br/>
**Iterated weighted least squares** - is used to solve certain optimization problems by an iterative method in which each step involves solving a weighted least squares problem

### Concept of code:

**First part:**
1. Random walk Metropolis-Hastings algorithm (rwmha) function
2. Using sample of size = 10000 to calculate rwmha variance
3. Phi value vs acceptance rate for rwmha where theoretically optimal value is 0.234
<img width="760" alt="Screen Shot 2020-06-24 at 20 31 36" src="https://user-images.githubusercontent.com/37827791/85619738-9ddc6500-b65a-11ea-8292-eef745b3f59c.png">
4. Using sample of size = 10000 with optimal phi value to calculate rwmha variance
<br/>
--------------------------------------------------------------------------------------------------------------------------------

**Second part:**
1. Iterated weighted least squares (iwls) function to return vector of binary responses and a model matrix will return the maximum likelihood estimates of — for a logistic regression
2. Using given dataset to use created iwls function and to compare with existing library of logistic regression
