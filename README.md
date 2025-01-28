# One-file-FFT

As much as this makes me sound like a cock, I've found myself to write code just as fast in C++ as python. In many of my data science projects, I've been forced to do part of my analysis in C++ due to sheer size of the datasets I'm working with. However, I've begun to find it incredibly annoying to constantly be porting between C++ and Python to answer a single question.

I coded this single header file to minimize the amount of time spent porting between the two, because I couldn't be arsed to come up with and remember a consistent file sharing strategy. So here is a FFT that can be included into any project by importing it as a header. I've decided to upload this since surprisingly neither ChatGPT nor copilot were able to successfully implement this for me at the time of writing, so hopefully this will save some time for somebody else who can't be bothered to install a proper fft library.

## FFT
The FFT is implemented using the Bluestein algorithm that calls the 2**n recursive algorithm. There are faster ways of doing this, but it is theoretically O(nlogn), which I have found to be good enough for my purposes. 

### Usage
There are two functions you should use - ```fft``` and ```ifft```. Each takes in ```std::vector<std::complex<T>>```, your choice of ```T```. The code technically will run using ```int```-types but you will find severe rounding errors so I guess stick with floating point. ```ifft``` contains the normalization factor. 

### Theoretical background
As stated previously, there are two algorithms under the hood, the original Cooley-Tukey for power of 2s, and the Bluestein to interpolate. 

#### Part 1: Powers of 2 (Cooley-Tukey)
The core idea is to note that the rows of the DFT matrix are all powers of some number:

$$w_{n,N} = \exp(-2\pi in/N)$$

Since we multiply by our input as a column on the right, we may recast the operation as evaluating a polynomial for each $w_n$:

$$X_n=P(w_{n,N})=\sum_k a_k w_{n,N}^k$$

Where the coefficients $a_k$ is our input vector. We next note that every polynomial can be expressed as a sum of its even and odd terms:

$$P(x) = \sum_k a_k x^k$$
$$=\sum_k a_{2k}x^{2k} + \sum_k a_{2k+1}x^{2k+1}$$
$$=P_e(x^2)+P_o(x^2)x$$

Then from some properties of roots of unity, we note that if $N$ is even we can very easily compute the two polynomials evaluated at $w_{n,N}^2$ from each evaluated at $w_{n,N/2}$. We can repeat this as long as our input is even, or in other words our original input was a power of 2. Analyzing the time complexity, this is $\mathcal O(n\log(n))$. 

#### Part 2: Interpolating to non-powers of 2 (Bluestein)
As previously noted, the previous method only works if each stage is even, or in other words the input length is a power of two.

This is where Bluestein comes in. For clarity, we set $w=w_{-1,N}$. For the hell of it, lets complete the square in the following formula:

$$X_n = \sum_k a_k w^{-kn}$$
$$= \sum_k a_k w^{k^2/2-kn+n^2/2}w^{-(k^2+n^2)/2}$$
$$= w^{-n^2/2}\sum_k (a_k w^{-k^2/2})w^{(n-k)^2/2}$$

Surprisingly a convolution falls out. Thus we have that:

$$X_n = w^{-n^2/2}(a_m w^{-m^2/2} \star w^{m^2/2})_n$$

If one recalls the convolution theorem for finite sequences/periodic sequences, that the convolution of two finite sequences of the same length is the same as the following:

$$(A_m \star B_m)_n = \mathcal F^{-1}[\mathcal F[A_m]\mathcal F[A_m]]$$

Thus, if we pad $a_m w^{-m^2/2}$ with zeros to be a power of two and compute the corresponding $B_m$, we can apply our Cooley-Turkey algorithm. 
