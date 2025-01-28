#define _USE_MATH_DEFINES
#include <vector>
#include <complex>

// #include <iostream>
// #include <fstream>
// #include <cstdlib>
// #include <ctime>

// recursive fft function
template <typename T>
void fft_(std::vector<std::complex<T>>& buf) {
    size_t n = buf.size();
    if (n <= 1) {
        return;
    }

    // Divide: Separate into even and odd indices
    std::vector<std::complex<T>> even(n / 2), odd(n / 2);
    for (size_t i = 0; i < n / 2; i++) {
        even[i] = buf[i * 2];
        odd[i] = buf[i * 2 + 1];
    }

    // Conquer: Recursive calls
    fft_(even);
    fft_(odd);

    // Combine
    for (size_t i = 0; i < n / 2; i++) {
        std::complex<T> t = std::polar((T)1.0, -(T)(2.0 * M_PI * i / n)) * odd[i];
        buf[i] = even[i] + t;           // First half of FFT output
        buf[i + n / 2] = even[i] - t;  // Second half of FFT output
    }
}

template <typename T>
void ifft_(std::vector<std::complex<T>> &buf) {
    size_t m = buf.size();
    for (auto &val : buf) {
        val = std::conj(val);
    }
    fft_(buf);
    for (auto &val : buf) {
        val = std::conj(val) / static_cast<T>(m);
    }
}


template <typename T>
void bluestein_fft(std::vector<std::complex<T>> &buf) {
    size_t n = buf.size();
    size_t m = 1;
    while (m <= 2 * n - 1) {
        m *= 2;
    }

    // std::complex<T> czero((T)0, (T)0);
    // std::vector<std::complex<T>> a(m, czero), b(m, czero), c(m);
    std::vector<std::complex<T>> a(m), b(m), pre(n);

    
    for (size_t i = 0; i < n; i++) {
        pre[i] = std::polar((T)1.0, -(T)M_PI * (T)(i*i) / (T)n);
    }
    
    for (size_t i = 0; i < n; i++) {
        a[i] = buf[i] * pre[i];
    }

    for (size_t i = 0; i < n ; i++) {
        b[i] = std::conj(pre[i]);
    }
    for (size_t i = 1; i < n ; i++){
        b[m - i ] = b[i];
    }
    fft_(a);
    fft_(b);

    for (size_t i = 0; i < m; i++) {
        b[i] = a[i] * b[i];
    }

    ifft_(b);

    for (size_t i = 0; i < n; i++) {
        buf[i] = b[i] * pre[i];
    }
}

template <typename T>
void bluestein_ifft(std::vector<std::complex<T>> &buf){
    for (size_t i = 0; i < buf.size(); i++) {
        buf[i] = std::conj(buf[i]);
    }
    bluestein_fft(buf);
    for (size_t i = 0; i < buf.size(); i++) {
        buf[i] = std::conj(buf[i]) / static_cast<T>(buf.size());
    }
}

template <typename T>
void fft(std::vector<std::complex<T>> &buf) {
    int n = buf.size();
    if (n & (n-1) == 0){
        fft_(buf);
    } else {
        bluestein_fft(buf);
    }
}

template <typename T>
void ifft(std::vector<std::complex<T>> &buf) {
    int n = buf.size();
    if (n & (n-1) == 0){
        ifft_(buf);
    } else {
        bluestein_ifft(buf);
    }
}

// template <typename T>
// void test_fft(int & in_, T &type){
    
//     // Creates int buffer with all 1s
//     std::vector<std::complex<T>> buf(in_);
//     for (int i = 0; i < buf.size(); i++) {
//         // buf[i] = std::complex<T>((T)i, (T)-i);
//         buf[i] = std::complex<T>((T)rand(), (T)rand());
//     }
//     std::vector<std::complex<T>> ans(in_);
//     for (size_t i = 0; i < buf.size(); i++) {
//         ans[i] = buf[i];
//     }

//     clock_t  tnow = clock();
//     asm volatile ("" ::: "memory");
//     fft(buf);
//     ifft(buf);
//     tnow = clock() - tnow;
//     asm volatile ("" ::: "memory");
//     std::cout<<"Time taken: "<<(double)tnow/CLOCKS_PER_SEC<<std::endl;

//     // prints error
//     int errors = 0;
//     for (size_t i = 0; i < buf.size(); i++) {
//         if (std::abs(buf[i] - ans[i])>0.001) {
//             errors++;
//             std::cout << "Error: " << i << " " << buf[i] <<" " << ans[i]<< std::endl;
//         }
//     }
//     if (errors == 0){
//         std::cout << "Correct" << std::endl;
//     }
// }

// // // Takes in input
// int main(int argc, char *argv[]) {
//     if (argc != 1) {
//         // takes in an integer
//         int user_input = std::stoi(argv[1]);
//         double test = 0;
//         test_fft(user_input, test);
//     } 
//     return 0;
// }