function m = matrix_vector_multi(x,eigen)
 %computes matrix-vector multiplication where matrix is Circulant
 
  m = real(fft(ifft(x).*eigen));