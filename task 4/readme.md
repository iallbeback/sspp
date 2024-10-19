Это 4-ое задание. Программа тестировалась на линуксе внтури виртуальной машине. 
Upd2: теперь идёт сравнение между 4-мя алгоритмами, 


Вот сами результаты:

N = 8, seq mull:  0.000000 secs

N = 8, diag mull:  0.000001 secs

N = 8, AVX mull:  0.000000 secs

N = 8, AVX diag mull:  0.000001 secs

equal

N = 512, seq mull:  0.249488 secs

N = 512, diag mull:  0.276892 secs

N = 512, AVX mull:  0.051912 secs

N = 512, AVX diag mull:  0.055954 secs

equal

N = 1024, seq mull:  2.987789 secs

N = 1024, diag mull:  3.281205 secs

N = 1024, AVX mull:  0.516477 secs

N = 1024, AVX diag mull:  0.780443 secs

equal

N = 2048, seq mull:  88.709071 secs

N = 2048, diag mull:  88.126881 secs

N = 2048, AVX mull:  4.357535 secs

N = 2048, AVX diag mull:  8.505144 secs

equal
