Это задание №6
Вот результаты:

qsort time:	0.146112 seconds

(это я запускал с p=3)
Threads:	3 
Parallel:	0.079873 seconds
Par/qsort:	54.67%
Arr equal

(а здесь стандартные тесты)
Threads:	1
Parallel:	0.143930 seconds
Par/qsort:	98.51%
Arr equal

Threads:	2
Parallel:	0.073517 seconds
Par/qsort:	50.32%
Arr equal

Threads:	4
Parallel:	0.056342 seconds
Par/qsort:	38.56%
Arr equal

Threads:	8
Parallel:	0.047070 seconds
Par/qsort:	32.22%
Arr equal

Threads:	16
Parallel:	0.036724 seconds
Par/qsort:	25.13%
Arr equal


Results Table:
Threads	Time (s)	Speedup	Efficiency
1	0.143930	1.015163	1.015163
2	0.073517	1.987473	0.993737
4	0.056342	2.593304	0.648326
8	0.047070	3.104144	0.388018
16	0.036724	3.978631	0.248664


![image](https://github.com/user-attachments/assets/82dd78a3-bae3-44ac-99af-f14baa213736)
![image](https://github.com/user-attachments/assets/0d659812-337a-4893-809c-de99bc1db638)
