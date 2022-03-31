Technical Caveats:
1.	Small module management: Because the computation time drastically increases with module number, small modules (voxel number < 10) are handled in the beginning: unite small and adjacent modules if their correlation between mean time series is high (>0.8). Then, merge small module with a larger neighboring module, using only the nearest voxels of that larger module to do correlation (by way of bwdist). Select the one with their mean correlation being the highest to merge.
2. Each cycle: permutation of the order of the modules to facilitate convergence.
-	Take modules A, B and C as an example. After partitioning, suppose that A is decomposed into A1 and A2, B into B1 and B2, and C into C1 and C2. A1 shows the highest similarity to B2, while B2 shows the highest similarity to C1. If changing the order does not occur, A1 will always be united with B2, and B2 will have no chance to unite with C1—this could not be the optimal solution. This permutation may facilitate the algorithm to converge.)
3. Each cycle: titrating up strategy; the initial partition is based on the previous partition (lower resolution)
-	 run MOSI from a lower gamma value to provide a good initial guess as 
	a starting point for each cycle.
-	In  other  words,  a  finer  FP  is  derived  from  a  coarser  FP  instead  of  breaking  the original  ROIs  into  many  small  modules  and  repeating  the splitting-unification  cycles  for  the  high 	gamma condition.
-	In this study, for example, the results of gamma 0.70 are based on the results of gamma 0.65, and  the results  of  gamma 0.75  are based  on  the  results  of gamma  0.70  and  so  forth.  This  simple strategy may substantially decrease computation time. 
4. Convergence: 
-	The modular structure does not change
-	Mutual information between the successive partition is higher than mu_cri (0.95)
5. Each hemisphere is handled separately in MOSI
6. Neighborhood consideration in a partition: 
-	Suppose Louvain algorithm divide a module to 2 sub-modules: sub-A and sub-B. If sub-A has 2 disconnected components, then sub-A is counted as two modules. 

Matlab code structure: AutoFP.m 
1.	The code (only 426 lines) has been uploaded to github: https://github.com/dwleeibru/MOSI. 
2.	Two small functions to facilitate the analytic flow: line 388-426
3.	Between line21 and line386 are the loops to process multiple subjects
4.	Between line27 and line383 are the loops to process multiple gamma values
5.	Input: from AFNI, line12 to line 17. The file content is organized as: coordinate (I,j,k) then time series.
- If a ROI has 50 voxels, and there are 250 fMRI scans, the size of the content is 50*253 (250+3=253) 
6.	Input setup: starting from previous gamma value, line28 to line 36 (Technical Caveat 3)
7.	Output: line 382, saved as “parti_' num2str(gamma) '_' side '.mat”, side is either left or right (Technical Caveat 5)
8.	For first cycle, initial guess is based on anatomy-based partition from AFNI: line42 to line71
9.	Parameters: line72 to line76
10. Between line84 and line380 are MOSI loops
11. Between line86 to line137 is the splitting stage of MOSI, in which Technical Caveat 6 is handled between line101 to line131
12. Between line140 and 280, small module condition is handled.
-	Unification of neighboring small and small modules, line152 to line207
-	Unification of neighboring small and larger modules, line210 to line285
-	Technical Caveat 1 is handled between line257 to line264
13. Between line288 to line371 is the unifying stage of MOSI, with relative distance calculation handled between line301 to line316
14. Line 377 handles the permutation of partitions, (Technical Caveat 2)
15. Between line351 and line375, convergence is examined (Technical Caveat 4)

Clarification: in general, unifying stage uses whole brain information (relative distance), but for small modules, correlation is applied.



