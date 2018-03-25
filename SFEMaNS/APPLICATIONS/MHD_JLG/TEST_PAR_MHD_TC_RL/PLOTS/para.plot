%TEST NAVIER STOKES
%April 23, 2008
%Zahir
% 16 modes, CYL40_TCM_SYM_Z.FEM
% Re=1, 100 pas de temps, dt=0.05
% Source term sur tous les modes.
%Speed up
1 1
2
4
8
16

$END

%Donnes brutes
1 7844.1 % 1 is Bogus since code does not run on one proc only.
2 5063.18999999999960
4 1945.46999999999753
8 1363.86999999999898
16 542.790000000000873


%Speed up
1 1 
2 1.549
4 4.032 
8 5.751
16 14.45

%eff
1 1
2 0.775
4 1.008
8 0.719
16 0.903


%Efficiency
