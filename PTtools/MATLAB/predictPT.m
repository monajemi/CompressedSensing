% Predict Phase Transition of Ell1 minimization for Gaussian matrices 
% Author: Hatef Monajemi 
% Date: June 12 2012
% delta : scalar determining the shape factor (undersampling)
% flag: 'Real', 'Cplex', 'Pos' , 'Bnd', 'Q', 'O'
% Citation: "Deterministic matrices matching the compressed sensing phase transitions of Gaussian random matrices.", 
% Monajemi et al., 2013 http://www.pnas.org/content/110/4/1181 

function eps_0 = predictPT(delta, flag)

PredPT = load('PredPT.mat');

switch flag 

case {'Real','R'}
eps_0 = interp1( PredPT.delta_Real, PredPT.eps_Real, delta,'PCHIP');

case {'Cplex','C'}
eps_0 = interp1( PredPT.delta_Cplex, PredPT.eps_Cplex, delta,'PCHIP');


case {'Pos', 'pos','R+'}
eps_0 = interp1( PredPT.delta_Pos, PredPT.eps_Pos, delta,'PCHIP');


case {'Bnd','bnd'}
eps_0 = interp1( PredPT.delta_Bnd, PredPT.eps_Bnd, delta,'PCHIP');


case 'Q'
eps_0 = interp1( PredPT.delta_Q, PredPT.eps_Q, delta,'PCHIP');
case 'O'
eps_0 = interp1( PredPT.delta_O, PredPT.eps_O, delta,'PCHIP');

otherwise
error('unknown flag')
end
