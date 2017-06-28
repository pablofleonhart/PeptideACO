from _acor import ACOR
from prediction import Prediction

prediction = Prediction()
params = ['psi1', 'phi2', 'psi2', 'phi3', 'psi3', 'phi4', 'psi4', 'phi5']
acor = ACOR( prediction.experimental, prediction.modified, params, False, 10 )
acor.evolve()