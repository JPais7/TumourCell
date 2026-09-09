from src.validation.selection_capabilities import phase43_selection_capabilities
from src.validation.temporal import HistoricalScaler
import numpy as np
def test_unimplemented_selection_apis_are_not_claimed_safe():
 c=phase43_selection_capabilities()
 assert c['feature_selection']=='NOT_IMPLEMENTED'; assert c['threshold_learning']=='NOT_IMPLEMENTED'; assert c['model_selection']=='NOT_IMPLEMENTED'
def test_future_only_feature_is_not_used_by_historical_scaler():
 historical=np.array([[0.,0.],[1.,0.],[2.,0.]])
 future=np.array([[0.,1e12]])
 a=HistoricalScaler().fit(historical); b=HistoricalScaler().fit(np.vstack([historical,future]))
 assert a.mean_[1]==0 and b.mean_[1]>0
