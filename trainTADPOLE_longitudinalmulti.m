longitudinal = 1;
delay_time = [24 24 24 24];
trainTADPOLE_longitudinal;
save stage1results_longitudinal24 mae cr bbb Yhat Ytrue foldid delay_time longitudinal target_predictions lls

longitudinal = 1;
delay_time = [12 12 24 12];
trainTADPOLE_longitudinal;
save stage1results_longitudinal12 mae cr bbb Yhat Ytrue foldid delay_time longitudinal target_predictions lls

longitudinal = 0;
delay_time = [24 24 24 24];
trainTADPOLE_longitudinal;
save stage1results_longitudinal0 mae cr bbb Yhat Ytrue foldid delay_time longitudinal target_predictions lls

