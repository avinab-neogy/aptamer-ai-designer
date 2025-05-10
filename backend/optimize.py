import pandas as pd
from sklearn.ensemble import RandomForestRegressor

def optimize_candidates(df):
    model = RandomForestRegressor()
    features = pd.get_dummies(df["structure"])
    model.fit(features, df["affinity"])
    df["score"] = model.predict(features)
    return df
