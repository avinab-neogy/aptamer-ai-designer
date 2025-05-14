import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import NotFittedError
import streamlit as st

def optimize_candidates(df):
    """Optimize aptamer candidates using machine learning"""
    try:
        # Check for required columns
        required_columns = {'sequence', 'structure', 'affinity'}
        missing_cols = required_columns - set(df.columns)
        if missing_cols:
            raise KeyError(f"Missing required columns: {missing_cols}")
        
        # Feature engineering
        df = df.copy()
        df['gc_content'] = df['sequence'].apply(
            lambda s: (s.count('G') + s.count('C')) / len(s) if len(s) > 0 else 0
        )
        df['structure_complexity'] = df['structure'].str.count(r'[()]')
        
        # Create features
        features = pd.get_dummies(df["structure"], prefix="struct")
        X = pd.concat([df[['gc_content', 'structure_complexity']], features], axis=1)
        y = df["affinity"]
        
        # Train model
        model = RandomForestRegressor(n_estimators=100, random_state=42)
        model.fit(X, y)
        
        # Predict scores
        df["score"] = model.predict(X)
        
        return df.sort_values("score", ascending=False)
    
    except KeyError as e:
        st.error(f"Data format error: {str(e)}")
        return pd.DataFrame()
    
    except ValueError as e:
        st.error(f"Invalid data: {str(e)}")
        return df
    
    except NotFittedError as e:
        st.error("Model failed to train. Check input data.")
        return df
    
    except Exception as e:
        st.error(f"Optimization failed: {str(e)}")
        return df
