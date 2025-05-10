from ViennaRNA import fold

def predict_secondary_structure(sequence):
    structure, mfe = fold(sequence)
    return structure, mfe
