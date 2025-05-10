import random

def generate_initial_candidates(num=20, length=20):
    bases = ['A', 'T', 'C', 'G']
    return [''.join(random.choices(bases, k=length)) for _ in range(num)]
